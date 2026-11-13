/*
 * This source file is part of the Alexandria Chemistry Toolkit.
 *
 * Copyright (C) 2026
 *
 * Developers:
 *             Mohammad Mehdi Ghahremanpour,
 *             Julian Marrades,
 *             Marie-Madeleine Walz,
 *             Paul J. van Maaren,
 *             David van der Spoel (Project leader)
 *
 * This program is free software; you can redistribute it and/or
 * modify it under the terms of the GNU General Public License
 * as published by the Free Software Foundation; either version 2
 * of the License, or (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this program; if not, write to the Free Software
 * Foundation, Inc., 51 Franklin Street, Fifth Floor,
 * Boston, MA  02110-1301, USA.
 */
/*! \internal \brief
 * Topologies and their conformations resident on the GPU, see
 * createSystem in gpu.h. They are copied once, after which only changed
 * parameters are copied to the GPU and energies back, such that the
 * communication between the processor and the GPU is small.
 *
 * The energy of each interaction in each slot is computed by one
 * work-item, then summed for each list of interactions in the order of the
 * interactions. The forces of each particle in each slot are computed by
 * one work-item from all its interactions, in the order of the lists and of
 * the interactions, separately for the interactions without shells, which
 * are the fixed forces when relaxing the shells. No atomic operations are
 * used, such that the results are the same in every run.
 * \author David van der Spoel <david.vanderspoel@icm.uu.se>
 */
#include "gpu.h"

#include <algorithm>
#include <cfloat>
#include <cstdint>
#include <cstdio>
#include <cstdlib>
#include <map>
#include <memory>
#include <stdexcept>
#include <string>
#include <vector>

#include <sycl/sycl.hpp>

#include "act/forcefield/setrecipe.h"
#include "act/gpu/pairrows.h"
#include "sycl/bondedkernels.h"
#include "sycl/context.h"
#include "sycl/pairkernels.h"

namespace alexandria
{

namespace gpu
{

namespace
{

//! The work items of a wavefront on AMD, a multiple of a warp on NVIDIA
constexpr uint32_t wavefrontSize = 64;

/*! \brief The number of interactions whose energies one work item sums
 * The energies of each chunk of a list in each slot are summed in order in
 * double precision, and then those of the chunks in order, such that the
 * energies do not depend on how the work items are scheduled, without
 * storing the energy of each interaction in each slot.
 */
constexpr uint32_t energyChunk = 32;

//! A block in the device arrays
struct DeviceSystemBlock
{
    //! Start of the coordinates and forces, rows of stride
    uint64_t xOffset;
    //! Start of the particles in the arrays of all particles of all blocks, such as the charges
    uint64_t atomOffset;
    //! Start of the particles in atomStart and shellStart, that of the block whose
    //! interactions this block shares, see SystemBlock::topologyOf
    uint64_t entryAtomOffset;
    //! What to add to the list of an interaction (pairList, bondList) for the list of this block
    uint32_t listShift;
    //! Start of the particles times slots
    uint64_t atomWorkOffset;
    //! Start of the slots of all blocks
    uint64_t slotOffset;
    //! Start of the shell rows times slots
    uint64_t shellRowWorkOffset;
    //! Start of the force and moving shells times slots
    uint64_t shellWorkOffset;
    //! Start of the force shells times slots
    uint64_t forceShellPartOffset;
    //! Start of the moving shells times slots
    uint64_t movingShellPartOffset;
    //! Start of the shell-parent pairs times slots
    uint64_t resetWorkOffset;
    //! The first force shell, moving shell, fragment and shell row
    uint32_t firstForceShell, nforceShell, firstMovingShell, nmovingShell;
    uint32_t firstFragment, nfragment, firstShellRow, nshellRow;
    //! The first shell-parent pair and their number
    uint32_t firstReset, nreset;
    //! Number of particles
    uint32_t natom;
    //! Distance between the rows
    uint32_t stride;
    //! Number of slots
    uint32_t nslot;
    //! Number of slots rounded up to a wavefront, for the forces on the shells in relaxBlocks
    uint32_t nslotPad;
    //! Number of energy arrays
    uint32_t nenergy;
};

//! A list of interactions in the device arrays
struct DeviceList
{
    //! Start of the chunks of energyChunk interactions times slots
    uint64_t workOffset;
    //! Start of the lists times slots, for the sums of the energies
    uint64_t listSlotOffset;
    //! Start of the interactions in the pair or bonded arrays
    uint64_t entryOffset;
    //! Start of the parameter table in the tables of all lists, see updateParameters
    uint64_t tableOffset;
    //! Start of the parameters of its interactions in params_, if ownParams,
    //! otherwise where those of the next list start, see entryParams
    uint64_t paramOffset;
    //! Number of interactions
    uint64_t n;
    //! The block
    uint32_t block;
    //! The potential
    uint32_t kernel;
    //! Whether bonded
    uint32_t bonded;
    //! Whether with shells
    uint32_t shells;
    //! The energy array
    uint32_t energyIndex;
    //! Number of energy terms, the others are zero
    uint32_t nterm;
    //! Whether its parameters are in params_, otherwise in its parameter table
    uint32_t ownParams;
};

//! A virtual site in the device arrays
struct DeviceVsite
{
    //! Start of the sites times slots
    uint64_t workOffset;
    //! The block
    uint32_t block;
    //! The site
    VsiteEntry entry;
};

//! The particles of a pair interaction, its parameters are elsewhere, see entryParams
struct DevicePair
{
    uint32_t i, j;
};

//! The particles of a bonded interaction, three for angles, see entryParams
struct DeviceBond
{
    uint32_t a[4];
};

/*! \brief The interactions of a system, which the blocks of the replicas share
 * The lists of each block have their own parameters: a row of their table,
 * or for Coulomb, which follows from the charges, their own for each pair.
 */
struct EntryData
{
    const DevicePair *pairs;
    const DeviceBond *bonds;
    const DeviceList *lists;
    //! The list of each interaction, of the block whose interactions are shared
    const uint32_t   *pairList;
    const uint32_t   *bondList;
    //! The row of each interaction in the parameter table of its list
    const uint32_t   *pairIndex;
    const uint32_t   *bondIndex;
    const float      *tables;
    const float      *params;
};

//! The eight parameters of interaction idx of list l, with the rows of index
inline const float *entryParams(const EntryData &d, const DeviceList &l, uint64_t idx, const uint32_t *index)
{
    return l.ownParams
        ? d.params + 8*(l.paramOffset + idx - l.entryOffset)
        : d.tables + 8*(l.tableOffset + index[idx]);
}

//! The element whose range holds index, by bisection
template <typename T>
inline uint32_t findRange(const T *a, uint32_t n, uint64_t index, uint64_t T::*start)
{
    uint32_t lo = 0, hi = n;
    while (hi - lo > 1)
    {
        uint32_t mid = (lo + hi)/2;
        if (a[mid].*start <= index)
        {
            lo = mid;
        }
        else
        {
            hi = mid;
        }
    }
    return lo;
}

/*! \brief As findRange, starting from the range of the first work item of the wavefront
 * Most work items of a wavefront are in the same range, so this needs one
 * or a few loads instead of a binary search, see makeRangeHints.
 * \param[in] a     The ranges
 * \param[in] n     Their number
 * \param[in] index The work item
 * \param[in] start The start of each range
 * \param[in] hint  The range of the first work item of each wavefront
 * \return The last range that starts at or before index
 */
template <typename T>
inline uint32_t findRangeFrom(const T *a, uint32_t n, uint64_t index, uint64_t T::*start, const uint32_t *hint)
{
    uint32_t lo = hint[index/wavefrontSize];
    while (lo + 1 < n && a[lo + 1].*start <= index)
    {
        lo++;
    }
    return lo;
}

//! A coordinate difference in nm, the exact difference of the integers rounded once
inline float differenceNm(int64_t a, int64_t b)
{
    return static_cast<float>(a - b)*static_cast<float>(coordinateUnitNm);
}

//! The interaction of a particle in the forces, as atomEntries holds them
enum : uint64_t
{
    //! The role of the particle in the interaction, in the lowest two bits
    roleMask   = 3,
    //! Bonded instead of pair, in the third bit
    bondedFlag = 4,
    //! The interaction, from the fourth bit on
    indexShift = 3
};

//! The state of the relaxation of the shells of a fragment in a slot, as in relaxation.cpp
struct FragmentState
{
    double step;
    double msPrev;
    double dxfPrev;
    int    moving;
    double fmaxRest;
    //! The mean square force against which progress is measured, see shellStallIterations
    double reference;
    //! Steps since the mean square force was halved
    int    stall;
};

//! The difference of two coordinates with wrap-around, see coordinateDifference
inline int64_t difference(int64_t a, int64_t b)
{
    return static_cast<int64_t>(static_cast<uint64_t>(a) - static_cast<uint64_t>(b));
}

//! The coordinate in nm, rounded to float as gmx::RVec stores it
inline float toNm(int64_t c)
{
    return static_cast<float>(static_cast<double>(c)*coordinateUnitNm);
}

//! A coordinate in nm in integer units, without the range check
inline int64_t toUnits(double x)
{
    return static_cast<int64_t>(sycl::rint(x*coordinateScaleNm));
}

//! Move a coordinate, see moveCoordinate in coordinate.h, valid is false for a displacement that is not finite
inline int64_t move(int64_t c, double d, bool *valid)
{
    if (!sycl::isfinite(d))
    {
        *valid = false;
        return c;
    }
    const double limit = 0.999*128.0;
    int64_t step = (sycl::fabs(d) < 1) ? toUnits(d) : toUnits(sycl::clamp(d, -limit, limit));
    int64_t big  = toUnits(limit);
    if (step > 0 && c > INT64_MAX - step)
    {
        return big;
    }
    if (step < 0 && c < INT64_MIN - step)
    {
        return -big;
    }
    return c + step;
}

//! The displacement factor along the force of a shell with a quartic spring, see calc_fscale
inline double fscale(double k1, double k2, double F0)
{
    if (F0 == 0)
    {
        return 0;
    }
    if (k2 == 0)
    {
        return -1/k1;
    }
    const double third = 1.0/3.0;
    const double three = sycl::pow(3.0, third);
    const double sqrt3 = sycl::sqrt(3.0);
    double k12    = k1*k2;
    double F0k22  = F0*k2*k2;
    double throot = sycl::pow(-9*F0k22 + sqrt3*sycl::sqrt(k12*k12*k12 + 27*F0k22*F0k22), third);
    double k      = (-three*three*k12 + three*throot*throot)/(6*k2*throot);
    return k/F0;
}

/*! \brief The force of one interaction on one of its particles
 * \param[in]  code  The interaction and the role of the particle, see atomEntries
 * \param[in]  x     The coordinates of the slot
 * \param[in]  st    The stride
 * \param[in]  d     The interactions
 * \param[in]  listShift The DeviceSystemBlock::listShift of the block
 * \param[out] fa    The force
 * \return whether the interaction has shells
 */
inline bool entryForce(uint64_t code, const int64_t *x, uint64_t st, const EntryData &d,
                       uint32_t listShift, float *fa)
{
    uint64_t idx  = code >> indexShift;
    uint32_t role = code & roleMask;
    if (code & bondedFlag)
    {
        const DeviceBond &be = d.bonds[idx];
        const DeviceList &l  = d.lists[d.bondList[idx] + listShift];
        float r[3][3], f[4][3];
        for (int k = 1; k < 4; k++)
        {
            for (int m = 0; m < 3; m++)
            {
                r[k-1][m] = differenceNm(x[(3*be.a[k]+m)*st], x[(3*be.a[0]+m)*st]);
            }
        }
        (void) bondedInteraction(static_cast<BondedKernel>(l.kernel), entryParams(d, l, idx, d.bondIndex), r, f);
        for (int m = 0; m < 3; m++)
        {
            fa[m] = f[role][m];
        }
        return l.shells != 0;
    }
    const DevicePair &pe = d.pairs[idx];
    const DeviceList &l  = d.lists[d.pairList[idx] + listShift];
    float dx[3], rsq = 0;
    for (int m = 0; m < 3; m++)
    {
        dx[m] = differenceNm(x[(3*pe.i+m)*st], x[(3*pe.j+m)*st]);
        rsq += dx[m]*dx[m];
    }
    float eu[pairEnergyTerms];
    float fscal = pairInteraction(static_cast<PairKernel>(l.kernel), entryParams(d, l, idx, d.pairIndex), rsq, eu);
    for (int m = 0; m < 3; m++)
    {
        // As the force of the pair on its first particle
        float fp = rsq > 0 ? fscal*dx[m] : 0.0f;
        fa[m] = role == 0 ? fp : -fp;
    }
    return l.shells != 0;
}

/*! \brief The force of an interaction on one of its particles, as entryForce,
 * with the coordinates of that particle given
 *
 * For pairs only the coordinates of the other particle are loaded. The
 * difference is exact before it is rounded once, so the result is that
 * of entryForce.
 */
inline void entryForceOn(uint64_t code, const int64_t own[3], const int64_t *x, uint64_t st,
                         const EntryData &d, uint32_t listShift, float *fa)
{
    if (code & bondedFlag)
    {
        (void) entryForce(code, x, st, d, listShift, fa);
        return;
    }
    uint64_t          idx     = code >> indexShift;
    uint32_t          role    = code & roleMask;
    const DevicePair &pe      = d.pairs[idx];
    const DeviceList &l       = d.lists[d.pairList[idx] + listShift];
    uint32_t          partner = role == 0 ? pe.j : pe.i;
    float dx[3], rsq = 0;
    for (int m = 0; m < 3; m++)
    {
        // The difference from the partner, which is that of entryForce for role 0
        dx[m] = differenceNm(own[m], x[(3*partner+m)*st]);
        rsq += dx[m]*dx[m];
    }
    float eu[pairEnergyTerms];
    float fscal = pairInteraction(static_cast<PairKernel>(l.kernel), entryParams(d, l, idx, d.pairIndex), rsq, eu);
    for (int m = 0; m < 3; m++)
    {
        // For role 1, entryForce gives -(fscal*(-d)), which is the same
        fa[m] = rsq > 0 ? fscal*dx[m] : 0.0f;
    }
}

} // namespace

class System
{
public:
    explicit System(sycl::queue *queue) : queue_(queue) {}

    ~System()
    {
        for (void *p : allocations_)
        {
            sycl::free(p, *queue_);
        }
    }

    //! The interactions and their parameters, for the kernels
    EntryData entryData() const
    {
        return { pairs_, bonds_, lists_, pairList_, bondList_, pairIndex_, bondIndex_, tables_, params_ };
    }

    /*! \brief Allocate device memory of n elements, freed with the system
     * \throws OutOfMemory if the GPU has no memory left, since the
     *         allocation then returns no memory, which the operations on
     *         the GPU would write to
     */
    template <typename T>
    T *allocate(size_t n)
    {
        size_t bytes = std::max<size_t>(n, 1)*sizeof(T);
        T     *p     = sycl::malloc_device<T>(std::max<size_t>(n, 1), *queue_);
        if (!p)
        {
            auto        device = queue_->get_device();
            std::string name   = device.get_info<sycl::info::device::name>();
            double      total  = device.get_info<sycl::info::device::global_mem_size>();
            throw OutOfMemory(name, total, allocatedBytes_, bytes);
        }
        allocations_.push_back(p);
        allocationBytes_.push_back(bytes);
        allocatedBytes_ += bytes;
        return p;
    }

    //! Free memory of allocate before destroySystem, and set the pointer to nullptr
    template <typename T>
    void release(T *&p)
    {
        auto a = std::find(allocations_.begin(), allocations_.end(), static_cast<void *>(p));
        if (p && a != allocations_.end())
        {
            sycl::free(p, *queue_);
            size_t i         = a - allocations_.begin();
            allocatedBytes_ -= allocationBytes_[i];
            allocationBytes_.erase(allocationBytes_.begin() + i);
            allocations_.erase(a);
        }
        p = nullptr;
    }

    //! Copy n elements to the device
    template <typename T>
    void upload(T *d, const T *h, size_t n)
    {
        if (n > 0)
        {
            queue_->memcpy(d, h, n*sizeof(T));
        }
    }

    /*! \brief A copy of a mask on the device, such as which blocks to compute
     * Each mask is copied once and then used again, since a system has
     * only a few, such that no copy is needed for each operation, and the
     * operations need not wait for a copy from host memory.
     * \param[in] h The mask on the host
     * \param[in] n Its length
     * \return the mask on the device
     */
    const char *deviceMask(const char *h, size_t n)
    {
        // The masks the GPU made, see chainMask
        if (chainMasks_ && h >= chainMasks_ && h < chainMasks_ + 3*nblock_)
        {
            return h;
        }
        std::string key(h, n);
        auto        found = masks_.find(key);
        if (found != masks_.end())
        {
            return found->second;
        }
        char *d = allocate<char>(std::max<size_t>(n, 1));
        if (n > 0)
        {
            queue_->memcpy(d, h, n).wait();
        }
        masks_.emplace(key, d);
        return d;
    }

    //! The masks on the device, see deviceMask
    std::map<std::string, char *>   masks_;

    sycl::queue                    *queue_;
    std::vector<void *>             allocations_;
    //! The size of each allocation in bytes, and their sum, see allocate
    std::vector<size_t>             allocationBytes_;
    size_t                          allocatedBytes_ = 0;
    std::vector<DeviceSystemBlock>  blocksHost_;
    std::vector<DeviceList>         listsHost_;
    //! The first list of each block
    std::vector<size_t>             firstList_;
    uint32_t                        nblock_ = 0, nlist_ = 0, nvsite_ = 0;
    uint64_t                        nx_ = 0, npair_ = 0, nbond_ = 0, natomWork_ = 0;
    uint64_t                        nwork_ = 0, nlistSlot_ = 0, nenergyAll_ = 0, nvsiteWork_ = 0;
    DeviceSystemBlock              *blocks_   = nullptr;
    DeviceList                     *lists_    = nullptr;
    DeviceVsite                    *vsites_   = nullptr;
    int64_t                        *x_        = nullptr;
    float                          *f_        = nullptr;
    float                          *fFixed_   = nullptr;
    DevicePair                     *pairs_    = nullptr;
    DeviceBond                     *bonds_    = nullptr;
    //! The parameters of the interactions of the lists without parameter table, see entryParams
    float                          *params_   = nullptr;
    uint64_t                        nparam_   = 0;
    //! The parameter tables of all lists and their row for each interaction, see updateParameters
    float                          *tables_     = nullptr;
    uint32_t                       *pairIndex_  = nullptr;
    uint32_t                       *bondIndex_  = nullptr;
    uint64_t                        ntable_     = 0;
    //! The charges and widths of all particles, four each
    double                         *charges_    = nullptr;
    uint64_t                        ncharge_    = 0;
    //! The recipes of the rows of the tables, see createRowRecipes
    SetRecipe                      *recipes_       = nullptr;
    RowRecipe                      *rowRecipes_    = nullptr;
    double                         *recipeFactors_ = nullptr;
    double                         *recipeValues_  = nullptr;
    //! The parameters of the sets of each replica, maxSetParameters each
    double                         *setValues_     = nullptr;
    uint64_t                        nset_ = 0, nrowRecipe_ = 0, nrecipeValue_ = 0, nrecipeReplica_ = 0;
    //! The recipes of the moving shells and the virtual sites, see createShellVsiteRecipes
    ShellRecipe                    *shellRecipes_  = nullptr;
    VsiteRecipe                    *vsiteRecipes_  = nullptr;
    uint64_t                        nshellRecipe_ = 0, nvsiteRecipe_ = 0;
    //! Whether updateRecipeValues gave values of which the rows were not made yet
    bool                            recipeValuesChanged_ = false;
    //! The virtual sites, as on the device
    std::vector<DeviceVsite>        vsitesHost_;
    //! The charge generation, see createAcm
    AcmMolecule                    *acmMolecules_ = nullptr;
    AcmFragment                    *acmFragments_ = nullptr;
    uint32_t                       *acmParticles_ = nullptr;
    uint32_t                       *acmBonds_     = nullptr;
    int32_t                        *acmSymmetry_  = nullptr;
    double                         *acmValues_    = nullptr;
    double                         *acmX_         = nullptr;
    uint32_t                       *acmParticleValues_ = nullptr;
    uint32_t                       *acmBondValues_     = nullptr;
    double                         *acmBondFactors_    = nullptr;
    uint32_t                       *acmFixedValues_    = nullptr;
    uint64_t                        nacmValue_ = 0;
    //! Whether the values or coordinates changed since the charges were generated
    bool                            acmChanged_ = false;
    //! The structures and virtual sites of the charge generation, see createAcmVsites
    int64_t                        *acmRefX_       = nullptr;
    uint32_t                       *acmVsiteStart_ = nullptr;
    AcmVsite                       *acmVsites_     = nullptr;
    uint64_t                        nacmVsite_     = 0;
    //! The recipes of the fixed parts of the charges, see createChargeRecipes
    ChargeRecipe                   *chargeRecipes_ = nullptr;
    uint64_t                        nchargeRecipe_ = 0;
    //! Whether the coordinates of the charge generation are to be constructed, see acmVsitesChanged
    bool                            acmVsitesChanged_ = false;
    //! Work space: the charges of each molecule and the matrices of each fragment
    double                         *acmCharges_   = nullptr;
    char                           *acmDone_      = nullptr;
    double                         *acmWork_      = nullptr;
    uint64_t                       *acmWorkStart_ = nullptr;
    uint32_t                        nacmMolecule_ = 0;
    //! Number of particles of the molecules
    uint64_t                        nacmAtom_ = 0;
    uint64_t                        nacm_ = 0, nacmBond_ = 0, nacmFragment_ = 0;
    bool                            acmGaussian_ = false;
    double                          acmCoulomb_ = 0, acmElectronVolt_ = 1;
    //! The list of each pair and each bonded interaction
    uint32_t                       *pairList_ = nullptr;
    uint32_t                       *bondList_ = nullptr;
    //! The first list of each block, and the number of lists at the end
    uint32_t                       *firstListDev_ = nullptr;
    //! Start of the energy arrays times terms times slots of each block
    uint64_t                       *energyWork_ = nullptr;
    uint64_t                        nenergyWork_ = 0;
    uint64_t                       *atomStart_ = nullptr;
    uint64_t                       *atomEntries_ = nullptr;
    //! The energy sums of each chunk of energyChunk interactions in each slot
    double                         *eChunk_   = nullptr;
    double                         *eList_    = nullptr;
    //! The energies of the blocks computed last, one after the other, see computeBlocks
    double                         *energies_ = nullptr;
    HostArray<double>               energiesHost_;
    //! Where the energies of each block computed start in energies_, see computeBlocks
    uint64_t                       *energyStart_ = nullptr;
    std::vector<uint64_t>           energyStartHost_;
    //! The shells, see SystemBlock
    ForceShell                     *forceShells_  = nullptr;
    MovingShell                    *movingShells_ = nullptr;
    ShellFragment                  *fragments_    = nullptr;
    //! The fragment of each moving shell in its block
    uint32_t                       *movingFrag_   = nullptr;
    //! The shells of each block, of which the forces are computed when relaxing
    uint32_t                       *shellRows_    = nullptr;
    //! The block or list of the first work item of each wavefront of the kernels, see findRangeFrom
    uint32_t                       *hintListWork_ = nullptr, *hintListSlot_ = nullptr, *hintAtomWork_ = nullptr;
    uint32_t                       *hintSlot_ = nullptr, *hintShellWork_ = nullptr, *hintMovingPart_ = nullptr;
    uint32_t                       *hintShellRow_ = nullptr;
    //! The shell-parent pairs
    uint32_t                       *resets_       = nullptr;
    //! The interactions with shells of each particle, as atomStart_ and atomEntries_
    uint64_t                       *shellStart_   = nullptr;
    uint64_t                       *shellEntries_ = nullptr;
    uint64_t                        nslotAll_ = 0, nshellRowWork_ = 0, nshellWork_ = 0;
    uint64_t                        nforcePart_ = 0, nmovingPart_ = 0, nresetWork_ = 0;
    //! The state of the relaxation of each fragment in each slot, two per slot
    void                           *state_ = nullptr;
    char                           *done_  = nullptr;
    char                           *doneHistory_ = nullptr;
    HostArray<char>                 doneHost_;
    float                          *dx_    = nullptr;
    double                         *partF2_ = nullptr, *partAfter_ = nullptr, *partDxF_ = nullptr;
    int                            *valid_ = nullptr;
    //! The number of slots whose shells did not converge, see relaxBlocks
    uint64_t                       *nfailed_ = nullptr;
    //! The statistics of the relaxations without polling, see relaxStatistics: computed, failed, not finite
    uint64_t                       *relaxStats_ = nullptr;
    //! Whether each particle of all blocks is a shell
    char                           *isShell_ = nullptr;
    //! The kept coordinates of each block, see saveCoordinates
    std::vector<int64_t *>          saved_;
    //! The kept shell coordinates of all blocks, in the layout of x_, see keepShells
    int64_t                        *keptShells_ = nullptr;
    //! Whether the shells of each slot converged in the last relaxation, see keepShells
    char                           *converged_ = nullptr;
    //! The energies of all blocks after and before relaxing, at fixedStart_, see keepEnergies
    double                         *energiesAfter_  = nullptr;
    double                         *energiesBefore_ = nullptr;
    uint64_t                       *fixedStart_     = nullptr;
    //! The deviations, see createDeviations
    DeviationSystem                 dev_;
    DeviationType                  *devTypes_          = nullptr;
    uint32_t                       *devBlockFirstType_ = nullptr;
    uint32_t                       *devSlotOfConf_     = nullptr;
    uint64_t                       *devBlockFirstConf_ = nullptr;
    DeviationDimer                 *devDimers_         = nullptr;
    int32_t                        *devMatches_        = nullptr;
    int32_t                        *devTargets_        = nullptr;
    double                         *devWeights_        = nullptr;
    double                         *devEqm_            = nullptr;
    char                           *devHave_           = nullptr;
    DeviationMonomer               *devMonomers_       = nullptr;
    double                         *devEpotWeights_    = nullptr;
    double                         *devEpotEqm_        = nullptr;
    char                           *devEpotHave_       = nullptr;
    //! The weight and loss of each target of each conformation of the dimers, and of the topologies
    double                         *devTerms_          = nullptr;
    double                         *devEpotTerms_      = nullptr;
    uint64_t                        devNDimerConf_     = 0, devNMonoConf_ = 0;
    //! The sums of each dimer, and the results of computeDeviations
    double                         *devDimerSums_      = nullptr;
    HostArray<double>               devSumsHost_;
    int                             maxIterAllocated_ = 0;
    //! The Monte Carlo chain, see createChain
    ChainMolecule                  *chainMols_      = nullptr;
    size_t                          nchainMol_      = 0;
    uint32_t                       *chainDeps_      = nullptr;
    ChainBound                     *chainBounds_    = nullptr;
    uint64_t                       *chainExtras_    = nullptr;
    ChainTerm                      *chainTerms_     = nullptr;
    uint32_t                        nchainTerm_     = 0;
    int                             chainLoss_      = 0;
    //! The first value of the first set kept of each molecule, and the size of a set
    uint64_t                       *chainEntry_     = nullptr;
    std::vector<uint64_t>           chainEntryHost_;
    std::vector<uint64_t>           chainSizeHost_;
    uint64_t                        chainCacheSize_ = 0;
    double                         *chainCache_     = nullptr;
    char                           *chainValid_     = nullptr;
    double                         *chainParams_    = nullptr;
    size_t                          nchainParam_    = 0;
    int32_t                        *chainHit_       = nullptr;
    char                           *chainActive_    = nullptr;
    double                         *chainBoundNow_  = nullptr;
    //! The masks of the blocks, for all, the monomers and the dimers, see chainMask
    char                           *chainMasks_     = nullptr;
    int32_t                        *chainBlockMol_  = nullptr;
    int32_t                        *chainBlockWhat_ = nullptr;
    double                         *chainChi2_      = nullptr;
};

/*! \brief The range of the first work item of each wavefront, for findRangeFrom
 * \param[in] sys   The system, which owns the device memory
 * \param[in] a     The ranges on the host
 * \param[in] total The number of work items
 * \param[in] start The start of each range
 * \return The device array
 */
template <typename T>
static uint32_t *makeRangeHints(System *sys, const std::vector<T> &a, uint64_t total, uint64_t T::*start)
{
    std::vector<uint32_t> hint(total/wavefrontSize + 1, 0);
    for (size_t g = 0; g < hint.size() && !a.empty(); g++)
    {
        hint[g] = findRange(a.data(), a.size(), g*wavefrontSize, start);
    }
    auto *d = sys->allocate<uint32_t>(hint.size());
    sys->upload(d, hint.data(), hint.size());
    sys->queue_->wait_and_throw();
    return d;
}

System *createSystem(const SystemBlock *blocks, size_t nblock)
{
    sycl::queue *queue = gpuQueue();
    if (!queue)
    {
        return nullptr;
    }
    // The system is freed if the GPU has no memory left for it
    std::unique_ptr<System> holder(new System(queue));
    auto *sys = holder.get();
    std::vector<int64_t>     x;
    std::vector<DevicePair>  pairs;
    std::vector<DeviceBond>  bonds;
    // The parameter tables and the parameters of lists without table, from
    // the interactions, which updateParameters changes
    std::vector<float>       tables, params;
    uint64_t                 natomAll = 0;
    std::vector<uint32_t>    pairList, bondList, firstList;
    std::vector<uint32_t>    pairIndex, bondIndex;
    std::vector<uint64_t>    energyWork;
    std::vector<uint64_t>    atomStart(1, 0);
    std::vector<uint64_t>    atomEntries;
    std::vector<DeviceVsite> vsites;
    std::vector<ForceShell>    forceShells;
    std::vector<MovingShell>   movingShells;
    std::vector<ShellFragment> fragments;
    std::vector<uint32_t>      movingFrag, shellRows, resets;
    std::vector<uint64_t>      shellStart(1, 0), shellEntries;
    std::vector<char>          isShell;
    for (size_t b = 0; b < nblock; b++)
    {
        const auto        &sb   = blocks[b];
        // The block whose interactions this one shares, itself if none
        size_t             base = sb.topologyOf >= 0 ? static_cast<size_t>(sb.topologyOf) : b;
        if (base != b)
        {
            const auto &bb = blocks[base];
            bool same = base < b && bb.topologyOf < 0 && bb.natom == sb.natom && bb.nlist == sb.nlist;
            for (size_t l = 0; same && l < sb.nlist; l++)
            {
                same = (bb.lists[l].n == sb.lists[l].n && bb.lists[l].bonded == sb.lists[l].bonded &&
                        bb.lists[l].kernel == sb.lists[l].kernel && bb.lists[l].shells == sb.lists[l].shells &&
                        bb.lists[l].ntable == sb.lists[l].ntable &&
                        (bb.lists[l].index == nullptr) == (sb.lists[l].index == nullptr));
            }
            if (!same)
            {
                throw std::runtime_error("createSystem: a block does not have the interactions of the block it shares them with");
            }
        }
        DeviceSystemBlock  db;
        db.xOffset         = sys->nx_;
        db.atomOffset      = natomAll;
        db.entryAtomOffset = base == b ? atomStart.size() - 1 : sys->blocksHost_[base].entryAtomOffset;
        db.listShift       = base == b ? 0 : sys->listsHost_.size() - sys->firstList_[base];
        natomAll          += sb.natom;
        db.atomWorkOffset = sys->natomWork_;
        db.natom          = sb.natom;
        db.stride         = sb.stride;
        db.nslot          = sb.nslot;
        // Such that the work items of a wavefront have the same shell
        db.nslotPad       = (sb.nslot + wavefrontSize - 1)/wavefrontSize*wavefrontSize;
        db.nenergy        = sb.nenergy;
        sys->blocksHost_.push_back(db);
        sys->firstList_.push_back(sys->listsHost_.size());
        firstList.push_back(sys->listsHost_.size());
        energyWork.push_back(sys->nenergyWork_);
        sys->nenergyWork_ += sb.nenergy*pairEnergyTerms*sb.nslot;
        x.insert(x.end(), sb.x, sb.x + 3*sb.natom*sb.stride);
        // The shells and fragments, with their indices in the arrays of all blocks
        auto &dbr = sys->blocksHost_.back();
        dbr.slotOffset            = sys->nslotAll_;
        dbr.firstForceShell       = forceShells.size();
        dbr.firstMovingShell      = movingShells.size();
        dbr.firstFragment         = fragments.size();
        dbr.nfragment             = sb.nfragment;
        dbr.firstReset            = resets.size()/2;
        dbr.nreset                = sb.nshellParent;
        dbr.resetWorkOffset       = sys->nresetWork_;
        std::vector<uint32_t> rows;
        for (size_t fi = 0; fi < sb.nfragment; fi++)
        {
            auto fr = sb.fragments[fi];
            for (uint32_t k = 0; k < fr.nforceShell; k++)
            {
                forceShells.push_back(sb.forceShells[fr.firstForceShell + k]);
                rows.push_back(sb.forceShells[fr.firstForceShell + k].shell);
            }
            for (uint32_t k = 0; k < fr.nmovingShell; k++)
            {
                movingShells.push_back(sb.movingShells[fr.firstMovingShell + k]);
                movingFrag.push_back(fi);
                rows.push_back(sb.movingShells[fr.firstMovingShell + k].shell);
            }
            fr.firstForceShell  = forceShells.size() - fr.nforceShell;
            fr.firstMovingShell = movingShells.size() - fr.nmovingShell;
            fragments.push_back(fr);
        }
        dbr.nforceShell  = forceShells.size() - dbr.firstForceShell;
        dbr.nmovingShell = movingShells.size() - dbr.firstMovingShell;
        std::sort(rows.begin(), rows.end());
        rows.erase(std::unique(rows.begin(), rows.end()), rows.end());
        dbr.firstShellRow         = shellRows.size();
        dbr.nshellRow             = rows.size();
        shellRows.insert(shellRows.end(), rows.begin(), rows.end());
        dbr.shellRowWorkOffset    = sys->nshellRowWork_;
        dbr.shellWorkOffset       = sys->nshellWork_;
        dbr.forceShellPartOffset  = sys->nforcePart_;
        dbr.movingShellPartOffset = sys->nmovingPart_;
        resets.insert(resets.end(), sb.shellParents, sb.shellParents + 2*sb.nshellParent);
        {
            std::vector<char> shellFlag(sb.natom, 0);
            for (auto r : rows)
            {
                shellFlag[r] = 1;
            }
            isShell.insert(isShell.end(), shellFlag.begin(), shellFlag.end());
        }
        sys->nslotAll_      += sb.nslot;
        sys->nshellRowWork_ += rows.size()*dbr.nslotPad;
        sys->nshellWork_    += (dbr.nforceShell + dbr.nmovingShell)*sb.nslot;
        sys->nforcePart_    += dbr.nforceShell*sb.nslot;
        sys->nmovingPart_   += dbr.nmovingShell*sb.nslot;
        sys->nresetWork_    += sb.nshellParent*sb.nslot;
        sys->nx_        += 3*sb.natom*sb.stride;
        sys->natomWork_ += sb.natom*sb.nslot;
        sys->nenergyAll_ += sb.nenergy*pairEnergyTerms*sb.stride;
        // The interactions of each particle, of all lists in their order,
        // and those of the lists with shells only
        std::vector<std::vector<uint64_t>> perAtom(base == b ? sb.natom : 0), perAtomShells(base == b ? sb.natom : 0);
        for (size_t l = 0; l < sb.nlist; l++)
        {
            const auto &il = sb.lists[l];
            DeviceList  dl;
            dl.workOffset     = sys->nwork_;
            dl.listSlotOffset = sys->nlistSlot_;
            dl.entryOffset    = base == b ? (il.bonded ? bonds.size() : pairs.size())
                                          : sys->listsHost_[sys->firstList_[base] + l].entryOffset;
            dl.tableOffset    = sys->ntable_;
            sys->ntable_     += il.ntable;
            dl.paramOffset    = sys->nparam_;
            dl.ownParams      = il.index ? 0 : 1;
            if (!il.index)
            {
                sys->nparam_ += il.n;
            }
            dl.n              = il.n;
            dl.block          = b;
            dl.kernel         = il.kernel;
            dl.bonded         = il.bonded ? 1 : 0;
            dl.shells         = il.shells ? 1 : 0;
            dl.energyIndex    = il.energyIndex;
            dl.nterm          = il.nterm;
            uint32_t listIndex = sys->listsHost_.size();
            sys->listsHost_.push_back(dl);
            sys->nwork_     += (il.n + energyChunk - 1)/energyChunk*sb.nslot;
            sys->nlistSlot_ += sb.nslot;
            tables.resize(8*sys->ntable_, 0.0f);
            for (size_t e = 0; e < il.n; e++)
            {
                const float *p = il.bonded ? il.bonds[e].p : il.pairs[e].p;
                if (il.index)
                {
                    std::copy(p, p + 8, tables.begin() + 8*(dl.tableOffset + il.index[e]));
                }
                else
                {
                    params.insert(params.end(), p, p + 8);
                }
                if (base != b)
                {
                    continue;
                }
                if (il.bonded)
                {
                    const auto &be  = il.bonds[e];
                    uint64_t    idx = bonds.size();
                    int         na  = (il.kernel == static_cast<uint32_t>(BondedKernel::FourierDihedral) ||
                                       il.kernel == static_cast<uint32_t>(BondedKernel::HarmonicDihedral)) ? 4 : 3;
                    for (int k = 0; k < na; k++)
                    {
                        perAtom[be.a[k]].push_back((idx << indexShift) | bondedFlag | k);
                        if (il.shells)
                        {
                            perAtomShells[be.a[k]].push_back((idx << indexShift) | bondedFlag | k);
                        }
                    }
                    bonds.push_back({ { be.a[0], be.a[1], be.a[2], be.a[3] } });
                    bondList.push_back(listIndex);
                    bondIndex.push_back(il.index ? il.index[e] : 0);
                }
                else
                {
                    const auto &pe  = il.pairs[e];
                    uint64_t    idx = pairs.size();
                    perAtom[pe.i].push_back((idx << indexShift) | 0);
                    perAtom[pe.j].push_back((idx << indexShift) | 1);
                    if (il.shells)
                    {
                        perAtomShells[pe.i].push_back((idx << indexShift) | 0);
                        perAtomShells[pe.j].push_back((idx << indexShift) | 1);
                    }
                    pairs.push_back({ pe.i, pe.j });
                    pairList.push_back(listIndex);
                    pairIndex.push_back(il.index ? il.index[e] : 0);
                }
            }
        }
        for (size_t a = 0; a < perAtom.size(); a++)
        {
            atomEntries.insert(atomEntries.end(), perAtom[a].begin(), perAtom[a].end());
            atomStart.push_back(atomEntries.size());
            shellEntries.insert(shellEntries.end(), perAtomShells[a].begin(), perAtomShells[a].end());
            shellStart.push_back(shellEntries.size());
        }
        for (size_t v = 0; v < sb.nvsite; v++)
        {
            DeviceVsite dv;
            dv.workOffset = sys->nvsiteWork_;
            dv.block      = b;
            dv.entry      = sb.vsites[v];
            vsites.push_back(dv);
            sys->nvsiteWork_ += sb.nslot;
        }
    }
    sys->nblock_ = nblock;
    sys->nlist_  = sys->listsHost_.size();
    sys->nvsite_ = vsites.size();
    sys->npair_  = pairs.size();
    sys->nbond_  = bonds.size();
    sys->blocks_      = sys->allocate<DeviceSystemBlock>(nblock);
    sys->lists_       = sys->allocate<DeviceList>(sys->nlist_);
    sys->vsites_      = sys->allocate<DeviceVsite>(sys->nvsite_);
    sys->x_           = sys->allocate<int64_t>(sys->nx_);
    sys->f_           = sys->allocate<float>(sys->nx_);
    sys->fFixed_      = sys->allocate<float>(sys->nx_);
    sys->pairs_       = sys->allocate<DevicePair>(sys->npair_);
    sys->bonds_       = sys->allocate<DeviceBond>(sys->nbond_);
    sys->params_      = sys->allocate<float>(8*sys->nparam_);
    firstList.push_back(sys->listsHost_.size());
    sys->pairList_     = sys->allocate<uint32_t>(pairList.size());
    sys->pairIndex_    = sys->allocate<uint32_t>(pairIndex.size());
    sys->bondIndex_    = sys->allocate<uint32_t>(bondIndex.size());
    sys->tables_       = sys->allocate<float>(8*sys->ntable_);
    sys->ncharge_      = natomAll;
    sys->charges_      = sys->allocate<double>(4*sys->ncharge_);
    sys->vsitesHost_   = vsites;
    sys->bondList_     = sys->allocate<uint32_t>(bondList.size());
    sys->firstListDev_ = sys->allocate<uint32_t>(firstList.size());
    sys->energyWork_   = sys->allocate<uint64_t>(energyWork.size());
    sys->atomStart_   = sys->allocate<uint64_t>(atomStart.size());
    sys->atomEntries_ = sys->allocate<uint64_t>(atomEntries.size());
    sys->eChunk_      = sys->allocate<double>(pairEnergyTerms*sys->nwork_);
    sys->eList_       = sys->allocate<double>(pairEnergyTerms*sys->nlistSlot_);
    sys->energies_    = sys->allocate<double>(sys->nenergyAll_);
    sys->energyStart_ = sys->allocate<uint64_t>(nblock);
    sys->forceShells_  = sys->allocate<ForceShell>(forceShells.size());
    sys->movingShells_ = sys->allocate<MovingShell>(movingShells.size());
    sys->fragments_    = sys->allocate<ShellFragment>(fragments.size());
    sys->movingFrag_   = sys->allocate<uint32_t>(movingFrag.size());
    sys->shellRows_    = sys->allocate<uint32_t>(shellRows.size());
    sys->resets_       = sys->allocate<uint32_t>(resets.size());
    sys->shellStart_   = sys->allocate<uint64_t>(shellStart.size());
    sys->shellEntries_ = sys->allocate<uint64_t>(shellEntries.size());
    sys->isShell_      = sys->allocate<char>(isShell.size());
    sys->upload(sys->isShell_, isShell.data(), isShell.size());
    sys->done_         = sys->allocate<char>(sys->nslotAll_);
    sys->dx_           = sys->allocate<float>(sys->nx_);
    sys->partF2_       = sys->allocate<double>(sys->nforcePart_);
    sys->partAfter_    = sys->allocate<double>(sys->nmovingPart_);
    sys->partDxF_      = sys->allocate<double>(sys->nmovingPart_);
    sys->valid_        = sys->allocate<int>(1);
    sys->nfailed_      = sys->allocate<uint64_t>(1);
    sys->relaxStats_   = sys->allocate<uint64_t>(3);
    sys->queue_->memset(sys->relaxStats_, 0, 3*sizeof(uint64_t));
    sys->converged_    = sys->allocate<char>(std::max<uint64_t>(sys->nslotAll_, 1));
    sys->queue_->memset(sys->converged_, 0, std::max<uint64_t>(sys->nslotAll_, 1));
    sys->saved_.assign(nblock, nullptr);
    sys->upload(sys->forceShells_, forceShells.data(), forceShells.size());
    sys->upload(sys->movingShells_, movingShells.data(), movingShells.size());
    sys->upload(sys->fragments_, fragments.data(), fragments.size());
    sys->upload(sys->movingFrag_, movingFrag.data(), movingFrag.size());
    sys->upload(sys->shellRows_, shellRows.data(), shellRows.size());
    sys->upload(sys->resets_, resets.data(), resets.size());
    sys->upload(sys->shellStart_, shellStart.data(), shellStart.size());
    sys->upload(sys->shellEntries_, shellEntries.data(), shellEntries.size());
    sys->upload(sys->blocks_, sys->blocksHost_.data(), nblock);
    sys->upload(sys->lists_, sys->listsHost_.data(), sys->nlist_);
    sys->upload(sys->vsites_, vsites.data(), vsites.size());
    sys->upload(sys->x_, x.data(), x.size());
    sys->upload(sys->pairs_, pairs.data(), pairs.size());
    sys->upload(sys->bonds_, bonds.data(), bonds.size());
    sys->upload(sys->params_, params.data(), params.size());
    sys->upload(sys->tables_, tables.data(), tables.size());
    sys->upload(sys->pairList_, pairList.data(), pairList.size());
    sys->upload(sys->pairIndex_, pairIndex.data(), pairIndex.size());
    sys->upload(sys->bondIndex_, bondIndex.data(), bondIndex.size());
    sys->upload(sys->bondList_, bondList.data(), bondList.size());
    sys->upload(sys->firstListDev_, firstList.data(), firstList.size());
    sys->upload(sys->energyWork_, energyWork.data(), energyWork.size());
    sys->upload(sys->atomStart_, atomStart.data(), atomStart.size());
    sys->upload(sys->atomEntries_, atomEntries.data(), atomEntries.size());
    queue->wait_and_throw();
    sys->hintListWork_   = makeRangeHints(sys, sys->listsHost_, sys->nwork_, &DeviceList::workOffset);
    sys->hintListSlot_   = makeRangeHints(sys, sys->listsHost_, sys->nlistSlot_, &DeviceList::listSlotOffset);
    sys->hintAtomWork_   = makeRangeHints(sys, sys->blocksHost_, sys->natomWork_, &DeviceSystemBlock::atomWorkOffset);
    sys->hintSlot_       = makeRangeHints(sys, sys->blocksHost_, sys->nslotAll_, &DeviceSystemBlock::slotOffset);
    sys->hintShellWork_  = makeRangeHints(sys, sys->blocksHost_, sys->nshellWork_, &DeviceSystemBlock::shellWorkOffset);
    sys->hintMovingPart_ = makeRangeHints(sys, sys->blocksHost_, sys->nmovingPart_, &DeviceSystemBlock::movingShellPartOffset);
    sys->hintShellRow_   = makeRangeHints(sys, sys->blocksHost_, sys->nshellRowWork_, &DeviceSystemBlock::shellRowWorkOffset);
    if (std::getenv("ACT_GPU_SYSTEM_INFO"))
    {
        // For developers, the amount of work of the kernels
        uint64_t shellEntryWork = 0, shellEntryPadded = 0, rowSlots = 0, padded = 0, maxEntries = 0;
        std::vector<uint64_t> rowEntries;
        for (const auto &b : sys->blocksHost_)
        {
            for (uint32_t r = 0; r < b.nshellRow; r++)
            {
                uint32_t a  = shellRows[b.firstShellRow + r];
                uint64_t ne = shellStart[b.entryAtomOffset + a + 1] - shellStart[b.entryAtomOffset + a];
                shellEntryWork   += ne*b.nslot;
                shellEntryPadded += ne*b.nslotPad;
                maxEntries        = std::max(maxEntries, ne);
                rowEntries.push_back(ne);
            }
            rowSlots += static_cast<uint64_t>(b.nshellRow)*b.nslot;
            padded   += static_cast<uint64_t>(b.nshellRow)*b.nslotPad;
        }
        fprintf(stderr, "GPU system: %.1f MB, %lu own parameters, %lu table rows\n",
                sys->allocatedBytes_/1e6, sys->nparam_, sys->ntable_);
        fprintf(stderr, "GPU system: %u blocks, %u lists, %lu slots, %lu pairs, %lu bonded, "
                "%lu energy chunk-slots, %lu shell-row-slots (%lu padded), %lu shell interaction-slots (%lu padded)\n",
                sys->nblock_, sys->nlist_, sys->nslotAll_, sys->npair_, sys->nbond_, sys->nwork_,
                rowSlots, padded, shellEntryWork, shellEntryPadded);
        if (!rowEntries.empty())
        {
            std::sort(rowEntries.begin(), rowEntries.end());
            auto q = [&rowEntries](double f) { return rowEntries[static_cast<size_t>(f*(rowEntries.size() - 1))]; };
            fprintf(stderr, "GPU system: %zu shell rows, interactions per row median %lu, 90%% %lu, 99%% %lu, max %lu\n",
                    rowEntries.size(), q(0.5), q(0.9), q(0.99), maxEntries);
        }
    }
    return holder.release();
}

void destroySystem(System *system)
{
    delete system;
}

void updateShells(System *sys, size_t block, const MovingShell *shells, size_t n)
{
    const auto &b = sys->blocksHost_[block];
    if (n != b.nmovingShell)
    {
        throw std::runtime_error("updateShells: the number of moving shells changed");
    }
    sys->upload(sys->movingShells_ + b.firstMovingShell, shells, n);
}

//! The combined width of two charge distributions, as in forces_gpu.cpp
inline double combinedWidth(double zi, double zj)
{
    if (zi == 0 || zj == 0)
    {
        return zi + zj;
    }
    return zi*zj/sycl::sqrt(zi*zi + zj*zj);
}

/*! \brief Solve a x = b with Gaussian elimination and complete pivoting
 *
 * As the column pivoted QR of MatrixWrapper::solve, the unknowns of
 * negligible pivots are zero if the matrix is singular, as for the bonds
 * of rings without bond hardness. The matrix and b are overwritten.
 */
inline void solveLinear(double *a, double *b, double *x, double *perm, uint32_t n)
{
    double amax = 0;
    for (uint32_t k = 0; k < n*n; k++)
    {
        amax = sycl::fmax(amax, sycl::fabs(a[k]));
    }
    for (uint32_t k = 0; k < n; k++)
    {
        perm[k] = k;
    }
    uint32_t rank = n;
    for (uint32_t k = 0; k < n; k++)
    {
        uint32_t p = k, q = k;
        double   best = 0;
        for (uint32_t i = k; i < n; i++)
        {
            for (uint32_t j = k; j < n; j++)
            {
                double v = sycl::fabs(a[i*n+j]);
                if (v > best)
                {
                    best = v;
                    p    = i;
                    q    = j;
                }
            }
        }
        if (best <= 1e-13*amax || best == 0)
        {
            rank = k;
            break;
        }
        if (p != k)
        {
            for (uint32_t j = 0; j < n; j++)
            {
                double t = a[k*n+j]; a[k*n+j] = a[p*n+j]; a[p*n+j] = t;
            }
            double t = b[k]; b[k] = b[p]; b[p] = t;
        }
        if (q != k)
        {
            for (uint32_t i = 0; i < n; i++)
            {
                double t = a[i*n+k]; a[i*n+k] = a[i*n+q]; a[i*n+q] = t;
            }
            double t = perm[k]; perm[k] = perm[q]; perm[q] = t;
        }
        for (uint32_t i = k+1; i < n; i++)
        {
            double f = a[i*n+k]/a[k*n+k];
            for (uint32_t j = k; j < n; j++)
            {
                a[i*n+j] -= f*a[k*n+j];
            }
            b[i] -= f*b[k];
        }
    }
    for (uint32_t k = rank; k < n; k++)
    {
        b[k] = 0;
    }
    for (uint32_t kk = rank; kk > 0; kk--)
    {
        uint32_t k = kk - 1;
        double   v = b[k];
        for (uint32_t j = k+1; j < rank; j++)
        {
            v -= a[k*n+j]*b[j];
        }
        b[k] = v/a[k*n+k];
    }
    for (uint32_t k = 0; k < n; k++)
    {
        x[static_cast<uint32_t>(perm[k])] = b[k];
    }
}

//! The Coulomb integral of two Gaussian charges without units, as Coulomb_GG
inline double coulombGG(double r, double zi, double zj)
{
    double zeff = 0;
    if (zi != 0 || zj != 0)
    {
        zeff = (zi == 0) ? zj : ((zj == 0) ? zi : zi*zj/sycl::sqrt(zi*zi + zj*zj));
    }
    if (zeff == 0)
    {
        return r == 0 ? 0.0 : 1.0/r;
    }
    if (r == 0)
    {
        // 1/sqrt(pi), as invsqrt_pi in Nuclear_GG
        return 2.0*zeff*0.56418958354775628695;
    }
    return sycl::erf(zeff*r)/r;
}

//! The work-group size of the charge generation of a fragment
constexpr uint32_t acmGroupSize = 64;

/*! \brief Generate the charges of the ACM particles with SQE, see QgenAcm::solveSQE
 *
 * First the charges of each molecule are those uploaded, then a work-group
 * per fragment generates those of its ACM particles, with the elimination
 * of solveLinear in parallel, and finally a work item per molecule makes
 * them symmetric as FragmentHandler::generateCharges does and stores them
 * in the blocks of the fragments and the molecule. Only if the values
 * changed, see updateAcm, otherwise the uploaded charges are the generated.
 */
void generateAcmCharges(System *sys)
{
    if (sys->nacmMolecule_ == 0 || !sys->acmChanged_)
    {
        return;
    }
    sys->acmChanged_ = false;
    auto     dM   = sys->acmMolecules_;
    auto     dF   = sys->acmFragments_;
    auto     dA   = sys->acmParticles_;
    auto     dBd  = sys->acmBonds_;
    auto     dS   = sys->acmSymmetry_;
    auto     dV   = sys->acmValues_;
    auto     dX   = sys->acmX_;
    auto     dPV  = sys->acmParticleValues_;
    auto     dBV  = sys->acmBondValues_;
    auto     dBF  = sys->acmBondFactors_;
    auto     dFV  = sys->acmFixedValues_;
    auto     dQm  = sys->acmCharges_;
    auto     dDn  = sys->acmDone_;
    auto     dW   = sys->acmWork_;
    auto     dWS  = sys->acmWorkStart_;
    auto     dB   = sys->blocks_;
    auto     dQ   = sys->charges_;
    bool     gauss = sys->acmGaussian_;
    double   coul  = sys->acmCoulomb_;
    double   ev    = sys->acmElectronVolt_;
    auto    *queue = sys->queue_;
    // The charges as uploaded
    queue->parallel_for(sycl::range<1>(sys->nacmMolecule_), [=](sycl::id<1> id)
    {
        const AcmMolecule &m  = dM[id[0]];
        for (uint32_t fi = 0; fi < m.nfragment; fi++)
        {
            const AcmFragment &f   = dF[m.firstFragment + fi];
            uint64_t           off = dB[f.block].atomOffset;
            for (uint32_t a = 0; a < f.natom; a++)
            {
                dQm[m.firstAtom + f.atomStart + a] = dQ[4*(off + a)];
            }
        }
    });
    queue->parallel_for(sycl::nd_range<1>(sys->nacmFragment_*acmGroupSize, acmGroupSize),
                        [=](sycl::nd_item<1> it)
    {
        auto               g   = it.get_group();
        uint32_t           fg  = it.get_group(0);
        uint32_t           lid = it.get_local_id(0);
        const uint32_t     G   = acmGroupSize;
        const AcmFragment &f   = dF[fg];
        uint32_t           n   = f.nacm, nb = f.nbond;
        if (n == 0)
        {
            return;
        }
        // The values of the parameters, zero if there are none
        auto value = [dV](uint32_t v) { return v == noValue ? 0.0 : dV[v]; };
        const double   *xa  = dX + 3*f.firstAcm;
        const uint32_t *pv  = dPV + 3*f.firstAcm;
        const uint32_t *bv  = dBV + 2*f.firstBond;
        const double   *bf  = dBF + f.firstBond;
        const uint32_t *bd  = dBd + 2*f.firstBond;
        double         *J   = dW + dWS[fg];
        double         *L   = J + n*n;
        double         *rhs = L + nb*nb;
        double         *sol = rhs + nb;
        double         *prm = sol + nb;
        double         *fac = prm + nb;
        double         *cc  = fac + nb;
        // For the reductions over the work-group, without local memory
        double         *red = cc + n;
        auto groupMax = [&](double v)
        {
            red[lid] = v;
            sycl::group_barrier(g);
            double r = red[0];
            for (uint32_t t = 1; t < G; t++)
            {
                r = sycl::fmax(r, red[t]);
            }
            sycl::group_barrier(g);
            return r;
        };
        auto groupMinIndex = [&](uint32_t v)
        {
            red[lid] = v;
            sycl::group_barrier(g);
            uint32_t r = static_cast<uint32_t>(red[0]);
            for (uint32_t t = 1; t < G; t++)
            {
                r = sycl::min(r, static_cast<uint32_t>(red[t]));
            }
            sycl::group_barrier(g);
            return r;
        };
        // As solveSQE
        double qfixed = f.qtotal;
        for (uint32_t k = 0; k < f.nfixed; k++)
        {
            qfixed -= value(dFV[f.firstFixed + k]);
        }
        double qshift = qfixed/n;
        for (uint32_t idx = lid; idx < n*n; idx += G)
        {
            uint32_t i = idx/n, j = idx % n;
            if (i == j)
            {
                double j0 = value(pv[3*i+1]);
                J[idx]    = j0 > 0 ? j0 : 0;
                continue;
            }
            // The distance as the processor computes it, in single precision
            float dx[3], r2 = 0;
            for (int d = 0; d < 3; d++)
            {
                dx[d] = static_cast<float>(xa[3*i+d]) - static_cast<float>(xa[3*j+d]);
                r2   += dx[d]*dx[d];
            }
            double r    = sycl::sqrt(r2);
            double eTot = gauss ? coulombGG(r, value(pv[3*i+2]), value(pv[3*j+2])) : 1.0/r;
            J[idx]      = 0.5*((coul*eTot)/ev);
        }
        sycl::group_barrier(g);
        for (uint32_t idx = lid; idx < nb*nb; idx += G)
        {
            uint32_t b  = idx/nb, c = idx % nb;
            uint32_t ai = bd[2*b], aj = bd[2*b+1], ak = bd[2*c], al = bd[2*c+1];
            double   v  = J[ai*n+ak] - J[ai*n+al] - J[aj*n+ak] + J[aj*n+al];
            if (b == c)
            {
                v += value(bv[2*b+1]);
            }
            L[idx] = v;
        }
        double qcorr = (1.0*f.qtotal)/n;
        for (uint32_t i = lid; i < n; i += G)
        {
            double c = value(pv[3*i]);
            for (uint32_t b = 0; b < nb; b++)
            {
                double dchi = value(bv[2*b])*bf[b];
                if (i == bd[2*b])
                {
                    c += dchi;
                }
                else if (i == bd[2*b+1])
                {
                    c -= dchi;
                }
            }
            for (uint32_t l = 0; l < n; l++)
            {
                c += J[i*n+l]*qcorr;
            }
            cc[i] = c;
        }
        sycl::group_barrier(g);
        for (uint32_t b = lid; b < nb; b += G)
        {
            rhs[b] = cc[bd[2*b+1]] - cc[bd[2*b]];
            prm[b] = b;
        }
        // Gaussian elimination with complete pivoting, as solveLinear
        double lmax = 0;
        for (uint32_t idx = lid; idx < nb*nb; idx += G)
        {
            lmax = sycl::fmax(lmax, sycl::fabs(L[idx]));
        }
        double   amax = groupMax(lmax);
        uint32_t rank = nb;
        for (uint32_t k = 0; k < nb; k++)
        {
            // The largest element, the first in the order of the rows if equal
            uint32_t m    = nb - k;
            double   best = 0;
            uint32_t bi   = 0xffffffff;
            for (uint32_t idx = lid; idx < m*m; idx += G)
            {
                uint32_t i = k + idx/m, j = k + idx % m;
                double   v = sycl::fabs(L[i*nb+j]);
                if (v > best)
                {
                    best = v;
                    bi   = i*nb + j;
                }
            }
            double   gbest = groupMax(best);
            uint32_t cand  = (best == gbest && gbest > 0) ? bi : 0xffffffff;
            uint32_t pick  = groupMinIndex(cand);
            if (gbest <= 1e-13*amax || gbest == 0)
            {
                rank = k;
                break;
            }
            uint32_t p = pick/nb, q = pick % nb;
            if (p != k)
            {
                for (uint32_t j = lid; j < nb; j += G)
                {
                    double t = L[k*nb+j]; L[k*nb+j] = L[p*nb+j]; L[p*nb+j] = t;
                }
                if (lid == 0)
                {
                    double t = rhs[k]; rhs[k] = rhs[p]; rhs[p] = t;
                }
            }
            sycl::group_barrier(g);
            if (q != k)
            {
                for (uint32_t i = lid; i < nb; i += G)
                {
                    double t = L[i*nb+k]; L[i*nb+k] = L[i*nb+q]; L[i*nb+q] = t;
                }
                if (lid == 0)
                {
                    double t = prm[k]; prm[k] = prm[q]; prm[q] = t;
                }
            }
            sycl::group_barrier(g);
            for (uint32_t i = k+1+lid; i < nb; i += G)
            {
                fac[i] = L[i*nb+k]/L[k*nb+k];
            }
            sycl::group_barrier(g);
            for (uint32_t idx = lid; idx < (m-1)*m; idx += G)
            {
                uint32_t i = k + 1 + idx/m, j = k + idx % m;
                L[i*nb+j] -= fac[i]*L[k*nb+j];
            }
            for (uint32_t i = k+1+lid; i < nb; i += G)
            {
                rhs[i] -= fac[i]*rhs[k];
            }
            sycl::group_barrier(g);
        }
        if (lid == 0)
        {
            for (uint32_t k = rank; k < nb; k++)
            {
                rhs[k] = 0;
            }
            for (uint32_t kk = rank; kk > 0; kk--)
            {
                uint32_t k = kk - 1;
                double   v = rhs[k];
                for (uint32_t j = k+1; j < rank; j++)
                {
                    v -= L[k*nb+j]*rhs[j];
                }
                rhs[k] = v/L[k*nb+k];
            }
            for (uint32_t k = 0; k < nb; k++)
            {
                sol[static_cast<uint32_t>(prm[k])] = rhs[k];
            }
        }
        sycl::group_barrier(g);
        for (uint32_t i = lid; i < n; i += G)
        {
            double q = 0;
            for (uint32_t b = 0; b < nb; b++)
            {
                if (bd[2*b] == i)
                {
                    q += sol[b];
                }
                else if (bd[2*b+1] == i)
                {
                    q -= sol[b];
                }
            }
            dQm[f.moleculeAtom + f.atomStart + dA[f.firstAcm + i]] = q + qshift;
        }
    });
    queue->parallel_for(sycl::range<1>(sys->nacmMolecule_), [=](sycl::id<1> id)
    {
        const AcmMolecule &m  = dM[id[0]];
        double            *qm = dQm + m.firstAtom;
        // As apply_symmetrized_charges
        const int32_t *sym  = dS + m.firstAtom;
        char          *done = dDn + m.firstAtom;
        for (uint32_t i = 0; i < m.natom; i++)
        {
            done[i] = 0;
        }
        for (uint32_t i = 0; i < m.natom; i++)
        {
            if (done[i])
            {
                continue;
            }
            double qsum = qm[i];
            int    nrq  = 1;
            for (uint32_t j = i+1; j < m.natom; j++)
            {
                if (sym[j] == sym[i])
                {
                    qsum += qm[j];
                    nrq++;
                }
            }
            if (nrq > 1)
            {
                double qaver = qsum/nrq;
                for (uint32_t j = 0; j < m.natom; j++)
                {
                    if (sym[j] == sym[i])
                    {
                        qm[j]   = qaver;
                        done[j] = 1;
                    }
                }
            }
            done[i] = 1;
        }
        for (uint32_t fi = 0; fi < m.nfragment; fi++)
        {
            const AcmFragment &f   = dF[m.firstFragment + fi];
            uint64_t           off = dB[f.block].atomOffset;
            for (uint32_t a = 0; a < f.natom; a++)
            {
                dQ[4*(off + a)] = qm[f.atomStart + a];
                if (m.dimerBlock != noBlock)
                {
                    dQ[4*(dB[m.dimerBlock].atomOffset + f.atomStart + a)] = qm[f.atomStart + a];
                }
            }
        }
    });
}

void createAcm(System *sys, const AcmSystem &acm)
{
    sys->nacmMolecule_   = acm.nmolecule;
    sys->nacm_           = acm.nacm;
    sys->nacmBond_       = acm.nbond;
    sys->nacmFragment_   = acm.nfragment;
    sys->acmGaussian_    = acm.gaussian;
    sys->acmCoulomb_     = acm.coulomb;
    sys->acmElectronVolt_ = acm.electronVolt;
    uint32_t natom = 0;
    for (size_t m = 0; m < acm.nmolecule; m++)
    {
        natom = std::max(natom, acm.molecules[m].firstAtom + acm.molecules[m].natom);
    }
    sys->nacmAtom_ = natom;
    std::vector<uint64_t> workStart;
    uint64_t              nwork = 0;
    for (size_t f = 0; f < acm.nfragment; f++)
    {
        const auto &fr = acm.fragments[f];
        workStart.push_back(nwork);
        nwork += fr.nacm*fr.nacm + fr.nbond*fr.nbond + 4*fr.nbond + fr.nacm + acmGroupSize;
    }
    sys->acmMolecules_ = sys->allocate<AcmMolecule>(acm.nmolecule);
    sys->acmFragments_ = sys->allocate<AcmFragment>(acm.nfragment);
    sys->acmParticles_ = sys->allocate<uint32_t>(acm.nacm);
    sys->acmBonds_     = sys->allocate<uint32_t>(2*acm.nbond);
    sys->acmSymmetry_  = sys->allocate<int32_t>(natom);
    sys->nacmValue_         = acm.nvalue;
    sys->acmValues_         = sys->allocate<double>(acm.nvalue);
    sys->acmX_              = sys->allocate<double>(3*acm.nacm);
    sys->acmParticleValues_ = sys->allocate<uint32_t>(3*acm.nacm);
    sys->acmBondValues_     = sys->allocate<uint32_t>(2*acm.nbond);
    sys->acmBondFactors_    = sys->allocate<double>(acm.nbond);
    sys->acmFixedValues_    = sys->allocate<uint32_t>(acm.nfixed);
    sys->upload(sys->acmX_, acm.x, 3*acm.nacm);
    sys->upload(sys->acmParticleValues_, acm.particleValues, 3*acm.nacm);
    sys->upload(sys->acmBondValues_, acm.bondValues, 2*acm.nbond);
    sys->upload(sys->acmBondFactors_, acm.bondFactors, acm.nbond);
    sys->upload(sys->acmFixedValues_, acm.fixedValues, acm.nfixed);
    sys->acmCharges_   = sys->allocate<double>(natom);
    sys->acmDone_      = sys->allocate<char>(natom);
    sys->acmWork_      = sys->allocate<double>(nwork);
    sys->acmWorkStart_ = sys->allocate<uint64_t>(workStart.size());
    sys->upload(sys->acmMolecules_, acm.molecules, acm.nmolecule);
    sys->upload(sys->acmFragments_, acm.fragments, acm.nfragment);
    sys->upload(sys->acmParticles_, acm.acmParticles, acm.nacm);
    sys->upload(sys->acmBonds_, acm.bonds, 2*acm.nbond);
    sys->upload(sys->acmSymmetry_, acm.symmetry, natom);
    sys->upload(sys->acmWorkStart_, workStart.data(), workStart.size());
    sys->queue_->wait_and_throw();
}

void downloadAcmCharges(System *sys, double *charges)
{
    sys->queue_->memcpy(charges, sys->acmCharges_, sys->nacmAtom_*sizeof(double));
    sys->queue_->wait_and_throw();
}

void updateAcm(System *sys, const double *values)
{
    sys->upload(sys->acmValues_, values, sys->nacmValue_);
    sys->acmChanged_ = true;
}

void createChargeRecipes(System *sys, const ChargeRecipe *charges, size_t n)
{
    sys->release(sys->chargeRecipes_);
    sys->nchargeRecipe_ = n;
    sys->chargeRecipes_ = sys->allocate<ChargeRecipe>(std::max<size_t>(n, 1));
    sys->upload(sys->chargeRecipes_, charges, n);
    sys->queue_->wait_and_throw();
}

/*! \brief Write the fixed parts of the charges from the values, see createChargeRecipes
 * \param[in] sys The system
 */
static void computeChargeRecipes(System *sys)
{
    uint64_t nv  = sys->nrecipeValue_;
    auto     dC  = sys->chargeRecipes_;
    auto     dF  = sys->recipeFactors_;
    auto     dV  = sys->recipeValues_;
    auto     dQ  = sys->charges_;
    sys->queue_->parallel_for(sycl::range<1>(sys->nchargeRecipe_), [=](sycl::id<1> id)
    {
        const ChargeRecipe &cr = dC[id[0]];
        double              v  = dV[cr.replica*nv + cr.slot];
        // As ForceFieldParameter::internalValue, see ParameterRecipes::slot
        dQ[cr.position] = cr.internal ? v*dF[cr.slot] : v;
    });
}

void createAcmVsites(System *sys, const int64_t *x, const uint32_t *vsiteStart,
                     const AcmVsite *vsites, size_t nvsite)
{
    sys->release(sys->acmRefX_);
    sys->release(sys->acmVsiteStart_);
    sys->release(sys->acmVsites_);
    sys->nacmVsite_     = nvsite;
    sys->acmRefX_       = sys->allocate<int64_t>(std::max<size_t>(3*sys->nacmAtom_, 1));
    sys->acmVsiteStart_ = sys->allocate<uint32_t>(sys->nacmMolecule_ + 1);
    sys->acmVsites_     = sys->allocate<AcmVsite>(std::max<size_t>(nvsite, 1));
    sys->upload(sys->acmRefX_, x, 3*sys->nacmAtom_);
    sys->upload(sys->acmVsiteStart_, vsiteStart, sys->nacmMolecule_ + 1);
    sys->upload(sys->acmVsites_, vsites, nvsite);
    sys->queue_->wait_and_throw();
    sys->acmVsitesChanged_ = true;
}

void acmVsitesChanged(System *sys)
{
    sys->acmVsitesChanged_ = true;
}

void acmChargesChanged(System *sys)
{
    sys->acmChanged_ = true;
}

void updateAcmCoordinates(System *sys, const double *x)
{
    sys->upload(sys->acmX_, x, 3*sys->nacm_);
    sys->acmChanged_ = true;
}

void downloadCharges(System *sys, double *charges)
{
    sys->queue_->memcpy(charges, sys->charges_, 4*sys->ncharge_*sizeof(double));
    sys->queue_->wait_and_throw();
}

//! The mathematical functions of SYCL, for the functions of combinationkernels.h
struct SyclMath
{
    //! \return the square root of x
    static double sqrt(double x) { return sycl::sqrt(x); }
    //! \return x to the power y
    static double pow(double x, double y) { return sycl::pow(x, y); }
    //! \return the absolute value of x
    static double abs(double x) { return sycl::fabs(x); }
};

void createRowRecipes(System *sys, const SetRecipe *recipes, size_t nset, const RowRecipe *rows, size_t nrow,
                      const double *factors, size_t nvalue, size_t nreplica)
{
    sys->release(sys->recipes_);
    sys->release(sys->rowRecipes_);
    sys->release(sys->recipeFactors_);
    sys->release(sys->recipeValues_);
    sys->release(sys->setValues_);
    sys->nset_           = nset;
    sys->nrowRecipe_     = nrow;
    sys->nrecipeValue_   = nvalue;
    sys->nrecipeReplica_ = nreplica;
    sys->recipes_        = sys->allocate<SetRecipe>(std::max<size_t>(nset, 1));
    sys->rowRecipes_     = sys->allocate<RowRecipe>(std::max<size_t>(nrow, 1));
    sys->recipeFactors_  = sys->allocate<double>(std::max<size_t>(nvalue, 1));
    sys->recipeValues_   = sys->allocate<double>(std::max<size_t>(nvalue*nreplica, 1));
    sys->setValues_      = sys->allocate<double>(std::max<size_t>(maxSetParameters*nset*nreplica, 1));
    sys->upload(sys->recipes_, recipes, nset);
    sys->upload(sys->rowRecipes_, rows, nrow);
    sys->upload(sys->recipeFactors_, factors, nvalue);
    sys->queue_->wait_and_throw();
    sys->recipeValuesChanged_ = false;
}

void createShellVsiteRecipes(System *sys, const ShellRecipe *shells, size_t nshell,
                             const VsiteRecipe *vsites, size_t nvsite)
{
    sys->release(sys->shellRecipes_);
    sys->release(sys->vsiteRecipes_);
    sys->nshellRecipe_ = nshell;
    sys->nvsiteRecipe_ = nvsite;
    sys->shellRecipes_ = sys->allocate<ShellRecipe>(std::max<size_t>(nshell, 1));
    sys->vsiteRecipes_ = sys->allocate<VsiteRecipe>(std::max<size_t>(nvsite, 1));
    sys->upload(sys->shellRecipes_, shells, nshell);
    sys->upload(sys->vsiteRecipes_, vsites, nvsite);
    sys->queue_->wait_and_throw();
}

void updateRecipeValues(System *sys, const double *values)
{
    if (sys->nrowRecipe_ == 0)
    {
        return;
    }
    sys->upload(sys->recipeValues_, values, sys->nrecipeValue_*sys->nrecipeReplica_);
    sys->recipeValuesChanged_ = true;
}

/*! \brief Make the sets of parameters of all replicas from the values, see createRowRecipes
 * \param[in] sys The system
 */
static void computeRecipeSets(System *sys)
{
    uint64_t nset = sys->nset_;
    uint64_t nv   = sys->nrecipeValue_;
    auto     dR   = sys->recipes_;
    auto     dF   = sys->recipeFactors_;
    auto     dV   = sys->recipeValues_;
    auto     dS   = sys->setValues_;
    if (nset > 0)
    {
        sys->queue_->parallel_for(sycl::range<1>(nset*sys->nrecipeReplica_), [=](sycl::id<1> id)
        {
            uint64_t w = id[0];
            uint64_t r = w / nset;
            (void) evaluateRecipe<SyclMath>(dR[w - r*nset], dV + r*nv, dF, dS + maxSetParameters*w);
        });
    }
}

/*! \brief Construct the coordinates of the charge generation, see createAcmVsites
 * One work item per molecule, which constructs its virtual sites in the
 * order of the processor, see constructStrided in forces/vsitehandler.cpp,
 * such that sites built from others are right too.
 * \param[in] sys The system
 */
static void computeAcmCoordinates(System *sys)
{
    auto     dM    = sys->acmMolecules_;
    auto     dF    = sys->acmFragments_;
    auto     dA    = sys->acmParticles_;
    auto     dX    = sys->acmX_;
    auto     dRef  = sys->acmRefX_;
    auto     dVS   = sys->acmVsiteStart_;
    auto     dVs   = sys->acmVsites_;
    auto     dS    = sys->setValues_;
    uint64_t nset  = sys->nset_;
    sys->queue_->parallel_for(sycl::range<1>(sys->nacmMolecule_), [=](sycl::id<1> id)
    {
        const AcmMolecule &m   = dM[id[0]];
        int64_t           *ref = dRef + 3*static_cast<uint64_t>(m.firstAtom);
        for (uint32_t v = dVS[id[0]]; v < dVS[id[0] + 1]; v++)
        {
            const AcmVsite   &av = dVs[v];
            const VsiteEntry &e  = av.entry;
            const double     *p  = dS + maxSetParameters*(av.replica*nset + av.set);
            float             pv[3];
            for (int k = 0; k < 3; k++)
            {
                // As VsiteRecipe
                pv[k] = av.index[k] >= 0 ? static_cast<float>(av.factor[k]*p[av.index[k]]) : e.p[k];
            }
            // As constructVsites
            uint32_t p0      = e.parent[0];
            float    r[3][3] = { { 0, 0, 0 }, { 0, 0, 0 }, { 0, 0, 0 } };
            for (uint32_t k = 1; k < e.nparent; k++)
            {
                for (int d = 0; d < 3; d++)
                {
                    r[k-1][d] = differenceNm(ref[3*e.parent[k]+d], ref[3*p0+d]);
                }
            }
            float xs[3];
            vsitePosition(e.kind, pv, r, xs);
            for (int d = 0; d < 3; d++)
            {
                ref[3*e.site+d] = ref[3*p0+d] + static_cast<int64_t>(sycl::rint(static_cast<double>(xs[d])*coordinateScaleNm));
            }
        }
        // The ACM particles in nm, as Coordinates::toRVec in single precision
        for (uint32_t f = m.firstFragment; f < m.firstFragment + m.nfragment; f++)
        {
            const AcmFragment &fr = dF[f];
            for (uint32_t i = fr.firstAcm; i < fr.firstAcm + fr.nacm; i++)
            {
                const int64_t *c = ref + 3*(fr.atomStart + dA[i]);
                for (int d = 0; d < 3; d++)
                {
                    dX[3*i+d] = static_cast<double>(static_cast<float>(c[d])*static_cast<float>(coordinateUnitNm));
                }
            }
        }
    });
}

/*! \brief Make the rows of the recipes from the sets, see computeRecipeSets
 * The rows as gpuPairEntry and gpuBondedEntry do on the processor, with the
 * charges of the shells as they are on the device, then the moving shells
 * and the virtual sites.
 * \param[in] sys     The system
 * \param[in] coulomb The Coulomb constant
 */
static void computeRecipeRows(System *sys, double coulomb)
{
    auto    *queue = sys->queue_;
    uint64_t nset  = sys->nset_;
    auto     dR    = sys->recipes_;
    auto     dS    = sys->setValues_;
    auto dRR = sys->rowRecipes_;
    auto dT  = sys->tables_;
    auto dQ  = sys->charges_;
    if (sys->nrowRecipe_ > 0)
    {
        queue->parallel_for(sycl::range<1>(sys->nrowRecipe_), [=](sycl::id<1> id)
        {
            const RowRecipe &rr = dRR[id[0]];
            const double    *p  = dS + maxSetParameters*(rr.replica*nset + rr.set);
            float           *t  = dT + 8*rr.row;
            if (rr.bonded)
            {
                rows::bondedRow(static_cast<BondedKernel>(rr.kernel), p, dR[rr.set].n, t);
            }
            else
            {
                double q = rr.charge >= 0 ? dQ[rr.charge] : 0.0;
                if (!rows::pairRow(static_cast<PairKernel>(rr.kernel), p, q, coulomb, t))
                {
                    // As stageRows
                    t[7] = 1;
                }
            }
        });
    }
    // The force constants of the moving shells, as addShellFragment in forces_simd.cpp
    auto dSR = sys->shellRecipes_;
    auto dMS = sys->movingShells_;
    if (sys->nshellRecipe_ > 0)
    {
        queue->parallel_for(sycl::range<1>(sys->nshellRecipe_), [=](sycl::id<1> id)
        {
            const ShellRecipe &sr    = dSR[id[0]];
            const double      *p     = dS + maxSetParameters*(sr.replica*nset + sr.set);
            double             alpha = p[polALPHA];
            double             q     = dQ[sr.charge];
            MovingShell       &ms    = dMS[sr.shell];
            ms.fcShell = (alpha > 0 && q != 0) ? q*q*coulomb/alpha : 0;
            ms.fcHyper = p[polFCHYPER];
            ms.rHyper  = p[polRHYPER];
        });
    }
    // The parameters of the virtual sites, as gpuVsiteEntry in forces_gpu.cpp
    auto dVR = sys->vsiteRecipes_;
    auto dVS = sys->vsites_;
    if (sys->nvsiteRecipe_ > 0)
    {
        queue->parallel_for(sycl::range<1>(sys->nvsiteRecipe_), [=](sycl::id<1> id)
        {
            const VsiteRecipe &vr = dVR[id[0]];
            const double      *p  = dS + maxSetParameters*(vr.replica*nset + vr.set);
            VsiteEntry        &ve = dVS[vr.vsite].entry;
            for (int k = 0; k < 3; k++)
            {
                if (vr.index[k] >= 0)
                {
                    ve.p[k] = static_cast<float>(vr.factor[k]*p[vr.index[k]]);
                }
            }
        });
    }
    sys->recipeValuesChanged_ = false;
}

void updateParameters(System *sys, const float *tables, const double *charges,
                      double coulomb, double epsilonr)
{
    auto *queue = sys->queue_;
    if (tables)
    {
        sys->upload(sys->tables_, tables, 8*sys->ntable_);
    }
    if (charges)
    {
        sys->upload(sys->charges_, charges, 4*sys->ncharge_);
    }
    // The sets of the recipes, of which the coordinates of the charge
    // generation may follow, and the rows after the charges of the shells
    bool recipes = sys->recipeValuesChanged_;
    if (recipes)
    {
        computeRecipeSets(sys);
    }
    if (recipes && sys->nchargeRecipe_ > 0)
    {
        // The charges are generated again if acmChargesChanged says so
        computeChargeRecipes(sys);
    }
    if (recipes && sys->acmVsitesChanged_ && sys->nacmMolecule_ > 0)
    {
        computeAcmCoordinates(sys);
        sys->acmVsitesChanged_ = false;
        sys->acmChanged_       = true;
    }
    generateAcmCharges(sys);
    if (recipes)
    {
        computeRecipeRows(sys, coulomb);
    }
    uint64_t nparam = sys->nparam_;
    if (nparam == 0)
    {
        return;
    }
    // The parameters of Coulomb, from the charges of the particles; the
    // other lists without parameter table keep those they were created with
    uint32_t nl  = sys->nlist_;
    auto     dL  = sys->lists_;
    auto     dB  = sys->blocks_;
    auto     dP  = sys->pairs_;
    auto     dPa = sys->params_;
    auto     dQ  = sys->charges_;
    queue->parallel_for(sycl::range<1>(nparam), [=](sycl::id<1> id)
    {
        uint64_t w = id[0];
        // The last list whose parameters start at or before w is that with
        // them, the lists without own parameters before it start there too
        const DeviceList &l = dL[findRange(dL, nl, w, &DeviceList::paramOffset)];
        if (l.bonded || l.kernel != static_cast<uint32_t>(PairKernel::Coulomb))
        {
            return;
        }
        // As gpuPairEntry
        const DevicePair &pe  = dP[l.entryOffset + w - l.paramOffset];
        uint64_t          off = dB[l.block].atomOffset;
        const double     *qi  = dQ + 4*(off + pe.i);
        const double     *qj  = dQ + 4*(off + pe.j);
        float            *p   = dPa + 8*w;
        for (int a = 0; a < 2; a++)
        {
            for (int c = 0; c < 2; c++)
            {
                int k     = 2*a + c;
                p[k]      = static_cast<float>(coulomb*qi[a]*qj[c]/epsilonr);
                p[4+k]    = static_cast<float>(combinedWidth(qi[2+a], qj[2+c]));
            }
        }
    });
}

void updateMovingShells(System *sys, const MovingShell *shells, size_t n)
{
    size_t total = 0;
    for (const auto &b : sys->blocksHost_)
    {
        total += b.nmovingShell;
    }
    if (n != total)
    {
        throw std::runtime_error("updateMovingShells: the number of moving shells changed");
    }
    sys->upload(sys->movingShells_, shells, n);
}

void updateVsites(System *sys, const VsiteEntry *vsites, size_t n)
{
    if (n != sys->vsitesHost_.size())
    {
        throw std::runtime_error("updateVsites: the number of virtual sites changed");
    }
    for (size_t v = 0; v < n; v++)
    {
        sys->vsitesHost_[v].entry = vsites[v];
    }
    sys->upload(sys->vsites_, sys->vsitesHost_.data(), n);
}

void constructVsites(System *sys)
{
    if (sys->nvsiteWork_ == 0)
    {
        return;
    }
    auto     dV   = sys->vsites_;
    auto     dB   = sys->blocks_;
    auto     dX   = sys->x_;
    uint32_t nv   = sys->nvsite_;
    sys->queue_->parallel_for(sycl::range<1>(sys->nvsiteWork_), [=](sycl::id<1> id)
    {
        uint64_t                 w = id[0];
        const DeviceVsite       &v = dV[findRange(dV, nv, w, &DeviceVsite::workOffset)];
        const DeviceSystemBlock &b = dB[v.block];
        uint64_t                 c = w - v.workOffset;
        uint64_t                 st = b.stride;
        const int64_t           *x  = dX + b.xOffset + c;
        uint32_t                 p0 = v.entry.parent[0];
        float r[3][3] = { { 0, 0, 0 }, { 0, 0, 0 }, { 0, 0, 0 } };
        for (uint32_t k = 1; k < v.entry.nparent; k++)
        {
            for (int d = 0; d < 3; d++)
            {
                r[k-1][d] = differenceNm(x[(3*v.entry.parent[k]+d)*st], x[(3*p0+d)*st]);
            }
        }
        float xs[3];
        vsitePosition(v.entry.kind, v.entry.p, r, xs);
        for (int d = 0; d < 3; d++)
        {
            dX[b.xOffset + (3*v.entry.site+d)*st + c] = x[(3*p0+d)*st] +
                static_cast<int64_t>(sycl::rint(static_cast<double>(xs[d])*coordinateScaleNm));
        }
    });
}

void computeAll(System *sys, double *const *energies)
{
    std::vector<char> use(sys->nblock_, 1);
    computeBlocks(sys, use.data(), energies);
}

void computeBlocks(System *sys, const char *use, double *const *energies, bool forces)
{
    if (energies == nullptr)
    {
        computeBlocksInPlace(sys, use, nullptr, forces);
        return;
    }
    std::vector<const double *> where(sys->nblock_, nullptr);
    computeBlocksInPlace(sys, use, where.data(), forces);
    for (size_t b = 0; b < sys->blocksHost_.size(); b++)
    {
        const auto &db = sys->blocksHost_[b];
        if (use[b] && energies[b] && where[b])
        {
            std::copy(where[b], where[b] + db.nenergy*pairEnergyTerms*db.stride, energies[b]);
        }
    }
}

static void computeBlocksLow(System *sys, const char *use, bool wantEnergies, const double **energies, bool forces,
                             double *keep = nullptr);

void computeBlocksInPlace(System *sys, const char *use, const double **energies, bool forces)
{
    computeBlocksLow(sys, use, energies != nullptr, energies, forces);
}

void computeBlocksKeep(System *sys, const char *use, bool forces, bool before)
{
    if (!sys->fixedStart_)
    {
        computeBlocksLow(sys, use, false, nullptr, forces);
        return;
    }
    // The energies go straight to where computeDeviations takes them
    computeBlocksLow(sys, use, true, nullptr, forces, before ? sys->energiesBefore_ : sys->energiesAfter_);
}

/*! \brief Compute the interactions of some blocks, as computeBlocksInPlace
 * \param[in]  sys          The system
 * \param[in]  use          Whether to compute each block
 * \param[in]  wantEnergies Whether to compute the energies
 * \param[out] energies     If not nullptr, where the energies of each block are on the host
 * \param[in]  forces       Whether to compute the forces too
 * \param[out] keep         If not nullptr, the energies are written here at fixedStart_
 *                          and stay on the GPU
 */
static void computeBlocksLow(System *sys, const char *use, bool wantEnergies, const double **energies, bool forces,
                             double *keep)
{
    auto *queue = sys->queue_;
    auto  dU    = sys->deviceMask(use, sys->nblock_);
    auto  dB    = sys->blocks_;
    auto  dL    = sys->lists_;
    auto  dX    = sys->x_;
    auto  dF    = sys->f_;
    auto  dFF   = sys->fFixed_;
    auto  dEnt  = sys->entryData();
    auto  dAS   = sys->atomStart_;
    auto  dAE   = sys->atomEntries_;
    auto  dEC   = sys->eChunk_;
    auto  dELs  = sys->eList_;
    auto  dE    = keep ? keep : sys->energies_;
    auto  dHLW  = sys->hintListWork_;
    auto  dHLS  = sys->hintListSlot_;
    auto  dHAW  = sys->hintAtomWork_;
    uint32_t nb = sys->nblock_, nl = sys->nlist_;
    // Without energies only the forces, which the relaxation needs, without
    // waiting for the GPU
    auto      &eStart       = sys->energyStartHost_;
    uint64_t   nenergy      = 0;
    if (wantEnergies)
    {
        // The energies of the blocks computed are stored one after the other,
        // such that only those are copied to the host
        eStart.resize(nb);
        for (uint32_t b = 0; b < nb; b++)
        {
            eStart[b] = nenergy;
            if (use[b])
            {
                const auto &db = sys->blocksHost_[b];
                nenergy += static_cast<uint64_t>(db.nenergy)*pairEnergyTerms*db.stride;
            }
        }
        if (!keep)
        {
            sys->upload(sys->energyStart_, eStart.data(), nb);
        }
        auto dES = keep ? sys->fixedStart_ : sys->energyStart_;
        uint64_t nwork = sys->nwork_, nlistSlot = sys->nlistSlot_;
        // The energy sums of each chunk of interactions in each slot
        if (nwork > 0)
        {
            queue->parallel_for(sycl::range<1>(nwork), [=](sycl::id<1> id)
            {
                uint64_t                 w  = id[0];
                const DeviceList        &l  = dL[findRangeFrom(dL, nl, w, &DeviceList::workOffset, dHLW)];
                const DeviceSystemBlock &b  = dB[l.block];
                if (!dU[l.block])
                {
                    return;
                }
                uint64_t                 local = w - l.workOffset;
                uint64_t                 first = local / b.nslot * energyChunk;
                uint64_t                 last  = sycl::min<uint64_t>(first + energyChunk, l.n);
                uint64_t                 c  = local % b.nslot;
                uint64_t                 st = b.stride;
                const int64_t           *x  = dX + b.xOffset + c;
                double sum[pairEnergyTerms] = { 0, 0, 0 };
                for (uint64_t e = first; e < last; e++)
                {
                    float en[pairEnergyTerms] = { 0, 0, 0 };
                    if (l.bonded)
                    {
                        const DeviceBond &be = dEnt.bonds[l.entryOffset + e];
                        float r[3][3], f[4][3];
                        for (int k = 1; k < 4; k++)
                        {
                            for (int d = 0; d < 3; d++)
                            {
                                r[k-1][d] = differenceNm(x[(3*be.a[k]+d)*st], x[(3*be.a[0]+d)*st]);
                            }
                        }
                        en[0] = bondedInteraction(static_cast<BondedKernel>(l.kernel),
                                                  entryParams(dEnt, l, l.entryOffset + e, dEnt.bondIndex), r, f);
                    }
                    else
                    {
                        const DevicePair &pe = dEnt.pairs[l.entryOffset + e];
                        float rsq = 0;
                        for (int d = 0; d < 3; d++)
                        {
                            float dd = differenceNm(x[(3*pe.i+d)*st], x[(3*pe.j+d)*st]);
                            rsq += dd*dd;
                        }
                        (void) pairInteraction(static_cast<PairKernel>(l.kernel),
                                               entryParams(dEnt, l, l.entryOffset + e, dEnt.pairIndex), rsq, en);
                    }
                    for (uint32_t t = 0; t < pairEnergyTerms; t++)
                    {
                        sum[t] += en[t];
                    }
                }
                for (uint32_t t = 0; t < l.nterm; t++)
                {
                    dEC[t*nwork + w] = sum[t];
                }
            });
        }
        // The sums of each list in each slot, in the order of the interactions
        if (nlistSlot > 0)
        {
            queue->parallel_for(sycl::range<1>(nlistSlot), [=](sycl::id<1> id)
            {
                uint64_t                 s = id[0];
                const DeviceList        &l = dL[findRangeFrom(dL, nl, s, &DeviceList::listSlotOffset, dHLS)];
                const DeviceSystemBlock &b = dB[l.block];
                uint64_t                 c = s - l.listSlotOffset;
                if (!dU[l.block])
                {
                    return;
                }
                for (uint32_t t = 0; t < pairEnergyTerms; t++)
                {
                    double sum = 0;
                    uint64_t nchunk = (l.n + energyChunk - 1)/energyChunk;
                    for (uint64_t k = 0; t < l.nterm && k < nchunk; k++)
                    {
                        sum += dEC[t*nwork + l.workOffset + k*b.nslot + c];
                    }
                    dELs[t*nlistSlot + s] = sum;
                }
            });
        }
        // The energy arrays of each block: the sums of its lists of the same
        // interaction type, in their order
        if (sys->nenergyWork_ > 0)
        {
            auto dEW = sys->energyWork_;
            auto dFL = sys->firstListDev_;
            queue->parallel_for(sycl::range<1>(sys->nenergyWork_), [=](sycl::id<1> id)
            {
                uint64_t w = id[0];
                uint32_t lo = 0, hi = nb;
                while (hi - lo > 1)
                {
                    uint32_t mid = (lo + hi)/2;
                    if (dEW[mid] <= w)
                    {
                        lo = mid;
                    }
                    else
                    {
                        hi = mid;
                    }
                }
                const DeviceSystemBlock &b = dB[lo];
                if (!dU[lo])
                {
                    return;
                }
                uint64_t local = w - dEW[lo];
                uint64_t c     = local % b.nslot;
                uint64_t et    = local / b.nslot;
                uint32_t eidx  = et / pairEnergyTerms;
                uint32_t t     = et % pairEnergyTerms;
                double   sum   = 0;
                for (uint32_t li = dFL[lo]; li < dFL[lo+1]; li++)
                {
                    if (dL[li].energyIndex == eidx)
                    {
                        sum += dELs[t*nlistSlot + dL[li].listSlotOffset + c];
                    }
                }
                dE[dES[lo] + (eidx*pairEnergyTerms + t)*b.stride + c] = sum;
            });
        }
    }
    // The forces on each particle in each slot, from all its interactions
    // in order, those of the interactions without shells also separately
    if (forces && sys->natomWork_ > 0)
    {
        queue->parallel_for(sycl::range<1>(sys->natomWork_), [=](sycl::id<1> id)
        {
            uint64_t                 w  = id[0];
            uint32_t                 bi = findRangeFrom(dB, nb, w, &DeviceSystemBlock::atomWorkOffset, dHAW);
            const DeviceSystemBlock &b  = dB[bi];
            uint64_t                 local = w - b.atomWorkOffset;
            uint64_t                 a  = local / b.nslot;
            uint64_t                 c  = local % b.nslot;
            uint64_t                 st = b.stride;
            const int64_t           *x  = dX + b.xOffset + c;
            if (!dU[bi])
            {
                return;
            }
            float fall[3] = { 0, 0, 0 }, ffix[3] = { 0, 0, 0 };
            for (uint64_t k = dAS[b.entryAtomOffset + a]; k < dAS[b.entryAtomOffset + a + 1]; k++)
            {
                float fa[3];
                bool  shells = entryForce(dAE[k], x, st, dEnt, b.listShift, fa);
                for (int m = 0; m < 3; m++)
                {
                    fall[m] += fa[m];
                    if (!shells)
                    {
                        ffix[m] += fa[m];
                    }
                }
            }
            for (int m = 0; m < 3; m++)
            {
                dF[b.xOffset + (3*a+m)*st + c]  = fall[m];
                dFF[b.xOffset + (3*a+m)*st + c] = ffix[m];
            }
        });
    }
    if (!wantEnergies || !energies)
    {
        return;
    }
    // The energies to the host
    double *hE = sys->energiesHost_.reserve(*queue, std::max<uint64_t>(sys->nenergyAll_, 1));
    if (nenergy > 0)
    {
        sys->queue_->memcpy(hE, dE, nenergy*sizeof(double));
    }
    queue->wait_and_throw();
    for (size_t b = 0; b < sys->blocksHost_.size(); b++)
    {
        energies[b] = use[b] ? hE + eStart[b] : nullptr;
    }
}

void resetShells(System *sys, const char *use)
{
    if (sys->nresetWork_ == 0)
    {
        return;
    }
    auto     dB = sys->blocks_;
    auto     dU = sys->deviceMask(use, sys->nblock_);
    auto     dX = sys->x_;
    auto     dR = sys->resets_;
    uint32_t nb = sys->nblock_;
    sys->queue_->parallel_for(sycl::range<1>(sys->nresetWork_), [=](sycl::id<1> id)
    {
        uint64_t                 w  = id[0];
        uint32_t                 bi = findRange(dB, nb, w, &DeviceSystemBlock::resetWorkOffset);
        const DeviceSystemBlock &b  = dB[bi];
        if (!dU[bi])
        {
            return;
        }
        uint64_t local = w - b.resetWorkOffset;
        uint64_t k     = local / b.nslot;
        uint64_t c     = local % b.nslot;
        uint32_t shell  = dR[2*(b.firstReset + k)];
        uint32_t parent = dR[2*(b.firstReset + k) + 1];
        int64_t *x = dX + b.xOffset + c;
        for (int m = 0; m < 3; m++)
        {
            x[(3*shell+m)*b.stride] = x[(3*parent+m)*b.stride];
        }
    });
}

void keepShells(System *sys, const char *use, bool restore)
{
    if (sys->nresetWork_ == 0)
    {
        return;
    }
    if (!sys->keptShells_)
    {
        sys->keptShells_ = sys->allocate<int64_t>(sys->nx_);
    }
    auto     dB = sys->blocks_;
    auto     dU = sys->deviceMask(use, sys->nblock_);
    auto     dX = sys->x_;
    auto     dK = sys->keptShells_;
    auto     dR = sys->resets_;
    auto     dC = sys->converged_;
    uint32_t nb = sys->nblock_;
    sys->queue_->parallel_for(sycl::range<1>(sys->nresetWork_), [=](sycl::id<1> id)
    {
        uint64_t                 w  = id[0];
        uint32_t                 bi = findRange(dB, nb, w, &DeviceSystemBlock::resetWorkOffset);
        const DeviceSystemBlock &b  = dB[bi];
        if (!dU[bi])
        {
            return;
        }
        uint64_t local = w - b.resetWorkOffset;
        uint64_t k     = local / b.nslot;
        uint64_t c     = local % b.nslot;
        uint32_t shell = dR[2*(b.firstReset + k)];
        uint64_t i     = b.xOffset + c + 3*shell*b.stride;
        // The shells that did not converge start from their parent next time
        uint64_t from  = dC[b.slotOffset + c] ? i : b.xOffset + c + 3*dR[2*(b.firstReset + k) + 1]*b.stride;
        for (int m = 0; m < 3; m++)
        {
            if (restore)
            {
                dX[i + m*b.stride] = dK[i + m*b.stride];
            }
            else
            {
                dK[i + m*b.stride] = dX[from + m*b.stride];
            }
        }
    });
}

void relaxBlocks(System *sys, const char *use, const char *relax,
                 const RelaxSettings &settings, size_t *nComputed, bool *notFinite,
                 size_t *nFailed)
{
    *nComputed = 0;
    *notFinite = false;
    *nFailed   = 0;
    auto *queue = sys->queue_;
    uint64_t nslotAll = sys->nslotAll_;
    if (nslotAll == 0)
    {
        return;
    }
    if (!sys->state_)
    {
        sys->state_ = sys->allocate<FragmentState>(2*nslotAll);
    }
    int maxIter = std::max(settings.maxIter, 1);
    if (maxIter > sys->maxIterAllocated_)
    {
        sys->doneHistory_      = sys->allocate<char>(maxIter*nslotAll);
        sys->maxIterAllocated_ = maxIter;
    }

    auto dB   = sys->blocks_;
    auto dU   = sys->deviceMask(use, sys->nblock_);
    auto dRel = sys->deviceMask(relax, 2*sys->nblock_);
    auto dX   = sys->x_;
    auto dF   = sys->f_;
    auto dFF  = sys->fFixed_;
    auto dDx  = sys->dx_;
    auto dSt  = static_cast<FragmentState *>(sys->state_);
    auto dDone = sys->done_;
    auto dHS   = sys->hintSlot_;
    auto dHSW  = sys->hintShellWork_;
    auto dHMP  = sys->hintMovingPart_;
    auto dHSR  = sys->hintShellRow_;
    auto dHist = sys->doneHistory_;
    auto dFS  = sys->forceShells_;
    auto dMS  = sys->movingShells_;
    auto dFr  = sys->fragments_;
    auto dMF  = sys->movingFrag_;
    auto dSR  = sys->shellRows_;
    auto dIS  = sys->isShell_;
    auto dPF2 = sys->partF2_;
    auto dPA  = sys->partAfter_;
    auto dPD  = sys->partDxF_;
    auto dVal = sys->valid_;
    auto dSS  = sys->shellStart_;
    auto dSE  = sys->shellEntries_;
    auto dEnt = sys->entryData();
    uint32_t nb = sys->nblock_;
    const double toler    = settings.msForceToler;
    const int    stallIterations = settings.stallIterations;
    const double maxDist2 = settings.maxShellDistance2;
    const double minStep  = settings.minStep;
    const double maxStep  = settings.maxStep;
    const double atWall   = 0.9999*maxDist2;
    const double wallAngle2 = settings.wallAngleToler*settings.wallAngleToler;
    queue->memset(dDx, 0, sys->nx_*sizeof(float));
    queue->memset(dVal, 0, sizeof(int));
    // The slots of the blocks not relaxed are done, the state of the others
    // starts with full steps and the largest force on the particles that
    // are not shells
    queue->parallel_for(sycl::range<1>(nslotAll), [=](sycl::id<1> id)
    {
        uint64_t                 s  = id[0];
        uint32_t                 bi = findRangeFrom(dB, nb, s, &DeviceSystemBlock::slotOffset, dHS);
        const DeviceSystemBlock &b  = dB[bi];
        uint64_t                 c  = s - b.slotOffset;
        dDone[s] = dU[bi] ? 0 : 1;
        for (uint32_t fi = 0; fi < 2; fi++)
        {
            double fmax = 0;
            if (dU[bi] && fi < b.nfragment)
            {
                const ShellFragment &fr = dFr[b.firstFragment + fi];
                const float         *f  = dF + b.xOffset + c;
                for (uint32_t a = fr.firstAtom; a < fr.firstAtom + fr.natom; a++)
                {
                    if (!dIS[b.atomOffset + a])
                    {
                        for (int m = 0; m < 3; m++)
                        {
                            fmax = sycl::fmax(fmax, static_cast<double>(sycl::fabs(f[(3*a+m)*b.stride])));
                        }
                    }
                }
            }
            dSt[2*s + fi] = { 1.0, DBL_MAX, 0.0, 0, fmax, DBL_MAX, 0 };
        }
    });
    int  iter       = 1;
    bool allDone    = false;
    int  nsubmitted = 0;
    const int checkEvery = 4;
    // Without polling, all iterations, those after all slots are done do nothing
    char *hDone = settings.poll ? sys->doneHost_.reserve(*queue, maxIter*nslotAll) : nullptr;
    while (iter < settings.maxIter && !allDone)
    {
        // The squared forces of the shells and the step products, see relaxation.cpp
        if (sys->nshellWork_ > 0)
        {
            queue->parallel_for(sycl::range<1>(sys->nshellWork_), [=](sycl::id<1> id)
            {
                uint64_t                 w  = id[0];
                uint32_t                 bi = findRangeFrom(dB, nb, w, &DeviceSystemBlock::shellWorkOffset, dHSW);
                const DeviceSystemBlock &b  = dB[bi];
                uint64_t local = w - b.shellWorkOffset;
                uint64_t c     = local % b.nslot;
                uint64_t k     = local / b.nslot;
                if (dDone[b.slotOffset + c])
                {
                    return;
                }
                uint64_t     st = b.stride;
                const float *f  = dF + b.xOffset + c;
                if (k < b.nforceShell)
                {
                    const ForceShell &fs = dFS[b.firstForceShell + k];
                    float  fv[3] = { f[(3*fs.shell)*st], f[(3*fs.shell+1)*st], f[(3*fs.shell+2)*st] };
                    double f2    = fv[0]*fv[0] + fv[1]*fv[1] + fv[2]*fv[2];
                    if (fs.wall)
                    {
                        const int64_t *x = dX + b.xOffset + c;
                        float d[3];
                        for (int m = 0; m < 3; m++)
                        {
                            d[m] = toNm(difference(x[(3*fs.shell+m)*st], x[(3*fs.parent+m)*st]));
                        }
                        double dx2 = d[0]*d[0] + d[1]*d[1] + d[2]*d[2];
                        double fdx = fv[0]*d[0] + fv[1]*d[1] + fv[2]*d[2];
                        if (dx2 >= atWall && fdx > 0)
                        {
                            // Converged if small compared to the force, from
                            // the cross product, see shellForce2
                            double c0 = static_cast<double>(fv[1])*d[2] - static_cast<double>(fv[2])*d[1];
                            double c1 = static_cast<double>(fv[2])*d[0] - static_cast<double>(fv[0])*d[2];
                            double c2 = static_cast<double>(fv[0])*d[1] - static_cast<double>(fv[1])*d[0];
                            double fperp2 = (c0*c0 + c1*c1 + c2*c2)/dx2;
                            f2 = (fperp2 <= wallAngle2*f2) ? 0.0 : fperp2;
                        }
                    }
                    dPF2[b.forceShellPartOffset + k*b.nslot + c] = f2;
                }
                else
                {
                    uint64_t     km = k - b.nforceShell;
                    uint32_t     a  = dMS[b.firstMovingShell + km].shell;
                    const float *dx = dDx + b.xOffset + c;
                    float        p  = 0;
                    for (int m = 0; m < 3; m++)
                    {
                        p += dx[(3*a+m)*st]*f[(3*a+m)*st];
                    }
                    dPA[b.movingShellPartOffset + km*b.nslot + c] = p;
                }
            });
        }
        // The convergence and the step of each fragment in each slot, and
        // which slots are done after this iteration in the history, which
        // the host reads, without a copy per iteration
        char *dH = dHist + nsubmitted*nslotAll;
        queue->parallel_for(sycl::range<1>(nslotAll), [=](sycl::id<1> id)
        {
            uint64_t s = id[0];
            if (dDone[s])
            {
                dH[s] = dDone[s];
                return;
            }
            uint32_t                 bi = findRangeFrom(dB, nb, s, &DeviceSystemBlock::slotOffset, dHS);
            const DeviceSystemBlock &b  = dB[bi];
            uint64_t                 c  = s - b.slotOffset;
            uint64_t                 st = b.stride;
            const float             *f  = dF + b.xOffset + c;
            bool converged = true;
            for (uint32_t fi = 0; fi < b.nfragment; fi++)
            {
                FragmentState &state = dSt[2*s + fi];
                if (!dRel[2*bi + fi])
                {
                    continue;
                }
                const ShellFragment &fr  = dFr[b.firstFragment + fi];
                uint64_t             fs0 = fr.firstForceShell - b.firstForceShell;
                uint64_t             ms0 = fr.firstMovingShell - b.firstMovingShell;
                if (state.moving)
                {
                    double dxForce = 0;
                    for (uint32_t k = 0; k < fr.nmovingShell; k++)
                    {
                        dxForce += dPD[b.movingShellPartOffset + (ms0 + k)*b.nslot + c];
                    }
                    state.dxfPrev = dxForce;
                    state.moving  = 0;
                }
                double msForce = 0;
                for (uint32_t k = 0; k < fr.nforceShell; k++)
                {
                    msForce += dPF2[b.forceShellPartOffset + (fs0 + k)*b.nslot + c];
                }
                if (fr.nShell > 0)
                {
                    msForce /= fr.nShell;
                }
                bool stalled = msForce > toler && msForce >= 0.5*state.msPrev;
                bool atFloor = false;
                if (stalled)
                {
                    double fmax = state.fmaxRest;
                    for (uint32_t k = 0; k < fr.nforceShell + fr.nmovingShell; k++)
                    {
                        uint32_t a = (k < fr.nforceShell) ? dFS[fr.firstForceShell + k].shell :
                            dMS[fr.firstMovingShell + k - fr.nforceShell].shell;
                        for (int m = 0; m < 3; m++)
                        {
                            fmax = sycl::fmax(fmax, static_cast<double>(sycl::fabs(f[(3*a+m)*st])));
                        }
                    }
                    double floor = 32*FLT_EPSILON*fmax;
                    atFloor = msForce <= floor*floor;
                }
                if (msForce <= toler || atFloor)
                {
                    continue;
                }
                // Without halving the mean square force for stallIterations
                // steps, the shells stop without having converged, see updateShells
                if (msForce <= 0.5*state.reference)
                {
                    state.reference = msForce;
                    state.stall     = 0;
                }
                else
                {
                    state.stall += 1;
                }
                if (stallIterations > 0 && state.stall >= stallIterations)
                {
                    continue;
                }
                converged = false;
                double after = 0;
                for (uint32_t k = 0; k < fr.nmovingShell; k++)
                {
                    after += dPA[b.movingShellPartOffset + (ms0 + k)*b.nslot + c];
                }
                if (state.dxfPrev > 0)
                {
                    double before = state.dxfPrev;
                    double scale  = (before > after) ? before/(before - after) : 1.0;
                    state.step    = sycl::clamp(state.step*scale, minStep, maxStep);
                }
                if (msForce > state.msPrev)
                {
                    state.step = sycl::fmin(state.step, 0.5);
                }
                state.msPrev = msForce;
                state.moving = 1;
            }
            dDone[s] = converged ? 1 : 0;
            dH[s]    = dDone[s];
        });
        nsubmitted++;
        // Move the shells of the fragments that did not converge
        if (sys->nmovingPart_ > 0)
        {
            queue->parallel_for(sycl::range<1>(sys->nmovingPart_), [=](sycl::id<1> id)
            {
                uint64_t                 w  = id[0];
                uint32_t                 bi = findRangeFrom(dB, nb, w, &DeviceSystemBlock::movingShellPartOffset, dHMP);
                const DeviceSystemBlock &b  = dB[bi];
                uint64_t local = w - b.movingShellPartOffset;
                uint64_t c     = local % b.nslot;
                uint64_t k     = local / b.nslot;
                uint64_t s     = b.slotOffset + c;
                const FragmentState &state = dSt[2*s + dMF[b.firstMovingShell + k]];
                if (dDone[s] || !state.moving)
                {
                    return;
                }
                uint64_t           st   = b.stride;
                int64_t           *x    = dX + b.xOffset + c;
                const float       *f    = dF + b.xOffset + c;
                float             *dx   = dDx + b.xOffset + c;
                const MovingShell &ms   = dMS[b.firstMovingShell + k];
                double             step = state.step;
                bool               valid = true;
                int64_t *xs = x + 3*ms.shell*st;
                float    fs[3];
                int64_t  xold[3];
                for (int m = 0; m < 3; m++)
                {
                    fs[m]   = f[(3*ms.shell+m)*st];
                    xold[m] = xs[m*st];
                }
                double disp[3];
                double rcs = 0;
                if (ms.fcHyper == 0)
                {
                    double scale = step*(1/ms.fcShell);
                    for (int m = 0; m < 3; m++)
                    {
                        disp[m] = scale*fs[m];
                    }
                }
                else
                {
                    float r[3];
                    for (int m = 0; m < 3; m++)
                    {
                        r[m] = toNm(difference(xs[m*st], x[(3*ms.core+m)*st]));
                    }
                    double d      = sycl::sqrt(r[0]*r[0] + r[1]*r[1] + r[2]*r[2]);
                    rcs           = d;
                    double invd   = 1/((d == 0) ? 1 : d);
                    double drh    = d - ms.rHyper;
                    bool   hyper  = d > ms.rHyper;
                    double kperp  = hyper ? ms.fcShell + 4*ms.fcHyper*drh*drh*drh*invd : ms.fcShell;
                    double kpar   = hyper ? ms.fcShell + 12*ms.fcHyper*drh*drh : ms.fcShell;
                    double ikperp = 1/kperp;
                    float  fr2    = fs[0]*r[0] + fs[1]*r[1] + fs[2]*r[2];
                    double fpar   = fr2*invd*invd*(1/kpar - ikperp);
                    for (int m = 0; m < 3; m++)
                    {
                        disp[m] = step*(fs[m]*ikperp + fpar*r[m]);
                    }
                }
                if (ms.fcHyper > 0 && rcs == 0)
                {
                    float  fnorm = sycl::sqrt(fs[0]*fs[0] + fs[1]*fs[1] + fs[2]*fs[2]);
                    double scale = step*fscale(ms.fcShell, ms.rHyper > 0 ? 0 : ms.fcHyper, fnorm);
                    for (int m = 0; m < 3; m++)
                    {
                        xs[m*st] = move(xs[m*st], -(fs[m]*scale), &valid);
                    }
                }
                else
                {
                    for (int m = 0; m < 3; m++)
                    {
                        xs[m*st] = move(xs[m*st], disp[m], &valid);
                    }
                }
                bool clipped = false;
                if (ms.wall)
                {
                    float d[3];
                    for (int m = 0; m < 3; m++)
                    {
                        d[m] = toNm(difference(xs[m*st], x[(3*ms.parent+m)*st]));
                    }
                    float dx2 = d[0]*d[0] + d[1]*d[1] + d[2]*d[2];
                    if (dx2 > maxDist2)
                    {
                        float scale = sycl::sqrt(static_cast<float>(maxDist2/dx2));
                        for (int m = 0; m < 3; m++)
                        {
                            double xp = static_cast<double>(x[(3*ms.parent+m)*st])*coordinateUnitNm;
                            xs[m*st] = toUnits(xp + scale*d[m]);
                        }
                        clipped = true;
                    }
                }
                float  dxj[3]  = { 0, 0, 0 };
                double dxForce = 0;
                if (!clipped)
                {
                    for (int m = 0; m < 3; m++)
                    {
                        dxj[m] = toNm(difference(xs[m*st], xold[m]));
                    }
                    dxForce = dxj[0]*fs[0] + dxj[1]*fs[1] + dxj[2]*fs[2];
                }
                for (int m = 0; m < 3; m++)
                {
                    dx[(3*ms.shell+m)*st] = dxj[m];
                }
                dPD[b.movingShellPartOffset + k*b.nslot + c] = dxForce;
                if (!valid)
                {
                    *dVal = 1;
                }
            });
        }
        // The forces on the shells: those of the interactions without shells
        // and those of their interactions with shells
        if (sys->nshellRowWork_ > 0)
        {
            // A work-group of a wavefront lies in a single shell row, as the
            // rows are padded to whole wavefronts. The block and the row
            // follow from the group alone, such that the compiler knows they
            // are the same for all its work items and can load the data of
            // the shell and its interactions once for the wavefront.
            queue->parallel_for(sycl::nd_range<1>(sys->nshellRowWork_, wavefrontSize),
                                [=](sycl::nd_item<1> item)
            {
                uint64_t                 w0 = item.get_group(0)*static_cast<uint64_t>(wavefrontSize);
                uint32_t                 bi = findRangeFrom(dB, nb, w0, &DeviceSystemBlock::shellRowWorkOffset, dHSR);
                const DeviceSystemBlock &b  = dB[bi];
                uint64_t local0 = w0 - b.shellRowWorkOffset;
                uint64_t row    = local0 / b.nslotPad;
                uint64_t c      = local0 % b.nslotPad + item.get_local_id(0);
                if (c >= b.nslot || dDone[b.slotOffset + c])
                {
                    return;
                }
                uint64_t       a    = dSR[b.firstShellRow + row];
                uint64_t       st   = b.stride;
                uint64_t       base = b.xOffset + 3*a*st + c;
                const int64_t *x    = dX + b.xOffset + c;
                float          fa[3];
                int64_t        own[3];
                for (int m = 0; m < 3; m++)
                {
                    fa[m]  = dFF[base + m*st];
                    own[m] = x[(3*a+m)*st];
                }
                for (uint64_t e = dSS[b.entryAtomOffset + a]; e < dSS[b.entryAtomOffset + a + 1]; e++)
                {
                    float fe[3];
                    entryForceOn(dSE[e], own, x, st, dEnt, b.listShift, fe);
                    for (int m = 0; m < 3; m++)
                    {
                        fa[m] += fe[m];
                    }
                }
                for (int m = 0; m < 3; m++)
                {
                    dF[base + m*st] = fa[m];
                }
            });
        }
        iter++;
        if (settings.poll && (nsubmitted % checkEvery == 0 || iter >= settings.maxIter))
        {
            queue->memcpy(hDone, dDone, nslotAll);
            queue->wait_and_throw();
            allDone = std::count(hDone, hDone + nslotAll, 0) == 0;
        }
    }
    // The slots whose shells did not converge: those left after the last
    // iteration and those that stopped without progress
    auto dNF = sys->nfailed_;
    auto dCv = sys->converged_;
    queue->memset(dNF, 0, sizeof(uint64_t));
    queue->parallel_for(sycl::range<1>(nslotAll), [=](sycl::id<1> id)
    {
        uint64_t s      = id[0];
        bool     failed = !dDone[s];
        for (uint32_t fi = 0; fi < 2; fi++)
        {
            failed = failed || (stallIterations > 0 && dSt[2*s + fi].stall >= stallIterations);
        }
        uint32_t bi = findRangeFrom(dB, nb, s, &DeviceSystemBlock::slotOffset, dHS);
        if (dU[bi])
        {
            dCv[s] = failed ? 0 : 1;
        }
        if (failed)
        {
            sycl::atomic_ref<uint64_t, sycl::memory_order::relaxed, sycl::memory_scope::device,
                             sycl::access::address_space::global_space> count(*dNF);
            count.fetch_add(1);
        }
    });
    if (!settings.poll)
    {
        // The statistics stay on the GPU, see relaxStatistics
        auto dRS = sys->relaxStats_;
        int  nit = nsubmitted;
        queue->parallel_for(sycl::range<1>(nslotAll), [=](sycl::id<1> id)
        {
            uint64_t s = id[0];
            uint64_t n = 0;
            for (int i = 0; i < nit; i++)
            {
                n += dHist[i*nslotAll + s] ? 0 : 1;
            }
            using Counter = sycl::atomic_ref<uint64_t, sycl::memory_order::relaxed, sycl::memory_scope::device,
                                             sycl::access::address_space::global_space>;
            if (n > 0)
            {
                Counter(dRS[0]).fetch_add(n);
            }
            if (s == 0)
            {
                Counter(dRS[1]).fetch_add(*dNF);
                if (*dVal)
                {
                    Counter(dRS[2]).fetch_or(1);
                }
            }
        });
        return;
    }
    // The slots computed in each iteration, those that did not converge
    int      validFlag = 0;
    uint64_t failed    = 0;
    queue->memcpy(hDone, dHist, nsubmitted*nslotAll);
    queue->memcpy(&validFlag, dVal, sizeof(int));
    queue->memcpy(&failed, dNF, sizeof(uint64_t));
    queue->wait_and_throw();
    *nFailed = failed;
    for (int i = 0; i < nsubmitted; i++)
    {
        *nComputed += std::count(hDone + i*nslotAll, hDone + (i+1)*nslotAll, 0);
    }
    *notFinite = validFlag != 0;
}

void relaxStatistics(System *sys, RelaxStatistics *stats)
{
    uint64_t h[3] = { 0, 0, 0 };
    if (sys->relaxStats_)
    {
        sys->queue_->memcpy(h, sys->relaxStats_, sizeof(h));
        sys->queue_->wait_and_throw();
        sys->queue_->memset(sys->relaxStats_, 0, sizeof(h)).wait();
    }
    stats->computed  = h[0];
    stats->failed    = h[1];
    stats->notFinite = h[2] != 0;
}

const int *deviceNotFinite(System *sys)
{
    return sys->valid_;
}

void copyCoordinates(System *sys, size_t from, size_t to, const uint32_t *pairs, size_t n,
                     const uint32_t *slots)
{
    const auto &bf = sys->blocksHost_[from];
    const auto &bt = sys->blocksHost_[to];
    if (n == 0 || bt.nslot == 0)
    {
        return;
    }
    // The pairs and the slots, freed at the end of the call
    auto dPairs = sycl::malloc_device<uint32_t>(2*n + bt.nslot, *sys->queue_);
    if (!dPairs)
    {
        auto device = sys->queue_->get_device();
        throw OutOfMemory(device.get_info<sycl::info::device::name>(),
                          device.get_info<sycl::info::device::global_mem_size>(),
                          sys->allocatedBytes_, (2*n + bt.nslot)*sizeof(uint32_t));
    }
    sys->queue_->memcpy(dPairs, pairs, 2*n*sizeof(uint32_t));
    sys->queue_->memcpy(dPairs + 2*n, slots, bt.nslot*sizeof(uint32_t));
    auto     dX    = sys->x_;
    uint64_t xf    = bf.xOffset, xt = bt.xOffset, sf = bf.stride, stt = bt.stride, nslot = bt.nslot;
    sys->queue_->parallel_for(sycl::range<1>(n*nslot), [=](sycl::id<1> id)
    {
        uint64_t k  = id[0] / nslot;
        uint64_t c  = id[0] % nslot;
        uint64_t cf = dPairs[2*n + c];
        for (int m = 0; m < 3; m++)
        {
            dX[xt + (3*dPairs[2*k+1]+m)*stt + c] = dX[xf + (3*dPairs[2*k]+m)*sf + cf];
        }
    });
    sys->queue_->wait_and_throw();
    sycl::free(dPairs, *sys->queue_);
}

void saveCoordinates(System *sys, size_t block)
{
    const auto &b = sys->blocksHost_[block];
    size_t      n = 3*b.natom*b.stride;
    if (!sys->saved_[block])
    {
        sys->saved_[block] = sys->allocate<int64_t>(n);
    }
    sys->queue_->memcpy(sys->saved_[block], sys->x_ + b.xOffset, n*sizeof(int64_t));
}

void restoreCoordinates(System *sys, size_t block)
{
    const auto &b = sys->blocksHost_[block];
    if (sys->saved_[block])
    {
        sys->queue_->memcpy(sys->x_ + b.xOffset, sys->saved_[block], 3*b.natom*b.stride*sizeof(int64_t));
    }
}

void downloadBlock(System *sys, size_t block, int64_t *x, float *f)
{
    const auto &b = sys->blocksHost_[block];
    size_t      n = 3*b.natom*b.stride;
    sys->queue_->memcpy(x, sys->x_ + b.xOffset, n*sizeof(int64_t));
    if (f)
    {
        sys->queue_->memcpy(f, sys->f_ + b.xOffset, n*sizeof(float));
    }
    sys->queue_->wait_and_throw();
}


//! The largest number of energy types of a block in computeDeviations
constexpr uint32_t maxDeviationTypes = 32;

/*! \brief The loss of a deviation, as loss in loss_function.h with a scaling factor of one
 * \param[in] lf     The loss function, see DeviationSystem::lossFunction
 * \param[in] deltaE The deviation
 * \return the loss
 */
static inline double deviationLoss(int lf, double deltaE)
{
#if defined(__clang__)
#pragma clang fp contract(off)
#endif
    switch (lf)
    {
    case 0:
        return deltaE*deltaE;
    case 1:
        return sycl::fabs(deltaE);
    case 2:
        return 1.0*(sycl::sqrt(deltaE*deltaE+1) - 1);
    default:
        {
            double sdea = sycl::sqrt(deltaE*deltaE+1);
            return 1.0*(1 - sdea + deltaE*sycl::log(deltaE + sdea));
        }
    }
}

/*! \brief The energy of a type of a block in a slot, as GpuTrainingSystem::typeEnergies makes it
 * \param[in] t     The type
 * \param[in] e     The energy arrays of the block
 * \param[in] st    The stride of the block
 * \param[in] s     The slot
 * \param[in] npot  The number of energy arrays of the block
 * \return the energy
 */
static inline double deviationTypeEnergy(const DeviationType &t, const double *e, uint64_t st,
                                         uint64_t s, uint32_t npot)
{
#if defined(__clang__)
#pragma clang fp contract(off)
#endif
    switch (t.kind)
    {
    case deviationTerm:
        return 0.0 + e[t.source*st + s];
    case deviationEpot:
        {
            double v = 0.0;
            for (uint32_t p = 0; p < npot; p++)
            {
                const double *ep = e + static_cast<uint64_t>(p)*pairEnergyTerms*st + s;
                v = v + (ep[0] + ep[st] + ep[2*st]);
            }
            return v;
        }
    case deviationExchange:
        {
            double x = t.source >= 0 ? 0.0 + e[t.source*st + s] : 0.0;
            return x + (0.0 + e[t.source2*st + s]);
        }
    default:
        return 0.0;
    }
}

void createDeviations(System *sys, const DeviationSystem &d)
{
    uint32_t nb    = sys->nblock_;
    uint64_t nconf = d.blockFirstConf[nb];
    // The energies at fixed places, as the blocks one after the other
    if (!sys->fixedStart_)
    {
        std::vector<uint64_t> start(nb + 1, 0);
        for (uint32_t b = 0; b < nb; b++)
        {
            const auto &db = sys->blocksHost_[b];
            start[b+1] = start[b] + static_cast<uint64_t>(db.nenergy)*pairEnergyTerms*db.stride;
        }
        sys->fixedStart_     = sys->allocate<uint64_t>(nb + 1);
        sys->energiesAfter_  = sys->allocate<double>(start[nb]);
        sys->energiesBefore_ = sys->allocate<double>(start[nb]);
        sys->upload(sys->fixedStart_, start.data(), nb + 1);
        sys->queue_->wait_and_throw();
    }
    // The buffers of an earlier call
    sys->queue_->wait_and_throw();
    sys->release(sys->devTypes_);
    sys->release(sys->devBlockFirstType_);
    sys->release(sys->devSlotOfConf_);
    sys->release(sys->devBlockFirstConf_);
    sys->release(sys->devDimers_);
    sys->release(sys->devMatches_);
    sys->release(sys->devTargets_);
    sys->release(sys->devWeights_);
    sys->release(sys->devEqm_);
    sys->release(sys->devHave_);
    sys->release(sys->devMonomers_);
    sys->release(sys->devEpotWeights_);
    sys->release(sys->devEpotEqm_);
    sys->release(sys->devEpotHave_);
    sys->release(sys->devDimerSums_);
    sys->release(sys->devTerms_);
    sys->release(sys->devEpotTerms_);
    sys->dev_               = d;
    sys->devTypes_          = sys->allocate<DeviationType>(d.blockFirstType[nb]);
    sys->devBlockFirstType_ = sys->allocate<uint32_t>(nb + 1);
    sys->devSlotOfConf_     = sys->allocate<uint32_t>(nconf);
    sys->devBlockFirstConf_ = sys->allocate<uint64_t>(nb + 1);
    sys->devDimers_         = sys->allocate<DeviationDimer>(d.ndimer);
    sys->devMatches_        = sys->allocate<int32_t>(d.nmatch);
    sys->devTargets_        = sys->allocate<int32_t>(d.ntarget);
    sys->devWeights_        = sys->allocate<double>(d.nweight);
    sys->devEqm_            = sys->allocate<double>(d.nweight*d.ntarget);
    sys->devHave_           = sys->allocate<char>(d.nweight*d.ntarget);
    sys->devMonomers_       = sys->allocate<DeviationMonomer>(d.nmonomer);
    sys->devEpotWeights_    = sys->allocate<double>(d.nepot);
    sys->devEpotEqm_        = sys->allocate<double>(d.nepot);
    sys->devEpotHave_       = sys->allocate<char>(d.nepot);
    sys->devDimerSums_      = sys->allocate<double>(2*d.ndimer*d.ntarget);
    sys->devNDimerConf_     = d.ndimer > 0 ? d.dimers[d.ndimer-1].confOffset + d.dimers[d.ndimer-1].nconf : 0;
    sys->devNMonoConf_      = d.nmonomer > 0 ? d.monomers[d.nmonomer-1].confOffset + d.monomers[d.nmonomer-1].nconf : 0;
    sys->devTerms_          = sys->allocate<double>(2*sys->devNDimerConf_*d.ntarget);
    sys->devEpotTerms_      = sys->allocate<double>(2*sys->devNMonoConf_);
    sys->upload(sys->devTypes_, d.types, d.blockFirstType[nb]);
    sys->upload(sys->devBlockFirstType_, d.blockFirstType, nb + 1);
    sys->upload(sys->devSlotOfConf_, d.slotOfConf, nconf);
    sys->upload(sys->devBlockFirstConf_, d.blockFirstConf, nb + 1);
    sys->upload(sys->devDimers_, d.dimers, d.ndimer);
    sys->upload(sys->devMatches_, d.matches, d.nmatch);
    sys->upload(sys->devTargets_, d.targets, d.ntarget);
    sys->upload(sys->devWeights_, d.weights, d.nweight);
    sys->upload(sys->devEqm_, d.eqm, d.nweight*d.ntarget);
    sys->upload(sys->devHave_, d.have, d.nweight*d.ntarget);
    sys->upload(sys->devMonomers_, d.monomers, d.nmonomer);
    sys->upload(sys->devEpotWeights_, d.epotWeights, d.nepot);
    sys->upload(sys->devEpotEqm_, d.epotEqm, d.nepot);
    sys->upload(sys->devEpotHave_, d.epotHave, d.nepot);
    sys->queue_->wait_and_throw();
}

void keepEnergies(System *sys, const char *use, bool before)
{
    if (!sys->fixedStart_ || sys->nenergyAll_ == 0)
    {
        return;
    }
    // The energies of the blocks computed last are one after the other, at
    // energyStart_, see computeBlocksInPlace
    auto     dU  = sys->deviceMask(use, sys->nblock_);
    auto     dE  = sys->energies_;
    auto     dES = sys->energyStart_;
    auto     dFS = sys->fixedStart_;
    auto     dK  = before ? sys->energiesBefore_ : sys->energiesAfter_;
    uint32_t nb  = sys->nblock_;
    sys->queue_->parallel_for(sycl::range<1>(sys->nenergyAll_), [=](sycl::id<1> id)
    {
        uint64_t w  = id[0];
        uint32_t lo = 0, hi = nb;
        while (hi - lo > 1)
        {
            uint32_t mid = (lo + hi)/2;
            if (dFS[mid] <= w)
            {
                lo = mid;
            }
            else
            {
                hi = mid;
            }
        }
        if (dU[lo])
        {
            dK[w] = dE[dES[lo] + (w - dFS[lo])];
        }
    });
}

void computeDeviations(System *sys, const char *use, double *dimerSums, double *epotTerms)
{
    const auto &d       = sys->dev_;
    auto        dU      = sys->deviceMask(use, sys->nblock_);
    auto       *queue   = sys->queue_;
    auto        dB      = sys->blocks_;
    auto        dT      = sys->devTypes_;
    auto        dBFT    = sys->devBlockFirstType_;
    auto        dSOC    = sys->devSlotOfConf_;
    auto        dBFC    = sys->devBlockFirstConf_;
    auto        dD      = sys->devDimers_;
    auto        dM      = sys->devMatches_;
    auto        dTg     = sys->devTargets_;
    auto        dW      = sys->devWeights_;
    auto        dQ      = sys->devEqm_;
    auto        dH      = sys->devHave_;
    auto        dMo     = sys->devMonomers_;
    auto        dEW     = sys->devEpotWeights_;
    auto        dEQ     = sys->devEpotEqm_;
    auto        dEH     = sys->devEpotHave_;
    auto        dTe     = sys->devTerms_;
    auto        dETe    = sys->devEpotTerms_;
    auto        dDS     = sys->devDimerSums_;
    auto        dFS     = sys->fixedStart_;
    auto        dA      = sys->energiesAfter_;
    auto        dBe     = sys->energiesBefore_;
    uint32_t    ntarget = d.ntarget;
    uint32_t    nmono   = d.nmonomer, ndimer = d.ndimer;
    uint64_t    ndc     = sys->devNDimerConf_, nmc = sys->devNMonoConf_;
    int         lf      = d.lossFunction;
    int32_t     allElecType = d.allElecType, exchIndType = d.exchIndType, inductionType = d.inductionType;
    // The weight and loss of each target of each conformation of each dimer,
    // with the energies made as compactInteraction makes them
    if (ndc > 0 && ntarget > 0)
    {
        queue->parallel_for(sycl::range<1>(ndc), [=](sycl::id<1> id)
        {
#if defined(__clang__)
#pragma clang fp contract(off)
#endif
            uint64_t gc = id[0];
            uint32_t lo = 0, hi = ndimer;
            while (hi - lo > 1)
            {
                uint32_t mid = (lo + hi)/2;
                if (dD[mid].confOffset <= gc)
                {
                    lo = mid;
                }
                else
                {
                    hi = mid;
                }
            }
            const DeviationDimer    &dm  = dD[lo];
            if (!dU[dm.dimer])
            {
                return;
            }
            uint64_t                 c   = gc - dm.confOffset;
            const DeviceSystemBlock &bd  = dB[dm.dimer];
            uint32_t                 t0  = dBFT[dm.dimer];
            uint32_t                 nt  = dBFT[dm.dimer + 1] - t0;
            uint32_t                 sD  = dSOC[dBFC[dm.dimer] + c];
            const double            *aD  = dA + dFS[dm.dimer];
            const double            *eD  = dBe + dFS[dm.dimer];
            int indEntry = -1, corrEntry = -1;
            for (uint32_t i = 0; i < nt; i++)
            {
                indEntry  = (dT[t0 + i].role & deviationInduction) ? static_cast<int>(i) : indEntry;
                corrEntry = (dT[t0 + i].role & deviationInductionCorrection) ? static_cast<int>(i) : corrEntry;
            }
            double act[maxDeviationTypes], bef[maxDeviationTypes];
            // subtractMonomers
            for (uint32_t i = 0; i < nt; i++)
            {
                const DeviationType &ty = dT[t0 + i];
                double edim = deviationTypeEnergy(ty, aD, bd.stride, sD, bd.nenergy);
                double edb4 = deviationTypeEnergy(ty, eD, bd.stride, sD, bd.nenergy);
                for (int j = 0; j < 2; j++)
                {
                    uint32_t m  = dm.monomer[j];
                    int32_t  mi = dM[dm.matchOffset + 2*i + j];
                    if (m != noBlock && mi >= 0 && c < dBFC[m + 1] - dBFC[m])
                    {
                        edim -= deviationTypeEnergy(dT[dBFT[m] + mi], dA + dFS[m], dB[m].stride,
                                                    dSOC[dBFC[m] + c], dB[m].nenergy);
                    }
                }
                for (int j = 0; j < 2; j++)
                {
                    uint32_t m  = dm.monomer[j];
                    int32_t  mi = dM[dm.matchOffset + 2*i + j];
                    if (m != noBlock && mi >= 0 && c < dBFC[m + 1] - dBFC[m])
                    {
                        edb4 -= deviationTypeEnergy(dT[dBFT[m] + mi], dBe + dFS[m], dB[m].stride,
                                                    dSOC[dBFC[m] + c], dB[m].nenergy);
                    }
                }
                act[i] = edim;
                bef[i] = edb4;
            }
            // calcInduction
            bool   haveInduction = false;
            double eInduction    = 0, eAllElec = 0, eExchInd = 0;
            for (uint32_t i = 0; i < nt; i++)
            {
                uint32_t role = dT[t0 + i].role;
                if (role & deviationElec)
                {
                    if (bef[i] != act[i])
                    {
                        eInduction   += act[i] - bef[i];
                        act[i]        = bef[i];
                        haveInduction = true;
                    }
                    eAllElec += act[i];
                }
                else if (role & deviationExchangeRole)
                {
                    eExchInd += act[i];
                }
            }
            double induction      = indEntry >= 0 ? act[indEntry] : 0.0;
            bool   haveInductionC = indEntry >= 0;
            if (haveInduction)
            {
                induction      = eInduction;
                haveInductionC = true;
                eAllElec      += eInduction;
                eExchInd      += eInduction;
            }
            else if (indEntry >= 0)
            {
                eExchInd += induction;
            }
            // potentialCorrection and addCorrection, to the induction
            if (corrEntry >= 0)
            {
                double correction = act[corrEntry];
                eAllElec      += correction;
                eExchInd      += correction;
                induction      = (haveInductionC ? induction : 0.0) + correction;
                haveInductionC = true;
            }
            double w = dW[dm.refOffset + c];
            for (uint32_t j = 0; j < ntarget; j++)
            {
                int32_t tt   = dTg[j];
                double  v    = 0;
                bool    have = false;
                if (tt == allElecType)
                {
                    v    = eAllElec;
                    have = true;
                }
                else if (tt == exchIndType)
                {
                    v    = eExchInd;
                    have = true;
                }
                else if (tt == inductionType)
                {
                    v    = induction;
                    have = haveInductionC;
                }
                else
                {
                    for (uint32_t i = 0; i < nt; i++)
                    {
                        if (dT[t0 + i].type == tt)
                        {
                            v    = act[i];
                            have = static_cast<int>(i) != corrEntry;
                        }
                    }
                }
                uint64_t r   = dm.refOffset*ntarget + static_cast<uint64_t>(j)*dm.nconf + c;
                bool     use = have && dH[r] && sycl::isfinite(v);
                double   de  = use ? dQ[r] - v : 0.0;
                double   wi  = use ? w : 0.0;
                dTe[2*(gc*ntarget + j)]     = wi;
                dTe[2*(gc*ntarget + j) + 1] = wi*deviationLoss(lf, de);
            }
        });
    }
    // The weight and loss of the potential energy of each conformation of each topology
    if (nmc > 0)
    {
        queue->parallel_for(sycl::range<1>(nmc), [=](sycl::id<1> id)
        {
#if defined(__clang__)
#pragma clang fp contract(off)
#endif
            uint64_t gc = id[0];
            uint32_t lo = 0, hi = nmono;
            while (hi - lo > 1)
            {
                uint32_t mid = (lo + hi)/2;
                if (dMo[mid].confOffset <= gc)
                {
                    lo = mid;
                }
                else
                {
                    hi = mid;
                }
            }
            const DeviationMonomer  &mo = dMo[lo];
            if (!dU[mo.block])
            {
                return;
            }
            uint64_t                 c  = gc - mo.confOffset;
            const DeviceSystemBlock &bm = dB[mo.block];
            int                      ep = -1;
            for (uint32_t i = dBFT[mo.block]; i < dBFT[mo.block + 1]; i++)
            {
                ep = dT[i].kind == deviationEpot ? static_cast<int>(i) : ep;
            }
            double   wq = 0, lq = 0;
            uint64_t q  = mo.refOffset + c;
            if (ep >= 0)
            {
                double e = deviationTypeEnergy(dT[ep], dA + dFS[mo.block], bm.stride,
                                               dSOC[dBFC[mo.block] + c], bm.nenergy);
                if (dEH[q] && sycl::isfinite(e))
                {
                    wq = dEW[q];
                    lq = dEW[q] * deviationLoss(lf, dEQ[q] - e);
                }
            }
            dETe[2*gc]     = wq;
            dETe[2*gc + 1] = lq;
        });
    }
    // The sum over the conformations of each dimer in order, as addLosses
    if (ndimer > 0 && ntarget > 0)
    {
        queue->parallel_for(sycl::range<1>(ndimer*ntarget), [=](sycl::id<1> id)
        {
#if defined(__clang__)
#pragma clang fp contract(off)
#endif
            uint32_t              k  = id[0]/ntarget;
            uint32_t              j  = id[0]%ntarget;
            const DeviationDimer &dm = dD[k];
            if (!dU[dm.dimer])
            {
                return;
            }
            double                sw = 0, sl = 0;
            for (uint64_t c = 0; c < dm.nconf; c++)
            {
                sw += dTe[2*((dm.confOffset + c)*ntarget + j)];
                sl += dTe[2*((dm.confOffset + c)*ntarget + j) + 1];
            }
            dDS[2*(k*ntarget + j)]     = sw;
            dDS[2*(k*ntarget + j) + 1] = sl;
        });
    }
    // The sums of the dimers and the terms of the potential energy to the host
    size_t  nd = 2*static_cast<size_t>(ndimer)*ntarget;
    size_t  ne = 2*nmc;
    double *hS = sys->devSumsHost_.reserve(*queue, std::max<size_t>(nd + ne, 1));
    if (nd > 0)
    {
        queue->memcpy(hS, dDS, nd*sizeof(double));
    }
    if (ne > 0)
    {
        queue->memcpy(hS + nd, dETe, ne*sizeof(double));
    }
    queue->wait_and_throw();
    std::copy(hS, hS + nd, dimerSums);
    std::copy(hS + nd, hS + nd + ne, epotTerms);
}

uint64_t blockAtomOffset(const System *sys, size_t block)
{
    return sys->blocksHost_[block].atomOffset;
}

void createChain(System *sys, const ChainSystem &c)
{
    for (auto p : { &sys->chainDeps_ })
    {
        sys->release(*p);
    }
    sys->release(sys->chainMols_);
    sys->release(sys->chainBounds_);
    sys->release(sys->chainExtras_);
    sys->release(sys->chainTerms_);
    sys->release(sys->chainEntry_);
    sys->release(sys->chainCache_);
    sys->release(sys->chainValid_);
    sys->release(sys->chainParams_);
    sys->release(sys->chainHit_);
    sys->release(sys->chainActive_);
    sys->release(sys->chainBoundNow_);
    sys->release(sys->chainMasks_);
    sys->release(sys->chainBlockMol_);
    sys->release(sys->chainBlockWhat_);
    sys->release(sys->chainChi2_);
    uint32_t ntarget = sys->dev_.ntarget;
    sys->nchainMol_   = c.nmol;
    sys->nchainTerm_  = c.nterm;
    sys->nchainParam_ = c.nparam;
    sys->chainLoss_   = c.lossFunction;
    // The sets kept of each molecule, see chainEntrySize
    sys->chainEntryHost_.assign(c.nmol, 0);
    sys->chainSizeHost_.assign(c.nmol, 0);
    uint64_t total = 0;
    for (size_t m = 0; m < c.nmol; m++)
    {
        const auto &cm = c.mols[m];
        uint64_t size  = cm.ndep + (cm.dimer >= 0 ? 2*ntarget : 0) + (cm.monomer >= 0 ? 2*cm.nepot : 0) + 1;
        sys->chainEntryHost_[m] = total;
        sys->chainSizeHost_[m]  = size;
        total += 2*size;
    }
    sys->chainCacheSize_ = total;
    sys->chainMols_      = sys->allocate<ChainMolecule>(c.nmol);
    sys->chainDeps_      = sys->allocate<uint32_t>(c.ndep);
    sys->chainBounds_    = sys->allocate<ChainBound>(c.nbound);
    sys->chainExtras_    = sys->allocate<uint64_t>(c.nextra);
    sys->chainTerms_     = sys->allocate<ChainTerm>(c.nterm);
    sys->chainEntry_     = sys->allocate<uint64_t>(c.nmol);
    sys->chainCache_     = sys->allocate<double>(total);
    sys->chainValid_     = sys->allocate<char>(2*c.nmol);
    sys->chainParams_    = sys->allocate<double>(c.nparam);
    sys->chainHit_       = sys->allocate<int32_t>(c.nmol);
    sys->chainActive_    = sys->allocate<char>(c.nmol);
    sys->chainBoundNow_  = sys->allocate<double>(c.nmol);
    sys->chainMasks_     = sys->allocate<char>(3*sys->nblock_);
    sys->chainBlockMol_  = sys->allocate<int32_t>(sys->nblock_);
    sys->chainBlockWhat_ = sys->allocate<int32_t>(sys->nblock_);
    sys->chainChi2_      = sys->allocate<double>(1);
    sys->upload(sys->chainMols_, c.mols, c.nmol);
    sys->upload(sys->chainDeps_, c.deps, c.ndep);
    sys->upload(sys->chainBounds_, c.bounds, c.nbound);
    sys->upload(sys->chainExtras_, c.extras, c.nextra);
    sys->upload(sys->chainTerms_, c.terms, c.nterm);
    sys->upload(sys->chainEntry_, sys->chainEntryHost_.data(), c.nmol);
    sys->upload(sys->chainBlockMol_, c.blockMol, sys->nblock_);
    sys->upload(sys->chainBlockWhat_, c.blockWhat, sys->nblock_);
    sys->queue_->memset(sys->chainValid_, 0, std::max<size_t>(2*c.nmol, 1));
    sys->queue_->memset(sys->chainCache_, 0, std::max<uint64_t>(total, 1)*sizeof(double));
    sys->queue_->memset(sys->chainActive_, 1, std::max<size_t>(c.nmol, 1));
    sys->queue_->memset(sys->chainHit_, 0xff, std::max<size_t>(c.nmol, 1)*sizeof(int32_t));
    sys->queue_->wait_and_throw();
}

size_t chainEntrySize(const System *sys, size_t mol)
{
    return sys->chainSizeHost_[mol];
}

void chainSetCache(System *sys, const double *values, const char *valid)
{
    if (sys->chainCacheSize_ > 0)
    {
        sys->queue_->memcpy(sys->chainCache_, values, sys->chainCacheSize_*sizeof(double));
    }
    if (sys->nchainMol_ > 0)
    {
        sys->queue_->memcpy(sys->chainValid_, valid, 2*sys->nchainMol_);
    }
    sys->queue_->wait_and_throw();
}

void chainGetCache(System *sys, double *values, char *valid)
{
    if (sys->chainCacheSize_ > 0)
    {
        sys->queue_->memcpy(values, sys->chainCache_, sys->chainCacheSize_*sizeof(double));
    }
    if (sys->nchainMol_ > 0)
    {
        sys->queue_->memcpy(valid, sys->chainValid_, 2*sys->nchainMol_);
    }
    sys->queue_->wait_and_throw();
}

void chainSetParameters(System *sys, const double *params)
{
    if (sys->nchainParam_ > 0)
    {
        sys->queue_->memcpy(sys->chainParams_, params, sys->nchainParam_*sizeof(double));
        sys->queue_->wait_and_throw();
    }
}

void chainSkip(System *sys, bool canSkip)
{
    auto    *queue = sys->queue_;
    auto     dM    = sys->chainMols_;
    auto     dDep  = sys->chainDeps_;
    auto     dEn   = sys->chainEntry_;
    auto     dC    = sys->chainCache_;
    auto     dV    = sys->chainValid_;
    auto     dP    = sys->chainParams_;
    auto     dHit  = sys->chainHit_;
    auto     dAct  = sys->chainActive_;
    size_t   nmol  = sys->nchainMol_;
    uint32_t ntarget = sys->dev_.ntarget;
    if (nmol > 0)
    {
        queue->parallel_for(sycl::range<1>(nmol), [=](sycl::id<1> id)
        {
            size_t               m  = id[0];
            const ChainMolecule &cm = dM[m];
            uint64_t size = cm.ndep + (cm.dimer >= 0 ? 2*ntarget : 0) + (cm.monomer >= 0 ? 2*cm.nepot : 0) + 1;
            int      hit  = -1;
            // The last set that matches, as prepareDeviation
            for (int e = 0; e < 2; e++)
            {
                const double *key = dC + dEn[m] + e*size;
                bool same = dV[2*m + e] != 0;
                for (uint32_t i = 0; i < cm.ndep && same; i++)
                {
                    same = key[i] == dP[dDep[cm.depOffset + i]];
                }
                hit = same ? e : hit;
            }
            dHit[m] = hit;
            dAct[m] = (canSkip && hit >= 0 && cm.skippable) ? 0 : 1;
        });
    }
    auto     dMask = sys->chainMasks_;
    auto     dBM   = sys->chainBlockMol_;
    auto     dBW   = sys->chainBlockWhat_;
    uint32_t nb    = sys->nblock_;
    if (nb > 0)
    {
        queue->parallel_for(sycl::range<1>(nb), [=](sycl::id<1> id)
        {
            uint32_t b   = id[0];
            int32_t  mol = dBM[b];
            char     act = (mol >= 0 && dAct[mol]) ? 1 : 0;
            dMask[b]        = act;
            dMask[nb + b]   = (act && dBW[b] < 2) ? 1 : 0;
            dMask[2*nb + b] = (act && dBW[b] == 2) ? 1 : 0;
        });
    }
}

const char *chainMask(System *sys, int what)
{
    int k = (what == -1) ? 0 : (what == -2 ? 1 : 2);
    return sys->chainMasks_ + k*sys->nblock_;
}

void chainActive(System *sys, char *active, int *hit)
{
    size_t nmol = sys->nchainMol_;
    if (nmol > 0)
    {
        std::vector<int32_t> h(nmol);
        sys->queue_->memcpy(active, sys->chainActive_, nmol);
        sys->queue_->memcpy(h.data(), sys->chainHit_, nmol*sizeof(int32_t));
        sys->queue_->wait_and_throw();
        std::copy(h.begin(), h.end(), hit);
    }
}

void chainComputeChi2(System *sys, bool keep)
{
    auto    *queue   = sys->queue_;
    auto     dM      = sys->chainMols_;
    auto     dDep    = sys->chainDeps_;
    auto     dBo     = sys->chainBounds_;
    auto     dEx     = sys->chainExtras_;
    auto     dT      = sys->chainTerms_;
    auto     dEn     = sys->chainEntry_;
    auto     dC      = sys->chainCache_;
    auto     dV      = sys->chainValid_;
    auto     dP      = sys->chainParams_;
    auto     dHit    = sys->chainHit_;
    auto     dAct    = sys->chainActive_;
    auto     dBN     = sys->chainBoundNow_;
    auto     dQa     = sys->acmCharges_;
    auto     dQ      = sys->charges_;
    auto     dDS     = sys->devDimerSums_;
    auto     dETe    = sys->devEpotTerms_;
    auto     dChi2   = sys->chainChi2_;
    size_t   nmol    = sys->nchainMol_;
    uint32_t nterm   = sys->nchainTerm_;
    uint32_t ntarget = sys->dev_.ntarget;
    int      lf      = sys->chainLoss_;
    if (nmol == 0)
    {
        queue->memset(dChi2, 0, sizeof(double));
        return;
    }
    // The bounds of the charges of the molecules computed, as BoundsDevComputer::moleculeBound
    queue->parallel_for(sycl::range<1>(nmol), [=](sycl::id<1> id)
    {
#if defined(__clang__)
#pragma clang fp contract(off)
#endif
        size_t               m     = id[0];
        const ChainMolecule &cm    = dM[m];
        double               bound = 0;
        if (dAct[m] && dQa)
        {
            for (uint32_t k = 0; k < cm.nbound; k++)
            {
                const ChainBound &cb = dBo[cm.boundOffset + k];
                double q = dQa[cb.charge] + (cb.secondCharge >= 0 ? dQ[cb.secondCharge] : 0.0);
                for (uint32_t e = 0; e < cb.nextra; e++)
                {
                    q += dQa[dEx[cb.extraOffset + e]];
                }
                if (q < cb.minimum)
                {
                    bound += deviationLoss(lf, cb.minimum - q);
                }
                else if (q > cb.maximum)
                {
                    bound += deviationLoss(lf, q - cb.maximum);
                }
            }
        }
        dBN[m] = bound;
    });
    // The sums of the targets in the order of the molecules and the total,
    // as finishDeviation and sumChiSquared
    queue->single_task([=]()
    {
#if defined(__clang__)
#pragma clang fp contract(off)
#endif
        double total = 0;
        for (uint32_t t = 0; t < nterm; t++)
        {
            const ChainTerm &term = dT[t];
            double sw = 0, sl = 0;
            for (size_t m = 0; m < nmol; m++)
            {
                const ChainMolecule &cm   = dM[m];
                uint64_t size  = cm.ndep + (cm.dimer >= 0 ? 2*ntarget : 0) + (cm.monomer >= 0 ? 2*cm.nepot : 0) + 1;
                const double *kept = dC + dEn[m] + (dHit[m] > 0 ? size : 0);
                bool fresh = dAct[m] != 0;
                if (term.kind == chainTermDimer && cm.dimer >= 0)
                {
                    const double *v = fresh ? dDS + 2*(static_cast<uint64_t>(cm.dimer)*ntarget + term.index)
                                            : kept + cm.ndep + 2*term.index;
                    sw += v[0];
                    sl += v[1];
                }
                else if (term.kind == chainTermEpot && cm.monomer >= 0)
                {
                    const double *v = fresh ? dETe + 2*cm.epotOffset
                                            : kept + cm.ndep + (cm.dimer >= 0 ? 2*ntarget : 0);
                    for (uint32_t c = 0; c < cm.nepot; c++)
                    {
                        sw += v[2*c];
                        sl += v[2*c + 1];
                    }
                }
                else if (term.kind == chainTermBounds)
                {
                    sw += 1;
                    sl += fresh ? dBN[m] : kept[size - 1];
                }
            }
            // As FittingTarget::chiSquaredWeighted
            if (sl > 0 && sw > 0)
            {
                total += sl*term.weight/sw;
            }
        }
        *dChi2 = total;
    });
    if (!keep)
    {
        return;
    }
    // The sets kept, as finishDeviation
    queue->parallel_for(sycl::range<1>(nmol), [=](sycl::id<1> id)
    {
        size_t               m    = id[0];
        const ChainMolecule &cm   = dM[m];
        uint64_t size  = cm.ndep + (cm.dimer >= 0 ? 2*ntarget : 0) + (cm.monomer >= 0 ? 2*cm.nepot : 0) + 1;
        double  *e0    = dC + dEn[m];
        double  *e1    = e0 + size;
        if (!dAct[m])
        {
            if (dHit[m] == 1)
            {
                for (uint64_t i = 0; i < size; i++)
                {
                    double t = e0[i];
                    e0[i]    = e1[i];
                    e1[i]    = t;
                }
                char v          = dV[2*m];
                dV[2*m]         = dV[2*m + 1];
                dV[2*m + 1]     = v;
            }
            return;
        }
        for (uint64_t i = 0; i < size; i++)
        {
            e1[i] = e0[i];
        }
        dV[2*m + 1] = dV[2*m];
        uint64_t k = 0;
        for (uint32_t i = 0; i < cm.ndep; i++)
        {
            e0[k++] = dP[dDep[cm.depOffset + i]];
        }
        if (cm.dimer >= 0)
        {
            for (uint32_t i = 0; i < 2*ntarget; i++)
            {
                e0[k++] = dDS[2*static_cast<uint64_t>(cm.dimer)*ntarget + i];
            }
        }
        if (cm.monomer >= 0)
        {
            for (uint32_t i = 0; i < 2*cm.nepot; i++)
            {
                e0[k++] = dETe[2*cm.epotOffset + i];
            }
        }
        e0[k]   = dBN[m];
        dV[2*m] = 1;
    });
}

double chainChi2(System *sys)
{
    double chi2 = 0;
    sys->queue_->memcpy(&chi2, sys->chainChi2_, sizeof(double));
    sys->queue_->wait_and_throw();
    return chi2;
}
} // namespace gpu

} // namespace alexandria

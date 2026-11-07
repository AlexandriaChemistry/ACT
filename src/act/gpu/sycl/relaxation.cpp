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
 * Relaxing the shells of the conformations of many molecules on the GPU
 * with SYCL, see relaxShells in gpu.h. Each iteration moves the shells of
 * the slots that did not converge as ForceComputer::updateShells moves
 * them, in three kernels with many work-items each, since a GPU needs many
 * to hide the time to read memory: the squared force and the step product
 * of each shell in each slot, the sums over the shells of each fragment
 * in their order, with the convergence and the step, and the move of each
 * shell. A fourth kernel computes the forces of each shell in each slot,
 * since the relaxation only needs those: the forces of the interactions
 * without shells, which do not change, and those of its interactions with
 * shells, which it computes, in their order, such that no forces of pairs
 * are stored and read again, which would take most of the bandwidth of the
 * memory. The forces on the other particles are computed once at the end,
 * if they are needed. For the rounding floor of the shell forces, see
 * updateShells, the largest force on these is taken from the start. Since the shells of a
 * slot do not move after it converged, the energies of its last
 * calculation are those of the final coordinates, which are computed once
 * at the end, in the order of the interactions. No atomic operations are
 * used, such that the results are the same in every run and do not depend
 * on the other molecules or slots. The flags of the slots that converged
 * are kept for each iteration and copied to the host every few iterations
 * only, since waiting for the GPU takes longer than an iteration of small
 * molecules.
 * \author David van der Spoel <david.vanderspoel@icm.uu.se>
 */
#include "gpu.h"

#include <algorithm>
#include <cfloat>
#include <chrono>
#include <cstdio>
#include <cstdlib>
#include <string>
#include <cstdint>
#include <limits>
#include <vector>

#include <sycl/sycl.hpp>

#include "sycl/context.h"
#include "sycl/pairkernels.h"

namespace alexandria
{

namespace gpu
{

namespace
{

//! A molecule in the device arrays
struct DeviceMolecule
{
    //! Start of the coordinates, forces and displacements, rows of stride
    uint64_t xOffset;
    //! Start of the slots
    uint64_t slotOffset;
    //! Start of the atoms times slots
    uint64_t atomWorkOffset;
    //! Start of the atoms in atomStart
    uint64_t atomOffset;
    //! Start of the force and moving shells times slots, see shellContributions
    uint64_t shellWorkOffset;
    //! Start of the force shells times slots, for the squared forces
    uint64_t forceShellPartOffset;
    //! Start of the moving shells times slots, for the moves and their sums
    uint64_t movingShellPartOffset;
    //! The first of the force shells
    uint32_t firstForceShell;
    //! Number of these
    uint32_t nforceShell;
    //! The first of the moving shells
    uint32_t firstMovingShell;
    //! Number of these
    uint32_t nmovingShell;
    //! Start of the shell rows times slots in the compact results
    uint64_t resultOffset;
    //! Start of the shells times slots, for their forces
    uint64_t shellRowWorkOffset;
    //! The first of the shells, in shellRows
    uint32_t firstShellRow;
    //! Number of shells
    uint32_t nshellRow;
    //! The first fragment
    uint32_t firstFragment;
    //! Number of fragments
    uint32_t nfragment;
    //! Number of particles
    uint32_t natom;
    //! Distance between the rows
    uint32_t stride;
    //! Number of slots
    uint32_t nslot;
};

//! A potential of a molecule in the device arrays
struct DeviceBlock
{
    //! Start of the interactions
    uint64_t pairOffset;
    //! Start of the interactions times slots
    uint64_t workOffset;
    //! Start of the slots of all blocks, for the energies
    uint64_t blockSlotOffset;
    //! Start of the energies, pairEnergyTerms rows of stride
    uint64_t energyOffset;
    //! Number of interactions
    uint64_t npair;
    //! The molecule
    uint32_t mol;
    //! The potential
    PairKernel kernel;
};

/*! \brief The element whose range holds index, by bisection
 * \param[in] a     The elements, sorted by start
 * \param[in] n     Number of elements
 * \param[in] index The index
 * \param[in] start Member with the start of the range of each element
 * \return the element
 */
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

//! The coordinate in nm, rounded to float as gmx::RVec stores it
inline float toNm(int64_t c)
{
    return static_cast<float>(static_cast<double>(c)*coordinateUnitNm);
}

//! The difference of two coordinates with wrap-around, see coordinateDifference
inline int64_t difference(int64_t a, int64_t b)
{
    return static_cast<int64_t>(static_cast<uint64_t>(a) - static_cast<uint64_t>(b));
}

//! A coordinate in nm in integer units, see toCoordinate, without the range check
inline int64_t toUnits(double x)
{
    return static_cast<int64_t>(sycl::rint(x*coordinateScaleNm));
}

/*! \brief Move a coordinate, see moveCoordinate in coordinate.h
 * \param[in]    c     The coordinate
 * \param[in]    d     The displacement in nm
 * \param[inout] valid Set to false if d is not finite, then c is returned
 * \return the moved coordinate
 */
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
    if (step > 0 && c > std::numeric_limits<int64_t>::max() - step)
    {
        return big;
    }
    if (step < 0 && c < std::numeric_limits<int64_t>::min() - step)
    {
        return -big;
    }
    return c + step;
}

/*! \brief The displacement factor of a shell with a quartic spring along the force, see calc_fscale
 * \param[in] k1 The harmonic force constant
 * \param[in] k2 The quartic force constant
 * \param[in] F0 The norm of the force
 * \return the factor
 */
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

//! The state of the relaxation of the shells of a fragment in a slot, see ShellRelaxation
struct FragmentState
{
    //! Fraction of the Newton step
    double step;
    //! Mean square force of the previous step
    double msPrev;
    //! Displacement times force of the previous step
    double dxfPrev;
    //! Whether the shells move in this iteration, then dxfPrev is their sum in the next
    int    moving;
    /*! The largest force component on the particles of the fragment that
     * are not shells, at the start, for the rounding floor of the shell
     * forces, see updateShells. During the relaxation only the forces on
     * the shells are computed, those on the other particles change little.
     */
    double fmaxRest;
    //! The mean square force against which progress is measured, see shellStallIterations
    double reference;
    //! Steps since the mean square force was halved
    int    stall;
};

//! The device memory for relaxShells
struct RelaxMemory
{
    DeviceArray<char>          in;
    HostArray<char>            inHost;
    DeviceArray<float>         dx;
    DeviceArray<FragmentState> state;
    DeviceArray<float>         epair;
    //! The squared force of each force shell in each slot
    DeviceArray<double>        partF2;
    //! The displacement times force of each moving shell in each slot, before and after it moves
    DeviceArray<double>        partAfter;
    DeviceArray<double>        partDxF;
    DeviceArray<double>        energies;
    HostArray<double>          energiesHost;
    DeviceArray<char>          doneHistory;
    HostArray<char>            doneHost;
    DeviceArray<int>           valid;
    //! The coordinates and forces of the shells, compact
    DeviceArray<char>          result;
    HostArray<char>            resultHost;
};

//! The device memory, never deleted, see gpuQueue
RelaxMemory &relaxMemory()
{
    static RelaxMemory *memory = new RelaxMemory;
    return *memory;
}

/*! \brief Times of the parts of relaxShells, for developers
 * With the environment variable ACT_GPU_PROFILE the host waits after each
 * part, which makes the calculation slower, and the times are printed to
 * standard error every 100 calls.
 */
class RelaxProfile
{
public:
    //! The parts
    enum Part { Upload, Update, Forces, Check, Energies, Download, NPART };

    //! \return whether the times are measured
    static bool on()
    {
        static const bool profile = std::getenv("ACT_GPU_PROFILE") != nullptr;
        return profile;
    }

    //! Start measuring a call
    void start(sycl::queue *queue)
    {
        queue->wait();
        t_ = std::chrono::steady_clock::now();
    }

    //! The time up to now belongs to a part
    void stop(sycl::queue *queue, Part part)
    {
        queue->wait();
        auto now   = std::chrono::steady_clock::now();
        time_[part] += std::chrono::duration<double>(now - t_).count();
        count_[part] += 1;
        t_ = now;
    }

    //! A call is done, with the work of its iterations
    void done(size_t iterations, size_t work, size_t slots)
    {
        ncall_ += 1;
        niter_ += iterations;
        work_  += work;
        slots_ += slots;
        if (ncall_ % 100 == 0)
        {
            static const char *names[NPART] = { "upload", "update", "forces", "check", "energies", "download" };
            std::string msg = "ACT_GPU_PROFILE relaxShells calls " + std::to_string(ncall_) +
                " iterations/call " + std::to_string(double(niter_)/ncall_) +
                " pairs*slots/call " + std::to_string(double(work_)/ncall_) +
                " slots/call " + std::to_string(double(slots_)/ncall_) + ", ms per call:";
            for (int p = 0; p < NPART; p++)
            {
                msg += std::string(" ") + names[p] + " " + std::to_string(1e3*time_[p]/ncall_);
            }
            std::fprintf(stderr, "%s\n", msg.c_str());
        }
    }

private:
    std::chrono::steady_clock::time_point t_;
    double time_[NPART]  = { 0 };
    long   count_[NPART] = { 0 };
    long   ncall_ = 0, niter_ = 0;
    double work_ = 0, slots_ = 0;
};

} // namespace

bool relaxShells(const RelaxMolecule *mols,
                 size_t               nmol,
                 const RelaxSettings &settings,
                 size_t              *nComputed,
                 bool                *notFinite)
{
    *nComputed = 0;
    *notFinite = false;
    sycl::queue *queue = gpuQueue();
    if (!queue)
    {
        return false;
    }
    auto &mem = relaxMemory();
    // Where everything is in the device arrays
    std::vector<DeviceMolecule> dmols(nmol);
    std::vector<DeviceBlock>    dblocks;
    uint64_t nx = 0, nslotAll = 0, natomWork = 0, natomAll = 0, nfrag = 0, nforceShell = 0, nmovingShell = 0;
    uint64_t npairAll = 0, nwork = 0, nblockSlot = 0, nenergy = 0, nentry = 0;
    uint64_t nshellWork = 0, nforcePart = 0, nmovingPart = 0, nshellRowAll = 0, nresult = 0;
    bool     allForces  = false;
    // The particles that are shells, of which the coordinates and forces
    // are copied back
    std::vector<uint32_t> shellRows;
    for (size_t k = 0; k < nmol; k++)
    {
        const auto &m = mols[k];
        auto       &d = dmols[k];
        d.xOffset        = nx;
        d.slotOffset     = nslotAll;
        d.atomWorkOffset = natomWork;
        d.atomOffset     = natomAll;
        d.firstFragment  = nfrag;
        d.nfragment      = m.nfragment;
        d.natom          = m.natom;
        d.stride         = m.stride;
        d.nslot          = m.nslot;
        nx              += 3*m.natom*m.stride;
        nslotAll        += m.nslot;
        natomWork       += m.natom*m.nslot;
        natomAll        += m.natom;
        nfrag           += m.nfragment;
        d.firstForceShell  = nforceShell;
        d.firstMovingShell = nmovingShell;
        for (size_t f = 0; f < m.nfragment; f++)
        {
            nforceShell  += m.fragments[f].nforceShell;
            nmovingShell += m.fragments[f].nmovingShell;
        }
        {
            std::vector<uint32_t> rows;
            for (size_t f = 0; f < m.nfragment; f++)
            {
                const auto &fr = m.fragments[f];
                for (size_t k = 0; k < fr.nforceShell; k++)
                {
                    rows.push_back(m.forceShells[fr.firstForceShell + k].shell);
                }
                for (size_t k = 0; k < fr.nmovingShell; k++)
                {
                    rows.push_back(m.movingShells[fr.firstMovingShell + k].shell);
                }
            }
            std::sort(rows.begin(), rows.end());
            rows.erase(std::unique(rows.begin(), rows.end()), rows.end());
            d.firstShellRow = nshellRowAll;
            d.nshellRow     = rows.size();
            d.resultOffset       = nresult;
            d.shellRowWorkOffset = nresult/3;
            shellRows.insert(shellRows.end(), rows.begin(), rows.end());
            nshellRowAll   += rows.size();
            nresult        += 3*rows.size()*m.nslot;
            allForces       = allForces || m.allForces;
        }
        d.nforceShell           = nforceShell - d.firstForceShell;
        d.nmovingShell          = nmovingShell - d.firstMovingShell;
        d.shellWorkOffset       = nshellWork;
        d.forceShellPartOffset  = nforcePart;
        d.movingShellPartOffset = nmovingPart;
        nshellWork  += (d.nforceShell + d.nmovingShell)*m.nslot;
        nforcePart  += d.nforceShell*m.nslot;
        nmovingPart += d.nmovingShell*m.nslot;
        for (size_t b = 0; b < m.nblock; b++)
        {
            DeviceBlock db;
            db.pairOffset      = npairAll;
            db.workOffset      = nwork;
            db.blockSlotOffset = nblockSlot;
            db.energyOffset    = nenergy;
            db.npair           = m.blocks[b].npair;
            db.mol             = k;
            db.kernel          = m.blocks[b].kernel;
            dblocks.push_back(db);
            npairAll   += m.blocks[b].npair;
            nwork      += m.blocks[b].npair*m.nslot;
            nblockSlot += m.nslot;
            nenergy    += pairEnergyTerms*m.stride;
            nentry     += 2*m.blocks[b].npair;
        }
    }
    if (nslotAll == 0)
    {
        return true;
    }
    uint32_t nm = nmol;
    uint32_t nb = dblocks.size();
    // All input in one buffer, also the arrays that change
    PackedArrays in;
    size_t oMols         = in.add<DeviceMolecule>(nmol);
    size_t oBlocks       = in.add<DeviceBlock>(dblocks.size());
    size_t oFragments    = in.add<ShellFragment>(nfrag);
    size_t oForceShells  = in.add<ForceShell>(nforceShell);
    size_t oMovingShells = in.add<MovingShell>(nmovingShell);
    size_t oPairs        = in.add<PairEntry>(npairAll);
    size_t oAtomStart    = in.add<uint64_t>(natomAll + nmol);
    size_t oEntries      = in.add<uint64_t>(nentry);
    size_t oPairKernel   = in.add<uint32_t>(npairAll);
    size_t oMovingFrag   = in.add<uint32_t>(nmovingShell);
    size_t oShellRows    = in.add<uint32_t>(nshellRowAll);
    size_t oIsShell      = in.add<char>(natomAll);
    size_t oFixed        = in.add<float>(nx);
    size_t oX            = in.add<int64_t>(nx);
    size_t oF            = in.add<float>(nx);
    size_t oDone         = in.add<char>(nslotAll);
    size_t oComputed     = in.add<char>(nslotAll);
    char  *h             = mem.inHost.reserve(*queue, in.size());
    std::copy(dmols.begin(), dmols.end(), reinterpret_cast<DeviceMolecule *>(h + oMols));
    std::copy(dblocks.begin(), dblocks.end(), reinterpret_cast<DeviceBlock *>(h + oBlocks));
    auto *hFrag   = reinterpret_cast<ShellFragment *>(h + oFragments);
    auto *hFS     = reinterpret_cast<ForceShell *>(h + oForceShells);
    auto *hMS     = reinterpret_cast<MovingShell *>(h + oMovingShells);
    auto *hPairs  = reinterpret_cast<PairEntry *>(h + oPairs);
    auto *hStart  = reinterpret_cast<uint64_t *>(h + oAtomStart);
    auto *hEntry  = reinterpret_cast<uint64_t *>(h + oEntries);
    auto *hPK     = reinterpret_cast<uint32_t *>(h + oPairKernel);
    auto *hMF     = reinterpret_cast<uint32_t *>(h + oMovingFrag);
    std::copy(shellRows.begin(), shellRows.end(), reinterpret_cast<uint32_t *>(h + oShellRows));
    {
        char *isShell = h + oIsShell;
        std::fill(isShell, isShell + natomAll, 0);
        for (size_t k = 0; k < nmol; k++)
        {
            const auto &d = dmols[k];
            for (uint32_t r = 0; r < d.nshellRow; r++)
            {
                isShell[d.atomOffset + shellRows[d.firstShellRow + r]] = 1;
            }
        }
    }
    uint64_t ifs = 0, ims = 0, ip = 0, ie = 0, ib = 0;
    for (size_t k = 0; k < nmol; k++)
    {
        const auto &m = mols[k];
        const auto &d = dmols[k];
        for (size_t f = 0; f < m.nfragment; f++)
        {
            auto fr = m.fragments[f];
            std::copy(m.forceShells + fr.firstForceShell, m.forceShells + fr.firstForceShell + fr.nforceShell, hFS + ifs);
            std::copy(m.movingShells + fr.firstMovingShell, m.movingShells + fr.firstMovingShell + fr.nmovingShell, hMS + ims);
            fr.firstForceShell  = ifs;
            fr.firstMovingShell = ims;
            std::fill(hMF + ims, hMF + ims + fr.nmovingShell, static_cast<uint32_t>(f));
            ifs += fr.nforceShell;
            ims += fr.nmovingShell;
            hFrag[d.firstFragment + f] = fr;
        }
        // The interactions of each atom, of all blocks in their order,
        // as the interaction times two, plus one for the second atom
        uint64_t *start = hStart + d.atomOffset + k;
        std::fill(start + 1, start + m.natom + 1, 0);
        start[0] = ie;
        for (size_t b = 0; b < m.nblock; b++)
        {
            for (size_t p = 0; p < m.blocks[b].npair; p++)
            {
                start[m.blocks[b].pairs[p].i + 1] += 1;
                start[m.blocks[b].pairs[p].j + 1] += 1;
            }
        }
        for (size_t a = 1; a <= m.natom; a++)
        {
            start[a] += start[a-1];
        }
        std::vector<uint64_t> next(start, start + m.natom);
        for (size_t b = 0; b < m.nblock; b++, ib++)
        {
            const auto &blk = m.blocks[b];
            std::copy(blk.pairs, blk.pairs + blk.npair, hPairs + ip);
            std::fill(hPK + ip, hPK + ip + blk.npair, static_cast<uint32_t>(blk.kernel));
            for (size_t p = 0; p < blk.npair; p++)
            {
                uint64_t gp = ip + p;
                hEntry[next[blk.pairs[p].i]++] = 2*gp;
                hEntry[next[blk.pairs[p].j]++] = 2*gp + 1;
            }
            ip += blk.npair;
        }
        ie = start[m.natom];
        size_t n = 3*m.natom*m.stride;
        std::copy(m.fixedForces, m.fixedForces + n, reinterpret_cast<float *>(h + oFixed) + d.xOffset);
        std::copy(m.x, m.x + n, reinterpret_cast<int64_t *>(h + oX) + d.xOffset);
        std::copy(m.f, m.f + n, reinterpret_cast<float *>(h + oF) + d.xOffset);
        std::copy(m.done, m.done + m.nslot, h + oDone + d.slotOffset);
    }
    std::fill(h + oComputed, h + oComputed + nslotAll, 0);

    char *din   = mem.in.reserve(*queue, in.size());
    auto  dMols = reinterpret_cast<const DeviceMolecule *>(din + oMols);
    auto  dBlk  = reinterpret_cast<const DeviceBlock *>(din + oBlocks);
    auto  dFrag = reinterpret_cast<const ShellFragment *>(din + oFragments);
    auto  dFS   = reinterpret_cast<const ForceShell *>(din + oForceShells);
    auto  dMS   = reinterpret_cast<const MovingShell *>(din + oMovingShells);
    auto  dPair = reinterpret_cast<const PairEntry *>(din + oPairs);
    auto  dStrt = reinterpret_cast<const uint64_t *>(din + oAtomStart);
    auto  dEnt  = reinterpret_cast<const uint64_t *>(din + oEntries);
    auto  dPK   = reinterpret_cast<const uint32_t *>(din + oPairKernel);
    auto  dMF   = reinterpret_cast<const uint32_t *>(din + oMovingFrag);
    auto  dSR   = reinterpret_cast<const uint32_t *>(din + oShellRows);
    auto  dIS   = din + oIsShell;
    auto  dFix  = reinterpret_cast<const float *>(din + oFixed);
    auto  dX    = reinterpret_cast<int64_t *>(din + oX);
    auto  dF    = reinterpret_cast<float *>(din + oF);
    auto  dDone = din + oDone;
    auto  dComp = din + oComputed;
    auto  dDx   = mem.dx.reserve(*queue, nx);
    auto  dSt   = mem.state.reserve(*queue, 2*nslotAll);
    auto  dPF2  = mem.partF2.reserve(*queue, std::max<uint64_t>(nforcePart, 1));
    auto  dPA   = mem.partAfter.reserve(*queue, std::max<uint64_t>(nmovingPart, 1));
    auto  dPD   = mem.partDxF.reserve(*queue, std::max<uint64_t>(nmovingPart, 1));
    auto  dEp   = mem.epair.reserve(*queue, pairEnergyTerms*std::max<uint64_t>(nwork, 1));
    auto  dE    = mem.energies.reserve(*queue, std::max<uint64_t>(nenergy, 1));
    auto  dVal  = mem.valid.reserve(*queue, 1);
    // The done flags after each iteration
    size_t maxIter = std::max(settings.maxIter, 1);
    auto   dHist   = mem.doneHistory.reserve(*queue, maxIter*nslotAll);
    char  *hDone   = mem.doneHost.reserve(*queue, maxIter*nslotAll);
    static RelaxProfile profile;
    const bool          prof = RelaxProfile::on();
    if (prof)
    {
        profile.start(queue);
    }
    queue->memcpy(din, h, in.size());
    queue->memset(dDx, 0, nx*sizeof(float));
    queue->memset(dVal, 0, sizeof(int));
    queue->parallel_for(sycl::range<1>(2*nslotAll), [=](sycl::id<1> id)
    {
        dSt[id[0]] = { 1.0, DBL_MAX, 0.0, 0, 0.0, DBL_MAX, 0 };
    });
    // The largest force on the particles that are not shells, from the
    // forces at the start
    queue->parallel_for(sycl::range<1>(nslotAll), [=](sycl::id<1> id)
    {
        uint64_t              s = id[0];
        const DeviceMolecule &m = dMols[findRange(dMols, nm, s, &DeviceMolecule::slotOffset)];
        uint64_t              c = s - m.slotOffset;
        const float          *f = dF + m.xOffset + c;
        for (uint32_t fi = 0; fi < m.nfragment; fi++)
        {
            const ShellFragment &fr   = dFrag[m.firstFragment + fi];
            double               fmax = 0;
            for (uint32_t a = fr.firstAtom; a < fr.firstAtom + fr.natom; a++)
            {
                if (!dIS[m.atomOffset + a])
                {
                    for (int mm = 0; mm < 3; mm++)
                    {
                        fmax = sycl::fmax(fmax, static_cast<double>(sycl::fabs(f[(3*a+mm)*m.stride])));
                    }
                }
            }
            dSt[2*s + fi].fmaxRest = fmax;
        }
    });
    const double toler    = settings.msForceToler;
    const int    stallIterations = settings.stallIterations;
    const double maxDist2 = settings.maxShellDistance2;
    const double minStep  = settings.minStep;
    const double maxStep  = settings.maxStep;
    if (prof)
    {
        profile.stop(queue, RelaxProfile::Upload);
    }
    int  iter    = 1;
    bool allDone = false;
    // The energy of each interaction with shells in each slot computed,
    // once at the end, since the shells of a slot do not move after it
    // converged
    auto energyKernel = [&]()
    {
        queue->parallel_for(sycl::range<1>(nwork), [=](sycl::id<1> id)
        {
            uint64_t              w = id[0];
            const DeviceBlock    &b = dBlk[findRange(dBlk, nb, w, &DeviceBlock::workOffset)];
            const DeviceMolecule &m = dMols[b.mol];
            uint64_t local = w - b.workOffset;
            uint64_t p     = local / m.nslot;
            uint64_t c     = local % m.nslot;
            if (!dComp[m.slotOffset + c])
            {
                return;
            }
            const PairEntry &pair = dPair[b.pairOffset + p];
            const int64_t   *x    = dX + m.xOffset + c;
            uint64_t         st   = m.stride;
            float rsq = 0;
            for (int mm = 0; mm < 3; mm++)
            {
                float d = static_cast<float>(x[(3*pair.i+mm)*st] - x[(3*pair.j+mm)*st])*
                    static_cast<float>(coordinateUnitNm);
                rsq += d*d;
            }
            float e[pairEnergyTerms];
            (void) pairInteraction(b.kernel, pair.p, rsq, e);
            for (int k = 0; k < pairEnergyTerms; k++)
            {
                dEp[k*nwork + w] = e[k];
            }
        });
    };
    /*! The forces of the shells of the slots that did not converge, or of
     * all particles of the slots computed: those of the interactions
     * without shells and those of each of its interactions with shells in
     * their order, computed here
     */
    auto forceKernel = [&](bool allAtoms)
    {
        size_t n = allAtoms ? natomWork : nresult/3;
        if (nwork == 0 || n == 0)
        {
            return;
        }
        queue->parallel_for(sycl::range<1>(n), [=](sycl::id<1> id)
        {
            uint64_t w = id[0];
            uint32_t k = allAtoms ? findRange(dMols, nm, w, &DeviceMolecule::atomWorkOffset) :
                findRange(dMols, nm, w, &DeviceMolecule::shellRowWorkOffset);
            const DeviceMolecule &m = dMols[k];
            uint64_t local = w - (allAtoms ? m.atomWorkOffset : m.shellRowWorkOffset);
            uint64_t row   = local / m.nslot;
            uint64_t c     = local % m.nslot;
            if (allAtoms ? !dComp[m.slotOffset + c] : dDone[m.slotOffset + c])
            {
                return;
            }
            uint64_t       a    = allAtoms ? row : dSR[m.firstShellRow + row];
            uint64_t       st   = m.stride;
            uint64_t       base = m.xOffset + 3*a*st + c;
            const int64_t *x    = dX + m.xOffset + c;
            float          fa[3];
            for (int mm = 0; mm < 3; mm++)
            {
                fa[mm] = dFix[base + mm*st];
            }
            const uint64_t *start = dStrt + m.atomOffset + k;
            for (uint64_t e = start[a]; e < start[a+1]; e++)
            {
                uint64_t         gp     = dEnt[e]/2;
                bool             second = (dEnt[e] % 2) == 1;
                const PairEntry &pair   = dPair[gp];
                float d[3], rsq = 0;
                for (int mm = 0; mm < 3; mm++)
                {
                    d[mm] = static_cast<float>(x[(3*pair.i+mm)*st] - x[(3*pair.j+mm)*st])*
                        static_cast<float>(coordinateUnitNm);
                    rsq  += d[mm]*d[mm];
                }
                float eunused[pairEnergyTerms];
                float fscal = pairInteraction(static_cast<PairKernel>(dPK[gp]), pair.p, rsq, eunused);
                for (int mm = 0; mm < 3; mm++)
                {
                    // As the force of the pair on its first particle
                    float fp = rsq > 0 ? fscal*d[mm] : 0.0f;
                    fa[mm]  += second ? -fp : fp;
                }
            }
            for (int mm = 0; mm < 3; mm++)
            {
                dF[base + mm*st] = fa[mm];
            }
            if (!allAtoms && row == 0)
            {
                dComp[m.slotOffset + c] = 1;
            }
        });
    };
    // Iterations are submitted without waiting, the host checks every few
    // whether all slots converged. Iterations after that do nothing.
    const int    checkEvery = 4;
    int          nsubmitted = 0;
    const double atWall     = 0.9999*maxDist2;
    while (iter < settings.maxIter && !allDone)
    {
        // The squared force of each shell in the convergence criterion and
        // the displacement of each moving shell in the previous step times
        // its force, see shellForces2 and ForceComputer::updateShells
        if (nshellWork > 0)
        {
            queue->parallel_for(sycl::range<1>(nshellWork), [=](sycl::id<1> id)
            {
                uint64_t              w = id[0];
                const DeviceMolecule &m = dMols[findRange(dMols, nm, w, &DeviceMolecule::shellWorkOffset)];
                uint64_t local = w - m.shellWorkOffset;
                uint64_t c     = local % m.nslot;
                uint64_t k     = local / m.nslot;
                if (dDone[m.slotOffset + c])
                {
                    return;
                }
                uint64_t     st = m.stride;
                const float *f  = dF + m.xOffset + c;
                if (k < m.nforceShell)
                {
                    const ForceShell &fs = dFS[m.firstForceShell + k];
                    float fv[3] = { f[(3*fs.shell)*st], f[(3*fs.shell+1)*st], f[(3*fs.shell+2)*st] };
                    double f2;
                    if (!fs.wall)
                    {
                        f2 = fv[0]*fv[0] + fv[1]*fv[1] + fv[2]*fv[2];
                    }
                    else
                    {
                        const int64_t *x = dX + m.xOffset + c;
                        float d[3];
                        for (int mm = 0; mm < 3; mm++)
                        {
                            d[mm] = toNm(difference(x[(3*fs.shell+mm)*st], x[(3*fs.parent+mm)*st]));
                        }
                        f2 = fv[0]*fv[0] + fv[1]*fv[1] + fv[2]*fv[2];
                        double dx2 = d[0]*d[0] + d[1]*d[1] + d[2]*d[2];
                        double fdx = fv[0]*d[0] + fv[1]*d[1] + fv[2]*d[2];
                        if (dx2 >= atWall && fdx > 0)
                        {
                            f2 -= fdx*fdx/dx2;
                        }
                    }
                    dPF2[m.forceShellPartOffset + k*m.nslot + c] = f2;
                }
                else
                {
                    uint64_t     km = k - m.nforceShell;
                    uint32_t     a  = dMS[m.firstMovingShell + km].shell;
                    const float *dx = dDx + m.xOffset + c;
                    float        p  = 0;
                    for (int mm = 0; mm < 3; mm++)
                    {
                        p += dx[(3*a+mm)*st]*f[(3*a+mm)*st];
                    }
                    dPA[m.movingShellPartOffset + km*m.nslot + c] = p;
                }
            });
        }
        // The convergence and the step of each fragment in each slot, from
        // the sums over its shells in their order
        queue->parallel_for(sycl::range<1>(nslotAll), [=](sycl::id<1> id)
        {
            uint64_t s = id[0];
            if (dDone[s])
            {
                return;
            }
            const DeviceMolecule &m  = dMols[findRange(dMols, nm, s, &DeviceMolecule::slotOffset)];
            uint64_t              c  = s - m.slotOffset;
            uint64_t              st = m.stride;
            const float          *f  = dF + m.xOffset + c;
            bool converged = true;
            for (uint32_t fi = 0; fi < m.nfragment; fi++)
            {
                const ShellFragment &fr    = dFrag[m.firstFragment + fi];
                FragmentState       &state = dSt[2*s + fi];
                uint64_t             fs0   = fr.firstForceShell - m.firstForceShell;
                uint64_t             ms0   = fr.firstMovingShell - m.firstMovingShell;
                if (state.moving)
                {
                    // The displacements of the previous step times the forces before it
                    double dxForce = 0;
                    for (uint32_t k = 0; k < fr.nmovingShell; k++)
                    {
                        dxForce += dPD[m.movingShellPartOffset + (ms0 + k)*m.nslot + c];
                    }
                    state.dxfPrev = dxForce;
                    state.moving  = 0;
                }
                double msForce = 0;
                for (uint32_t k = 0; k < fr.nforceShell; k++)
                {
                    msForce += dPF2[m.forceShellPartOffset + (fs0 + k)*m.nslot + c];
                }
                if (fr.nShell > 0)
                {
                    msForce /= fr.nShell;
                }
                // Within the rounding error of the forces, see updateShells
                bool stalled = msForce > toler && msForce >= 0.5*state.msPrev;
                bool atFloor = false;
                if (stalled)
                {
                    // The largest force of the shells now and of the others at the start
                    double fmax = state.fmaxRest;
                    for (uint32_t k = 0; k < fr.nforceShell + fr.nmovingShell; k++)
                    {
                        uint32_t a = (k < fr.nforceShell) ? dFS[fr.firstForceShell + k].shell :
                            dMS[fr.firstMovingShell + k - fr.nforceShell].shell;
                        for (int mm = 0; mm < 3; mm++)
                        {
                            fmax = sycl::fmax(fmax, static_cast<double>(sycl::fabs(f[(3*a+mm)*st])));
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
                if (state.stall >= stallIterations)
                {
                    continue;
                }
                converged = false;
                // The step from the forces before and after the previous one
                double after = 0;
                for (uint32_t k = 0; k < fr.nmovingShell; k++)
                {
                    after += dPA[m.movingShellPartOffset + (ms0 + k)*m.nslot + c];
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
        });
        queue->memcpy(dHist + nsubmitted*nslotAll, dDone, nslotAll);
        nsubmitted++;
        // Move the shells of the fragments that did not converge
        if (nmovingPart > 0)
        {
            queue->parallel_for(sycl::range<1>(nmovingPart), [=](sycl::id<1> id)
            {
                uint64_t              w = id[0];
                const DeviceMolecule &m = dMols[findRange(dMols, nm, w, &DeviceMolecule::movingShellPartOffset)];
                uint64_t local = w - m.movingShellPartOffset;
                uint64_t c     = local % m.nslot;
                uint64_t k     = local / m.nslot;
                uint64_t s     = m.slotOffset + c;
                const FragmentState &state = dSt[2*s + dMF[m.firstMovingShell + k]];
                if (dDone[s] || !state.moving)
                {
                    return;
                }
                uint64_t           st   = m.stride;
                int64_t           *x    = dX + m.xOffset + c;
                const float       *f    = dF + m.xOffset + c;
                float             *dx   = dDx + m.xOffset + c;
                const MovingShell &ms   = dMS[m.firstMovingShell + k];
                double             step = state.step;
                bool               valid = true;
                int64_t *xs = x + 3*ms.shell*st;
                float    fs[3];
                int64_t  xold[3];
                for (int mm = 0; mm < 3; mm++)
                {
                    fs[mm]   = f[(3*ms.shell+mm)*st];
                    xold[mm] = xs[mm*st];
                }
                double disp[3];
                double rcs = 0;
                if (ms.fcHyper == 0)
                {
                    double scale = step*(1/ms.fcShell);
                    for (int mm = 0; mm < 3; mm++)
                    {
                        disp[mm] = scale*fs[mm];
                    }
                }
                else
                {
                    float r[3];
                    for (int mm = 0; mm < 3; mm++)
                    {
                        r[mm] = toNm(difference(xs[mm*st], x[(3*ms.core+mm)*st]));
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
                    for (int mm = 0; mm < 3; mm++)
                    {
                        disp[mm] = step*(fs[mm]*ikperp + fpar*r[mm]);
                    }
                }
                if (ms.fcHyper > 0 && rcs == 0)
                {
                    // On top of the core, along the force, see update_shell_pos
                    float  fnorm = sycl::sqrt(fs[0]*fs[0] + fs[1]*fs[1] + fs[2]*fs[2]);
                    double scale = step*fscale(ms.fcShell, ms.rHyper > 0 ? 0 : ms.fcHyper, fnorm);
                    for (int mm = 0; mm < 3; mm++)
                    {
                        xs[mm*st] = move(xs[mm*st], -(fs[mm]*scale), &valid);
                    }
                }
                else
                {
                    for (int mm = 0; mm < 3; mm++)
                    {
                        xs[mm*st] = move(xs[mm*st], disp[mm], &valid);
                    }
                }
                bool clipped = false;
                if (ms.wall)
                {
                    float d[3];
                    for (int mm = 0; mm < 3; mm++)
                    {
                        d[mm] = toNm(difference(xs[mm*st], x[(3*ms.parent+mm)*st]));
                    }
                    float dx2 = d[0]*d[0] + d[1]*d[1] + d[2]*d[2];
                    if (dx2 > maxDist2)
                    {
                        float scale = sycl::sqrt(static_cast<float>(maxDist2/dx2));
                        for (int mm = 0; mm < 3; mm++)
                        {
                            double xp = static_cast<double>(x[(3*ms.parent+mm)*st])*coordinateUnitNm;
                            xs[mm*st] = toUnits(xp + scale*d[mm]);
                        }
                        clipped = true;
                    }
                }
                // Shells held at the maximum distance do not tell how far to step
                float  dxj[3]  = { 0, 0, 0 };
                double dxForce = 0;
                if (!clipped)
                {
                    for (int mm = 0; mm < 3; mm++)
                    {
                        dxj[mm] = toNm(difference(xs[mm*st], xold[mm]));
                    }
                    dxForce = dxj[0]*fs[0] + dxj[1]*fs[1] + dxj[2]*fs[2];
                }
                for (int mm = 0; mm < 3; mm++)
                {
                    dx[(3*ms.shell+mm)*st] = dxj[mm];
                }
                dPD[m.movingShellPartOffset + k*m.nslot + c] = dxForce;
                if (!valid)
                {
                    *dVal = 1;
                }
            });
        }
        if (prof)
        {
            profile.stop(queue, RelaxProfile::Update);
        }
        forceKernel(false);
        if (prof)
        {
            profile.stop(queue, RelaxProfile::Forces);
        }
        iter++;
        if (nsubmitted % checkEvery == 0 || iter >= settings.maxIter)
        {
            queue->memcpy(hDone, dDone, nslotAll);
            queue->wait_and_throw();
            allDone = std::count(hDone, hDone + nslotAll, 0) == 0;
            if (prof)
            {
                profile.stop(queue, RelaxProfile::Check);
            }
        }
    }
    // The slots computed in each iteration, those that did not converge
    queue->memcpy(hDone, dHist, nsubmitted*nslotAll);
    queue->wait_and_throw();
    for (int i = 0; i < nsubmitted && nwork > 0; i++)
    {
        *nComputed += std::count(hDone + i*nslotAll, hDone + (i+1)*nslotAll, 0);
    }
    // The forces on the other particles, from the coordinates of the last
    // calculation of each slot, if they are needed
    if (allForces)
    {
        forceKernel(true);
    }
    // The energies of the last calculation of each slot
    if (nwork > 0)
    {
        energyKernel();
        // The energies of each potential in the order of the interactions
        queue->parallel_for(sycl::range<1>(nblockSlot), [=](sycl::id<1> id)
        {
            uint64_t           s = id[0];
            const DeviceBlock &b = dBlk[findRange(dBlk, nb, s, &DeviceBlock::blockSlotOffset)];
            const DeviceMolecule &m = dMols[b.mol];
            uint64_t c = s - b.blockSlotOffset;
            if (!dComp[m.slotOffset + c])
            {
                return;
            }
            for (int k = 0; k < pairEnergyTerms; k++)
            {
                double sum = 0;
                for (uint64_t p = 0; p < b.npair; p++)
                {
                    sum += dEp[k*nwork + b.workOffset + p*m.nslot + c];
                }
                dE[b.energyOffset + k*m.stride + c] = sum;
            }
        });
    }
    if (prof)
    {
        profile.stop(queue, RelaxProfile::Energies);
    }
    // The results: only the shells move, and only their forces are needed
    // unless the forces on all particles are, which saves most of the copying
    PackedArrays out;
    size_t oRX   = out.add<int64_t>(nresult);
    size_t oRF   = out.add<float>(allForces ? 0 : nresult);
    char  *dRes  = mem.result.reserve(*queue, std::max<size_t>(out.size(), 16));
    char  *hRes  = mem.resultHost.reserve(*queue, std::max<size_t>(out.size(), 16));
    if (nresult > 0)
    {
        auto dRX   = reinterpret_cast<int64_t *>(dRes + oRX);
        auto dRF   = reinterpret_cast<float *>(dRes + oRF);
        bool onlyX = allForces;
        queue->parallel_for(sycl::range<1>(nresult), [=](sycl::id<1> id)
        {
            uint64_t              w = id[0];
            const DeviceMolecule &m = dMols[findRange(dMols, nm, w, &DeviceMolecule::resultOffset)];
            uint64_t local = w - m.resultOffset;
            uint64_t c     = local % m.nslot;
            uint64_t rm    = local / m.nslot;
            uint64_t a     = dSR[m.firstShellRow + rm/3];
            uint64_t i     = m.xOffset + (3*a + rm%3)*m.stride + c;
            dRX[w] = dX[i];
            if (!onlyX)
            {
                dRF[w] = dF[i];
            }
        });
        queue->memcpy(hRes, dRes, out.size());
    }
    if (allForces)
    {
        queue->memcpy(h + oF, dF, nx*sizeof(float));
    }
    // The slots whose shells stopped without progress did not converge,
    // as those left after the last iteration
    queue->parallel_for(sycl::range<1>(nslotAll), [=](sycl::id<1> id)
    {
        uint64_t              s = id[0];
        const DeviceMolecule &m = dMols[findRange(dMols, nm, s, &DeviceMolecule::slotOffset)];
        for (uint32_t fi = 0; fi < m.nfragment; fi++)
        {
            if (dSt[2*s + fi].stall >= stallIterations)
            {
                dDone[s] = 0;
            }
        }
    });
    // The done and computed flags, which are next to each other
    queue->memcpy(h + oDone, dDone, (oComputed - oDone) + nslotAll);
    double *hE = mem.energiesHost.reserve(*queue, std::max<uint64_t>(nenergy, 1));
    if (nenergy > 0)
    {
        queue->memcpy(hE, dE, nenergy*sizeof(double));
    }
    int validFlag = 0;
    queue->memcpy(&validFlag, dVal, sizeof(int));
    queue->wait_and_throw();
    *notFinite = validFlag != 0;
    if (prof)
    {
        profile.stop(queue, RelaxProfile::Download);
        profile.done(nsubmitted, nwork, nslotAll);
    }
    ib = 0;
    for (size_t k = 0; k < nmol; k++)
    {
        const auto &m = mols[k];
        const auto &d = dmols[k];
        const char *comp = h + oComputed + d.slotOffset;
        const auto *rx   = reinterpret_cast<const int64_t *>(hRes + oRX) + d.resultOffset;
        const auto *rf   = reinterpret_cast<const float *>(hRes + oRF) + d.resultOffset;
        for (size_t rm = 0; rm < 3*d.nshellRow; rm++)
        {
            size_t row = 3*shellRows[d.firstShellRow + rm/3] + rm%3;
            for (size_t c = 0; c < m.nslot; c++)
            {
                m.x[row*m.stride + c] = rx[rm*m.nslot + c];
                if (!allForces && comp[c])
                {
                    m.f[row*m.stride + c] = rf[rm*m.nslot + c];
                }
            }
        }
        if (allForces)
        {
            const auto *hf = reinterpret_cast<const float *>(h + oF) + d.xOffset;
            for (size_t row = 0; row < 3*m.natom; row++)
            {
                for (size_t c = 0; c < m.nslot; c++)
                {
                    if (comp[c])
                    {
                        m.f[row*m.stride + c] = hf[row*m.stride + c];
                    }
                }
            }
        }
        std::copy(h + oDone + d.slotOffset, h + oDone + d.slotOffset + m.nslot, m.done);
        std::copy(comp, comp + m.nslot, m.computed);
        for (size_t b = 0; b < m.nblock; b++, ib++)
        {
            const double *eb = hE + dblocks[ib].energyOffset;
            for (int t = 0; t < pairEnergyTerms; t++)
            {
                for (size_t c = 0; c < m.nslot; c++)
                {
                    if (comp[c])
                    {
                        m.blocks[b].energies[t*m.stride + c] = eb[t*m.stride + c];
                    }
                }
            }
        }
    }
    return true;
}

} // namespace gpu

} // namespace alexandria

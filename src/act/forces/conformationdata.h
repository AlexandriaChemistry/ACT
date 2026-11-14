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
 * Coordinates, forces and energies of all conformations of a topology,
 * stored such that loops over conformations access memory contiguously.
 *
 * The ACT energies of the conformations are kept here rather than in the
 * energy maps of their CoordsEnergy, which hold the QM energies, since
 * visiting the energy maps of thousands of conformations repeatedly is
 * limited by memory access. The coordinates and forces are kept here as
 * well, for all conformations together.
 * \author David van der Spoel <david.vanderspoel@icm.uu.se>
 */
#ifndef ACT_FORCES_CONFORMATIONDATA_H
#define ACT_FORCES_CONFORMATIONDATA_H

#include <algorithm>
#include <array>
#include <cstddef>
#include <initializer_list>
#include <limits>
#include <map>
#include <type_traits>
#include <utility>
#include <vector>

#include "act/basics/act_particle.h"
#include "act/basics/interactiontype.h"
#include "act/forcefield/forcefield_parameterlist.h"
#include "act/forcefield/potential.h"
#include "act/forces/forces_explicit_simd.h"
#include "act/topology/actatom.h"
#include "act/topology/coords_energy.h"
#include "act/utility/timings.h"
#include "gromacs/math/vectypes.h"
#include "gromacs/utility/gmxassert.h"

namespace alexandria
{

/*! \brief Swap two columns of arrays in the layout of ConformationArrays
 * \param[inout] v      The arrays
 * \param[in]    stride The distance between the rows
 * \param[in]    a      The first column
 * \param[in]    b      The second column
 */
template <typename U>
void swapColumns(std::vector<U> *v, size_t stride, size_t a, size_t b)
{
    size_t nrow = stride > 0 ? v->size()/stride : 0;
    U     *p    = v->data();
    for(size_t row = 0; row < nrow; row++)
    {
        std::swap(p[row*stride + a], p[row*stride + b]);
    }
}

/*! \brief Coordinates and forces of the conformations of a topology.
 *
 * Component m of particle a in slot s is stored at (3*a+m)*stride + s,
 * such that loops over conformations access memory contiguously and can be
 * vectorized. Each conformation has its own slot, see slot. While the
 * shells are relaxed, the conformations for which they did not converge
 * yet are put in the first nslot slots, see arrange, such that only those
 * are computed. Otherwise slot c holds conformation c.
 * The computations use nslot slots, a multiple of the explicit SIMD
 * width. Those beyond the conformations are copies of the last one, see
 * pad, they are computed but not used. The stride, the distance between
 * the rows, is at least the number of conformations rounded up to the
 * width, and one width more if that is a multiple of 256. It is larger
 * when conformations were added one at a time, such that the rows do not
 * have to be moved each time.
 * The coordinates are integers, see coordinate.h, in both precisions.
 * \tparam T The floating point type of the forces, float or double, such
 *           that the explicit SIMD of that precision computes them, see
 *           computeExplicitSimd
 */
template <typename T>
class ConformationArrays
{
public:
    //! Number of conformations
    size_t              nconf   = 0;
    //! Number of conformations to compute, in the first nslot slots
    size_t              nactive = 0;
    //! Number of slots that are computed, nactive rounded up to the width
    size_t              nslot   = 0;
    //! Distance between the rows
    size_t              stride  = 0;
    //! Number of particles
    size_t              natom   = 0;
    //! Coordinates in units of coordinateUnit, see coordinate.h
    std::vector<Coordinate> x;
    //! Forces
    std::vector<T>      f;
    /*! Energies of one potential, kept here to reuse the memory, in double
     * precision since sums of many interactions are taken
     */
    std::vector<double> energies;

    /*! \return the slot of a conformation
     * \param[in] c The conformation
     */
    size_t slot(size_t c) const { return slot_.empty() ? c : slot_[c]; }

    /*! \return the conformation in a slot, of the first nconf
     * \param[in] s The slot
     */
    size_t conformation(size_t s) const { return conf_.empty() ? s : conf_[s]; }

    //! \return Whether slot c holds conformation c for all conformations
    bool inOrder() const { return slot_.empty(); }

    //! \return Whether the last calculation computed conformation c
    bool computed(size_t c) const { return slot(c) < nslot; }

    //! \return the number of conformations that the last calculation computed
    size_t ncomputed() const { return std::min(nslot, nconf); }

    /*! \brief Change the number of conformations and particles
     * The coordinates and forces of the conformations and particles that
     * are in both the old and the new arrays are kept, the others are zero.
     * The slots are put in order first.
     * \param[in] nc The number of conformations
     * \param[in] na The number of particles
     */
    void resize(size_t nc, size_t na)
    {
        restoreOrder();
        size_t width = explicitSimdWidth<T>();
        size_t ns    = ((nc + width - 1)/width)*width;
        size_t st    = stride;
        // Grow by half of the stride if conformations are added one at a
        // time, otherwise use as few slots as possible
        if (ns > stride || na != natom || 2*ns < stride)
        {
            st = (na == natom && ns > stride && stride > 0) ? std::max(ns, stride + stride/2) : ns;
            st = ((st + width - 1)/width)*width;
            // Rows a multiple of 2 kB apart map to the same sets of the caches,
            // which then hold few of them
            if (st % 256 == 0)
            {
                st += width;
            }
        }
        if (st != stride || na != natom)
        {
            std::vector<Coordinate> xn(DIM*na*st, 0);
            std::vector<T>          fn(DIM*na*st, 0);
            size_t         nrow = DIM*std::min(na, natom);
            size_t         nc0  = std::min(nc, nconf);
            for(size_t row = 0; row < nrow; row++)
            {
                std::copy(x.begin() + row*stride, x.begin() + row*stride + nc0, xn.begin() + row*st);
                std::copy(f.begin() + row*stride, f.begin() + row*stride + nc0, fn.begin() + row*st);
            }
            x.swap(xn);
            f.swap(fn);
        }
        nconf   = nc;
        nactive = nc;
        nslot   = ns;
        stride  = st;
        natom   = na;
    }

    /*! \brief Put the conformations to compute in the first slots
     * These are the nslot computed slots, the active conformations rounded
     * up to the explicit SIMD width, the conformations that are not active
     * in them are computed as well. Active conformations in other slots
     * are swapped with those that are not active in the computed slots,
     * which saves swaps when few of them converged. The slots of other
     * conformations are swapped with theirs, together
     * with the columns of other arrays in the same order, such as the forces
     * of the entries without shells.
     * \param[in]    active The conformations to compute, in increasing order
     * \param[inout] otherU Arrays in the same order and their strides, may be empty
     * \param[inout] otherV Likewise, of another type
     */
    template <typename U, typename V>
    void arrange(const std::vector<size_t>                            &active,
                 std::initializer_list<std::pair<std::vector<U> *, size_t>> otherU,
                 std::initializer_list<std::pair<std::vector<V> *, size_t>> otherV)
    {
        if (active.size() == nconf)
        {
            restoreOrder(otherU, otherV);
        }
        else
        {
            if (slot_.empty())
            {
                slot_.resize(nconf);
                conf_.resize(nconf);
                for(size_t c = 0; c < nconf; c++)
                {
                    slot_[c] = c;
                    conf_[c] = c;
                }
            }
            // The computed slots, those of conformations that are not
            // active are computed as well, with the same result as before
            size_t width = explicitSimdWidth<T>();
            size_t na    = ((active.size() + width - 1)/width)*width;
            // The active conformations in slots from na on swap with the
            // conformations that are not active in slots below na. Active
            // conformations below na stay where they are.
            std::vector<char> isActive(nconf, 0);
            for(auto c : active)
            {
                isActive[c] = 1;
            }
            size_t free = 0;
            for(auto c : active)
            {
                if (slot_[c] < na)
                {
                    continue;
                }
                while (isActive[conf_[free]])
                {
                    free++;
                }
                swapSlots(slot_[c], free, otherU, otherV);
            }
        }
        size_t width = explicitSimdWidth<T>();
        nactive      = active.size();
        nslot        = ((nactive + width - 1)/width)*width;
    }

    /*! \brief Put slot c back to conformation c for all conformations
     * \param[inout] otherU Arrays in the same order and their strides
     * \param[inout] otherV Likewise
     */
    template <typename U = T, typename V = T>
    void restoreOrder(std::initializer_list<std::pair<std::vector<U> *, size_t>> otherU = {},
                      std::initializer_list<std::pair<std::vector<V> *, size_t>> otherV = {})
    {
        if (!slot_.empty())
        {
            for(size_t s = 0; s < nconf; s++)
            {
                // Bring conformation s into slot s
                if (conf_[s] != s)
                {
                    swapSlots(s, slot_[s], otherU, otherV);
                }
            }
            slot_.clear();
            conf_.clear();
        }
        size_t width = explicitSimdWidth<T>();
        nactive      = nconf;
        nslot        = ((nconf + width - 1)/width)*width;
    }

    /*! \brief Copy the coordinates of the last conformation to the padding
     * Has to be called after changing the coordinates of the conformations
     * and before computing, such that the computed slots beyond the
     * conformations hold a meaningful conformation.
     */
    void pad()
    {
        if (nconf == 0)
        {
            return;
        }
        for(size_t row = 0; row < DIM*natom; row++)
        {
            Coordinate *xr = x.data() + row*stride;
            std::fill(xr + std::min(nconf, nslot), xr + nslot, xr[nconf - 1]);
        }
    }

    //! Set the forces of the computed slots to zero
    void clearForces()
    {
        for(size_t row = 0; row < DIM*natom; row++)
        {
            std::fill(f.begin() + row*stride, f.begin() + row*stride + nslot, 0);
        }
    }

    /*! \brief Copy the coordinates of the computed slots of other arrays
     * The slots are those of src, the stride and the number of computed
     * slots are those of this precision, such that the arrays can be
     * computed in another precision.
     * \param[in] src The arrays to copy
     */
    template <typename U>
    void copyCoordinates(const ConformationArrays<U> &src)
    {
        static const size_t timer = Timings::instance().index("Copy coordinates to arrays");
        const ScopedTimer   scoped(timer);
        if (src.nconf != nconf || src.natom != natom)
        {
            slot_.clear();
            conf_.clear();
            resize(src.nconf, src.natom);
        }
        slot_ = src.slot_;
        conf_ = src.conf_;
        size_t width = explicitSimdWidth<T>();
        nactive      = src.nactive;
        // The conformations to compute may be in any of the computed slots
        // of src, which may be more than this precision has
        nslot        = ((std::min(src.nslot, nconf) + width - 1)/width)*width;
        for(size_t row = 0; row < DIM*natom; row++)
        {
            const Coordinate *s = src.x.data() + row*src.stride;
            std::copy(s, s + ncomputed(), x.begin() + row*stride);
        }
        pad();
    }

    /*! \brief Copy the forces of the computed slots to other arrays in the same order
     * \param[inout] dst The arrays
     */
    template <typename U>
    void storeForces(ConformationArrays<U> *dst) const
    {
        static const size_t timer = Timings::instance().index("Copy forces from arrays");
        const ScopedTimer   scoped(timer);
        GMX_RELEASE_ASSERT(dst->nconf == nconf && dst->natom == natom, "Different conformations");
        for(size_t row = 0; row < DIM*natom; row++)
        {
            const T *s = f.data() + row*stride;
            std::copy(s, s + ncomputed(), dst->f.begin() + row*dst->stride);
        }
    }

    /*! \brief Copy the forces of the computed slots of other arrays in the same order
     * \param[in] src The arrays
     */
    template <typename U>
    void copyForces(const ConformationArrays<U> &src)
    {
        static const size_t timer = Timings::instance().index("Copy forces to arrays");
        const ScopedTimer   scoped(timer);
        GMX_RELEASE_ASSERT(src.nconf == nconf && src.natom == natom, "Different conformations");
        for(size_t row = 0; row < DIM*natom; row++)
        {
            const U *s = src.f.data() + row*src.stride;
            std::copy(s, s + ncomputed(), f.begin() + row*stride);
        }
    }

    /*! \brief Set the forces of the computed slots from forces kept in the same layout
     * \param[in] src The forces, with the same stride and order
     */
    void loadForces(const std::vector<T> &src)
    {
        static const size_t timer = Timings::instance().index("Copy fixed forces to arrays");
        const ScopedTimer   scoped(timer);
        GMX_RELEASE_ASSERT(src.size() == f.size(), "Forces kept with another layout");
        for(size_t row = 0; row < DIM*natom; row++)
        {
            std::copy(src.begin() + row*stride, src.begin() + row*stride + nslot, f.begin() + row*stride);
        }
    }

private:
    template <typename U> friend class ConformationArrays;
    //! The slot of each conformation, empty if slot c holds conformation c
    std::vector<size_t> slot_;
    //! The conformation in each slot, empty likewise
    std::vector<size_t> conf_;

    /*! \brief Swap the conformations of two slots
     * \param[in]    a      A slot
     * \param[in]    b      Another slot
     * \param[inout] otherU Arrays in the same order and their strides
     * \param[inout] otherV Likewise
     */
    template <typename U, typename V>
    void swapSlots(size_t a, size_t b,
                   std::initializer_list<std::pair<std::vector<U> *, size_t>> otherU,
                   std::initializer_list<std::pair<std::vector<V> *, size_t>> otherV)
    {
        swapColumns(&x, stride, a, b);
        swapColumns(&f, stride, a, b);
        for(auto &o : otherU)
        {
            if (!o.first->empty())
            {
                swapColumns(o.first, o.second, a, b);
            }
        }
        for(auto &o : otherV)
        {
            if (!o.first->empty())
            {
                swapColumns(o.first, o.second, a, b);
            }
        }
        std::swap(conf_[a], conf_[b]);
        slot_[conf_[a]] = a;
        slot_[conf_[b]] = b;
    }
};

/*! \brief The ACT energies of all conformations, per interaction type
 * For each interaction type and conformation whether there is an energy,
 * the energy and its backup. The energies of one interaction type are
 * stored contiguously, see energies() for those of one conformation.
 */
class ConformationEnergies
{
public:
    //! Pointers to the data of one interaction type, indexed by conformation
    struct Term
    {
        //! The energies
        double *act    = nullptr;
        //! The backups of the energies
        double *before = nullptr;
        //! Whether there is an energy
        char   *have   = nullptr;

        /*! \brief Whether there is an energy for conformation c
         * \param[in] c The conformation index
         * \return Whether there is an energy for conformation c
         */
        bool has(size_t c) const { return have && have[c]; }
        /*! \brief Set the energy of conformation c, a new energy has no backup
         * \param[in] c     The conformation index
         * \param[in] value The energy
         */
        void set(size_t c, double value)
        {
            if (!have[c])
            {
                before[c] = 0;
                have[c]   = 1;
            }
            act[c] = value;
        }
        /*! \brief Restore the backup of conformation c, as ACTEnergy::restoreBackup
         * \param[in] c The conformation index
         */
        void restoreBackup(size_t c)
        {
            act[c]    = before[c];
            before[c] = 0;
        }
    };

    //! Constructor, without conformations or interaction types
    ConformationEnergies() { index_.fill(-1); }

    /*! \brief Set the number of conformations
     * All energies are removed if the number changes.
     * \param[in] nconf The number of conformations
     */
    void setNumberOfConformations(size_t nconf)
    {
        if (nconf != nconf_)
        {
            index_.fill(-1);
            types_.clear();
            data_.clear();
            nconf_ = nconf;
        }
    }

    //! \return The number of conformations
    size_t nconf() const { return nconf_; }

    /*! \brief Find the data for an interaction type, add it if not present
     * \param[in] t The interaction type
     * \return The data for interaction type t
     */
    Term add(InteractionType t)
    {
        auto &i = index_[static_cast<size_t>(t)];
        if (i < 0)
        {
            i = static_cast<int>(data_.size());
            data_.push_back({ std::vector<double>(nconf_, 0.0), std::vector<double>(nconf_, 0.0),
                              std::vector<char>(nconf_, 0) });
            types_.insert(std::lower_bound(types_.begin(), types_.end(), t), t);
        }
        return term(i);
    }

    /*! \brief Find the data for an interaction type
     * \param[in] t The interaction type
     * \return The data for interaction type t, with null pointers if not present
     */
    Term find(InteractionType t)
    {
        auto i = index_[static_cast<size_t>(t)];
        return i < 0 ? Term() : term(i);
    }

    //! \return The interaction types present, in increasing order like in an energy map
    const std::vector<InteractionType> &types() const { return types_; }

    /*! \brief Set whether the backups are the energies before relaxing shells
     * Without shells the energies do not change, and the backups are not
     * made. They are then taken to be the same as the energies.
     * \param[in] backups Whether there are backups
     */
    void setBackups(bool backups) { backups_ = backups; }

    //! \return Whether the backups are the energies before relaxing shells, see setBackups
    bool haveBackups() const { return backups_; }

    /*! \brief Set whether new energies are added to the present ones
     * Normally the energies of an interaction type are replaced, but when the
     * entries of a type are computed in parts they have to be added.
     * \param[in] accumulate Whether to add the energies
     */
    void setAccumulate(bool accumulate) { accumulate_ = accumulate; }

    //! \return Whether new energies are added to the present ones, see setAccumulate
    bool accumulate() const { return accumulate_; }

    /*! \brief Set the energy terms of the next potential
     * Some potentials give the energy of different terms, depending on
     * the interaction type they are used for.
     * \param[in] types The terms, or nullptr for those of the potential
     */
    void setTermTypes(const std::vector<InteractionType> *types) { termTypes_ = types; }

    //! \return The energy terms of the next potential, see setTermTypes
    const std::vector<InteractionType> *termTypes() const { return termTypes_; }

    /*! \brief The energies of one conformation
     * \param[in] c The conformation
     * \return The energy of each interaction type that conformation c has
     */
    std::map<InteractionType, double> energies(size_t c) const
    {
        std::map<InteractionType, double> emap;
        for(auto t : types_)
        {
            const auto &d = data_[index_[static_cast<size_t>(t)]];
            if (c < nconf_ && d.have[c])
            {
                emap.insert({ t, d.act[c] });
            }
        }
        return emap;
    }

private:
    //! Number of interaction types
    static constexpr size_t ntype_ = static_cast<size_t>(InteractionType::CHARGE) + 1;
    //! The data of one interaction type
    struct Data
    {
        //! The energies
        std::vector<double> act;
        //! The backups
        std::vector<double> before;
        //! Whether there is an energy
        std::vector<char>   have;
    };
    //! Number of conformations
    size_t                      nconf_ = 0;
    //! Whether the backups are the energies before relaxing shells, see setBackups
    bool                        backups_ = false;
    //! Whether new energies are added to the present ones, see setAccumulate
    bool                        accumulate_ = false;
    //! The energy terms of the next potential, see setTermTypes
    const std::vector<InteractionType> *termTypes_ = nullptr;
    //! Position of each interaction type in data_, -1 if not present
    std::array<int, ntype_>     index_;
    //! The interaction types present, sorted
    std::vector<InteractionType> types_;
    //! The data, the vectors keep their memory when data_ grows
    std::vector<Data>           data_;

    /*! \brief Get pointers to the data
     * \param[in] i Position in data_
     * \return Pointers to the data at position i
     */
    Term term(int i)
    {
        auto &d = data_[i];
        return { d.act.data(), d.before.data(), d.have.data() };
    }
};

/*! \brief Energies of the entries of one interaction type in all conformations
 * They can be used again as long as the parameters, the coordinates and,
 * if the potential uses them, the charges did not change.
 */
struct CachedEnergies
{
    //! Whether energies have been stored
    bool                valid             = false;
    //! The potential
    Potential           potential         = Potential::NONE;
    //! Topology::parameterVersion of the interaction type
    size_t              parameterVersion  = 0;
    //! ConformationData::coordinates
    size_t              coordinateVersion = 0;
    //! Number of entries
    size_t              nentries          = 0;
    //! Number of conformations
    size_t              nconf             = 0;
    //! The charges of the atoms, if the potential uses them
    std::vector<double> charges;
    //! The energies, as passed to storeEnergies
    std::vector<double> energies;
    //! Distance between the energies of conformations
    size_t              confStride        = 0;
    //! Distance between the energy terms of a conformation
    size_t              termStride        = 0;
    /*! For Coulomb, the energy of each term of each interaction in each
     * slot of the arrays for charges of one, see coulombKernelsExplicitSimd.
     * The energies for other charges follow from it for the interactions
     * whose widths and coordinates do not change.
     */
    std::vector<real>   kernel;
    //! Whether kernel has been computed
    bool                kernelValid             = false;
    //! ConformationData::coordinates of kernel
    size_t              kernelCoordinateVersion = 0;
    //! Number of interactions in kernel
    size_t              kernelEntries           = 0;
    //! Number of conformations in kernel
    size_t              kernelConformations     = 0;
    //! Distance between the interactions in kernel
    size_t              kernelStride            = 0;
    /*! For each interaction, the widths kernel was computed with: those of
     * the charges of the first particle, then of the second, such that only
     * the interactions with changed widths are computed again
     */
    std::vector<double> kernelWidths;
    //! For each interaction, bit k tells whether term k is in kernel
    std::vector<unsigned char> kernelTerms;
    /*! The coordinates of kernel, as ConformationArrays<real>::x, such that only
     * the interactions of particles that moved are computed again
     */
    std::vector<Coordinate> kernelX;
};

/*! \brief Energies and forces of the entries without shells in all conformations
 * While relaxing the shells, only the entries with shells change, see
 * Topology::separateShells. The others are computed once, in the first
 * calculation of the conformations, and kept here.
 */
struct FixedEnergies
{
    //! Whether energies and forces have been stored
    bool                                valid = false;
    //! The number of conformations
    size_t                              nconf = 0;
    //! The forces in the layout of ConformationArrays<real>
    std::vector<real>                   arrayForces;
    //! The same for ConformationArrays<double>, in a single precision build
    std::vector<double>                 arrayForcesDouble;
    //! The interaction types, including EPOT
    std::vector<InteractionType>        types;
    //! For each type the energy of each conformation
    std::vector<std::vector<double>>    energies;
    //! For each type whether each conformation has an energy
    std::vector<std::vector<char>>      have;

    /*! \brief The forces in the layout of ConformationArrays<T>
     * \tparam T The floating point type of the arrays
     * \return arrayForces or arrayForcesDouble
     */
    template <typename T>
    std::vector<T> &forcesOfArrays()
    {
        if constexpr (std::is_same_v<T, real>)
        {
            return arrayForces;
        }
        else
        {
            return arrayForcesDouble;
        }
    }
};

/*! \brief Reference data of all conformations for computing deviations
 * The quantum chemistry energies and the Boltzmann weights do not change
 * when the parameters change, they are kept here by the code computing the
 * deviations, see ForceEnergyDevComputer.
 */
struct DeviationReference
{
    //! ConformationData::conformations of the conformations
    size_t                           conformationVersion = 0;
    //! The inverse Boltzmann temperature
    double                           beta              = 0;
    //! Whether the induction correction is separate from the induction
    bool                             separate          = false;
    //! The interaction types
    std::vector<InteractionType>     itypes;
    //! Weight of each conformation
    std::vector<double>              weight;
    //! Reference energy of each interaction type and conformation
    std::vector<std::vector<double>> energy;
    //! Whether there is a reference energy for each interaction type and conformation
    std::vector<std::vector<char>>   have;
};

/*! \brief State of the shell relaxation of the conformations, see ForceComputer::updateShells
 * In the layout of ConformationArrays and with the slots of the arrays, such
 * that the loops over the conformations access it in order. Row i of step,
 * force and dxForce holds the values of the shells of fragment i of a dimer.
 */
struct ShellRelaxation
{
    //! Distance between the rows
    size_t              stride = 0;
    //! Fraction of the Newton step by which the shells move
    std::vector<double> step;
    //! Mean square shell force in the previous relaxation step
    std::vector<double> force;
    /*! The sum over the shells of their displacement in the previous
     * relaxation step times the force before it, zero if there was none
     */
    std::vector<double> dxForce;
    /*! The mean square shell force against which progress is measured,
     * the last that was at most half the one before, see ForceComputer::setStallIterations
     */
    std::vector<double> reference;
    //! The number of relaxation steps since the mean square force was halved, see reference
    std::vector<double> stall;
    //! Displacements of the shells in the previous step, component m of particle a in row DIM*a+m
    std::vector<real>   dx;
    //! Work vectors of ForceComputer::updateShells, kept to use the memory again
    std::vector<double> msForce, fmax, after, dxStep, disp, rcs;
    //! Likewise, whether the shells in each slot are relaxed
    std::vector<char>   relax;
    //! Likewise, whether a shell in each slot moved in the last step
    std::vector<char>   moved;

    /*! \brief Start the shell relaxation with full steps
     * \param[in] st    The stride of the arrays
     * \param[in] natom The number of particles
     */
    void reset(size_t st, size_t natom)
    {
        stride = st;
        step.assign(2*st, 1.0);
        force.assign(2*st, std::numeric_limits<double>::max());
        dxForce.assign(2*st, 0.0);
        reference.assign(2*st, std::numeric_limits<double>::max());
        stall.assign(2*st, 0.0);
        dx.assign(DIM*natom*st, 0);
    }
};

//! Coordinates, forces and energies of all conformations of a topology
struct ConformationData
{
    //! The conformations, see x and f for their coordinates and forces
    std::vector<CoordsEnergy> confs;

    //! \return the number of conformations
    size_t nconf() const { return confs.size(); }

    /*! \brief Return the coordinates of a conformation
     * \param[in] c The conformation
     * \return the coordinates
     */
    CoordinateVectors<Coordinate> x(size_t c)
    {
        return { arrays.x.data() + arrays.slot(c), arrays.stride, arrays.natom };
    }

    //! \copydoc x
    CoordinateVectors<const Coordinate> x(size_t c) const
    {
        return { arrays.x.data() + arrays.slot(c), arrays.stride, arrays.natom };
    }

    /*! \brief Return the computed forces of a conformation
     * \param[in] c The conformation
     * \return the forces, as many as there are particles
     */
    ParticleVectors<real> f(size_t c)
    {
        return { arrays.f.data() + arrays.slot(c), arrays.stride, arrays.natom };
    }

    //! \copydoc f
    ParticleVectors<const real> f(size_t c) const
    {
        return { arrays.f.data() + arrays.slot(c), arrays.stride, arrays.natom };
    }

    /*! \brief Set the coordinates of a conformation, with zero forces
     * The arrays are resized to the conformations in confs, all
     * conformations have to have the same number of particles. With
     * another number than before, the other conformations have to be set
     * again as well.
     * \param[in] c     The conformation
     * \param[in] begin The coordinates of the first particle
     * \param[in] end   Past the coordinates of the last particle
     */
    void setX(size_t                                 c,
              std::vector<gmx::RVec>::const_iterator begin,
              std::vector<gmx::RVec>::const_iterator end)
    {
        size_t natom = end - begin;
        if (arrays.nconf != confs.size() || arrays.natom != natom)
        {
            arrays.resize(confs.size(), natom);
        }
        const auto xc = x(c);
        const auto fc = f(c);
        for(size_t a = 0; a < natom; a++)
        {
            for(int m = 0; m < DIM; m++)
            {
                xc[a][m] = begin[a][m];
                fc[a][m] = 0;
            }
        }
    }

    /*! \brief Set the coordinates of a conformation, with zero forces
     * \param[in] c The conformation
     * \param[in] x The coordinates of all particles
     */
    void setX(size_t c, const std::vector<gmx::RVec> &x) { setX(c, x.begin(), x.end()); }

    /*! \brief Set the coordinates of a conformation from integers, with zero forces
     * As the other setX, the integers are copied, without rounding.
     * \param[in] c The conformation
     * \param[in] x The coordinates of all particles
     */
    void setX(size_t c, const CoordinateVectors<const Coordinate> &x)
    {
        size_t natom = x.size();
        if (arrays.nconf != confs.size() || arrays.natom != natom)
        {
            arrays.resize(confs.size(), natom);
        }
        const auto xc = this->x(c);
        const auto fc = f(c);
        for(size_t a = 0; a < natom; a++)
        {
            for(int m = 0; m < DIM; m++)
            {
                xc[a][m] = x[a][m];
                fc[a][m] = 0;
            }
        }
    }

    //! Coordinates and forces of the conformations, see x and f
    ConformationArrays<real>   arrays;
    //! Copies of the coordinates and forces in double precision, used if doublePrecision is set
    ConformationArrays<double> arraysDouble;
    /*! Whether to compute the potentials in double precision in a single
     * precision build, for methods that need it such as the normal mode
     * analysis. The coordinates are then copied to arraysDouble and the
     * forces back to arrays. Potentials without explicit SIMD, and single
     * conformations, are still computed in real precision.
     */
    bool                 doublePrecision = false;

    /*! \brief The arrays of precision T
     * \tparam T The floating point type
     * \return arrays or arraysDouble
     */
    template <typename T>
    ConformationArrays<T> &arraysOf()
    {
        if constexpr (std::is_same_v<T, real>)
        {
            return arrays;
        }
        else
        {
            return arraysDouble;
        }
    }

    /*! \brief Put the conformations to compute in the first slots of the arrays
     * The forces of the entries without shells and the state of the shell
     * relaxation are kept in the same order.
     * With all conformations, slot c holds conformation c again.
     * \param[in] active The conformations to compute, in increasing order
     */
    void arrangeSlots(const std::vector<size_t> &active)
    {
        arrays.template arrange<real, double>(active,
                                              { { &fixed.arrayForces, arrays.stride },
                                                { &shells.dx, shells.stride } },
                                              { { &fixed.arrayForcesDouble, arraysDouble.stride },
                                                { &shells.step, shells.stride },
                                                { &shells.force, shells.stride },
                                                { &shells.dxForce, shells.stride },
                                                { &shells.reference, shells.stride },
                                                { &shells.stall, shells.stride } });
    }

    //! \return Whether arraysDouble is used, see doublePrecision
    bool useArraysDouble() const { return doublePrecision && !std::is_same_v<real, double>; }
    //! Energies
    ConformationEnergies energies;
    //! Whether the virtual sites in the conformations were constructed
    bool                 haveVsites   = false;
    //! The Topology::vsiteVersion the virtual sites were constructed for
    size_t               vsiteVersion = 0;
    /*! Whether the next energy calculation has to construct the virtual
     * sites of the conformations first.
     */
    bool                 constructVsites   = false;
    /*! Whether the energy calculations spread the forces on the virtual
     * sites over their parents, in the arrays if these are used. This is
     * done in each calculation, which does not change the forces on the
     * shells, but saves spreading them afterwards.
     */
    bool                 spreadVsites      = false;
    //! Whether the forces of all conformations were spread, see spreadVsites
    bool                 vsiteForcesSpread = false;
    //! Changes when the coordinates of the conformations change
    StructureVersion     coordinates;
    /*! Changes when the conformations are made anew, but not when only
     * the virtual sites move, unlike coordinates
     */
    StructureVersion     conformations;
    //! Energies of each interaction type, to use again if nothing changed
    std::map<InteractionType, CachedEnergies> cache;
    //! Reference data for computing deviations
    DeviationReference   deviation;
    //! Energies and forces of the entries without shells
    FixedEnergies        fixed;
    //! State of the shell relaxation, in the order of the slots
    ShellRelaxation      shells;
    //! The conformations computed in the last energy calculation, kept to use the memory again
    std::vector<size_t>  active;
    /*! The number of conformations if all that are not in active are done,
     * such that the next calculation leaving out those that are done only
     * has to visit the conformations in active, see CoordsEnergy::done.
     * Zero if that is not known.
     */
    size_t               activeKnown = 0;

    //! Mark the coordinates and virtual sites as outdated
    void invalidate()
    {
        fixed.valid = false;
        activeKnown = 0;
        haveVsites = false;
        constructVsites   = false;
        vsiteForcesSpread = false;
        coordinates.update();
        conformations.update();
    }
};

} // namespace alexandria

#endif

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
 * The energies of each conformation are reported in the energy map of its
 * CoordsEnergy, but during a calculation they are kept here, since visiting
 * the energy maps of thousands of conformations repeatedly is limited by
 * memory access. The coordinates are kept between calculations as well,
 * since only those of virtual sites and shells change.
 * \author David van der Spoel <david.vanderspoel@icm.uu.se>
 */
#ifndef ACT_FORCES_CONFORMATIONDATA_H
#define ACT_FORCES_CONFORMATIONDATA_H

#include <algorithm>
#include <array>
#include <cstddef>
#include <map>
#include <vector>

#include "act/basics/act_particle.h"
#include "act/basics/interactiontype.h"
#include "act/forcefield/forcefield_parameterlist.h"
#include "act/forcefield/potential.h"
#include "act/forces/forces_explicit_simd.h"
#include "act/topology/actatom.h"
#include "act/topology/coords_energy.h"
#include "gromacs/math/vectypes.h"
#include "gromacs/utility/gmxassert.h"

namespace alexandria
{

/*! \brief Coordinates and forces of the conformations of a topology.
 *
 * Component m of particle a in slot c is stored at (3*a+m)*stride + c,
 * such that loops over conformations access memory contiguously and can be
 * vectorized. The slots hold all conformations, or only those for which the
 * shells did not converge yet, see CoordsEnergy::done, and conformation(c)
 * is the conformation in slot c. The stride is the number of slots rounded
 * up to a multiple of the explicit SIMD width. The extra slots are copies of
 * the last one, they are computed but not used.
 */
class ConformationArrays
{
public:
    //! Number of slots, that is conformations in the arrays
    size_t              nconf  = 0;
    //! Number of slots including the padding
    size_t              stride = 0;
    //! Coordinates
    std::vector<double> x;
    //! Forces
    std::vector<double> f;
    //! Energies of one potential, kept here to reuse the memory
    std::vector<double> energies;

    /*! \brief Mark the coordinates as outdated
     * Has to be called when the coordinates of the atoms in the
     * conformations change, the next call of setCoordinates then
     * copies all of them.
     */
    void invalidate() { valid_ = false; }

    /*! \brief Mark the coordinates of the virtual sites as outdated
     * Has to be called when the positions of the virtual sites in the
     * conformations change, the next call of setCoordinates then
     * copies them.
     */
    void invalidateVsites() { vsitesValid_ = false; }

    /*! \brief Copy the coordinates of the conformations and clear the forces
     * Only the coordinates of shells, and of virtual sites if they changed,
     * are copied if those of all particles of the same conformations were
     * copied before, since the atoms do not move.
     * \param[in] coords   The conformations, all with the same number of particles
     * \param[in] atoms    The particles
     * \param[in] forces   Whether forces will be computed, if not they are not cleared
     * \param[in] skipDone Whether to leave out the conformations that are done
     * \return false if the conformations do not have the same size or
     *         there is no conformation left
     */
    bool setCoordinates(const std::vector<CoordsEnergy> &coords,
                        const std::vector<ActAtom>      &atoms,
                        bool                             forces   = true,
                        bool                             skipDone = false)
    {
        size_t natom = coords.empty() ? 0 : coords[0].x.size();
        for(const auto &ce : coords)
        {
            if (ce.x.size() != natom)
            {
                valid_ = false;
                return false;
            }
        }
        // The same slots as before unless conformations are left out
        bool sameSlots = (valid_ && confs_.size() == coords.size() && !skipDone);
        if (!sameSlots)
        {
            std::vector<size_t> confs;
            for(size_t c = 0; c < coords.size(); c++)
            {
                if (!skipDone || !coords[c].done)
                {
                    confs.push_back(c);
                }
            }
            sameSlots = (valid_ && confs == confs_);
            confs_.swap(confs);
        }
        size_t nc = confs_.size();
        if (nc == 0)
        {
            valid_ = false;
            return false;
        }
        size_t width = explicitSimdWidth();
        size_t st    = ((nc + width - 1)/width)*width;
        bool   reuse = sameSlots && nc == nconf && st == stride && x.size() == DIM*natom*st;
        nconf  = nc;
        stride = st;
        if (forces)
        {
            f.assign(DIM*natom*stride, 0.0);
        }
        if (reuse)
        {
            for(auto a : shells_)
            {
                copyParticle(coords, a);
            }
            if (!vsitesValid_)
            {
                for(auto a : vsites_)
                {
                    copyParticle(coords, a);
                }
            }
        }
        else
        {
            x.resize(DIM*natom*stride);
            shells_.clear();
            vsites_.clear();
            for(size_t a = 0; a < natom; a++)
            {
                copyParticle(coords, a);
                if (atoms.size() != natom || atoms[a].pType() == ActParticle::Shell)
                {
                    shells_.push_back(a);
                }
                else if (atoms[a].pType() != ActParticle::Atom)
                {
                    vsites_.push_back(a);
                }
            }
            valid_ = true;
        }
        vsitesValid_ = true;
        return true;
    }

    /*! \brief Copy the forces of the conformations into this structure
     * \param[in] coords The conformations
     */
    void loadForces(const std::vector<CoordsEnergy> &coords)
    {
        size_t natom = f.size()/(DIM*stride);
        for(size_t a = 0; a < natom; a++)
        {
            for(int m = 0; m < DIM; m++)
            {
                for(size_t c = 0; c < nconf; c++)
                {
                    f[(DIM*a+m)*stride + c] = coords[confs_[c]].fcalc[a][m];
                }
            }
        }
    }

    /*! \brief Copy the forces back to the conformations
     * \param[inout] coords The conformations
     */
    void storeForces(std::vector<CoordsEnergy> *coords) const
    {
        size_t natom = f.size()/(DIM*stride);
        for(size_t a = 0; a < natom; a++)
        {
            for(int m = 0; m < DIM; m++)
            {
                for(size_t c = 0; c < nconf; c++)
                {
                    (*coords)[confs_[c]].fcalc[a][m] = f[(DIM*a+m)*stride + c];
                }
            }
        }
    }

    /*! \brief Copy the forces on the shells back to the conformations
     * The forces on the other particles are not changed.
     * \param[inout] coords The conformations
     */
    void storeShellForces(std::vector<CoordsEnergy> *coords) const
    {
        for(auto a : shells_)
        {
            for(int m = 0; m < DIM; m++)
            {
                for(size_t c = 0; c < nconf; c++)
                {
                    (*coords)[confs_[c]].fcalc[a][m] = f[(DIM*a+m)*stride + c];
                }
            }
        }
    }

    /*! \brief Set the forces from forces kept in the same layout
     * The forces may have been kept for more conformations than are in
     * the slots now, see conformations.
     * \param[in] src       The forces
     * \param[in] srcStride The stride of src
     * \param[in] srcConfs  The conformation in each slot of src, it has to
     *                      contain all conformations in the slots now
     */
    void loadForces(const std::vector<double> &src,
                    size_t                     srcStride,
                    const std::vector<size_t> &srcConfs)
    {
        size_t nrow = x.size()/stride;
        f.resize(nrow*stride);
        if (srcConfs == confs_ && srcStride == stride)
        {
            std::copy(src.begin(), src.begin() + nrow*stride, f.begin());
            return;
        }
        // The slot of each conformation in src
        std::vector<size_t> srcSlot(stride);
        size_t              k = 0;
        for(size_t c = 0; c < stride; c++)
        {
            size_t conf = confs_[std::min(c, nconf-1)];
            while (k < srcConfs.size() && srcConfs[k] < conf)
            {
                k++;
            }
            GMX_RELEASE_ASSERT(k < srcConfs.size() && srcConfs[k] == conf,
                               "Forces of a conformation were not kept");
            srcSlot[c] = k;
        }
        for(size_t row = 0; row < nrow; row++)
        {
            const double *s = src.data() + row*srcStride;
            double       *d = f.data() + row*stride;
            for(size_t c = 0; c < stride; c++)
            {
                d[c] = s[srcSlot[c]];
            }
        }
    }

    //! \return The conformations in the slots
    const std::vector<size_t> &conformations() const { return confs_; }

private:
    //! Whether x holds the coordinates of all particles
    bool                valid_       = false;
    //! Whether x holds the current coordinates of the virtual sites
    bool                vsitesValid_ = false;
    //! The conformation in each slot
    std::vector<size_t> confs_;
    //! The shells, or all particles if their types are not known, their coordinates change
    std::vector<size_t> shells_;
    //! The other particles that are not atoms, their coordinates change with their parameters
    std::vector<size_t> vsites_;

    /*! \brief Copy the coordinates of particle a of all slots, including the padding
     * \param[in] coords The conformations
     * \param[in] a      The particle index
     */
    void copyParticle(const std::vector<CoordsEnergy> &coords, size_t a)
    {
        for(size_t c = 0; c < stride; c++)
        {
            const auto &xc = coords[confs_[std::min(c, nconf-1)]].x[a];
            for(int m = 0; m < DIM; m++)
            {
                x[(DIM*a+m)*stride + c] = xc[m];
            }
        }
    }
};

/*! \brief The ACT energies of all conformations, per interaction type
 * This holds the same information as the ACT energies in the energy maps
 * of the conformations: for each interaction type and conformation whether
 * there is an energy, the energy and its backup. The energies of one
 * interaction type are stored contiguously, and are copied to the energy
 * maps by store().
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

    /*! \brief Copy the energies to the energy maps of the conformations
     * \param[inout] coords The conformations
     */
    void store(std::vector<CoordsEnergy> *coords) const
    {
        for(size_t c = 0; c < coords->size() && c < nconf_; c++)
        {
            auto &eMap = (*coords)[c].eMap;
            for(auto t : types_)
            {
                const auto &d = data_[index_[static_cast<size_t>(t)]];
                if (d.have[c])
                {
                    auto &e = eMap[t];
                    e.setACT(backups_ ? d.before[c] : d.act[c]);
                    e.backup();
                    e.setACT(d.act[c]);
                }
            }
        }
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
    //! Whether the forces are in arrayForces rather than forces
    bool                                inArrays = false;
    //! The forces in the layout of ConformationArrays
    std::vector<double>                 arrayForces;
    //! The stride of arrayForces
    size_t                              stride = 0;
    //! The conformation in each slot of arrayForces
    std::vector<size_t>                 slots;
    //! The interaction types, including EPOT
    std::vector<InteractionType>        types;
    //! For each type the energy of each conformation
    std::vector<std::vector<double>>    energies;
    //! For each type whether each conformation has an energy
    std::vector<std::vector<char>>      have;
    //! The forces of each conformation, unless inArrays
    std::vector<std::vector<gmx::RVec>> forces;
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

//! Coordinates, forces and energies of all conformations of a topology
struct ConformationData
{
    //! Coordinates and forces
    ConformationArrays   arrays;
    //! Energies
    ConformationEnergies energies;
    //! Whether the virtual sites in the conformations were constructed
    bool                 haveVsites   = false;
    //! The Topology::vsiteVersion the virtual sites were constructed for
    size_t               vsiteVersion = 0;
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

    //! Mark the coordinates and virtual sites as outdated
    void invalidate()
    {
        arrays.invalidate();
        fixed.valid = false;
        haveVsites = false;
        coordinates.update();
        conformations.update();
    }
};

} // namespace alexandria

#endif

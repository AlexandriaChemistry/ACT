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
#include <vector>

#include "act/basics/act_particle.h"
#include "act/basics/interactiontype.h"
#include "act/forces/forces_explicit_simd.h"
#include "act/topology/actatom.h"
#include "act/topology/coords_energy.h"
#include "gromacs/math/vectypes.h"

namespace alexandria
{

/*! \brief Coordinates and forces of all conformations of a topology.
 *
 * Component m of particle a in conformation c is stored at (3*a+m)*stride + c,
 * such that loops over conformations access memory contiguously and can be
 * vectorized. The stride is the number of conformations rounded up to a
 * multiple of the explicit SIMD width. The extra conformations are copies of
 * the last one, they are computed but not used.
 */
class ConformationArrays
{
public:
    //! Number of conformations
    size_t              nconf  = 0;
    //! Number of conformations including the padding
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

    /*! \brief Copy the coordinates of all conformations and clear the forces
     * Only the coordinates of shells, and of virtual sites if they changed,
     * are copied if those of all particles were copied before, since the
     * atoms do not move.
     * \param[in] coords The conformations, all with the same number of particles
     * \param[in] atoms  The particles
     * \param[in] forces Whether forces will be computed, if not they are not cleared
     * \return false if the conformations do not have the same size
     */
    bool setCoordinates(const std::vector<CoordsEnergy> &coords,
                        const std::vector<ActAtom>      &atoms,
                        bool                             forces = true)
    {
        size_t nc = coords.size();
        if (nc == 0)
        {
            valid_ = false;
            return false;
        }
        size_t natom = coords[0].x.size();
        for(const auto &ce : coords)
        {
            if (ce.x.size() != natom)
            {
                valid_ = false;
                return false;
            }
        }
        size_t width = explicitSimdWidth();
        size_t st    = ((nc + width - 1)/width)*width;
        bool   reuse = valid_ && nc == nconf && st == stride && x.size() == DIM*natom*st;
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

    //! Copy the forces of all conformations into this structure
    void loadForces(const std::vector<CoordsEnergy> &coords)
    {
        size_t natom = f.size()/(DIM*stride);
        for(size_t a = 0; a < natom; a++)
        {
            for(int m = 0; m < DIM; m++)
            {
                for(size_t c = 0; c < nconf; c++)
                {
                    f[(DIM*a+m)*stride + c] = coords[c].fcalc[a][m];
                }
            }
        }
    }

    //! Copy the forces back to the conformations
    void storeForces(std::vector<CoordsEnergy> *coords) const
    {
        size_t natom = f.size()/(DIM*stride);
        for(size_t a = 0; a < natom; a++)
        {
            for(int m = 0; m < DIM; m++)
            {
                for(size_t c = 0; c < nconf; c++)
                {
                    (*coords)[c].fcalc[a][m] = f[(DIM*a+m)*stride + c];
                }
            }
        }
    }

private:
    //! Whether x holds the coordinates of all particles
    bool                valid_       = false;
    //! Whether x holds the current coordinates of the virtual sites
    bool                vsitesValid_ = false;
    //! The shells, or all particles if their types are not known, their coordinates change
    std::vector<size_t> shells_;
    //! The other particles that are not atoms, their coordinates change with their parameters
    std::vector<size_t> vsites_;

    //! Copy the coordinates of particle a of all conformations, including the padding
    void copyParticle(const std::vector<CoordsEnergy> &coords, size_t a)
    {
        for(size_t c = 0; c < stride; c++)
        {
            const auto &xc = coords[std::min(c, nconf-1)].x[a];
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

        //! \return Whether there is an energy for conformation c
        bool has(size_t c) const { return have && have[c]; }
        //! Set the energy of conformation c, a new energy has no backup
        void set(size_t c, double value)
        {
            if (!have[c])
            {
                before[c] = 0;
                have[c]   = 1;
            }
            act[c] = value;
        }
        //! Restore the backup of conformation c, as ACTEnergy::restoreBackup
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

    //! \return The data for interaction type t, added if not present
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

    //! \return The data for interaction type t, with null pointers if not present
    Term find(InteractionType t)
    {
        auto i = index_[static_cast<size_t>(t)];
        return i < 0 ? Term() : term(i);
    }

    //! \return The interaction types present, in increasing order like in an energy map
    const std::vector<InteractionType> &types() const { return types_; }

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
                    e.setACT(d.before[c]);
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
    //! Position of each interaction type in data_, -1 if not present
    std::array<int, ntype_>     index_;
    //! The interaction types present, sorted
    std::vector<InteractionType> types_;
    //! The data, the vectors keep their memory when data_ grows
    std::vector<Data>           data_;

    //! \return Pointers to the data at position i
    Term term(int i)
    {
        auto &d = data_[i];
        return { d.act.data(), d.before.data(), d.have.data() };
    }
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

    //! Mark the coordinates and virtual sites as outdated
    void invalidate()
    {
        arrays.invalidate();
        haveVsites = false;
    }
};

} // namespace alexandria

#endif

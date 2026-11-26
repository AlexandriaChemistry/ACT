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
 * Implements part of the alexandria program.
 * \author David van der Spoel <david.vanderspoel@icm.uu.se>
 */

#ifndef ACT_FRAGMENTWALL_H
#define ACT_FRAGMENTWALL_H

#include <cstddef>
#include <vector>

#include "gromacs/math/vectypes.h"

namespace alexandria
{

/*! \brief The distance between the compounds of a dimer or cluster, and a wall that keeps them together
 *
 * In vacuum, compounds that move apart do not come back. For two
 * compounds (fragments) with the distance R between their centers of mass
 * the wall is the potential
 * \f[ U = \frac{k}{2}(R - r_w)^2 \quad\mbox{for}\quad R > r_w \f]
 * and zero otherwise. For more than two compounds R is the distance of
 * the center of mass of each compound to that of all compounds, and U is
 * the sum over the compounds. The force is divided over the atoms of a
 * compound in proportion to their masses, such that it does not change
 * the structure or the orientation of the compound. The forces add up to
 * zero.
 */
class FragmentWall
{
public:
    /*! \brief Constructor
     * \param[in] atomStart     The first particle of each fragment, and one
     *                          more than the last particle of the last fragment
     * \param[in] mass          The mass of each particle
     * \param[in] radius        The radius of the wall in nm
     * \param[in] forceConstant The force constant in kJ/mol nm^2. There is
     *                          no wall if the radius or the force constant
     *                          is zero, or with one fragment.
     * \throws gmx::InvalidInputError for fragments without mass or without
     * particles, or a negative radius or force constant
     */
    FragmentWall(const std::vector<size_t> &atomStart,
                 const std::vector<double> &mass,
                 double                     radius,
                 double                     forceConstant);

    //! \return the number of fragments
    size_t nfragments() const { return fragmentMass_.size(); }

    //! \return whether there is a wall
    bool active() const { return nfragments() > 1 && radius_ > 0 && forceConstant_ > 0; }

    //! \return the radius of the wall in nm
    double radius() const { return radius_; }

    /*! \brief The distance between the fragments
     * \param[in] x The coordinates of all particles in nm
     * \return the distance between the centers of mass of two fragments,
     *         or the largest distance of the center of mass of a fragment
     *         to that of all of them, in nm. Zero for one fragment.
     */
    double distance(const std::vector<gmx::DVec> &x) const;

    /*! \brief The energy and the forces of the wall
     * \param[in]  x The coordinates of all particles in nm
     * \param[out] f The forces on all particles in kJ/mol nm
     * \return the energy in kJ/mol
     */
    double compute(const std::vector<gmx::DVec> &x,
                   std::vector<gmx::DVec>       *f) const;

private:
    //! The first particle of each fragment and one more than the last one
    std::vector<size_t> atomStart_;
    //! The mass of each particle
    std::vector<double> mass_;
    //! The mass of each fragment
    std::vector<double> fragmentMass_;
    //! The radius of the wall
    double              radius_        = 0;
    //! The force constant of the wall
    double              forceConstant_ = 0;

    /*! \brief The centers of mass of the fragments
     * \param[in] x The coordinates of all particles
     * \return the centers of mass
     */
    std::vector<gmx::DVec> centers(const std::vector<gmx::DVec> &x) const;
};

} // namespace alexandria

#endif

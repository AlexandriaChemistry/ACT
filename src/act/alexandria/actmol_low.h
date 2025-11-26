/*
 * This source file is part of the Alexandria Chemistry Toolkit.
 *
 * Copyright (C) 2014-2025
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
 * \author Mohammad Mehdi Ghahremanpour <mohammad.ghahremanpour@icm.uu.se>
 * \author David van der Spoel <david.vanderspoel@icm.uu.se>
 */

#ifndef ACTMOL_LOW_H
#define ACTMOL_LOW_H

#include <cassert>
#include <cstdio>
#include <cstring>
#include <vector>

#include "gromacs/pbcutil/pbc.h"

#include "act/basics/atomization_energy.h"
#include "act/basics/chargemodel.h"
#include "act/basics/identifier.h"

namespace alexandria
{

class ActAtom;
class ForceField;
class QgenAcm;

/*! \brief
 * Class to determine how to treat missing parameters.
 */
enum class missingParameters 
{ 
    //! Will generate an error if a parameter is missing
    Error,
    //! Missing parameters will be ignored
    Ignore, 
    //! Missing parameters will be generated
    Generate 
}; 

//! \return true if the four atoms are in a plane
bool is_planar(const rvec   xi,  const rvec xj,
               const rvec   xk,  const rvec xl,
               const t_pbc *pbc, real phi_toler);

//! \return true if the three atoms are on a line
bool is_linear(const rvec xi, const rvec xj,
               const rvec xk, const t_pbc *pbc,
               real th_toler);

/*! \brief relative position of an atom on a bond
 * This needs at least three atoms and two bond orders or it will throw
 * \param[in] pd  Force Field to look up bond lengths
 * \param[in] atoms Should be exactly three atoms
 * \param[in] bondOrders Should be exactly two
 * \return relative position
 */
real calc_relposition(const ForceField               *pd,
                      const std::vector<std::string> &atoms,
                      const std::vector<double>      &bondOrders);

//! \brief Convert string to an array of doubles
std::vector<double> getDoubles(const std::string &s);

//! \brief Make sure atoms are within a box
void put_in_box(int natom, matrix box, rvec x[], real dbox);

/*! \brief Compute a rotation matrix
 * \param[in]  target_vec The vector to rotate to
 * \param[in]  ref_vec    The vector to rotate from
 * \param[out] rotmatrix  The matrix that accomplishes the rotation
 */
void calc_rotmatrix(rvec target_vec, rvec ref_vec, matrix rotmatrix);

/*! \brief Compute atomization energy for a molecule
 * \param[in] atoms The atom elements etc.
 * \param[in] atomenergy Table of components of the atomization energy
 * \param[in] temperature The temperature at which to compute the energy
 * \return The total atomization energy
 */
double computeAtomizationEnergy(const std::vector<ActAtom> &atoms,
                                const AtomizationEnergy    &atomenergy,
                                double                      temperature);

} // namespace alexandria
#endif

/*
 * This source file is part of the Alexandria Chemistry Toolkit.
 *
 * Copyright (C) 2026
 *
 * Developers:
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
 * Computes the rotational symmetry number of a structure.
 * \author David van der Spoel <david.vanderspoel@icm.uu.se>
 */
#ifndef ACT_PROPERTIES_SYMMETRYNUMBER_H
#define ACT_PROPERTIES_SYMMETRYNUMBER_H

#include <cstddef>
#include <vector>

#include "gromacs/math/vectypes.h"

namespace alexandria
{

/*! \brief Compute the rotational symmetry number of a structure
 *
 * The rotational symmetry number is the number of proper rotations
 * (including the identity) that map the structure onto itself, where
 * atoms can only be mapped onto atoms of the same element. It is the
 * order of the rotational subgroup of the point group, e.g. 2 for water,
 * 3 for ammonia, 12 for methane and 2 for carbon dioxide.
 *
 * The rotations are found by enumerating permutations of atoms of the
 * same element that preserve all interatomic distances, after which
 * only permutations that can be achieved by a proper rotation, rather
 * than a reflection, are counted.
 *
 * \param[in] atomicNumbers The atomic number of each atom
 * \param[in] masses        The mass of each atom, for the center of mass
 * \param[in] coords        The coordinates (nm)
 * \param[in] tolerance     Tolerance for distances and the fit (nm)
 * \param[in] maxNodes      Maximum number of steps in the search
 * \return The symmetry number, or 0 if the search was not completed
 *         within maxNodes steps
 */
int rotationalSymmetryNumber(const std::vector<int>       &atomicNumbers,
                             const std::vector<double>    &masses,
                             const std::vector<gmx::RVec> &coords,
                             double                        tolerance = 0.01,
                             size_t                        maxNodes  = 1000000);

} // namespace alexandria

#endif

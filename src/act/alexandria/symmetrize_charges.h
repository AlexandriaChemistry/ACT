/*
 * This source file is part of the Alexandria Chemistry Toolkit.
 *
 * Copyright (C) 2014-2023
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
 
 
#ifndef ACT_SYMMETRIZE_CHARGES_H
#define ACT_SYMMETRIZE_CHARGES_H

#include <vector>

#include "act/forcefield/forcefield.h"
#include "topology.h"

struct gmx_atomprop;

namespace alexandria
{

/*! Generate a list of symmetry-related charges.
 * \param[in]  topology    The molecular topology
 * \param[in]  pd          The force field, containing symmetrization info
 * \param[in]  symm_string Optional (may be nullptr) user provided list
 * \param[out] sym_charges The final list
 */
void get_symmetrized_charges(Topology         *topology,
                             const ForceField *pd,
                             const char       *symm_string,
                             std::vector<int> *sym_charges);

/*! Symmetrize the charges in an array.
 * \param[in]  q           array of charges
 * \param[out] sym_charges The list of charges identities
 * \throws if arrays are not equally long.
 */
void apply_symmetrized_charges(std::vector<double>    *q,
                               const std::vector<int> &sym_charges);

} // namespace alexandria

#endif

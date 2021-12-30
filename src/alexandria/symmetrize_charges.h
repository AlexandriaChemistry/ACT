/*
 * This source file is part of the Alexandria Chemistry Toolkit.
 *
 * Copyright (C) 2014-2020 
 *
 * Developers:
 *             Mohammad Mehdi Ghahremanpour, 
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

#include "gromacs/topology/atoms.h"

#include "poldata.h"
#include "topology.h"

struct gmx_atomprop;

namespace alexandria
{

/*! Generate a list of symmetry-related charges
 * \param[in]  bQsym       If false, make an identity list
 * \param[in]  atoms       The atoms in the molecules
 * \param[in]  pd          The force field, containing symmetrization info
 * \param[in]  symm_string Optional (may be nullptr) user provided list
 * \param[out] sym_charges The final list
 */
void symmetrize_charges(bool                                bQsym,
                        const t_atoms                      *atoms,
                        const std::vector<TopologyEntry *> &bonds,
                        const Poldata                      *pd,
                        const char                         *symm_string,
                        std::vector<int>                   *sym_charges);

} // namespace alexandria

#endif

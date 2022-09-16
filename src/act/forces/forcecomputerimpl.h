/*
 * This source file is part of the Alexandria Chemistry Toolkit.
 *
 * Copyright (C) 2014-2022
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
#include <vector>

#include "act/basics/interactiontype.h"
#include "alexandria/topology.h"
#include "gromacs/math/vectypes.h"

namespace alexandria
{

typedef void (*bondForceComputer)(const std::vector<TopologyEntry *> &bonds,
                                  const std::vector<ActAtom>         &atoms,
                                  const std::vector<gmx::RVec>       *coordinates,
                                  std::vector<gmx::RVec>             *forces,
                                  std::map<InteractionType, double>  *energies);

/*! \brief Return a bonded force computer according to typedef.
 * \param[in] gromacs_index Number corresponding to GROMACS list of energy terms
 * \return pointer to the appropriate function
 * \throws with gmx::InternalError if no such function is implemented.
 */
bondForceComputer getBondForceComputer(int gromacs_index);

} // namespace


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

#ifndef ACT_MINIMIZE_H
#define ACT_MINIMIZE_H

#include <vector>

#include "act/alexandria/molhandler.h"
#include "gromacs/math/vectypes.h"

namespace alexandria
{

class ACTMol;
class ForceComputer;
class ForceField;
class JsonTree;
class MsgHandler;
class SimulationConfigHandler;

/*! \brief Minimize the energy of a molecule and report the results
 *
 * The energies and multipoles are written before and after the
 * minimization, and for a dimer the interaction energies after the
 * minimization. This is used by alexandria minimize and by
 * alexandria simulate -minimize.
 * \param[in]    msghandler For output and errors
 * \param[in]    pd         The force field
 * \param[in]    forceComp  The force computer
 * \param[in]    sch        The minimization settings
 * \param[in]    freezeFile Index file with atoms to freeze, may be nullptr
 * \param[inout] actmol     The molecule
 * \param[inout] xmin       The starting coordinates, on return the minimized ones
 * \param[inout] jtree      The energies before and after are added here
 * \return the status of the minimization
 */
eMinimizeStatus minimizeMolecule(MsgHandler                    *msghandler,
                                 const ForceField              *pd,
                                 const ForceComputer           *forceComp,
                                 const SimulationConfigHandler &sch,
                                 const char                    *freezeFile,
                                 ACTMol                        *actmol,
                                 std::vector<gmx::RVec>        *xmin,
                                 JsonTree                      *jtree);

} // namespace alexandria

#endif

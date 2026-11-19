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
/*! \internal \file
 * \brief
 * Declares functions that describe a molecule in a JsonTree, for programs
 * that read the output of the alexandria tools, such as the Python
 * interface and the graphical user interface.
 *
 * \author David van der Spoel <david.vanderspoel@icm.uu.se>
 */
#ifndef ACT_TOPOLOGY_ACTMOL_JSON_H
#define ACT_TOPOLOGY_ACTMOL_JSON_H

#include <string>
#include <vector>

#include "act/utility/jsontree.h"
#include "gromacs/math/vectypes.h"

namespace alexandria
{

class ACTMol;

/*! \brief Describe a molecule with its particles and bonds
 *
 * The tree has the name, the formula, the total charge and the mass of the
 * molecule, an array "particles" and an array "bonds".
 * Each particle has its name, element, atomic number, force field type,
 * particle type (Atom, Shell, VSite, ...), mass and charge. Particles that
 * belong to atoms, that is shells and virtual sites, have an array
 * "parents" with the indices of these atoms. The second charge of a
 * particle, if it has one, is given as charge_vs.
 * Each bond has an array "atoms" with two indices and the bond order.
 * All indices count from zero in the array of particles, which has the
 * same order as the coordinates of the molecule.
 *
 * \param[in] mol The molecule, of which the topology must have been made
 * \param[in] key The key of the tree
 * \return the tree
 */
JsonTree moleculeToJson(const ACTMol *mol, const std::string &key = "molecule");

/*! \brief Make an array with the coordinates of all particles
 * \param[in] key    The key of the array
 * \param[in] coords The coordinates in nm
 * \return an array with for each particle an array of x, y and z in nm
 */
JsonTree coordinatesToJson(const std::string            &key,
                           const std::vector<gmx::RVec> &coords);

} // namespace alexandria

#endif

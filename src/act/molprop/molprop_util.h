/*
 * This source file is part of the Alexandria Chemistry Toolkit.
 *
 * Copyright (C) 2014-2026
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

#ifndef MOLPROP_UTIL_H
#define MOLPROP_UTIL_H

#include <string>
#include <vector>

#include "act/molprop/molprop.h"
#include "gromacs/utility/arrayref.h"

namespace alexandria
{
/*! \brief Reset the molprop index number
 * \param[in] mp The molprops
 */
void generate_index(std::vector<MolProp> *mp);

/*! \brief Will merge multiple molprop files.
 * \param[in] msg_handler For info and debugging
 * \param[in] infiles     MolProp filenames
 * \param[in] mp          The molprop vector
 */
void merge_xml(MsgHandler                             *msg_handler,
               const gmx::ArrayRef<const std::string> &infiles,
               std::vector<alexandria::MolProp>       *mp);

/*! \brief Check that a set of molecules is consistent, as edit_mp -validate
 * First each molecule is checked with MolProp::validate, which removes the
 * experiments with the atoms in another order than the reference
 * experiment. Then all molecules are checked together:
 * - each compound, a fragment identified by its InChI, has its atoms in the
 *   same order (names and types) in all molecules, the first one it occurs in
 *   is the reference, such that it can be treated as one compound, and
 * - the fragments of a dimer A#B are in the order of A and B, if their
 *   names are those of the dimer. Both A#B and B#A may occur.
 * Each problem is an error of msghandler.
 * \param[in]    msghandler For the problems
 * \param[inout] mps        The molecules
 * \return The number of problems
 */
size_t validateMolProps(MsgHandler *msghandler, std::vector<MolProp> *mps);

} // namespace alexandria

#endif

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

} // namespace alexandria

#endif

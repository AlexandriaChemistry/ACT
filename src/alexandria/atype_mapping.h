/*
 * This source file is part of the Alexandria Chemistry Toolkit.
 *
 * Copyright (C) 2021
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
 * \author David van der Spoel <david.vanderspoel@icm.uu.se>
 */

#ifndef ATYPE_MAPPING_H
#define ATYPE_MAPPING_H

#include <map>

#include "molprop/molprop.h"

/*! \brief
 * Read a mapping file
 * \param[in]  filenm The name of the file to read from. If empty a
 *                    library file will be read
 * \param[out] g2a    Map for the atom names going from Gaff to Alexandria
 */
void gaffToAlexandria(const std::string                  &filenm,
                      std::map<std::string, std::string> *g2a);

/*! \brief
 * Rename atoms in molprop structure according to the table
 * \param[inout] mp  The molecule properties
 * \param[in]    g2a The mapping table
 * \return true if successful, false otherwise.
 */
bool renameAtomTypes(alexandria::MolProp                      *mp,
                     const std::map<std::string, std::string> &g2a);

#endif

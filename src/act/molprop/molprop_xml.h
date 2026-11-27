/*
 * This source file is part of the Alexandria Chemistry Toolkit.
 *
 * Copyright (C) 2014-2020 
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
 
 
#ifndef MOLPROP_XML_H
#define MOLPROP_XML_H

#include <vector>

#include "act/molprop/molprop.h"
#include "gromacs/utility/basedefinitions.h"

namespace alexandria
{

    class MsgHandler;
    
/*! \brief
 * Write a vector of molprops to an XML file
 *
 * Function that uses the libxml2 library to store an array of molprops
 * to file. The file is written one piece at a time, an experiment of a
 * compound or the rest of it, such that there is no XML tree of the whole
 * file in memory, which would take many times the size of the file.
 * \param[in] fn         The file name to write to
 * \param[in] mpt        The vector of MolProp
 * \param[in] bCompress  Determines whether zlib compression is used when writing
 * \ingroup group_molprop
 */
void MolPropWrite(const char                 *fn,
                  const std::vector<MolProp> &mpt,
                  gmx_bool                    bCompress);

/*! \brief
 * Reads a vector of molprops from an XML file
 *
 * Function that uses the libxml2 library to read an array of molprops
 * from file. The file is read as a stream, one experiment of a compound
 * or other part of it at a time, such that there is no XML tree of the
 * whole file in memory, which would take many times the size of the file.
 * The compounds are added to those that are there. If the file cannot be
 * read to the end, none of its compounds are added.
 * \param[in]    msg_handler For info and debugging
 * \param[in]    fn          The file name to read from
 * \param[inout] mpt         The vector of MolProp
 * \throws gmx::FileIOError if the file is not there or is not XML
 * \ingroup group_molprop
 */
void MolPropRead(MsgHandler           *msg_handler,
                 const char           *fn,
                 std::vector<MolProp> *mpt);

} // namespace alexandria
#endif

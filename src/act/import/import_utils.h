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

#ifndef ACT_IMPORT_IMPORT_UTILS_H
#define ACT_IMPORT_IMPORT_UTILS_H

#include <map>
#include <string>
#include <vector>

#include "act/basics/allmols.h"
#include "act/molprop/fragment.h"
#include "act/molprop/topologyentry.h"

namespace alexandria
{
class MsgHandler;

typedef struct
{
    //! The atom indices starting from zero indicating the bond
    int    ai, aj;
    //! The bond order
    double order;
} SimpleBond;

typedef struct
{
    //! The name of this group
    std::string              name;
    //! The Smarts key for this entry
    std::string              smarts;
    //! Total charge
    int                      charge;
    //! Total multiplicity
    int                      multiplicity;
    //! Strings corresponding to atom types
    std::vector<std::string> atomtypes;
    //! SimpleBond structures, containing atom numbers and bondorders
    std::vector<SimpleBond>  bonds;
} AtomBondtypeEntry;

/*! \brief Fetch database of special cases.
 * \return a vector of AtomBondtype structures.
 */
std::vector<AtomBondtypeEntry> getAtomBondtypeDB();

/*! \brief Fetch database of special cases.
 * \param[in]  msghandler For warnings and errors
 * \param[in]  dbname     Name of the database to read. If empty the default atom_bond.xml will be used.
 * \param[out] abdb       The atom and bond type database
 */
void readAtomBondtypeDB(MsgHandler                     *msghandler,
                        const std::string              &dbname,
                        std::vector<AtomBondtypeEntry> *abdb);
/*! \brief Write the database to a json file
 * \param[in] msghandler For warnings and errors
 * \param[in] filenm     The filename to write to
 * \param[in] db         The AtomBondtype database
 */
void writeAtomBondtypeDB(MsgHandler                           *msghandler,
                         const std::string                    &filenm,
                         const std::vector<AtomBondtypeEntry> &db);
} // namespace

#endif

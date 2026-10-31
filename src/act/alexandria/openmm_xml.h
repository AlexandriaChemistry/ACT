/*
 * This source file is part of the Alexandria Chemistry Toolkit.
 *
 * Copyright (C) 2021-2025
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
 * \author Marie-Madeleine Walz <marie-madeleine.walz@icm.uu.se>
 * \author David van der Spoel <david.vanderspoel@icm.uu.se>
 */
 
#ifndef OPENMM_XML_H
#define OPENMM_XML_H

#include "act/forcefield/forcefield.h"
#include "act/topology/actmol.h"

namespace alexandria
{
    class MsgHandler;

    /*! \brief Determine the residue names of the compounds for OpenMM
     *
     * The same names should be used in the force field file and the PDB
     * file, such that OpenMM can find the bonds from the residue and atom
     * names. A name has at most four letters or digits and cannot be a
     * residue that OpenMM knows itself, such as HOH. Compounds that do not
     * get a name from the user get one from the first letters of their
     * IUPAC name, or a number if that does not work.
     *
     * \param[in] msghandler For errors in the user input
     * \param[in] actmols    The ACT molecule structures
     * \param[in] userNames  Names from the user separated by spaces or
     *                       commas. An entry compound=NAME is for the
     *                       compound with that name or InChi, an entry
     *                       NAME is for the next compound without a name,
     *                       in the order in which they appear in actmols.
     * \return A map from the InChi of each compound to its residue name
     */
    std::map<std::string, std::string> openMMResidueNames(MsgHandler                *msghandler,
                                                          const std::vector<ACTMol> &actmols,
                                                          const std::string         &userNames);

    /*! \brief Store the ForceField force field to an OpenMM XML file
     *
     * \param[in] msghandler For status and messages
     * \param[in] fileName   The filename to save the force field to
     * \param[in] pd         Pointer to a ForceField class instance
     * \param[in] actmols    The ACT molecule structures
     * \param[in] mDrude     Mass to use for the drude particle if any
     * \param[in] ntrain     Minimum number of training data to store a parameter to OpenMM
     * \param[in] addNumbersToAtoms Will add integer index to all atom types
     * \param[in] residueNames From the InChi of each compound to its
     *                         residue name, see openMMResidueNames. If
     *                         empty the names are generated.
     */
    void writeOpenMM(MsgHandler                               *msghandler,
                     const std::string                        &fileName,
                     const ForceField                         *pd,
                     const std::vector<ACTMol>                &actmols,
                     double                                    mDrude,
                     int                                       ntrain,
                     bool                                      addNumbersToAtoms = true,
                     const std::map<std::string, std::string> &residueNames = {});

} // namespace alexandria

#endif

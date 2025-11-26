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

#include "act/alexandria/actmol.h"
#include "act/forcefield/forcefield.h"

namespace alexandria
{
    class MsgHandler;
    
    /*! \brief Store the ForceField force field to an OpenMM XML file
     *
     * \param[in] msghandler For status and messages
     * \param[in] fileName   The filename to save the force field to
     * \param[in] pd         Pointer to a ForceField class instance
     * \param[in] actmols    The ACT molecule structures
     * \param[in] mDrude     Mass to use for the drude particle if any
     * \param[in] ntrain     Minimum number of training data to store a parameter to OpenMM
     * \param[in] addNumbersToAtoms Will add integer index to all atom types
     */
    void writeOpenMM(MsgHandler                *msghandler,
                     const std::string         &fileName,
                     const ForceField          *pd,
                     const std::vector<ACTMol> &actmols,
                     double                     mDrude,
                     int                        ntrain,
                     bool                       addNumbersToAtoms = true);

} // namespace alexandria

#endif

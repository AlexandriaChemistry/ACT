/*
 * This source file is part of the Alexandria Chemistry Toolkit.
 *
 * Copyright (C) 2021 
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
    /*! \brief Store the ForceField force field to an OpenMM XML file
     *
     * \param[in] fileName  The filename to save the force field to
     * \param[in] simParams The filename to save the simulation params to
     * \param[in] pd        Pointer to a ForceField class instance
     * \param[in] actmols   The ACT molecule structures
     * \param[in] mDrude    Mass to use for the drude particle if any
     * \param[in] compress  Whether or not to write a compressed file
     * \param[in] addNumbersToAtoms Will add integer index to all atom types
     */
    void writeOpenMM(const std::string         &fileName,
                     const std::string         &simParams,
                     const ForceField          *pd,
                     const std::vector<ACTMol> &actmols,
                     double                     mDrude,
                     bool                       addNumbersToAtoms = true);

} // namespace alexandria

#endif

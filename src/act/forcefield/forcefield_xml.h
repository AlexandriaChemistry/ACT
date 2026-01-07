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

#ifndef FORCEFIELD_XML_H
#define FORCEFIELD_XML_H

#include <string>

namespace alexandria
{
    class ForceField;

    /*! \brief Store the ForceField class to an XML file
     *
     * \param[in] fileName The filename to save to
     * \param[in] pd       Pointer to a ForceField class instance
     * \param[in] compress Whether or not to write a compressed file
     */
    void writeForceField(const std::string &fileName,
                      const ForceField     *pd,
                      bool compress = true);

    /*! \brief Read a ForceField class from an XML file
     *
     * \param[in]  fileName The filename to read from
     * \param[out] pd       The ForceField class instance
     */
    void readForceField(const std::string &fileName,
                     ForceField           *pd);

} // namespace alexandria

#endif

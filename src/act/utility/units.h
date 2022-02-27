/*
 * This source file is part of the Alexandria Chemistry Toolkit.
 *
 * Copyright (C) 2020 
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

#ifndef UNITS_H
#define UNITS_H

#include <string>

namespace alexandria
{
/*! Routine for unit conversion to GROMACS internal units
 * \param[in] x    The value to convert
 * \param[in] unit Normal string, e.g. pm or Angstrom
 * \return x*factor
 */
double convertToGromacs(double x, const std::string &unit);

/*! Routine to determine the corresponding GROMACS unit
 * after conversion.
 * \param[in] unit The external unit
 * \return the GROMACS unit
 * \throws if not found
 */
std::string gromacsUnit(const std::string &unit);

/*! Routine for unit conversion from GROMACS internal units
 * \param[in] x    The value to convert
 * \param[in] unit Normal string, e.g. pm or Angstrom
 * \return x/factor
 */
double convertFromGromacs(double x, const std::string &unit);

} // namespace gmx

#endif

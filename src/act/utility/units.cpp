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

#include "units.h"

#include <cmath>

#include <map>
#include <string>

#include "gromacs/math/units.h"
#include "gromacs/utility/exceptions.h"
#include "gromacs/utility/fatalerror.h"
#include "gromacs/utility/stringutil.h"

namespace alexandria
{

const static std::map<const std::string, double> unitConversion =
    {
        { "Angstrom",   A2NM    },
        { "nm",         1       },
        { "pm",         0.001   },
        { "Bohr",       BOHR2NM },
        { "Bohr3",      BOHR2NM*BOHR2NM*BOHR2NM },
        { "kcal/mol",   CAL2JOULE },
        { "kJ/mol",     1 },
        { "J/mol K",    1 },
        { "cal/mol K",  CAL2JOULE },
        { "Hartree",    ONE_4PI_EPS0/BOHR2NM },
        { "Hartree/e",  ONE_4PI_EPS0/BOHR2NM },
        { "Angstrom3",  A2NM*A2NM*A2NM },
        { "Coulomb",    1.0/E_CHARGE },
        { "Debye",      DEBYE2ENM },
        { "Electron",   1 },
        { "Buckingham", A2NM*DEBYE2ENM },
        { "DAngstrom2", A2NM*A2NM*DEBYE2ENM },
        { "DAngstrom3", A2NM*A2NM*A2NM*DEBYE2ENM },
        { "1/nm",       1 },
        { "degree",     M_PI/180.0 },
        { "kJ/mol/rad2",1 },
        { "kJ/mol/nm2", 1 },
        { "e",          1 },
        { "", 1 }
    };

double convertToGromacs(double x, const std::string &unit)
{
    auto uc = unitConversion.find(unit);
    if (uc  == unitConversion.end())
    {
#if CMAKE_BUILD_TYPE == CMAKE_BUILD_TYPE_DEBUG
        GMX_THROW(gmx::InternalError(gmx::formatString("Unknown unit '%s'\n", unit.c_str()).c_str()));
#else
        fprintf(stderr, "Unknown unit %s\n", unit.c_str());
#endif
        return 1;
    }
    return x*uc->second;
}

double convertFromGromacs(double x, const std::string &unit)
{
    auto uc = unitConversion.find(unit);
    if (uc  == unitConversion.end())
    {
#if CMAKE_BUILD_TYPE == CMAKE_BUILD_TYPE_DEBUG
        GMX_THROW(gmx::InternalError(gmx::formatString("Unknown unit %s\n", unit.c_str()).c_str()));
#else
        fprintf(stderr, "Unknown unit %s\n", unit.c_str());
#endif
        return 1;
    }
    return x/uc->second;
}

} // namespace alexandria

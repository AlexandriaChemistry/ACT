/*
 * This source file is part of the Alexandria Chemistry Toolkit.
 *
 * Copyright (C) 2020-2022
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

const static std::map<const std::string, std::pair<const std::string, double> > unitConversion =
    {
        { "Angstrom",    { "nm", A2NM    } },
        { "nm",          { "nm", 1       } },
        { "pm",          { "nm", 0.001   } },
        { "Bohr",        { "nm", BOHR2NM } },
        { "Bohr3",       { "nm3", BOHR2NM*BOHR2NM*BOHR2NM } },
        { "nm3",         { "nm3", 1 } },
        { "kcal/mol",    { "kJ/mol", CAL2JOULE } },
        { "kJ/mol",      { "kJ/mol", 1 } },
        { "kJ/mol e",    { "kJ/mol e", 1 } },
        { "J/mol K",     { "J/mol K", 1 } },
        { "cal/mol K",   { "J/mol K", CAL2JOULE } },
        { "Hartree",     { "kJ/mol", ONE_4PI_EPS0/BOHR2NM } },
        { "Hartree/e",   { "kJ/mol e", ONE_4PI_EPS0/BOHR2NM } },
        { "Hartree/Bohr",{ "kJ/mol nm", ONE_4PI_EPS0/(BOHR2NM*BOHR2NM) } },
        { "kJ/mol nm",   { "kJ/mol nm", 1 } },
        { "kJ/mol nm2",  { "kJ/mol nm2", 1 } },
        { "kJ/mol nm3",  { "kJ/mol nm3", 1 } },
        { "kJ/mol nm4",  { "kJ/mol nm4", 1 } },
        { "kJ/mol nm6",  { "kJ/mol nm6", 1 } },
        { "kJ/mol nm8",  { "kJ/mol nm8", 1 } },
        { "kJ/mol nm10", { "kJ/mol nm10", 1 } },
        { "Angstrom3",   { "nm3", A2NM*A2NM*A2NM } },
        { "A^3",         { "nm3", A2NM*A2NM*A2NM } },
        { "Coulomb",     { "e", 1.0/E_CHARGE } },
        { "e nm",        { "e nm", 1 } },
        { "Debye",       { "e nm", DEBYE2ENM } },
        { "D",           { "e nm", DEBYE2ENM } },
        { "Electron",    { "e", 1 } },
        { "Buckingham",  { "e nm2", A2NM*DEBYE2ENM } },
        { "B",           { "e nm2", A2NM*DEBYE2ENM } },
        { "e nm2",       { "e nm2", 1 } },
        { "e nm3",       { "e nm3", 1 } },
        { "DAngstrom2",  { "e nm2", A2NM*A2NM*DEBYE2ENM } },
        { "DAngstrom3",  { "e nm3", A2NM*A2NM*A2NM*DEBYE2ENM } },
        { "D.Angstrom",  { "e nm", A2NM*DEBYE2ENM } },
        { "D.Angstrom2", { "e nm2", A2NM*A2NM*DEBYE2ENM } },
        { "D.Angstrom3", { "e nm3", A2NM*A2NM*A2NM*DEBYE2ENM } },
        { "D*Angstrom",  { "e nm", A2NM*DEBYE2ENM } },
        { "D*Angstrom2", { "e nm2", A2NM*A2NM*DEBYE2ENM } },
        { "D*Angstrom3", { "e nm3", A2NM*A2NM*A2NM*DEBYE2ENM } },
        { "1/nm",        { "1/nm", 1 } },
        { "Da",          { "Da", 1 } },
        { "eV",          { "kJ/mol", ELECTRONVOLT } },
        { "eV/e",        { "kJ/mol e", ELECTRONVOLT } },
        { "ps",          { "ps", 1 } },
        { "fs",          { "ps", 1/KILO } },
        { "nm2",         { "nm2", 1 } },
        { "degree",      { "", M_PI/180.0 } },
        { "kJ/mol/rad2", { "kJ/mol/rad2", 1 } },
        { "kJ/mol/nm2",  { "kJ/mol/nm2", 1 } },
        { "e",           { "e", 1 } },
        { "cm^-1",       { "THz", SPEED_OF_LIGHT/1e7 } },
        { "km/mol",      { "nm/mol", 1e-12 } },
        { "nm/mol",      { "nm/mol", 1 } },
        { "",            { "", 1 } }
    };

std::string gromacsUnit(const std::string &unit)
{
    auto uc = unitConversion.find(unit);
    if (uc  == unitConversion.end())
    {
#if CMAKE_BUILD_TYPE == CMAKE_BUILD_TYPE_DEBUG
        GMX_THROW(gmx::InternalError(gmx::formatString("Unknown unit '%s'\n", unit.c_str()).c_str()));
#else
        fprintf(stderr, "Unknown unit getting GROMACS unit %s\n", unit.c_str());
#endif
        return "";
    }
    return uc->second.first;
    
}

double convertToGromacs(double x, const std::string &unit)
{
    auto uc = unitConversion.find(unit);
    if (uc  == unitConversion.end())
    {
#if CMAKE_BUILD_TYPE == CMAKE_BUILD_TYPE_DEBUG
        GMX_THROW(gmx::InternalError(gmx::formatString("Unknown unit '%s'\n", unit.c_str()).c_str()));
#else
        fprintf(stderr, "Unknown unit converting to GROMACS %s\n", unit.c_str());
#endif
        return 1;
    }
    return x*uc->second.second;
}

double convertFromGromacs(double x, const std::string &unit)
{
    auto uc = unitConversion.find(unit);
    if (uc  == unitConversion.end())
    {
#if CMAKE_BUILD_TYPE == CMAKE_BUILD_TYPE_DEBUG
        GMX_THROW(gmx::InternalError(gmx::formatString("Unknown unit %s\n", unit.c_str()).c_str()));
#else
        fprintf(stderr, "Unknown unit converting from GROMACS %s\n", unit.c_str());
#endif
        return 1;
    }
    return x/uc->second.second;
}

} // namespace alexandria

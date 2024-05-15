/*
 * This source file is part of the Alexandria Chemistry Toolkit.
 *
 * Copyright (C) 2024
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
#include "potential.h"

#include <cstring>
#include <map>

#include "gromacs/topology/ifunc.h"
#include "gromacs/utility/cstringutil.h"

namespace alexandria
{

static const std::map<Potential, const std::string> pot2str =
    {
        { Potential::NONE,                   "" },
        { Potential::LJ8_6,                  "LJ8_6" },
        { Potential::LJ12_6,                 "LJ12_6" },
        { Potential::LJ14_7,                 "LJ14_7" },
        { Potential::GENERALIZED_BUCKINGHAM, "GENERALIZED_BUCKINGHAM" },
        { Potential::WANG_BUCKINGHAM,        "WANG_BUCKINGHAM" },
        { Potential::EXPONENTIAL,            "EXPONENTIAL" },
        { Potential::COULOMB_POINT,          "COULOMB_POINT" },
        { Potential::COULOMB_GAUSSIAN,       "COULOMB_GAUSSIAN" },
        { Potential::COULOMB_SLATER,         "COULOMB_SLATER" },
        { Potential::HARMONIC_BONDS,         "HARMONIC_BONDS" },
        { Potential::MORSE_BONDS,            "MORSE_BONDS" },
        { Potential::CUBIC_BONDS,            "CUBIC_BONDS" },
        { Potential::HARMONIC_ANGLES,        "HARMONIC_ANGLES" },
        { Potential::LINEAR_ANGLES,          "LINEAR_ANGLES" },
        { Potential::UREY_BRADLEY_ANGLES,    "UREY_BRADLEY_ANGLES" },
        { Potential::HARMONIC_DIHEDRALS,     "HARMONIC_DIHEDRALS" },
        { Potential::FOURIER_DIHEDRALS,      "FOURIER_DIHEDRALS" },
        { Potential::PROPER_DIHEDRALS,       "PROPER_DIHEDRALS" },
        { Potential::POLARIZATION,           "POLARIZATION" },
        { Potential::VSITE1,                 "VSITE1" },
        { Potential::VSITE2,                 "VSITE2" },
        { Potential::VSITE2FD,               "VSITE2FD" },
        { Potential::VSITE3,                 "VSITE3" },
        { Potential::VSITE3FD,               "VSITE3FD" },
        { Potential::VSITE3FAD,              "VSITE3FAD" },
        { Potential::VSITE3OUT,              "VSITE3OUT" },
        { Potential::VSITE3OUTS,             "VSITE3OUTS" }
    };

static const std::map<Potential, int> act2gmx = {
    { Potential::LJ12_6,                 F_LJ },
    { Potential::WANG_BUCKINGHAM,        F_WBHAM },
    { Potential::GENERALIZED_BUCKINGHAM, F_GBHAM },
    { Potential::COULOMB_POINT,          F_COUL_SR },
    { Potential::HARMONIC_BONDS,         F_BONDS },
    { Potential::MORSE_BONDS,            F_MORSE },
    { Potential::CUBIC_BONDS,            F_CUBICBONDS },
    { Potential::HARMONIC_ANGLES,        F_ANGLES },
    { Potential::LINEAR_ANGLES,          F_LINEAR_ANGLES },
    { Potential::UREY_BRADLEY_ANGLES,    F_UREY_BRADLEY },
    { Potential::HARMONIC_DIHEDRALS,     F_IDIHS },
    { Potential::FOURIER_DIHEDRALS,      F_FOURDIHS },
    { Potential::PROPER_DIHEDRALS,       F_PDIHS },
    { Potential::POLARIZATION,           F_POLARIZATION },
    { Potential::VSITE1,                 F_VSITEN },
    { Potential::VSITE2,                 F_VSITE2 },
    { Potential::VSITE2FD,               F_VSITE2FD },
    { Potential::VSITE3,                 F_VSITE3 },
    { Potential::VSITE3FD,               F_VSITE3FD },
    { Potential::VSITE3FAD,              F_VSITE3FAD },
    { Potential::VSITE3OUT,              F_VSITE3OUT },
    { Potential::VSITE3OUTS,             F_VSITE3OUTS }
};

const std::string &potentialToString(Potential p)
{
    return pot2str.find(p)->second;
}

bool stringToPotential(const std::string &pname, Potential *p)
{
    // Search for native ACT names first.
    for (const auto &ps : pot2str)
    {
        if (pname == ps.second)
        {
            *p = ps.first;
            return true;
        }
    }
    // Now search for old-style GROMACS names
    for (const auto &ps : act2gmx)
    {
        size_t n1 = strlen(interaction_function[ps.second].name);
        if (n1 == pname.size() &&
            (gmx_strncasecmp(pname.c_str(), interaction_function[ps.second].name, n1) == 0))
        {
            *p = ps.first;
            return true;
        }
    }
    return false;
}

int potentialToGromacsType(Potential p)
{
    auto ppp = act2gmx.find(p);
    if (act2gmx.end() != ppp)
    {
        return ppp->second;
    }
    return -1;
}

const char *potentialToGromacsString(Potential p)
{
    int ftype = potentialToGromacsType(p);
    if (-1 == ftype)
    {
        return nullptr;
    }
    else
    {
        return interaction_function[ftype].name;
    }
}

static const std::map<ChargeType, Potential> cp = {
    { ChargeType::Point, Potential::COULOMB_POINT },
    { ChargeType::Gaussian, Potential::COULOMB_GAUSSIAN },
    { ChargeType::Slater, Potential::COULOMB_SLATER }
};

Potential chargeTypeToPotential(ChargeType c)
{
    return cp.find(c)->second;
}

ChargeType potentialToChargeType(Potential p)
{
    ChargeType ct = ChargeType::Point;
    for (const auto &cccp : cp)
    {
        if (cccp.second == p)
        {
            ct = cccp.first;
        }
    }
    return ct;
}

} // namespace alexandria

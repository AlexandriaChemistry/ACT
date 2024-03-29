/*
 * This source file is part of the Alexandria Chemistry Toolkit.
 *
 * Copyright (C) 2014-2024
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

#include "interactiontype.h"

#include <map>
#include <set>

#include "gromacs/utility/exceptions.h"
#include "gromacs/utility/stringutil.h"

namespace alexandria
{

typedef struct
{
    std::string name, description;
} NameDescr;

std::map<InteractionType, NameDescr> eitNames = {
    { InteractionType::BONDS,              { "BONDS", "bonded forces" } },
    { InteractionType::ANGLES,             { "ANGLES", "angle forces" } },
    { InteractionType::LINEAR_ANGLES,      { "LINEAR_ANGLES", "linear angle forces" } },
    { InteractionType::PROPER_DIHEDRALS,   { "PROPER_DIHEDRALS", "proper dihedrals" } },
    { InteractionType::IMPROPER_DIHEDRALS, { "IMPROPER_DIHEDRALS", "improper dihedrals" } },
    { InteractionType::VDW,                { "VANDERWAALS", "Van der Waals interactions" } },
    { InteractionType::DISPERSION,         { "DISPERSION", "London dispersion" } },
    { InteractionType::EXCHANGE,           { "EXCHANGE", "Exchange repulsion due to Pauli principle" } },
    { InteractionType::CHARGETRANSFER,     { "CHARGETRANSFER", "Charge transfer energy correction" } },
    { InteractionType::EPOT,               { "EPOT", "Potential energy" } },
    { InteractionType::POLARIZATION,       { "POLARIZATION", "polarization" } },
    { InteractionType::CONSTR,             { "CONSTR", "constraints" } },
    { InteractionType::VSITE2,             { "VSITE2", "virtual sites with two constructing atoms" } },
    { InteractionType::VSITE2FD,           { "VSITE2FD", "virtual sites with two constructing atoms with fixed distance to atom" } },
    { InteractionType::VSITE3,             { "VSITE3", "virtual sites with three constructing atoms" } },
    { InteractionType::VSITE3FD,           { "VSITE3FD", "virtual site with three constructing atoms and fixed distance" } },
    { InteractionType::VSITE3FAD,          { "VSITE3FAD", "virtual sites with three constructing atoms and fixed angle and distance" } },
    { InteractionType::VSITE3OUT,          { "VSITE3OUT", "virtual sites with three contructing atoms, out of plane" } },
    { InteractionType::VSITE3OUTS,         { "VSITE3OUTS", "virtual sites with three contructing atoms, out of plane, symmetric with respect to the plane" } },
    { InteractionType::COULOMB,            { "COULOMB", "Coulomb interactions" } },
    { InteractionType::BONDCORRECTIONS,    { "BONDCORRECTIONS", "bond charge corrections" } },
    { InteractionType::ELECTRONEGATIVITYEQUALIZATION, { "ELECTRONEGATIVITYEQUALIZATION", "electronegativity equalization" } },
    { InteractionType::CHARGE, { "CHARGE", "charge" } }
};

const std::string &interactionTypeToString(InteractionType iType)
{
    auto en = eitNames.find(iType);
    if (en == eitNames.end())
    {
        GMX_THROW(gmx::InternalError(gmx::formatString("No string corresponding to interaction type %d",
                                                       static_cast<int>(iType)).c_str()));
    }
    return en->second.name;
}

const std::string &interactionTypeToDescription(InteractionType iType)
{
    auto en = eitNames.find(iType);
    if (en == eitNames.end())
    {
        GMX_THROW(gmx::InternalError(gmx::formatString("No string corresponding to interaction type %d",
                                                       static_cast<int>(iType)).c_str()));
    }
    return en->second.description;
}

InteractionType stringToInteractionType(const std::string &name)
{
    for (auto &eit : eitNames)
    {
        if (name == eit.second.name)
        {
            return eit.first;
        }
    }
    GMX_THROW(gmx::InvalidInputError(gmx::formatString("No such interaction type '%s'",
                                                       name.c_str()).c_str()));
    return InteractionType::BONDS;
}

int interactionTypeToNatoms(InteractionType iType)
{
    switch (iType)
    {
    case InteractionType::PROPER_DIHEDRALS:
    case InteractionType::IMPROPER_DIHEDRALS:
    case InteractionType::VSITE3FD:
    case InteractionType::VSITE3FAD:
    case InteractionType::VSITE3OUT:
    case InteractionType::VSITE3OUTS:
    case InteractionType::VSITE3:
        return 4;
    case InteractionType::ANGLES:
    case InteractionType::LINEAR_ANGLES:
    case InteractionType::VSITE2:
    case InteractionType::VSITE2FD:
        return 3;
    case InteractionType::VDW:
    case InteractionType::DISPERSION:
    case InteractionType::EXCHANGE:
    case InteractionType::CHARGETRANSFER:
    case InteractionType::POLARIZATION:
    case InteractionType::COULOMB:
    case InteractionType::ELECTRONEGATIVITYEQUALIZATION:
        return 1;
    case InteractionType::BONDS:
    case InteractionType::CONSTR:
    case InteractionType::BONDCORRECTIONS:
        return 2;
    case InteractionType::CHARGE:
    case InteractionType::EPOT:
        return 0;
    }
    return 0;
}

bool isVsite(InteractionType iType)
{
    std::set myvs = { InteractionType::VSITE2, InteractionType::VSITE2FD,
                      InteractionType::VSITE3, InteractionType::VSITE3FD, InteractionType::VSITE3FAD,
                      InteractionType::VSITE3OUT, InteractionType::VSITE3OUTS };
    return (myvs.end() != myvs.find(iType));
}

} // namespace alexandria

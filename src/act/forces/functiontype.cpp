/*
 * This source file is part of the Alexandria Chemistry Toolkit.
 *
 * Copyright (C) 2014-2022
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

#include "functiontype.h"

#include <map>

#include "gromacs/utility/exceptions.h"
#include "gromacs/utility/stringutil.h"

namespace alexandria
{


typedef struct
{
    std::string name, description;
} FuncDescr;

std::map<FunctionType, FuncDescr> eftNames = {
    { FunctionType::MORSE,              { "Morse", "Morse potential function" } },
    { FunctionType::HARMONIC,           { "Harmonic", "Harmonic potential function" } },
    { FunctionType::UREY_BRADLEY,       { "Urey_Bradley", "UREY_BRADLEY for angles between three atoms" } },
    { FunctionType::COSINE,   			 { "Cosine", "Cosine function for angles between three atoms" } },
    { FunctionType::LINEAR_ANGLE, 		 { "LINEAR_ANGLE", "Linear angles" } },
    { FunctionType::IDIHS,              { "IDIHS", "Out of plane dihedral angles" } },
    { FunctionType::PDIHS, 		       { "PDIHS", "Proper dihedrals" } },
    { FunctionType::FOURDIHS,           { "FOURDIHS", "Proper dihedral" } },
    { FunctionType::DIRAC_DELTA,        { "Dirac", "Coulomb between two point charges" } },
    { FunctionType::GAUSSIAN,           { "Gaussian", "Coulomb between two Gaussian-type charges" } },
    { FunctionType::SLATER,	          { "Slater", "Coulomb between two Slater-type charges" } },
    { FunctionType::LENNARD_JONES, 		 { "Lennard Jones", "VDW functions" } },
    { FunctionType::BUCKINGHAM,    		 { "Buckingham", "VDW functions" } },
    { FunctionType::WANG_BUCKINGHAM, 	 { "Wang_Buckingham", "Buffered Buckingham potential" } }
};

const std::string &functionTypeToString(FunctionType fType)
{
    auto en = eftNames.find(fType);
    if (en == eftNames.end())
    {
        GMX_THROW(gmx::InternalError(gmx::formatString("No string corresponding to function type %d",
                                                       static_cast<int>(fType)).c_str()));
    }
    return en->second.name;
}

const std::string &functionTypeToDescription(FunctionType fType)
{
    auto en = eftNames.find(fType);
    if (en == eftNames.end())
    {
        GMX_THROW(gmx::InternalError(gmx::formatString("No string corresponding to function type %d",
                                                       static_cast<int>(fType)).c_str()));
    }
    return en->second.description;
}

FunctionType stringToFunctionType(const std::string &name)
{
    for (auto &eft : eftNames)
    {
        if (name == eft.second.name)
        {
            return eft.first;
        }
    }
    GMX_THROW(gmx::InvalidInputError(gmx::formatString("No such function type %s",
                                                       name.c_str()).c_str()));
    return FunctionType::MORSE;
}

int functionTypeToNparams(FunctionType fType)
{
    switch (fType)
    {
    case FunctionType::COSINE:
        return 5;
    case FunctionType::FOURDIHS:
    case FunctionType::UREY_BRADLEY:
        return 4;
    case FunctionType::MORSE:
    case FunctionType::BUCKINGHAM:
    case FunctionType::WANG_BUCKINGHAM:
        return 3;
    case FunctionType::HARMONIC:
    case FunctionType::LINEAR_ANGLE:
    case FunctionType::IDIHS:
    case FunctionType::LENNARD_JONES:
    case FunctionType::COUL_SR:
    case FunctionType::PDIHS:
        return 2;
    case FunctionType::GAUSSIAN:
    case FunctionType::SLATER:
        return 1;
    case FunctionType::DIRAC_DELTA:
        return 0;
    }
    return 0;
}
} // namespace

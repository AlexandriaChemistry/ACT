/*
 * This source file is part of the Alexandria Chemistry Toolkit.
 *
 * Copyright (C) 2019
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
 * \author Mohammad Mehdi Ghahremanpour <mohammad.ghahremanpour@icm.uu.se>
 * \author David van der Spoel <david.vanderspoel@icm.uu.se>
 */

#include "chargemodel.h"

#include <algorithm>
#include <map>
#include <string>

#include "gromacs/utility/exceptions.h"
#include "gromacs/utility/stringutil.h"

namespace alexandria
{

std::map<ChargeType, const std::string> ct2Name =
    {
        { ChargeType::Point,    "Point"    },
        { ChargeType::Gaussian, "Gaussian" },
        { ChargeType::Slater,   "Slater"   }
    };

std::map<const std::string, ChargeType> name2CT;

ChargeType name2ChargeType(const std::string &name)
{
    if (name2CT.empty())
    {
        for(auto &k : ct2Name)
        {
            name2CT.emplace(k.second, k.first);
        }
    }
    auto cc = name2CT.find(name);
    if (cc == name2CT.end())
    {
        GMX_THROW(gmx::InvalidInputError(gmx::formatString("Unknown charge type %s. Note that input is case-sensitive.\n", name.c_str()).c_str()));
    }
    return cc->second;
}

const std::string &chargeTypeName(ChargeType ct)
{
    auto cc = ct2Name.find(ct);

    if (cc == ct2Name.end())
    {
        GMX_THROW(gmx::InternalError("Invalid ChargeType."));
    }
    
    return cc->second;
}

std::map<ChargeGenerationAlgorithm, const std::string> cg2Name =
    {
        { ChargeGenerationAlgorithm::NONE,      "None"      },
        { ChargeGenerationAlgorithm::EEM,       "EEM"       },
        { ChargeGenerationAlgorithm::SQE,       "SQE"       },
        { ChargeGenerationAlgorithm::ESP,       "ESP"       },
        { ChargeGenerationAlgorithm::Custom,    "Custom"    },
        { ChargeGenerationAlgorithm::CM5,       "CM5"       },
        { ChargeGenerationAlgorithm::Hirshfeld, "Hirshfeld" },
        { ChargeGenerationAlgorithm::Mulliken,  "Mulliken"  }
    };

std::map<const std::string, ChargeGenerationAlgorithm> name2CG;

ChargeGenerationAlgorithm nameToChargeGenerationAlgorithm(const std::string &name)
{
    if (name2CG.empty())
    {
        for(auto &k : cg2Name)
        {
            name2CG.emplace(k.second, k.first);
        }
    }
    auto cc = name2CG.find(name);
    if (cc == name2CG.end())
    {
        GMX_THROW(gmx::InvalidInputError(gmx::formatString("Unknown charge generation algorithm %s. Note that input is case-sensitive.\n",
                                                           name.c_str()).c_str()));
    }
    return cc->second;
}

const std::string &chargeGenerationAlgorithmName(ChargeGenerationAlgorithm cg)
{
    auto cc = cg2Name.find(cg);

    GMX_RELEASE_ASSERT(cc != cg2Name.end(), "Internal error.");
    
    return cc->second;
}

} // namespace alexandria

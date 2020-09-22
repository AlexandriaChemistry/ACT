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

#include "gromacs/utility/gmxassert.h"

namespace alexandria
{

std::map<ChargeType, const std::string> ct2Name =
    {
        { eqtPoint,    "Point"    },
        { eqtGaussian, "Gaussian" },
        { eqtSlater,   "Slater"   }
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
    if (cc != name2CT.end())
    {
        return cc->second;
    }
    fprintf(stderr, "Unknown charge type %s. Note that input is case-sensitive.\n",
            name.c_str());
    return eqtNR;
}

const std::string &chargeTypeName(ChargeType ct)
{
    auto cc = ct2Name.find(ct);

    GMX_RELEASE_ASSERT(cc != ct2Name.end(), "Internal error.");
    
    return cc->second;
}

std::map<ChargeGenerationAlgorithm, const std::string> cg2Name =
    {
        { eqgNONE, "None" },
        { eqgEEM,  "EEM"  },
        { eqgSQE,  "SQE"  },
        { eqgESP,  "ESP"  }
    };

std::map<const std::string, ChargeGenerationAlgorithm> name2CG;

ChargeGenerationAlgorithm name2ChargeGenerationAlgorithm(const std::string &name)
{
    if (name2CG.empty())
    {
        for(auto &k : cg2Name)
        {
            name2CG.emplace(k.second, k.first);
        }
    }
    auto cc = name2CG.find(name);
    if (cc != name2CG.end())
    {
        return cc->second;
    }
    fprintf(stderr, "Unknown charge generation algorithm %s. Note that input is case-sensitive.\n",
            name.c_str());
    return eqgNR;
}

const std::string &chargeGenerationAlgorithmName(ChargeGenerationAlgorithm cg)
{
    auto cc = cg2Name.find(cg);

    GMX_RELEASE_ASSERT(cc != cg2Name.end(), "Internal error.");
    
    return cc->second;
}

} // namespace alexandria

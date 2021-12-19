/*
 * This source file is part of the Alexandria Chemistry Toolkit.
 *
 * Copyright (C) 2014-2020 
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

#include "phase.h"

#include <map>

#include "gromacs/utility/exceptions.h"

namespace alexandria
{

static std::map<const std::string, ePhase> stringToPhase = 
    {
        { "gas",    ePhase::GAS },
        { "liquid", ePhase::LIQUID }, 
        { "solid",  ePhase::SOLID },
        { "plasma", ePhase::PLASMA }
    };
    
const std::string &phase2string(ePhase ep)
{
    for(auto &i : stringToPhase)
    {
        if (i.second == ep)
        {
            return i.first;
        }
    }
    GMX_THROW(gmx::InvalidInputError("Invalid phase"));
    return stringToPhase.begin()->first;
}

ePhase string2phase(const std::string &phase)
{
    auto ss = stringToPhase.find(phase);
    if (ss != stringToPhase.end())
    {
        return ss->second;
    }
    GMX_THROW(gmx::InvalidInputError("Invalid phase"));
    return ePhase::GAS;
}

} // namespace alexandria

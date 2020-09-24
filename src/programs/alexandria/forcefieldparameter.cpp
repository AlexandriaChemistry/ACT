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

#include "forcefieldparameter.h"

#include <map>

#include "gromacs/utility/exceptions.h"
#include "gromacs/utility/stringutil.h"

namespace alexandria
{

static std::map<Mutability, const std::string> mut2string =
    {
        { Mutability::Fixed,   "Fixed"   },
        { Mutability::Bounded, "Bounded" },
        { Mutability::Free,    "Free"    }
    };

static std::map<const std::string, Mutability> string2mut;
    
const std::string &mutabilityName(Mutability mutability)
{
    auto m2s = mut2string.find(mutability);
    
    return m2s->second;
}

bool nameToMutability(const std::string &name, Mutability *mutability)
{
    if (string2mut.empty())
    {
        for (auto iter = mut2string.begin(); iter != mut2string.end(); ++iter)
        {
            string2mut.insert({iter->second, iter->first});
        }
    }
    auto s2m = string2mut.find(name);
    if (s2m != string2mut.end())
    {
        *mutability = s2m->second;
        return true;
    }
    else
    {
        return false;
    }
}

void ForceFieldParameter::setValue(double value)
{
    switch (mutability_)
    {
    case Mutability::Free:
        value_ = value;
        break;
    case Mutability::Bounded:
        if (value >= minimum_ && value <= maximum_)
        {
            value_ = value;
        }
        else
        {
            double newval = std::min(maximum_, std::max(minimum_, value));
            if (strict_)
            {
                auto buf = gmx::formatString("Can not modify value %s outside its bounds of %g-%g. Setting it to %g.",
                                             identifier_.c_str(), minimum_, maximum_, newval);
                    GMX_THROW(gmx::InvalidInputError(buf));
            }
            value_ = newval;
        }
        break;
    case Mutability::Fixed:
        if (strict_)
        {
            auto buf = gmx::formatString("Cannot modify parameter %s since it is fixed", identifier_.c_str());
            GMX_THROW(gmx::InvalidInputError(buf));
        }
        break;
    }
}

void ForceFieldParameter::setUncertainty(double uncertainty)
{ 
    if (mutability_ != Mutability::Fixed)
    {
        uncertainty_ = std::max(0.0, uncertainty);
    }
    else if (strict_)
    {
        auto buf = gmx::formatString("Cannot modify uncertainty %s since the parameter is fixed", identifier_.c_str());
        GMX_THROW(gmx::InternalError(buf));
    }
}

} // namespace alexandria

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

#include "forcefieldparameterlist.h"

#include "gromacs/utility/exceptions.h"
#include "gromacs/utility/stringutil.h"

namespace alexandria
{

void ForceFieldParameterList::addParameter(const std::string         &identifier,
                                           const ForceFieldParameter &param)
{
    auto params = parameters_.find(identifier);
    
    if (params == parameters_.end())
    {
        // New parameter!
        std::vector<ForceFieldParameter> entry;
        entry.push_back(param);
        parameters_.insert({identifier, entry});
    }
    else
    {
        // TODO: Check whether a parameter of this type already exists
        params->second.push_back(param);
    }
}

const std::vector<ForceFieldParameter> &ForceFieldParameterList::searchParameter(const std::string &identifier) const
{
    auto params = parameters_.find(identifier);
    
    if (params == parameters_.end())
    {
        GMX_THROW(gmx::InvalidInputError(gmx::formatString("No such identifier %s in parameter list for %s", identifier.c_str(), function_.c_str()).c_str()));
    }
    
    return params->second;
}

std::vector<ForceFieldParameter> &ForceFieldParameterList::findParameter(const std::string &identifier)
{
    auto params = parameters_.find(identifier);
    
    GMX_RELEASE_ASSERT(params != parameters_.end(),
                       gmx::formatString("No such identifier %s in parameter list for %s", identifier.c_str(), function_.c_str()).c_str());
    
    return params->second;
}

} // namespace

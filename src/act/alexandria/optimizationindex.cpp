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
 * \author Julian Ramon Marrades Furquet <julian.marrades@hotmail.es>
 */

#include "optimizationindex.h"
#include "act/forcefield/forcefield.h"
#include "act/forcefield/forcefield_parameter.h"
#include "gromacs/utility/stringutil.h"

namespace alexandria
{

void OptimizationIndex::findForceFieldParameter(ForceField *pd)
{
    if (iType_ != InteractionType::CHARGE)
    {
        ffp_ = pd->findForces(iType_)->findParameterType(parameterId_, parameterType_);
    }
    else if (pd->hasParticleType(particleType_))
    {
        auto pt = pd->findParticleType(particleType_);
        ffp_    = pt->parameter(parameterType_);
    }
    GMX_RELEASE_ASSERT(ffp_, gmx::formatString("Could not find parameter %s",
                                               parameterId_.id().c_str()).c_str());
}

std::string OptimizationIndex::name() const
{
    if (InteractionType::CHARGE == iType_)
    {
        return gmx::formatString("%s-%s",
                                 particleType_.c_str(),
                                 parameterType_.c_str());
    }
    else
    {
        return gmx::formatString("%s-%s",
                                 parameterId_.id().c_str(),
                                 parameterType_.c_str());
    }
}

}

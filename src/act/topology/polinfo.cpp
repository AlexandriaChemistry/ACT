/*
 * This source file is part of the Alexandria Chemistry Toolkit.
 *
 * Copyright (C) 2014-2026
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
 
#include "polinfo.h"

#include "act/basics/identifier.h"
#include "act/forcefield/forcefield.h"
#include "act/forcefield/forcefield_parametername.h"
#include "act/topology/topology.h"

namespace alexandria
{

void PolInfo::init(const ForceField *pd,
                   const Topology   *top)
{
    isShell_.clear();
    nShell_ = 0;
    if (!pd->interactionPresent(InteractionType::POLARIZATION))
    {
        return;
    }
    for(auto &aa : top->atoms())
    {
        bool bIS = aa.pType() == ActParticle::Shell;
        isShell_.push_back(bIS);
        if (bIS)
        {
            nShell_ += 1;
        }
    }
}

} // namespace

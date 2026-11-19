/*
 * This source file is part of the Alexandria Chemistry Toolkit.
 *
 * Copyright (C) 2026
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
/*! \internal \file
 * \brief
 * Implements functions that describe a molecule in a JsonTree.
 *
 * \author David van der Spoel <david.vanderspoel@icm.uu.se>
 */
#include "actpre.h"

#include <cmath>
#include <string>
#include <vector>

#include "act/basics/act_particle.h"
#include "act/basics/interactiontype.h"
#include "act/topology/actatom.h"
#include "act/topology/actmol.h"
#include "act/topology/topology.h"
#include "act/utility/jsontree.h"
#include "actmol_json.h"
#include "gromacs/math/vectypes.h"

namespace alexandria
{

JsonTree moleculeToJson(const ACTMol *mol, const std::string &key)
{
    JsonTree    jmol(key);
    const auto &atoms = mol->atomsConst();
    double      qtot  = 0;
    for (const auto &atom : atoms)
    {
        qtot += atom.totalCharge();
    }
    jmol.addObject("name", mol->getMolname());
    jmol.addObject("formula", mol->formula());
    jmol.addObject("charge", static_cast<int>(std::lround(qtot)));
    jmol.addValueUnit("mass", mol->totalMass(), "g/mol");

    JsonTree jparticles = JsonTree::array("particles");
    for (const auto &atom : atoms)
    {
        JsonTree jp("particle");
        jp.addObject("name", atom.name());
        jp.addObject("element", atom.element());
        jp.addObject("atomicnumber", atom.atomicNumber());
        jp.addObject("type", atom.ffType());
        jp.addObject("particletype", actParticleToString(atom.pType()));
        jp.addObject("mass", atom.mass());
        jp.addObject("charge", atom.charge());
        if (atom.chargeVS() != 0)
        {
            jp.addObject("charge_vs", atom.chargeVS());
        }
        if (atom.pType() != ActParticle::Atom)
        {
            jp.addArray("parents", atom.parents());
        }
        jparticles.addObject(jp);
    }
    jmol.addObject(jparticles);

    JsonTree jbonds = JsonTree::array("bonds");
    if (mol->topology()->hasEntry(InteractionType::BONDS))
    {
        for (const auto &bond : mol->topology()->entry(InteractionType::BONDS))
        {
            JsonTree jb("bond");
            jb.addArray("atoms", std::vector<int>({ bond->atomIndex(0), bond->atomIndex(1) }));
            jb.addObject("order", bond->bondOrder());
            jbonds.addObject(jb);
        }
    }
    jmol.addObject(jbonds);
    return jmol;
}

JsonTree coordinatesToJson(const std::string            &key,
                           const std::vector<gmx::RVec> &coords)
{
    JsonTree jx = JsonTree::array(key);
    for (const auto &x : coords)
    {
        jx.addArray("", std::vector<real>({ x[XX], x[YY], x[ZZ] }));
    }
    return jx;
}

} // namespace alexandria

/*
 * This source file is part of the Alexandria Chemistry Toolkit.
 *
 * Copyright (C) 2014-2025
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


#include "act/forcefield/combinationrules.h"

    
    
/*! \brief Generation parameters using combination rules
 * \param[in] pd    The force field
 * \param[in] itype The interaction type
 * \param[in] force If set, all interaction parameters will be recomputed, if not only the ones for which 
 *                  at least one of the two constituting parameters changed will be reset.
 */
static void generateParameterPairs(ForceField      *pd,
                                   InteractionType  itype,
                                   bool             force)
{
    // Do not crash if e.g. there is no VDWCORRECTION.
    if (!pd->interactionPresent(itype))
    {
        return;
    }
    auto forcesVdw = pd->findForces(itype);
    auto comb_rule = getCombinationRule(*forcesVdw);

    // We temporarily store the new parameters here
    ForceFieldParameterListMap *parm = forcesVdw->parameters();;
    
    // Now do the double loop
    int nid = 0;
    for (auto &ivdw : *forcesVdw->parameters())
    {
        auto &iid    = ivdw.first;
        // Check whether this is a single atom parameter
        if (iid.atoms().size() > 1)
        {
            continue;
        }
        auto &iparam = ivdw.second;
        bool iupdated = false;
        for(const auto &ip : iparam)
        {
            iupdated = iupdated || ip.second.updated();
        }
        for (auto &jvdw : *forcesVdw->parameters())
        {
            auto &jid    = jvdw.first;
            // Check whether this is a single atom parameter and
            // whether this is is larger or equal to iid.
            if (jid.atoms().size() > 1 || jid.id() < iid.id())
            {
                continue;
            }
            // Test whether or not to include this pair.
            // This will only be used with the Kronecker combination rule.
            bool includePair = true;
            if (Potential::BORN_MAYER == forcesVdw->potential() || 
                Potential::MORSE_BONDS == forcesVdw->potential())
            {
                auto ai = iid.atoms()[0];
                auto aj = jid.atoms()[0];
                if (pd->hasParticleType(ai) && pd->hasParticleType(aj))
                {
                    auto pti = pd->findParticleType(ai)->apType();
                    auto ptj = pd->findParticleType(aj)->apType();
                    // The pair will be included only if one particle is a vsite
                    // and the other an atom.
                    includePair = ((ActParticle::Vsite == pti && ActParticle::Atom == ptj) ||
                                   (ActParticle::Vsite == ptj && ActParticle::Atom == pti));
                }
            }
            auto &jparam = jvdw.second;
            bool jupdated = false;
            for(const auto &jp : jparam)
            {
                jupdated = jupdated || jp.second.updated();
            }
            if (!(iupdated || jupdated || force))
            {
                continue;
            }
            // Fill the parameters, potential dependent
            ForceFieldParameterMap pmap;
            evalCombinationRule(forcesVdw->potential(),
                                comb_rule, ivdw.second, jvdw.second, includePair, &pmap);

            parm->insert_or_assign(Identifier({ iid.id(), jid.id() }, { 1 }, CanSwap::Yes),
                                   std::move(pmap));
            nid += 1;
        }
    }
    if (debug)
    {
        int np = forcesVdw->parameters()->size();
        fprintf(debug, "Made %d/%d identifiers for %s in generateParameterPairs\n", nid, np,
                interactionTypeToString(itype).c_str());
    }
    // Phew, we're done!
}

/*! \brief Generate coulomb pair parameters
 * \param[in] pd    The force field
 * \param[in] force If set, all interaction parameters will be recomputed, if not only the ones for which 
 *                  at least one of the two constituting parameters changed will be reset.
 */
static void generateCoulombParameterPairs(ForceField *pd, bool force)
{
    auto forcesCoul = pd->findForces(InteractionType::ELECTROSTATICS);
    
    // Fudge unit
    std::string unit("kJ/mol");
    
    // We use dependent mutability to show these are not independent params
    auto mutd = Mutability::Dependent;
    auto cname = potentialToParameterName(Potential::COULOMB_GAUSSIAN);
    auto zeta = cname[coulZETA];
    // Finally add the new parameters to the exisiting list
    auto fold = forcesCoul->parameters();
    // Now do the double loop
    int nid = 0;
    for (auto &icoul : *forcesCoul->parameters())
    {
        auto &iid    = icoul.first;
        // Check whether this is a single atom parameter
        if (iid.atoms().size() > 1)
        {
            continue;
        }
        double izeta = icoul.second[zeta].internalValue();
        for (auto &jcoul : *forcesCoul->parameters())
        {
            auto &jid    = jcoul.first;
            // Check whether this is a single atom parameter and
            // whether this is is larger or equal to iid.
            if (jid.atoms().size() > 1 || jid.id() < iid.id() ||
                (!icoul.second[zeta].updated() && !jcoul.second[zeta].updated() && !force))
            {
                continue;
            }
            double     jzeta  = jcoul.second[zeta].internalValue();
            Identifier pairID({ iid.id(), jid.id() }, { 1 }, CanSwap::Yes);
            nid += 1;
            auto       oldfp  = fold->find(pairID);
            if (oldfp == fold->end())
            {
                ForceFieldParameterMap ffpm = {
                    { cname[coulZETA], 
                      ForceFieldParameter(unit, izeta, 0, 0, izeta, izeta, mutd, false, true) },
                    { cname[coulZETA2],
                      ForceFieldParameter(unit, jzeta, 0, 0, jzeta, jzeta, mutd, false, true) }
                };
                fold->insert({pairID, ffpm});
            }
            else
            {
                auto &pi = oldfp->second.find(cname[coulZETA])->second;
                pi.forceSetValue(izeta);
                auto &pj = oldfp->second.find(cname[coulZETA2])->second;
                pj.forceSetValue(jzeta);
            }
        }
    }
    if (debug)
    {
        int np = forcesCoul->parameters()->size();
        fprintf(debug, "Made %d/%d identifiers in generateCoulombParameterPairs\n", nid, np);
    }
    // Phew, we're done!
}

    
void generateDependentParameter(ForceField *pd, bool force)
{
    generateParameterPairs(pd, InteractionType::VDW, force);
    generateParameterPairs(pd, InteractionType::VDWCORRECTION, force);
    generateParameterPairs(pd, InteractionType::INDUCTIONCORRECTION, force);
    generateCoulombParameterPairs(pd, force);
}


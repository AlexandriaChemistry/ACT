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
#include "generate_dependent.h"

#include "act/forcefield/forcefield.h"
#include "act/forcefield/forcefield_parametername.h"

namespace alexandria
{

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
    auto comb_rule = forcesVdw->combinationRules();
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
        // Check whether the parameter for atom i has been updated
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
            // Check whether the parameter for atom j has been updated
            auto &jparam = jvdw.second;
            bool jupdated = false;
            for(const auto &jp : jparam)
            {
                jupdated = jupdated || jp.second.updated();
            }
            // Check whether the combination rule was updated
            bool crupdated = false;
            for(const auto &cr : comb_rule)
            {
                crupdated = crupdated || cr.second.ffplConst().updated();
            }
            if (!(iupdated || jupdated || crupdated || force))
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
    // Finally add the new parameters to the existing list
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

/*! \brief Generate quadrupole polarization parameters
 * \param[in] pd    The force field
 * \param[in] force If set, all interaction parameters will be recomputed, if not only the ones for which 
 *                  at least one of the two constituting parameters changed will be reset.
 */
static void generateQpolParameterPairs(ForceField *pd, bool force)
{
    auto itQPol = InteractionType::QUADRUPOLE_POLARIZATION;
    if (!pd->interactionPresent(itQPol))
    {
        return;
    }
    auto forcesQpol = pd->findForces(itQPol);
    
    // Fudge unit
    std::string unit("kJ/(mol e) nm6");
    
    // We use dependent mutability to show these are not independent params
    auto mutd = Mutability::Dependent;
    auto cname = potentialToParameterName(Potential::QUADRUPOLE_POLARIZATION);
    auto c6    = cname[qpolC6];
    auto b     = cname[qpolB];
    // Finally add the new parameters to the existing list
    auto fold = forcesQpol->parameters();
    // Now do the double loop
    int nid = 0;
    for (auto &iqpol : *forcesQpol->parameters())
    {
        auto &iid    = iqpol.first;
        // Check whether this is a single atom parameter
        if (iid.atoms().size() > 1)
        {
            continue;
        }
        double c6i = iqpol.second[c6].internalValue();
        double bi  = iqpol.second[b].internalValue();
        for (auto &jqpol : *forcesQpol->parameters())
        {
            auto &jid    = jqpol.first;
            // Check whether this is a single atom parameter and
            // whether this is is larger or equal to iid.
            if (jid.atoms().size() > 1 || jid.id() < iid.id() ||
                (!iqpol.second[c6].updated() && !jqpol.second[c6].updated() &&
                 !iqpol.second[b].updated() && !jqpol.second[b].updated() && !force))
            {
                continue;
            }
            double     c6j = jqpol.second[c6].internalValue();
            double     bj  = jqpol.second[b].internalValue();
            Identifier pairID({ iid.id(), jid.id() }, { 1 }, CanSwap::Yes);
            nid += 1;
            auto       oldfp  = fold->find(pairID);
            if (oldfp == fold->end())
            {
                ForceFieldParameterMap ffpm = {
                    { cname[qpolC6i], 
                      ForceFieldParameter(unit, c6i, 0, 0, c6i, c6i, mutd, false, true) },
                    { cname[qpolC6j], 
                      ForceFieldParameter(unit, c6j, 0, 0, c6j, c6j, mutd, false, true) },
                    { cname[qpolBi],
                      ForceFieldParameter(unit, bi, 0, 0, bi, bi, mutd, false, true) },
                    { cname[qpolBj],
                      ForceFieldParameter(unit, bj, 0, 0, bj, bj, mutd, false, true) }
                };
                fold->insert({pairID, ffpm});
            }
            else
            {
                auto &pc6i = oldfp->second.find(cname[qpolC6i])->second;
                pc6i.forceSetValue(c6i);
                auto &pc6j = oldfp->second.find(cname[qpolC6j])->second;
                pc6j.forceSetValue(c6j);
                auto &pbi = oldfp->second.find(cname[qpolBi])->second;
                pbi.forceSetValue(bi);
                auto &pbj = oldfp->second.find(cname[qpolBj])->second;
                pbj.forceSetValue(bj);
            }
        }
    }
    if (debug)
    {
        int np = forcesQpol->parameters()->size();
        fprintf(debug, "Made %d/%d identifiers in generateQpolParameterPairs\n", nid, np);
    }
    // Phew, we're done!
}

void generateDependentParameter(ForceField *pd, bool force)
{
    generateParameterPairs(pd, InteractionType::VDW, force);
    generateParameterPairs(pd, InteractionType::VDWCORRECTION, force);
    generateParameterPairs(pd, InteractionType::INDUCTIONCORRECTION, force);
    generateCoulombParameterPairs(pd, force);
    generateQpolParameterPairs(pd, force);
}

} // namespace

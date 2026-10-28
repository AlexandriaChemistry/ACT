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

#include <algorithm>
#include <vector>

#include "act/forcefield/forcefield.h"
#include "act/forcefield/forcefield_parametername.h"

namespace alexandria
{

//! The parameters of a single particle type, from which pair parameters are made
struct SingleParameters
{
    //! The identifier
    const Identifier       *id;
    //! The parameters
    ForceFieldParameterMap *params;
};

/*! \brief Collect the parameters of single particle types
 * The list also contains the pair parameters, which are skipped. Since
 * the list is a std::map, the pointers stay valid when pairs are added.
 * \param[in] parameters All parameters of an interaction type
 * \return The single particle parameters, sorted by identifier, such that
 *         for pairs (i, j) with i <= j the identifier of i is not larger
 */
static std::vector<SingleParameters> singleParameters(ForceFieldParameterListMap *parameters)
{
    std::vector<SingleParameters> singles;
    for (auto &p : *parameters)
    {
        if (p.first.atoms().size() <= 1)
        {
            singles.push_back({ &p.first, &p.second });
        }
    }
    std::sort(singles.begin(), singles.end(),
              [](const SingleParameters &a, const SingleParameters &b)
              { return a.id->id() < b.id->id(); });
    return singles;
}

/*! \brief Generation parameters using combination rules
 * \param[in]    msghandler For debugging info
 * \param[inout] pd         The force field
 * \param[in]    itype      The interaction type
 * \param[in]    force      If set, all interaction parameters will be
 *                          recomputed, if not only the ones for which
 *                          at least one of the two constituting parameters
 *                          changed will be reset.
 */
static void generateParameterPairs(MsgHandler      *msghandler,
                                   ForceField      *pd,
                                   InteractionType  itype,
                                   bool             force)
{
    // Do not crash if e.g. there is no VDWCORRECTION.
    if (!pd->interactionPresent(itype))
    {
        return;
    }
    auto        forcesVdw = pd->findForces(itype);
    const auto &comb_rule = forcesVdw->combinationRules();
    // The new parameters are added to the same list
    ForceFieldParameterListMap *parm = forcesVdw->parametersNoRemoval();

    // If the combination rule was updated all pairs change
    bool crupdated = false;
    for(const auto &cr : comb_rule)
    {
        crupdated = crupdated || cr.second.ffplConst().updated();
    }
    auto singles = singleParameters(parm);
    // Check whether the parameters for each particle type have been updated
    std::vector<bool> updated;
    for(const auto &single : singles)
    {
        bool supdated = false;
        for(const auto &p : *single.params)
        {
            supdated = supdated || p.second.updated();
        }
        updated.push_back(supdated);
    }
    // With the Kronecker combination rule pairs are only included if one
    // particle is a vsite and the other an atom.
    bool kronecker = (Potential::BORN_MAYER == forcesVdw->potential() ||
                      Potential::MORSE_BONDS == forcesVdw->potential());
    std::vector<const ParticleType *> ptypes;
    if (kronecker)
    {
        for(const auto &single : singles)
        {
            const auto &atom = single.id->atoms()[0];
            ptypes.push_back(pd->hasParticleType(atom) ? pd->findParticleType(atom) : nullptr);
        }
    }

    // Now do the double loop
    int nid = 0;
    for (size_t i = 0; i < singles.size(); i++)
    {
        for (size_t j = i; j < singles.size(); j++)
        {
            if (!(updated[i] || updated[j] || crupdated || force))
            {
                continue;
            }
            bool includePair = true;
            if (kronecker && ptypes[i] && ptypes[j])
            {
                auto pti = ptypes[i]->apType();
                auto ptj = ptypes[j]->apType();
                includePair = ((ActParticle::Vsite == pti && ActParticle::Atom == ptj) ||
                               (ActParticle::Vsite == ptj && ActParticle::Atom == pti));
            }
            // Fill the parameters, potential dependent
            ForceFieldParameterMap pmap;
            evalCombinationRule(forcesVdw->potential(),
                                comb_rule, *singles[i].params, *singles[j].params, includePair, &pmap);
            if (!pmap.empty())
            {
                auto myId = Identifier({ singles[i].id->id(), singles[j].id->id() }, { 1 }, CanSwap::Yes);
                if (msghandler && msghandler->debug())
                {
                    auto mystr = gmx::formatString("Generated pair parameters for %s / %s",
                                                   potentialToString(forcesVdw->potential()).c_str(),
                                                   myId.id().c_str());
                    for(const auto &p : pmap)
                    {
                        mystr += gmx::formatString(" %s: %g", p.first.c_str(), p.second.value());
                    }
                    msghandler->writeDebug(mystr);
                }
                parm->insert_or_assign(std::move(myId), std::move(pmap));
                nid += 1;
            }
        }
    }
    if (msghandler && msghandler->debug())
    {
        int np = forcesVdw->numberOfParameters();
        msghandler->writeDebug(gmx::formatString("Made %d/%d identifiers for %s in generateParameterPairs\n",
                                                 nid, np, interactionTypeToString(itype).c_str()));
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
    auto fold = forcesCoul->parametersNoRemoval();
    // The zeta of each particle type and whether it was updated
    auto singles = singleParameters(fold);
    std::vector<double> zetas;
    std::vector<bool>   updated;
    for(const auto &single : singles)
    {
        const auto &p = (*single.params)[zeta];
        zetas.push_back(p.internalValue());
        updated.push_back(p.updated());
    }
    // Now do the double loop
    int nid = 0;
    for (size_t i = 0; i < singles.size(); i++)
    {
        double izeta = zetas[i];
        for (size_t j = i; j < singles.size(); j++)
        {
            if (!updated[i] && !updated[j] && !force)
            {
                continue;
            }
            double     jzeta  = zetas[j];
            Identifier pairID({ singles[i].id->id(), singles[j].id->id() }, { 1 }, CanSwap::Yes);
            nid += 1;
            auto       oldfp  = fold->find(pairID);
            if (oldfp == fold->end())
            {
                ForceFieldParameterMap ffpm = {
                    { cname[coulZETA], 
                      ForceFieldParameter(unit, izeta, 0, 1, izeta, izeta, mutd, false, true) },
                    { cname[coulZETA2],
                      ForceFieldParameter(unit, jzeta, 0, 1, jzeta, jzeta, mutd, false, true) }
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
        int np = forcesCoul->numberOfParameters();
        fprintf(debug, "Made %d/%d identifiers in generateCoulombParameterPairs\n", nid, np);
    }
    // Phew, we're done!
}
    
void generateDependentParameter(MsgHandler *msghandler,
                                ForceField *pd,
                                bool        force)
{
    generateParameterPairs(msghandler, pd, InteractionType::VDW, force);
    generateParameterPairs(msghandler, pd, InteractionType::VDWCORRECTION, force);
    generateParameterPairs(msghandler, pd, InteractionType::INDUCTIONCORRECTION, force);
    generateCoulombParameterPairs(pd, force);
}

} // namespace

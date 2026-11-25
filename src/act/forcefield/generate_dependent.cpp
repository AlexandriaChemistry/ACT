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
#include <algorithm>
#include <cmath>
#include <cstddef>
#include <memory>
#include <set>
#include <type_traits>
#include <utility>
#include <vector>

#include "act/basics/act_particle.h"
#include "act/basics/identifier.h"
#include "act/basics/interactiontype.h"
#include "act/forcefield/combinationrules.h"
#include "act/forcefield/forcefield.h"
#include "act/forcefield/forcefield_parameter.h"
#include "act/forcefield/forcefield_parameterlist.h"
#include "act/forcefield/forcefield_parametername.h"
#include "act/forcefield/particletype.h"
#include "act/forcefield/potential.h"
#include "act/utility/units.h"
#include "generate_dependent.h"
#include "gromacs/math/units.h"
#include "gromacs/utility/exceptions.h"
#include "gromacs/utility/stringutil.h"

namespace alexandria
{

const char *kshellName = "kshell";

//! The parameters of a single particle type, from which pair parameters are made
struct SingleParameters
{
    //! The identifier
    const Identifier       *id;
    //! The parameters
    ForceFieldParameterMap *params;
};

/*! \brief The single particle types of a parameter list and their pairs
 * Finding the single particle types means going through all entries of
 * the list, which contains all pairs too, and sorting them. Since this is
 * done after each change of a parameter during training, the result is
 * kept, together with pointers to the pair entries. They stay valid as
 * long as entries are only added, see ForceFieldParameterList::structureVersion,
 * and the number of entries does not change.
 */
struct SingleParameterCache
{
    //! The list for which this was made
    const ForceFieldParameterList *list      = nullptr;
    //! Structure version of the list
    size_t                         version   = 0;
    //! Number of entries in the list
    size_t                         nentries  = 0;
    //! The single particle types, sorted by identifier
    std::vector<SingleParameters>  singles;
    //! The pair entries for i <= j at index i*singles.size() + j, nullptr if absent
    std::vector<ForceFieldParameterMap *> pairs;
    //! The particle types of the singles, if needed
    std::vector<const ParticleType *>    ptypes;
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

/*! \brief Find the single particle types of a list and their pairs
 * \param[in]    pd     The force field, for the particle types
 * \param[inout] list   The parameter list
 * \param[in]    ptypes Whether the particle types are needed
 * \return The singles and pairs, which stay valid until the next call
 *         for a different list or after the list changed
 */
static SingleParameterCache *singleParameterCache(const ForceField        *pd,
                                                  ForceFieldParameterList *list,
                                                  bool                     ptypes)
{
    // One for each interaction type with pair parameters
    static std::vector<SingleParameterCache> caches;
    auto parm = list->parametersNoRemoval();
    SingleParameterCache *cache = nullptr;
    for(auto &c : caches)
    {
        if (c.list == list)
        {
            cache = &c;
        }
    }
    if (nullptr == cache)
    {
        caches.emplace_back();
        cache       = &caches.back();
        cache->list = list;
    }
    if (cache->version == list->structureVersion() && cache->nentries == parm->size() &&
        !(ptypes && cache->ptypes.size() != cache->singles.size()))
    {
        return cache;
    }
    cache->version  = list->structureVersion();
    cache->nentries = parm->size();
    cache->singles  = singleParameters(parm);
    size_t n        = cache->singles.size();
    cache->pairs.assign(n*n, nullptr);
    for(size_t i = 0; i < n; i++)
    {
        for(size_t j = i; j < n; j++)
        {
            const auto myId = Identifier({ cache->singles[i].id->id(), cache->singles[j].id->id() }, { 1 }, CanSwap::Yes);
            const auto pp   = parm->find(myId);
            if (parm->end() != pp)
            {
                cache->pairs[i*n+j] = &pp->second;
            }
        }
    }
    cache->ptypes.clear();
    if (ptypes)
    {
        for(const auto &single : cache->singles)
        {
            const auto &atom = single.id->atoms()[0];
            cache->ptypes.push_back(pd->hasParticleType(atom) ? pd->findParticleType(atom) : nullptr);
        }
    }
    return cache;
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
    // Without a function the parameters are not used, for instance those
    // of a correction that gen_ff did not turn on.
    if (Potential::NONE == forcesVdw->potential())
    {
        return;
    }
    const auto &comb_rule = forcesVdw->combinationRules();
    // The new parameters are added to the same list
    ForceFieldParameterListMap *parm = forcesVdw->parametersNoRemoval();
    // With the Kronecker combination rule pairs are only included if one
    // particle is a vsite and the other an atom.
    bool kronecker = (Potential::BORN_MAYER == forcesVdw->potential() ||
                      Potential::MORSE_BONDS == forcesVdw->potential());

    // If the combination rule was updated all pairs change
    bool crupdated = false;
    for(const auto &cr : comb_rule)
    {
        crupdated = crupdated || cr.second.ffplConst().updated();
    }
    auto        cache   = singleParameterCache(pd, forcesVdw, kronecker);
    const auto &singles = cache->singles;
    size_t      nsingle = singles.size();
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
    const auto &ptypes = cache->ptypes;
    // The rules are looked up once, for the particle types that have the
    // parameters of the first that needs pairs, see CombinationPlan.
    std::unique_ptr<CombinationPlan> plan;
    std::vector<bool>                planned;

    // Now do the double loop
    int nid = 0;
    for (size_t i = 0; i < nsingle; i++)
    {
        for (size_t j = i; j < nsingle; j++)
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
            if (!plan)
            {
                plan = std::make_unique<CombinationPlan>(forcesVdw->potential(), comb_rule,
                                                         *singles[i].params);
                for(const auto &single : singles)
                {
                    planned.push_back(plan->size() > 0 && plan->matches(*single.params));
                }
            }
            auto pair = cache->pairs[i*nsingle+j];
            if (pair && planned[i] && planned[j] && !(msghandler && msghandler->debug()))
            {
                // Change the parameters of the pair in place
                plan->combine(*singles[i].params, *singles[j].params, includePair, pair);
                nid += 1;
                continue;
            }
            // Fill the parameters, potential dependent
            ForceFieldParameterMap pmap;
            evalCombinationRule(forcesVdw->potential(),
                                comb_rule, *singles[i].params, *singles[j].params, includePair, &pmap);
            if (!pmap.empty())
            {
                if (msghandler && msghandler->debug())
                {
                    const auto myId  = Identifier({ singles[i].id->id(), singles[j].id->id() }, { 1 }, CanSwap::Yes);
                    auto mystr = gmx::formatString("Generated pair parameters for %s / %s",
                                                   potentialToString(forcesVdw->potential()).c_str(),
                                                   myId.id().c_str());
                    for(const auto &p : pmap)
                    {
                        mystr += gmx::formatString(" %s: %g", p.first.c_str(), p.second.value());
                    }
                    msghandler->writeDebug(mystr);
                }
                if (pair)
                {
                    // The same as insert_or_assign, without looking up the identifier
                    *pair = std::move(pmap);
                }
                else
                {
                    auto myId = Identifier({ singles[i].id->id(), singles[j].id->id() }, { 1 }, CanSwap::Yes);
                    parm->insert_or_assign(std::move(myId), std::move(pmap));
                }
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

//! The name of the charge of a particle type
static const char *chargeName = "charge";

/*! \brief The parameters of a shell type with a fixed force constant
 * \param[in] pd    The force field
 * \param[in] id    The identifier of the polarization parameters
 * \param[in] pmap  The polarization parameters
 * \param[out] kshell The force constant
 * \param[out] alpha  The polarizability
 * \param[out] charge The charge of the shell
 * \return Whether the force constant is fixed
 * \throws gmx::InvalidInputError if the parameters are not consistent
 */
template <class Map, class FF, class Param>
static bool shellParameters(FF               *pd,
                            const Identifier &id,
                            Map              &pmap,
                            Param           **kshell,
                            Param           **alpha,
                            Param           **charge)
{
    const auto kp = pmap.find(kshellName);
    if (pmap.end() == kp || Mutability::Dependent == kp->second.mutability())
    {
        return false;
    }
    const auto ap = pmap.find(potentialToParameterName(Potential::POLARIZATION)[polALPHA]);
    if (pmap.end() == ap)
    {
        GMX_THROW(gmx::InvalidInputError(gmx::formatString("The polarization parameters of %s have a fixed %s but no alpha",
                                                           id.id().c_str(), kshellName)));
    }
    if (!pd->hasParticleType(id) || !pd->findParticleType(id)->hasParameter(chargeName))
    {
        GMX_THROW(gmx::InvalidInputError(gmx::formatString("The polarization parameters of %s have a fixed %s but there is no shell particle type %s with a charge",
                                                           id.id().c_str(), kshellName, id.id().c_str())));
    }
    *kshell = &kp->second;
    *alpha  = &ap->second;
    if constexpr (std::is_const_v<Param>)
    {
        *charge = &pd->findParticleType(id)->parameterConst(chargeName);
    }
    else
    {
        *charge = pd->findParticleType(id)->parameter(chargeName);
    }
    bool alphaDep  = Mutability::Dependent == (*alpha)->mutability();
    bool chargeDep = Mutability::Dependent == (*charge)->mutability();
    if (alphaDep == chargeDep)
    {
        GMX_THROW(gmx::InvalidInputError(gmx::formatString("The force constant %s of shell type %s is fixed (mutability %s), so either alpha or the charge of %s must have mutability Dependent, not %s",
                                                           kshellName, id.id().c_str(),
                                                           mutabilityName((*kshell)->mutability()).c_str(),
                                                           id.id().c_str(), alphaDep ? "both" : "neither")));
    }
    if ((*kshell)->internalValue() <= 0)
    {
        GMX_THROW(gmx::InvalidInputError(gmx::formatString("The force constant %s of shell type %s is %g %s, it must be positive",
                                                           kshellName, id.id().c_str(),
                                                           (*kshell)->value(), (*kshell)->unit().c_str())));
    }
    return true;
}

std::vector<ShellForceConstant> fixedShellForceConstants(const ForceField *pd)
{
    std::vector<ShellForceConstant> result;
    if (!pd->interactionPresent(InteractionType::POLARIZATION))
    {
        return result;
    }
    for (const auto &p : pd->findForcesConst(InteractionType::POLARIZATION).parametersConst())
    {
        const ForceFieldParameter *kshell = nullptr;
        const ForceFieldParameter *alpha  = nullptr;
        const ForceFieldParameter *charge = nullptr;
        if (shellParameters(pd, p.first, p.second, &kshell, &alpha, &charge))
        {
            ShellForceConstant s;
            s.shellType      = p.first.id();
            s.kshell         = kshell->value();
            s.kshellUnit     = kshell->unit();
            s.alpha          = alpha->value();
            s.alphaUnit      = alpha->unit();
            s.charge         = charge->value();
            s.alphaDependent = Mutability::Dependent == alpha->mutability();
            result.push_back(std::move(s));
        }
    }
    std::sort(result.begin(), result.end(),
              [](const ShellForceConstant &a, const ShellForceConstant &b)
              { return a.shellType < b.shellType; });
    return result;
}

/*! \brief Give a dependent parameter a new value, keeping its other properties
 * \param[inout] p     The parameter
 * \param[in]    value The value in the unit of the parameter
 * \return Whether the value changed
 */
static bool setDependentValue(ForceFieldParameter *p, double value)
{
    if (p->value() == value && p->minimum() == value && p->maximum() == value)
    {
        return false;
    }
    // Without bounds that the new value would violate
    p->setMaximum(value);
    (void) p->setMinimum(value);
    p->forceSetValue(value);
    return true;
}

ShellParametersChanged generateShellParameters(ForceField *pd)
{
    ShellParametersChanged changed;
    if (!pd->interactionPresent(InteractionType::POLARIZATION))
    {
        return changed;
    }
    for (auto &p : *pd->findForces(InteractionType::POLARIZATION)->parametersNoRemoval())
    {
        ForceFieldParameter *kshell = nullptr;
        ForceFieldParameter *alpha  = nullptr;
        ForceFieldParameter *charge = nullptr;
        if (!shellParameters(pd, p.first, p.second, &kshell, &alpha, &charge))
        {
            continue;
        }
        // In the units of the program, as in the force computation,
        // where the force constant is ONE_4PI_EPS0*q*q/alpha
        double k = kshell->internalValue();
        if (Mutability::Dependent == alpha->mutability())
        {
            double q        = convertToGromacs(charge->value(), charge->unit());
            double newAlpha = ONE_4PI_EPS0*q*q/k;
            changed.alpha   = setDependentValue(alpha, convertFromGromacs(newAlpha, alpha->unit())) || changed.alpha;
        }
        else
        {
            double a       = alpha->internalValue();
            double q       = std::sqrt(k*a/ONE_4PI_EPS0);
            double newQ    = charge->value() > 0 ? q : -q;
            changed.charge = setDependentValue(charge, convertFromGromacs(newQ, charge->unit())) || changed.charge;
        }
    }
    return changed;
}

std::vector<Identifier> shellTypes(const ForceField *pd)
{
    std::vector<Identifier> types;
    if (!pd->interactionPresent(InteractionType::POLARIZATION))
    {
        return types;
    }
    const auto alphaName = potentialToParameterName(Potential::POLARIZATION)[polALPHA];
    for (const auto &p : pd->findForcesConst(InteractionType::POLARIZATION).parametersConst())
    {
        if (p.second.count(alphaName) > 0 && pd->hasParticleType(p.first) &&
            ActParticle::Shell == pd->findParticleType(p.first)->apType() &&
            pd->findParticleType(p.first)->hasParameter(chargeName))
        {
            types.push_back(p.first);
        }
    }
    std::sort(types.begin(), types.end(),
              [](const Identifier &a, const Identifier &b) { return a.id() < b.id(); });
    return types;
}

/*! \brief The polarization parameters of a shell type
 * \param[in] pd        The force field
 * \param[in] shellType The shell type
 * \return The parameters
 * \throws gmx::InvalidInputError if it is not a shell type
 */
static ForceFieldParameterMap *shellPolarization(ForceField *pd, const Identifier &shellType)
{
    const auto types = shellTypes(pd);
    if (std::find(types.begin(), types.end(), shellType) == types.end())
    {
        GMX_THROW(gmx::InvalidInputError(gmx::formatString("%s is not a shell type with a polarizability",
                                                           shellType.id().c_str())));
    }
    return &pd->findForces(InteractionType::POLARIZATION)->parametersNoRemoval()->find(shellType)->second;
}

/*! \brief The force constant of the current polarizability and charge of a shell type
 * \param[in] alpha     The polarizability
 * \param[in] charge    The charge of the shell
 * \param[in] shellType The shell type, for the message
 * \return The force constant in kJ/mol nm2
 * \throws gmx::InvalidInputError if alpha or the charge is zero
 */
static double currentForceConstant(const ForceFieldParameter &alpha,
                                   const ForceFieldParameter &charge,
                                   const Identifier          &shellType)
{
    const double a = alpha.internalValue();
    const double q = convertToGromacs(charge.value(), charge.unit());
    if (a <= 0 || q == 0)
    {
        GMX_THROW(gmx::InvalidInputError(gmx::formatString("The force constant of shell type %s is not defined, since its polarizability is %g %s and its charge %g %s",
                                                           shellType.id().c_str(), alpha.value(), alpha.unit().c_str(),
                                                           charge.value(), charge.unit().c_str())));
    }
    return ONE_4PI_EPS0*q*q/a;
}

void fixShellForceConstant(ForceField       *pd,
                           const Identifier &shellType,
                           ShellFollows      follows,
                           double            kshell)
{
    auto *pmap   = shellPolarization(pd, shellType);
    auto *alpha  = &pmap->find(potentialToParameterName(Potential::POLARIZATION)[polALPHA])->second;
    auto *charge = pd->findParticleType(shellType)->parameter(chargeName);
    const double k = kshell > 0 ? kshell : currentForceConstant(*alpha, *charge, shellType);
    if (!std::isfinite(k))
    {
        GMX_THROW(gmx::InvalidInputError(gmx::formatString("The force constant %g of shell type %s is not a number",
                                                           k, shellType.id().c_str())));
    }
    // The unit of the program, such that no conversion is needed
    pmap->insert_or_assign(kshellName, ForceFieldParameter("kJ/mol nm2", k, 0, 1, k, k,
                                                           Mutability::Fixed, true, true));
    auto *follower = ShellFollows::Alpha == follows ? alpha : charge;
    auto *other    = ShellFollows::Alpha == follows ? charge : alpha;
    if (Mutability::Dependent == other->mutability())
    {
        other->setMutability(Mutability::Fixed);
    }
    follower->setMutability(Mutability::Dependent);
    (void) generateShellParameters(pd);
}

bool releaseShellForceConstant(ForceField       *pd,
                               const Identifier &shellType)
{
    auto *pmap   = shellPolarization(pd, shellType);
    const auto kp     = pmap->find(kshellName);
    const bool fixed = pmap->end() != kp && Mutability::Dependent != kp->second.mutability();
    auto *alpha  = &pmap->find(potentialToParameterName(Potential::POLARIZATION)[polALPHA])->second;
    auto *charge = pd->findParticleType(shellType)->parameter(chargeName);
    for (auto *p : { alpha, charge })
    {
        if (Mutability::Dependent == p->mutability())
        {
            p->setMutability(Mutability::Fixed);
        }
    }
    if (pmap->end() != kp)
    {
        // As it would be computed, if it could
        const double k = (alpha->internalValue() > 0 && charge->value() != 0) ?
            currentForceConstant(*alpha, *charge, shellType) : 0;
        kp->second.setDependent(convertFromGromacs(k, kp->second.unit()), k);
    }
    return fixed;
}

void generateDependentParameter(MsgHandler                      *msghandler,
                                ForceField                      *pd,
                                bool                             force,
                                const std::set<InteractionType> *changed)
{
    const auto todo = [force, changed](InteractionType itype)
          {
              return force || nullptr == changed || changed->count(itype) > 0;
          };
    // The shells first, since a charge may change. Particle types may
    // have changed if changed is nullptr.
    if (todo(InteractionType::POLARIZATION))
    {
        (void) generateShellParameters(pd);
    }
    for(auto itype : { InteractionType::VDW, InteractionType::VDWCORRECTION,
                       InteractionType::INDUCTIONCORRECTION })
    {
        if (todo(itype))
        {
            generateParameterPairs(msghandler, pd, itype, force);
        }
    }
}

} // namespace

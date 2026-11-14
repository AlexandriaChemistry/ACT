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
/*! \internal \brief
 * Compact storage of the non-bonded pairs of a topology, see nonbondedpairs.h.
 * \author David van der Spoel <david.vanderspoel@icm.uu.se>
 */
#include "nonbondedpairs.h"


#include <algorithm>
#include <cmath>
#include <numeric>
#include <string>

#include "act/basics/msg_handler.h"
#include "act/forcefield/forcefield.h"
#include "act/forcefield/forcefield_parametername.h"
#include "act/forcefield/potential.h"
#include "act/topology/actatom.h"
#include "gromacs/utility/exceptions.h"
#include "gromacs/utility/stringutil.h"

namespace alexandria
{

const std::vector<double> PairRange::noParams_;
const Identifier          PairRange::noId_;

void NonBondedPairs::make(size_t natom, const std::vector<std::set<size_t>> &exclusions)
{
    pairs_.clear();
    types_.clear();
    intermolecular_ = false;
    shells_         = false;
    // The same pairs in the same order as Topology::makePairs
    for (size_t i = 0; i < natom; i++)
    {
        for (size_t j = i+1; j < natom; j++)
        {
            if (exclusions[i].find(j) == exclusions[i].end())
            {
                pairs_.push_back({ static_cast<int>(i), static_cast<int>(j) });
            }
        }
    }
}

void NonBondedPairs::addType(const ForceField           *pd,
                             InteractionType             itype,
                             const std::vector<ActAtom> &atoms)
{
    if (pairs_.empty())
    {
        return;
    }
    TypePairs tp;
    tp.pos.resize(pairs_.size());
    std::iota(tp.pos.begin(), tp.pos.end(), 0);
    // As Topology::setEntryIdentifiers, which only gives identifiers if
    // the force field has the interaction type, the Coulomb interactions
    // never since they take the parameters from the particles
    tp.hasTypes = (pd && pd->interactionPresent(itype));
    if (tp.hasTypes)
    {
        auto       canSwap = pd->findForcesConst(itype).canSwap();
        // The type pair of each ordered pair of particle types
        std::map<std::pair<std::string, std::string>, uint32_t> typeIndex;
        tp.type.reserve(pairs_.size());
        for (const auto &p : pairs_)
        {
            std::pair<std::string, std::string> key(atoms[p[0]].ffType(), atoms[p[1]].ffType());
            auto found = typeIndex.find(key);
            if (found == typeIndex.end())
            {
                std::vector<std::string> btype;
                for (const auto &ffType : { key.first, key.second })
                {
                    auto atype = pd->findParticleType(ffType);
                    if (atype->hasInteractionType(itype))
                    {
                        btype.push_back(atype->interactionTypeToIdentifier(itype).id());
                    }
                }
                Identifier id;
                if (btype.size() == 1)
                {
                    id = Identifier(btype[0]);
                }
                else if (btype.size() == 2)
                {
                    // The fake bond order of AtomPair
                    id = Identifier(btype, { 1.0 }, canSwap);
                }
                found = typeIndex.insert({ key, tp.types.ids.size() }).first;
                tp.types.ids.push_back(id);
                tp.types.params.emplace_back();
            }
            tp.type.push_back(found->second);
        }
    }
    types_[itype] = std::move(tp);
}

/*! \brief The position of epsilon among the parameters of a potential whose energy is proportional to it
 * \param[in] pot The potential
 * \return The position, -1 if the energy does not vanish with epsilon, such
 *         as for LJ12_6_4 with its term in gamma, or the potential has none
 */
static int vanishingEpsilon(Potential pot)
{
    switch (pot)
    {
    case Potential::LJ12_6:                 return lj12_6EPSILON;
    case Potential::LJ8_6:                  return lj8_6EPSILON;
    case Potential::LJ14_7:                 return lj14_7EPSILON;
    case Potential::WANG_BUCKINGHAM:        return wbhEPSILON;
    case Potential::GENERALIZED_BUCKINGHAM: return gbhEPSILON;
    default:                                return -1;
    }
}

/*! \brief Whether a pair of particle types does not interact, and cannot during training
 * Its epsilon is zero, or not a number, as 0/0 of HogervorstEpsilon, and
 * the epsilons of both particle types are fixed, and so is the exponent of
 * GeneralizedMean, which could make it differ from zero otherwise.
 * \param[in] fs      The parameters of the interaction type
 * \param[in] pair    The identifier of the pair
 * \param[in] name    The name of epsilon
 * \param[in] epsilon The epsilon of the pair
 * \return whether the pair can be left out
 */
static bool fixedZeroEpsilon(const ForceFieldParameterList &fs,
                             const Identifier              &pair,
                             const char                    *name,
                             double                         epsilon)
{
    if (!(epsilon == 0 || std::isnan(epsilon)) || pair.atoms().size() != 2)
    {
        return false;
    }
    const auto &params = fs.parametersConst();
    for (const auto &atom : pair.atoms())
    {
        const auto single = params.find(Identifier(atom));
        if (params.end() == single)
        {
            return false;
        }
        const auto eps = single->second.find(name);
        if (single->second.end() == eps || eps->second.mutability() != Mutability::Fixed)
        {
            return false;
        }
    }
    const auto &rules = fs.combinationRules();
    const auto  rule  = rules.find(name);
    return !(rules.end() != rule && CombRule::GeneralizedMean == rule->second.rule() &&
             rule->second.ffplConst().mutability() != Mutability::Fixed);
}

bool NonBondedPairs::fillParameters(MsgHandler       *msghandler,
                                    const ForceField *pd,
                                    InteractionType   itype,
                                    missingParameters missing)
{
    const auto it = types_.find(itype);
    if (it == types_.end() || !it->second.hasTypes || !pd->interactionPresent(itype))
    {
        return false;
    }
    auto       &tp  = it->second;
    const auto &fs  = pd->findForcesConst(itype);
    const auto  pp  = potprops.find(fs.potential());
    if (potprops.end() == pp)
    {
        GMX_THROW(gmx::InternalError(gmx::formatString("Missing case %s when filling the non-bonded pairs.",
                                                       potentialToString(fs.potential()).c_str()).c_str()));
    }
    const auto &names = pp->second.param;
    bool        changed = false;
    size_t      ntype   = tp.types.ids.size();
    // The entries of the type pairs in the force field, found again only
    // when pointers into its parameter map may have become invalid. Entries
    // may be added without that, so those not found are looked for again.
    const auto &ffparams = fs.parametersConst();
    bool        refind   = tp.entriesVersion != fs.structureVersion() || tp.entries.size() != ntype;
    if (refind)
    {
        tp.entries.assign(ntype, nullptr);
        tp.entriesVersion = fs.structureVersion();
    }
    // Which type pairs can be left out, see fixedZeroEpsilon, decided again
    // when the entries are looked up
    int  epsilonIndex = vanishingEpsilon(fs.potential());
    bool decideDrop   = refind || tp.dropped.size() != ntype;
    if (decideDrop)
    {
        tp.dropped.assign(ntype, 0);
    }
    bool anyMissing = false, anyDropped = false;
    for (size_t t = 0; t < ntype; t++)
    {
        // As fillParams in topology.cpp, for each type pair rather than each pair
        auto &entry = tp.entries[t];
        if (nullptr == entry)
        {
            const auto ffiter = ffparams.find(tp.types.ids[t]);
            if (ffparams.end() != ffiter)
            {
                entry = &ffiter->second;
            }
        }
        auto &param = tp.types.params[t];
        if (nullptr == entry)
        {
            if (missingParameters::Ignore != missing)
            {
                msghandler->msg(ACTStatus::Error, ACTMessage::MissingFFParameter,
                                potentialToString(fs.potential()) + " " + tp.types.ids[t].id());
            }
            changed    = changed || !param.empty();
            param.clear();
            anyMissing = true;
            continue;
        }
        if (param.size() != names.size())
        {
            param.assign(names.size(), 0);
            changed = true;
        }
        size_t found = 0;
        for (size_t i = 0; i < names.size(); i++)
        {
            double value = 0;
            const auto fp    = entry->find(names[i]);
            if (entry->end() != fp)
            {
                value  = fp->second.internalValue();
                found += 1;
            }
            if (value != param[i])
            {
                param[i] = value;
                changed  = true;
            }
        }
        if (found != names.size())
        {
            msghandler->msg(ACTStatus::Warning, ACTMessage::MissingFFParameter,
                            gmx::formatString("%s %s found %zu, expected %zu independent parameters",
                                              potentialToString(fs.potential()).c_str(),
                                              tp.types.ids[t].id().c_str(), found, names.size()));
        }
        if (decideDrop && epsilonIndex >= 0)
        {
            tp.dropped[t] = fixedZeroEpsilon(fs, tp.types.ids[t], names[epsilonIndex], param[epsilonIndex]) ? 1 : 0;
        }
        if (tp.dropped[t])
        {
            // Without parameters, such that its pairs are left out below
            param.clear();
            anyDropped = true;
        }
    }
    // Pairs without parameters are left out, as are such TopologyEntry
    // objects, which also keeps the positions of the separations right.
    // Once all type pairs had parameters no pairs are left out any more.
    if (!anyMissing && !tp.mayRemove && !(anyDropped && !tp.droppedApplied))
    {
        return changed;
    }
    size_t n = 0, intermolecularStart = tp.intermolecularStart, shellStart = tp.shellStart;
    for (size_t k = 0; k < tp.pos.size(); k++)
    {
        if (!tp.types.params[tp.type[k]].empty())
        {
            tp.pos[n]  = tp.pos[k];
            tp.type[n] = tp.type[k];
            n         += 1;
        }
        else
        {
            changed              = true;
            intermolecularStart -= (k < tp.intermolecularStart) ? 1 : 0;
            shellStart          -= (k < tp.shellStart) ? 1 : 0;
        }
    }
    tp.pos.resize(n);
    tp.type.resize(n);
    tp.intermolecularStart = intermolecularStart;
    tp.shellStart          = shellStart;
    tp.mayRemove           = anyMissing;
    tp.droppedApplied      = tp.droppedApplied || anyDropped;
    return changed;
}

void NonBondedPairs::leaveOut(InteractionType itype, const std::vector<char> &particles)
{
    const auto it = types_.find(itype);
    if (it == types_.end())
    {
        return;
    }
    // As the pairs without parameters in fillParameters
    auto  &tp = it->second;
    size_t n  = 0, intermolecularStart = tp.intermolecularStart, shellStart = tp.shellStart;
    for (size_t k = 0; k < tp.pos.size(); k++)
    {
        const auto &p = pairs_[tp.pos[k]];
        if (particles[p[0]] || particles[p[1]])
        {
            intermolecularStart -= (k < tp.intermolecularStart) ? 1 : 0;
            shellStart          -= (k < tp.shellStart) ? 1 : 0;
            continue;
        }
        tp.pos[n] = tp.pos[k];
        if (!tp.type.empty())
        {
            tp.type[n] = tp.type[k];
        }
        n += 1;
    }
    tp.pos.resize(n);
    if (!tp.type.empty())
    {
        tp.type.resize(n);
    }
    tp.intermolecularStart = intermolecularStart;
    tp.shellStart          = shellStart;
}

void NonBondedPairs::separate(const std::function<bool(const std::array<int, 2> &)> &keepFirst,
                              bool                                                  shells)
{
    // The new position of each pair
    std::vector<uint32_t> order(pairs_.size());
    std::iota(order.begin(), order.end(), 0);
    const auto last   = std::stable_partition(order.begin(), order.end(),
                                              [this, &keepFirst](uint32_t k) { return keepFirst(pairs_[k]); });
    size_t nfirst = last - order.begin();
    std::vector<uint32_t>           newPos(pairs_.size());
    std::vector<std::array<int, 2>> pairs(pairs_.size());
    for (size_t k = 0; k < order.size(); k++)
    {
        newPos[order[k]] = k;
        pairs[k]         = pairs_[order[k]];
    }
    pairs_.swap(pairs);
    for (auto &t : types_)
    {
        auto &tp = t.second;
        // The pairs of the type in the new order, which is that of a
        // stable partition of them as well, since their positions increase
        std::vector<std::pair<uint32_t, uint32_t>> posType(tp.pos.size());
        for (size_t k = 0; k < tp.pos.size(); k++)
        {
            posType[k] = { newPos[tp.pos[k]], tp.hasTypes ? tp.type[k] : 0 };
        }
        std::sort(posType.begin(), posType.end());
        size_t start = 0;
        for (size_t k = 0; k < posType.size(); k++)
        {
            tp.pos[k] = posType[k].first;
            if (tp.hasTypes)
            {
                tp.type[k] = posType[k].second;
            }
            start += (posType[k].first < nfirst) ? 1 : 0;
        }
        if (shells)
        {
            tp.shellStart = start;
        }
        else
        {
            tp.intermolecularStart = start;
        }
    }
    if (shells)
    {
        // As Topology::separateShells, the order of the intermolecular pairs is lost
        shells_         = true;
        intermolecular_ = false;
    }
    else
    {
        intermolecular_ = true;
    }
}

void NonBondedPairs::setParameters(InteractionType itype, const std::vector<double> &params)
{
    const auto it = types_.find(itype);
    if (it == types_.end() || !it->second.hasTypes)
    {
        GMX_THROW(gmx::InternalError(gmx::formatString("No type pairs of %s",
                                                       interactionTypeToString(itype).c_str()).c_str()));
    }
    for (auto &p : it->second.types.params)
    {
        p = params;
    }
}

std::vector<InteractionType> NonBondedPairs::types() const
{
    std::vector<InteractionType> itypes;
    for (const auto &t : types_)
    {
        itypes.push_back(t.first);
    }
    return itypes;
}

PairRange NonBondedPairs::makeRange(const TypePairs &tp, size_t begin, size_t end) const
{
    return PairRange(pairs_.data(), tp.pos.data(), tp.hasTypes ? tp.type.data() : nullptr,
                     tp.hasTypes ? &tp.types : nullptr, begin, end);
}

PairRange NonBondedPairs::range(InteractionType itype) const
{
    const auto it = types_.find(itype);
    if (it == types_.end())
    {
        return PairRange();
    }
    return makeRange(it->second, 0, it->second.pos.size());
}

PairRange NonBondedPairs::intermolecular(InteractionType itype) const
{
    const auto it = types_.find(itype);
    if (it == types_.end())
    {
        return PairRange();
    }
    const auto &tp = it->second;
    return makeRange(tp, intermolecular_ ? tp.intermolecularStart : 0, tp.pos.size());
}

PairRange NonBondedPairs::shellEntries(InteractionType itype, bool shells) const
{
    const auto it = types_.find(itype);
    if (it == types_.end())
    {
        return PairRange();
    }
    const auto &tp    = it->second;
    // Without separation all pairs count as pairs with shells, as in Topology::shellEntries
    size_t      start = shells_ ? tp.shellStart : 0;
    return shells ? makeRange(tp, start, tp.pos.size()) : makeRange(tp, 0, start);
}

const PairTypes *NonBondedPairs::pairTypes(InteractionType itype) const
{
    const auto it = types_.find(itype);
    if (it == types_.end() || !it->second.hasTypes)
    {
        return nullptr;
    }
    return &it->second.types;
}

size_t NonBondedPairs::memory() const
{
    size_t bytes = pairs_.capacity()*sizeof(pairs_[0]);
    for (const auto &t : types_)
    {
        const auto &tp = t.second;
        bytes += (tp.pos.capacity() + tp.type.capacity())*sizeof(uint32_t);
        for (size_t k = 0; k < tp.types.ids.size(); k++)
        {
            bytes += sizeof(Identifier) + tp.types.params[k].capacity()*sizeof(double);
        }
    }
    return bytes;
}

} // namespace alexandria

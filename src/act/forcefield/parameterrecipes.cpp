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
#include "parameterrecipes.h"

#include <string>

#include "act/forcefield/combinationrules.h"
#include "act/forcefield/forcefield.h"
#include "act/forcefield/potential.h"
#include "act/utility/units.h"
#include "gromacs/utility/exceptions.h"
#include "gromacs/utility/stringutil.h"

namespace alexandria
{

int32_t ParameterRecipes::slot(const ForceFieldParameter *p)
{
    auto s = slots_.find(p);
    if (s != slots_.end())
    {
        return s->second;
    }
    auto index = static_cast<int32_t>(values_.size());
    slots_.insert({ p, index });
    params_.push_back(p);
    values_.push_back(p->value());
    // As ForceFieldParameter::internalValue
    factors_.push_back(convertToGromacs(1.0, p->unit()));
    return index;
}

int32_t ParameterRecipes::findSlot(const ForceFieldParameter *p) const
{
    auto s = slots_.find(p);
    return s == slots_.end() ? -1 : s->second;
}

int ParameterRecipes::directSet(const ForceField *pd, InteractionType itype, const Identifier &id)
{
    auto key = std::make_pair(itype, id);
    auto d   = direct_.find(key);
    if (d != direct_.end())
    {
        return d->second;
    }
    if (!pd->interactionPresent(itype))
    {
        return -1;
    }
    const auto &fs     = pd->findForcesConst(itype);
    const auto &params = fs.parametersConst();
    auto        entry  = params.find(id);
    if (params.end() == entry)
    {
        return -1;
    }
    // In the order of the potential, as Topology::fillParameters
    auto names = potentialToParameterName(fs.potential());
    GMX_RELEASE_ASSERT(names.size() <= maxSetParameters, "Too many parameters for a recipe");
    SetRecipe recipe;
    recipe.pair = 0;
    recipe.n    = static_cast<uint8_t>(names.size());
    for (size_t k = 0; k < names.size(); k++)
    {
        auto p           = entry->second.find(names[k]);
        recipe.slotI[k]  = entry->second.end() == p ? -1 : slot(&p->second);
    }
    recipes_.push_back(recipe);
    int index = static_cast<int>(recipes_.size()) - 1;
    direct_.insert({ key, index });
    return index;
}

int ParameterRecipes::pairSet(const ForceField *pd, InteractionType itype, const Identifier &pair)
{
    auto key = std::make_pair(itype, pair);
    auto d   = pairs_.find(key);
    if (d != pairs_.end())
    {
        return d->second;
    }
    const auto &atoms = pair.atoms();
    if (!pd->interactionPresent(itype) || atoms.size() != 2)
    {
        return -1;
    }
    const auto &fs     = pd->findForcesConst(itype);
    const auto &params = fs.parametersConst();
    auto        ti     = params.find(Identifier(atoms[0]));
    auto        tj     = params.find(Identifier(atoms[1]));
    if (params.end() == ti || params.end() == tj)
    {
        return -1;
    }
    SetRecipe recipe;
    recipe.pair = 1;
    // The values of both particle types in the order of the first, with
    // their rules, as evalCombinationRule
    std::vector<std::string> names;
    const auto              &rules = fs.combinationRules();
    for (const auto &value : ti->second)
    {
        auto vj   = tj->second.find(value.first);
        auto rule = rules.find(value.first);
        if (tj->second.end() == vj || rules.end() == rule)
        {
            GMX_THROW(gmx::InvalidInputError(gmx::formatString("Cannot combine %s of %s and %s",
                                                               value.first.c_str(), atoms[0].c_str(),
                                                               atoms[1].c_str()).c_str()));
        }
        int q = static_cast<int>(names.size());
        GMX_RELEASE_ASSERT(q < maxSetParameters, "Too many values for a recipe");
        names.push_back(value.first);
        recipe.slotI[q]    = slot(&value.second);
        recipe.slotJ[q]    = slot(&vj->second);
        recipe.rules[q]    = rule->second.rule();
        recipe.exponent[q] = slot(&rule->second.ffplConst());
    }
    recipe.nin = static_cast<uint8_t>(names.size());
    // The positions of sigma, epsilon and gamma among the values
    auto pc        = pairCombination(fs.potential());
    auto potnames  = potentialToParameterName(fs.potential());
    auto where     = [&names, &potnames](int index)
    {
        if (index < 0)
        {
            return -1;
        }
        for (size_t q = 0; q < names.size(); q++)
        {
            if (names[q] == potnames[index])
            {
                return static_cast<int>(q);
            }
        }
        return -1;
    };
    if (!pc.simple)
    {
        pc.dist    = where(pc.dist);
        pc.epsilon = where(pc.epsilon);
        pc.gamma   = where(pc.gamma);
    }
    recipe.pc = pc;
    // For Kronecker, only pairs of a virtual site and an atom interact, see
    // generateParameterPairs in generate_dependent.cpp
    bool kronecker = (Potential::BORN_MAYER == fs.potential() ||
                      Potential::MORSE_BONDS == fs.potential());
    recipe.include = 1;
    if (kronecker && pd->hasParticleType(atoms[0]) && pd->hasParticleType(atoms[1]))
    {
        auto pti = pd->findParticleType(atoms[0])->apType();
        auto ptj = pd->findParticleType(atoms[1])->apType();
        recipe.include = ((ActParticle::Vsite == pti && ActParticle::Atom == ptj) ||
                          (ActParticle::Vsite == ptj && ActParticle::Atom == pti)) ? 1 : 0;
    }
    // The parameters in the order of the potential, as NonBondedPairs::fillParameters
    GMX_RELEASE_ASSERT(potnames.size() <= maxSetParameters, "Too many parameters for a recipe");
    recipe.n = static_cast<uint8_t>(potnames.size());
    for (size_t k = 0; k < potnames.size(); k++)
    {
        recipe.order[k] = -1;
        for (size_t q = 0; q < names.size(); q++)
        {
            if (names[q] == potnames[k])
            {
                recipe.order[k] = static_cast<int8_t>(q);
            }
        }
    }
    recipes_.push_back(recipe);
    int index = static_cast<int>(recipes_.size()) - 1;
    pairs_.insert({ key, index });
    return index;
}

void ParameterRecipes::update()
{
    for (size_t s = 0; s < params_.size(); s++)
    {
        values_[s] = params_[s]->value();
    }
}

std::vector<double> ParameterRecipes::evaluate(int set) const
{
    const auto         &recipe = recipes_[set];
    std::vector<double> out(recipe.n, 0.0);
    auto status = evaluateRecipe<combination::StdMath>(recipe, values_.data(), factors_.data(), out.data());
    if (combination::PairStatus::Ok != status)
    {
        GMX_THROW(gmx::InvalidInputError("The combination rules cannot be applied to these parameters"));
    }
    return out;
}

} // namespace alexandria

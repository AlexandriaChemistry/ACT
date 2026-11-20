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
 * Implements functions that describe a force field in a JsonTree.
 *
 * \author David van der Spoel <david.vanderspoel@icm.uu.se>
 */
#include <algorithm>
#include <map>
#include <string>
#include <vector>

#include "act/basics/act_particle.h"
#include "act/basics/chargemodel.h"
#include "act/basics/identifier.h"
#include "act/basics/interactiontype.h"
#include "act/basics/mutability.h"
#include "act/forcefield/combinationrules.h"
#include "act/forcefield/forcefield.h"
#include "act/forcefield/forcefield_parameter.h"
#include "act/forcefield/forcefield_parameterlist.h"
#include "act/forcefield/forcefield_validation.h"
#include "act/forcefield/particletype.h"
#include "act/forcefield/potential.h"
#include "act/utility/jsontree.h"
#include "forcefield_json.h"

namespace alexandria
{

/*! \brief Make an object with options
 * \param[in] options The options
 * \return the tree, with the key options
 */
static JsonTree optionsToJson(const std::map<std::string, std::string> &options)
{
    // An object also without options
    JsonTree jopt = JsonTree::object("options");
    for (const auto &opt : options)
    {
        jopt.addObject(opt.first, opt.second);
    }
    return jopt;
}

JsonTree parameterToJson(const std::string         &type,
                         const ForceFieldParameter &param,
                         const std::string         &key)
{
    JsonTree jp(key);
    jp.addObject("type", type);
    jp.addObject("unit", param.unit());
    jp.addObject("value", param.value());
    jp.addObject("uncertainty", param.uncertainty());
    jp.addObject("minimum", param.minimum());
    jp.addObject("maximum", param.maximum());
    jp.addObject("ntrain", static_cast<int>(param.ntrain()));
    jp.addObject("mutability", mutabilityName(param.mutability()));
    jp.addBool("nonnegative", param.nonNegative());
    return jp;
}

JsonTree forceFieldToJson(const ForceField &pd, const std::string &key)
{
    JsonTree jff(key);
    jff.addObject("file", pd.filename());
    jff.addObject("checksum", pd.checkSum());
    jff.addObject("timestamp", pd.timeStamp());
    jff.addBool("polarizable", pd.polarizable());
    jff.addObject("charge_generation",
                  chargeGenerationAlgorithmName(pd.chargeGenerationAlgorithm()));
    jff.addObject("nexcl", std::max(0, pd.nexcl()));
    jff.addObject("epsilonr", pd.epsilonr());

    JsonTree jptypes("particletypes");
    jptypes.addObject(optionsToJson(pd.particleTypesOptions()));
    JsonTree jtypes = JsonTree::array("types");
    for (const auto &aType : pd.particleTypesConst())
    {
        JsonTree jt("particletype");
        jt.addObject("identifier", aType.second.id().id());
        jt.addObject("type", actParticleToString(aType.second.apType()));
        jt.addObject("description", aType.second.description());
        jt.addObject(optionsToJson(aType.second.optionsConst()));
        JsonTree jparams = JsonTree::array("parameters");
        for (const auto &param : aType.second.parametersConst())
        {
            jparams.addObject(parameterToJson(param.first, param.second));
        }
        jt.addObject(jparams);
        jtypes.addObject(jt);
    }
    jptypes.addObject(jtypes);
    jff.addObject(jptypes);

    JsonTree jinters = JsonTree::array("interactions");
    for (const auto &fs : pd.forcesConst())
    {
        JsonTree ji("interaction");
        ji.addObject("type", interactionTypeToString(fs.first));
        ji.addObject("function", potentialToString(fs.second.potential()));
        ji.addObject("canswap", canSwapToString(fs.second.canSwap()));
        ji.addObject(optionsToJson(fs.second.option()));
        JsonTree jrules = JsonTree::object("combination_rules");
        for (const auto &cr : fs.second.combinationRules())
        {
            JsonTree jr(cr.first);
            jr.addObject("rule", combinationRuleName(cr.second.rule()));
            jr.addObject(parameterToJson("exponent", cr.second.ffplConst(), "exponent"));
            jrules.addObject(jr);
        }
        ji.addObject(jrules);
        JsonTree jlists = JsonTree::array("parameterlists");
        for (const auto &params : fs.second.parametersConst())
        {
            // As in the force field file, entries that are derived from
            // others are not part of the description
            const bool dependent = std::all_of(params.second.begin(), params.second.end(),
                                               [](const auto &p)
                                               { return Mutability::Dependent == p.second.mutability(); });
            if (dependent)
            {
                continue;
            }
            JsonTree jl("parameterlist");
            jl.addObject("identifier", params.first.id());
            JsonTree jparams = JsonTree::array("parameters");
            for (const auto &param : params.second)
            {
                jparams.addObject(parameterToJson(param.first, param.second));
            }
            jl.addObject(jparams);
            jlists.addObject(jl);
        }
        ji.addObject(jlists);
        jinters.addObject(ji);
    }
    jff.addObject(jinters);

    JsonTree jsym = JsonTree::array("symmetric_charges");
    for (auto sc = pd.getSymchargesBegin(); sc != pd.getSymchargesEnd(); ++sc)
    {
        JsonTree js("symmetric_charge");
        js.addObject("central", sc->getCentral());
        js.addObject("attached", sc->getAttached());
        js.addObject("numattach", sc->getNumattach());
        jsym.addObject(js);
    }
    jff.addObject(jsym);
    return jff;
}

const char *forceFieldIssueKindName(ForceFieldIssue::Kind kind)
{
    switch (kind)
    {
    case ForceFieldIssue::Kind::ResetToBounds:
        return "ResetToBounds";
    case ForceFieldIssue::Kind::OutOfBounds:
        return "OutOfBounds";
    case ForceFieldIssue::Kind::MinimumAboveMaximum:
        return "MinimumAboveMaximum";
    case ForceFieldIssue::Kind::Negative:
        return "Negative";
    case ForceFieldIssue::Kind::Duplicate:
        return "Duplicate";
    }
    return "Unknown";
}

JsonTree forceFieldIssuesToJson(const std::vector<ForceFieldIssue> &issues,
                                const std::string                  &key)
{
    JsonTree jval(key);
    jval.addBool("ok", issues.empty());
    const auto nrepaired = std::count_if(issues.begin(), issues.end(),
                                         [](const ForceFieldIssue &i) { return i.repaired; });
    jval.addObject("repaired", static_cast<int>(nrepaired));
    JsonTree jissues = JsonTree::array("issues");
    for (const auto &issue : issues)
    {
        JsonTree ji("issue");
        ji.addObject("kind", forceFieldIssueKindName(issue.kind));
        ji.addObject("where", issue.where);
        ji.addObject("identifier", issue.identifier);
        ji.addObject("parameter", issue.parameter);
        ji.addObject("text", issue.text);
        ji.addBool("repaired", issue.repaired);
        jissues.addObject(ji);
    }
    jval.addObject(jissues);
    return jval;
}

} // namespace alexandria

/*
 * This source file is part of the Alexandria Chemistry Toolkit.
 *
 * Copyright (C) 2023-2025
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
 * \author David van der Spoel <david.vanderspoel@icm.uu.se>
 */

#include "combruleutil.h"

#include <cstring>

#include <map>
#include <string>
#include <vector>

#include "act/basics/interactiontype.h"
#include "act/forcefield/forcefield_parametername.h"
#include "act/forces/combinationrules.h"
#include "act/utility/stringutil.h"
#include "gromacs/commandline/pargs.h"
#include "gromacs/utility/fatalerror.h"

namespace alexandria
{

typedef struct {
    const char *flag;
    const char *var;
    const std::string deft;
    size_t index;
} cr_param;

const std::map<InteractionType, std::vector<cr_param> > mycr =
    {
        {
            InteractionType::VDW, {
                { "-cr_eps",  "epsilon", combinationRuleName(CombRule::Geometric), 0 },
                { "-cr_sig",  "sigma",   combinationRuleName(CombRule::Arithmetic), 1 },
                { "-cr_rmin", "rmin",    combinationRuleName(CombRule::Arithmetic), 2 },
                { "-cr_gam",  "gamma",   combinationRuleName(CombRule::Geometric), 3 },
                { "-cr_del",  "delta",   combinationRuleName(CombRule::Geometric), 4 },
                { "-cr_abh",  "Abh",     combinationRuleName(CombRule::Geometric), 5 },
                { "-cr_bbh",  "bbh",     combinationRuleName(CombRule::Arithmetic), 6 },
                { "-cr_c6bh", "c6bh",    combinationRuleName(CombRule::Geometric), 7 },
                { "-cr_ttA",  "Att",     combinationRuleName(CombRule::Geometric), 8 },
                { "-cr_ttB",  "btt",     combinationRuleName(CombRule::Arithmetic), 9 },
                { "-cr_ttC6", "c6tt",    combinationRuleName(CombRule::Geometric), 10 },
                { "-cr_ttC8", "c8tt",    combinationRuleName(CombRule::Geometric), 11 },
                { "-cr_ttC10","c10tt",   combinationRuleName(CombRule::Geometric), 12 },
                { "-cr_tt2bA",  "Att2b", combinationRuleName(CombRule::Geometric), 13 },
                { "-cr_tt2bBexch", "bExchtt2b", combinationRuleName(CombRule::Arithmetic), 14 },
                { "-cr_tt2bBdisp", "bDisptt2b", combinationRuleName(CombRule::Arithmetic), 15 },
                { "-cr_tt2bC6", "c6tt2b", combinationRuleName(CombRule::Geometric), 16 },
                { "-cr_tt2bC8", "c8tt2b", combinationRuleName(CombRule::Geometric), 17 },
                { "-cr_tt2bC10","c10tt2b", combinationRuleName(CombRule::Geometric), 18 }
            },
        },
        {
            InteractionType::VDWCORRECTION, {
                { "-cr_aexp",  "aexp", combinationRuleName(CombRule::Geometric), 19 },
                { "-cr_bexp",  "bexp", combinationRuleName(CombRule::Arithmetic), 20 }
            }
        },
        {
            InteractionType::INDUCTIONCORRECTION, {
                { "-cr_a1dexp",  "a1dexp", combinationRuleName(CombRule::Geometric), 21 },
                { "-cr_a2dexp",  "a2dexp", combinationRuleName(CombRule::Geometric), 22 },
                { "-cr_bdexp",  "bdexp",  combinationRuleName(CombRule::Arithmetic), 23 },
                { "-cr_De", "De", combinationRuleName(CombRule::Geometric), 24 },
                { "-cr_D0", "D0", combinationRuleName(CombRule::Geometric), 25 },
                { "-cr_beta", "beta", combinationRuleName(CombRule::Arithmetic), 26 },
                { "-cr_length", "bondlength", combinationRuleName(CombRule::Arithmetic), 27 },
            }
        }
    };

void CombRuleUtil::addInfo(std::vector<const char *> *crinfo)
{
    crinfo->push_back("[PAR]Combination rules can be specified for all Van der Waals parameters separately, depending");
    crinfo->push_back("on the potential chosen. You can choose from:[BR]");
    for (auto cr : combRuleName)
    {
        std::string newstr = gmx::formatString("%s ", cr.second.c_str());
        crinfo->push_back(strdup(newstr.c_str()));
    }
    crinfo->push_back("[PAR]Make sure to use the exact strings above including capitalization.");
    crinfo->push_back("Some of the rules that include parameter names should only be used for that parameter.");
    std::string last = gmx::formatString("If not specified, the %s rule will be selected.", combinationRuleName(CombRule::Geometric).c_str());
    crinfo->push_back(strdup(last.c_str()));
}

void CombRuleUtil::addPargs(std::vector<t_pargs> *pa)
{
    for(const auto &cr : mycr)
    {
        cr_flag_.resize(cr_flag_.size() + cr.second.size(), nullptr);
    }
    for(const auto &cr : mycr)
    {
        for (const auto &mm : cr.second)
        {
            auto desc = gmx::formatString("Combination rule to use for %s parameter %s",
                                          interactionTypeToString(cr.first).c_str(),
                                          mm.var);
            pa->push_back({ mm.flag, FALSE, etSTR, {&cr_flag_[mm.index]},
                    strdup(desc.c_str()) });
        }
    }
}

static bool hasParameter(const ForceFieldParameterList *ffpl,
                         const char                    *param)
{
    for (auto &ffp : ffpl->parametersConst())
    {
        for (auto &pp : ffp.second)
        {
            if (pp.first.compare(param) == 0)
            {
                return true;
            }
        }
    }
    return false;
}

int CombRuleUtil::extract(const std::vector<t_pargs> &pa,
                          ForceFieldParameterList    *vdw,
                          ForceFieldParameterList    *vdwcorr,
                          ForceFieldParameterList    *induccorr)
{
    // This will crash if there is no geometric combrule in the map...
    const char *defval = combRuleName.find(CombRule::Geometric)->second.c_str();
    int changed = 0;
    for(const auto &mcr : mycr)
    {
        for(const auto &mm : mcr.second)
        {
            const char *value;
            // If this has not been touched on the command line, do nothing.
            if (!opt2parg_bSet(mm.flag, pa.size(), pa.data()))
            {
                value = mm.deft.c_str();
            }
            else
            {
                value = cr_flag_[mm.index];
            }
            if (!value || strlen(value) == 0)
            {
                // Check whether we have existing values
                if ((mcr.first == InteractionType::VDW && vdw && 
                     vdw->combinationRuleExists(mm.var)) ||
                    (mcr.first == InteractionType::VDWCORRECTION && vdwcorr && 
                     vdwcorr->combinationRuleExists(mm.var)) ||
                    (mcr.first == InteractionType::INDUCTIONCORRECTION && induccorr && 
                     induccorr->combinationRuleExists(mm.var)))
                {
                    continue;
                }
                value = defval;
            }
            // Will throw if incorrect string
            CombRule cr;
            if (!combinationRuleRule(value, &cr))
            {
                GMX_THROW(gmx::InvalidInputError(gmx::formatString("Invalid combination rule name %s for parameter %s", value, mm.var).c_str()));
            }
            bool doChange = true;
            if (mcr.first == InteractionType::VDW)
            {
                if (vdw)
                {
                    if (vdw->combinationRuleExists(mm.var))
                    {
                        auto oldRule = vdw->combinationRule(mm.var);
                        if (oldRule != value)
                        {
                            printf("Changing combination rule for %s from %s to %s\n",
                                   mm.var, oldRule.c_str(), value);
                        }
                        else
                        {
                            doChange = false;
                        }
                    }
                    if (doChange)
                    {
                        if (hasParameter(vdw, mm.var))
                        {
                            vdw->addCombinationRule(mm.var, value);
                            changed += 1;
                        }
                    }
                }
            }
            else if (mcr.first == InteractionType::VDWCORRECTION)
            {
                if (vdwcorr && hasParameter(vdwcorr, mm.var))
                {
                    vdwcorr->addCombinationRule(mm.var, value);
                }
            }
            else if (mcr.first == InteractionType::INDUCTIONCORRECTION)
            {
                if (induccorr && hasParameter(induccorr, mm.var))
                {
                    induccorr->addCombinationRule(mm.var, value);
                }
            }
        }
    }
    return changed;
}

int CombRuleUtil::convert(ForceFieldParameterList *vdw)
{
    int changed = 0;
    if (vdw)
    {
        std::string crule("combination_rule");
        if (vdw->optionExists(crule))
        {
            for(auto &crule : getCombinationRule(*vdw))
            {
                vdw->addCombinationRule(crule.first, combinationRuleName(crule.second));
                changed += 1;
            }
        }
    }
        
    return changed;
}

} // namespace alexandria

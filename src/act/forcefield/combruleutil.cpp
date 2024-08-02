/*
 * This source file is part of the Alexandria Chemistry Toolkit.
 *
 * Copyright (C) 2023,2024
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
#include "act/forces/combinationrules.h"
#include "act/utility/stringutil.h"
#include "gromacs/commandline/pargs.h"
#include "gromacs/utility/fatalerror.h"

namespace alexandria
{

typedef struct {
    const char *flag;
    const char *var;
    size_t index;
} cr_param;

const std::map<InteractionType, std::vector<cr_param> > mycr =
    {
        {
            InteractionType::VDW, {
                { "-cr_eps",  "epsilon", 0 },
                { "-cr_sig",  "sigma",   1 },
                { "-cr_rmin", "rmin",    2 },
                { "-cr_gam",  "gamma",   3 },
                { "-cr_del",  "delta",   4 }
            },
        },
        {
            InteractionType::VDWCORRECTION, {
                { "-cr_aexp",  "aexp", 5 },
                { "-cr_bexp",  "bexp", 6 }
            }
        },
        {
            InteractionType::INDUCTIONCORRECTION, {
                { "-cr_a1dexp",  "a1dexp", 7 },
                { "-cr_a2dexp",  "a2dexp", 8 },
                { "-cr_bdexp",  "bdexp",  9 },
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
        cr_flag_.resize(cr_flag_.size() + cr.second.size());
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
            // If this has not been touched on the command line, do nothing.
            if (!opt2parg_bSet(mm.flag, pa.size(), pa.data()))
            {
                continue;
            }
            auto value = cr_flag_[mm.index];
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
                        vdw->addCombinationRule(mm.var, value);
                        changed += 1;
                    }
                }
            }
            else if (mcr.first == InteractionType::VDWCORRECTION)
            {
                if (vdwcorr)
                {
                    vdwcorr->addCombinationRule(mm.var, value);
                }
            }
            else if (mcr.first == InteractionType::INDUCTIONCORRECTION)
            {
                if (induccorr)
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

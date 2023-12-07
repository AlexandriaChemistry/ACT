/*
 * This source file is part of the Alexandria Chemistry Toolkit.
 *
 * Copyright (C) 2023
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
#include <vector>

#include "act/forces/combinationrules.h"
#include "act/utility/stringutil.h"
#include "gromacs/commandline/pargs.h"
#include "gromacs/utility/fatalerror.h"

namespace alexandria
{

const std::map<const char *, const char *> mycr = {
    { "-cr_eps",  "epsilon" },
    { "-cr_sig",  "sigma"   },
    { "-cr_rmin", "rmin"    },
    { "-cr_gam",  "gamma"   },
    { "-cr_del",  "delta"   }
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
}

void CombRuleUtil::addPargs(std::vector<t_pargs> *pa)
{
    cr_flag_.resize(mycr.size());
    desc_.resize(mycr.size());
    size_t i = 0;
    for (const auto &mm : mycr)
    {
        desc_[i] = gmx::formatString("Combination rule to use for Van der Waals interaction parameter %s",
                                     mm.second);
        t_pargs mp = { mm.first, FALSE, etSTR, {&cr_flag_[i]}, desc_[i].c_str() };
        pa->push_back(mp);
        i += 1;
    }
}

void CombRuleUtil::extract(ForceFieldParameterList *vdw)
{
    size_t i = 0;
    for(const auto &mm : mycr)
    {
        if (strlen(cr_flag_[i]) > 0)
        {
            // Will throw if incorrect string
            CombRule cr;
            if (!combinationRuleRule(cr_flag_[i], &cr))
            {
                GMX_THROW(gmx::InvalidInputError(gmx::formatString("Invalid combination rule name %s for parameter %s", cr_flag_[i], mm.second).c_str()));
            }
            vdw->addCombinationRule(mm.first, cr_flag_[i]);
        }
        i += 1;
    }
}

} // namespace alexandria

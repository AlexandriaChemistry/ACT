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
/*! \internal \brief
 * Implements part of the alexandria program.
 * \author Mohammad Mehdi Ghahremanpour <mohammad.ghahremanpour@icm.uu.se>
 * \author David van der Spoel <david.vanderspoel@icm.uu.se>
 * \author Julian Ramon Marrades Furquet <julian@marrad.es>
 */
#include "sensitivity.h"

#include "act/basics/msg_handler.h"
#include "act/utility/jsontree.h"
#include "act/alexandria/acmfitnesscomputer.h"

namespace alexandria
{

bool SensitivityAnalysis::run(MsgHandler           *msghandler,
                              StaticIndividualInfo *sii,
                              ACMFitnessComputer   *fitComp,
                              ga::Genome           *genome,
                              iMolSelect            ims,
                              JsonTree             *jtree,
                              bool                  quiet)
{
    std::vector<double> *param = genome->basesPtr();
    const auto upperBound      = sii->upperBound();
    const auto lowerBound      = sii->lowerBound();
    const auto paramNames      = sii->paramNames();

    if (param->size() == 0)
    {
        return true;
    }
    std::set<int> changed;
    sii->updateForceField(msghandler, changed, genome->bases());
    auto cdc    = CalcDev::Compute;
    fitComp->distributeTasks(cdc);
    auto chi2_0 = fitComp->calcDeviation(msghandler, cdc, ims);
    auto tw = msghandler->tw();
    if (quiet)
    {
        tw = nullptr;
    }
    if (tw)
    {
        tw->writeStringFormatted("\nStarting sensitivity analysis. chi2_0 = %g nParam = %zu\n",
                                 chi2_0, param->size());
    }
    JsonTree sens("sensitivity");
    bool minimum = true;
    // Reset force constants every run
    forceConstant_.clear();
    for (size_t i = 0; i < param->size(); ++i)
    {
        Sensitivity s;
        double pstore = (*param)[i];
        double deltap = (upperBound[i]-lowerBound[i])/200;
        double pmin   = std::max((*param)[i]-deltap, lowerBound[i]);
        double pmax   = std::min((*param)[i]+deltap, upperBound[i]);
        double p_0    = 0.5*(pmin+pmax);
        std::set<int> changed;
        changed.insert(i);
        (*param)[i]     = pmin;
        sii->updateForceField(msghandler, changed, *param);
        fitComp->distributeTasks(cdc);
        s.add((*param)[i], fitComp->calcDeviation(msghandler, cdc, ims));
        (*param)[i]     = p_0;
        sii->updateForceField(msghandler, changed, *param);
        fitComp->distributeTasks(cdc);
        s.add((*param)[i], fitComp->calcDeviation(msghandler, cdc, ims));
        (*param)[i]     = pmax;
        sii->updateForceField(msghandler, changed, *param);
        fitComp->distributeTasks(cdc);
        s.add((*param)[i],  fitComp->calcDeviation(msghandler, cdc, ims));
        (*param)[i]     = pstore;
        sii->updateForceField(msghandler, changed, *param);
        s.computeForceConstants(tw);
        minimum = minimum && s.a() >= 0;
        // Compute dimensionless force constant, or zero if undefined
        double fc = 0;
        if (pstore != 0)
        {
            fc = s.a() / (pstore*pstore);
        }
        forceConstant_.push_back(fc);
        s.print(tw, &sens, std::to_string(i), paramNames[i]);
    }
    if (tw)
    {
        tw->writeString("Sensitivity analysis done.");
    }
    if (jtree)
    {
        jtree->addObject(sens);
    }
    return minimum;
}                                      

}

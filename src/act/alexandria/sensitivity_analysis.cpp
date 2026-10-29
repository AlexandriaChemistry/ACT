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
#include "sensitivity_analysis.h"

#include "act/basics/msg_handler.h"
#include "act/utility/communicationrecord.h"
#include "act/utility/jsontree.h"
#include "act/alexandria/acmfitnesscomputer.h"
#include "act/alexandria/staticindividualinfo.h"
#include "act/ga/genome.h"
#include "act/utility/timings.h"

namespace alexandria
{

bool SensitivityAnalysis::run(MsgHandler           *msghandler,
                              StaticIndividualInfo *sii,
                              ACMFitnessComputer   *fitComp,
                              ga::Genome           *genome,
                              iMolSelect            ims,
                              JsonTree             *jtree,
                              bool                  quiet,
                              bool                  useMiddlemen)
{
    static const size_t timer = Timings::instance().index("Sensitivity analysis");
    ScopedTimer         scoped(timer);
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
    // Each parameter is changed to three values, the chi2 of all of them
    // is computed first. Each value is computed on its own, by the master
    // or by one of the middlemen, which only compute the training set.
    auto                cr       = sii->commRec();
    std::vector<int>    workers;
    if (useMiddlemen && iMolSelect::Train == ims && cr->isMaster())
    {
        workers = cr->middlemen();
    }
    size_t              np       = param->size();
    std::vector<double> pvalue(3*np), chi2(3*np);
    for (size_t i = 0; i < np; ++i)
    {
        double deltap   = (upperBound[i]-lowerBound[i])/200;
        double pmin     = std::max((*param)[i]-deltap, lowerBound[i]);
        double pmax     = std::min((*param)[i]+deltap, upperBound[i]);
        pvalue[3*i]     = pmin;
        pvalue[3*i + 1] = 0.5*(pmin+pmax);
        pvalue[3*i + 2] = pmax;
    }
    for (size_t first = 0; first < pvalue.size(); first += workers.size() + 1)
    {
        // The middlemen get the values after the first one
        size_t nsent = std::min(workers.size(), pvalue.size() - first - 1);
        for (size_t w = 0; w < nsent; w++)
        {
            size_t     k = first + 1 + w;
            ga::Genome g = *genome;
            (*g.basesPtr())[k/3] = pvalue[k];
            cr->send_data(workers[w]);
            g.Send(cr, workers[w]);
            cr->send(workers[w], TrainFFMiddlemanMode::FITNESS);
        }
        // The master computes the first one
        size_t        i      = first/3;
        double        pstore = (*param)[i];
        std::set<int> changed = { static_cast<int>(i) };
        (*param)[i] = pvalue[first];
        sii->updateForceField(msghandler, changed, *param);
        fitComp->distributeTasks(cdc);
        chi2[first] = fitComp->calcDeviation(msghandler, cdc, ims);
        (*param)[i] = pstore;
        sii->updateForceField(msghandler, changed, *param);
        for (size_t w = 0; w < nsent; w++)
        {
            cr->recv(workers[w], &chi2[first + 1 + w]);
        }
    }
    JsonTree sens("sensitivity");
    bool minimum = true;
    // Reset force constants every run
    forceConstant_.clear();
    for (size_t i = 0; i < np; ++i)
    {
        Sensitivity s;
        for (size_t k = 3*i; k < 3*i + 3; k++)
        {
            s.add(pvalue[k], chi2[k]);
        }
        double pstore = (*param)[i];
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

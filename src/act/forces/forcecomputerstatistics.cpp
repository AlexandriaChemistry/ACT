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
#include <string>

#include "forcecomputerstatistics.h"

#include "act/forces/forcecomputer.h"
#include "act/utility/communicationrecord.h"

namespace alexandria
{

static void makeTheString(size_t       totalEval,
                          size_t       totalShellIter,
                          size_t       totalShellFailed,
                          bool         aggregate,
                          std::string *stats)
{
    double ratio = 0;
    if (totalEval > 0)
    {
        ratio = (totalShellIter*1.0)/totalEval;
    }
    stats->assign(gmx::formatString("ForceComputer %s#eval %zu #shell iter %zu ratio %.2f #failed shell convergence %zu",
                                    aggregate ? "total " : "",
                                    totalEval, totalShellIter, ratio, totalShellFailed));
}

std::string ForceComputerStatistics::statistics(const CommunicationRecord *cr,
                                                const ForceComputer       *forceComp,
                                                int                        nElites,
                                                bool                       aggregate)
{
    collect(cr, forceComp, nElites);
    return report(aggregate);
}

void ForceComputerStatistics::collect(const CommunicationRecord *cr,
                                      const ForceComputer       *forceComp,
                                      int                        nElites)
{
    size_t totalEval        = forceComp->numEvaluations();
    size_t totalShellIter   = forceComp->numShellIterations();
    size_t totalShellFailed = forceComp->numShellConvergenceFailed();
    forceComp->resetStatistics();
    master_ = cr->isMaster();
    // Now communicate
    if (master_)
    {
        // Middlemen does not include the master, hence <=
        for(size_t i = std::max(1, nElites); i <= cr->middlemen().size(); i++)
        {
            auto src = cr->middlemen()[i-1];
            size_t nEval = 0, nShellIter = 0, nShellFailed = 0;
            cr->recv(src, &nEval);
            cr->recv(src, &nShellIter);
            cr->recv(src, &nShellFailed);
            totalEval        += nEval;
            totalShellIter   += nShellIter;
            totalShellFailed += nShellFailed;
        }
        // Update totals
        totalEval_                += totalEval;
        totalShellIter_           += totalShellIter;
        totalShellConvergeFailed_ += totalShellFailed;
        lastEval_                  = totalEval;
        lastShellIter_             = totalShellIter;
        lastShellConvergeFailed_   = totalShellFailed;
    }
    else
    {
        cr->send(0, totalEval);
        cr->send(cr->superior(), totalShellIter);
        cr->send(cr->superior(), totalShellFailed);
    }
}

std::string ForceComputerStatistics::report(bool aggregate) const
{
    std::string stats;
    if (master_)
    {
        if (aggregate)
        {
            makeTheString(totalEval_, totalShellIter_, totalShellConvergeFailed_, aggregate, &stats);
        }
        else
        {
            makeTheString(lastEval_, lastShellIter_, lastShellConvergeFailed_, aggregate, &stats);
        }
    }
    return stats;
}

}

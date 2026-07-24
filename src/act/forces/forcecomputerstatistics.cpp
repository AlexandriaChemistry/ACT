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

static std::string makeTheString(size_t totalEval,
                                 size_t totalShellIter,
                                 size_t totalShellFailed,
                                 bool   aggregate)
{
    double ratio = 0;
    if (totalEval > 0)
    {
        ratio = (totalShellIter*1.0)/totalEval;
    }
    return gmx::formatString("ForceComputer %s#eval %zu #shell iter %zu ratio %.2f #failed shell convergence %zu",
                             aggregate ? "total " : "",
                             totalEval, totalShellIter, ratio, totalShellFailed);
}

std::string ForceComputerStatistics::statistics(const CommunicationRecord *cr,
                                                const ForceComputer       *forceComp,
                                                int                        nElites,
                                                bool                       aggregate)
{
    size_t totalEval        = forceComp->numEvaluations();
    size_t totalShellIter   = forceComp->numShellIterations();
    size_t totalShellFailed = forceComp->numShellConvergenceFailed();
    forceComp->resetStatistics();
    std::string stats;
    // Now communicate
    if (cr->isMaster())
    {
        for(size_t i = std::max(1, nElites); i < cr->middlemen().size(); i++)
        {
            auto src = cr->middlemen()[i-1];
            size_t nEval, nShellIter, nShellFailed;
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
        if (aggregate)
        {
            return makeTheString(totalEval_, totalShellIter_, totalShellConvergeFailed_, aggregate);
        }
        else
        {
            return makeTheString(totalEval, totalShellIter, totalShellFailed, aggregate);
        }
    }
    else
    {
        cr->send(cr->superior(), totalEval);
        cr->send(cr->superior(), totalShellIter);
        cr->send(cr->superior(), totalShellFailed);
    }
    return stats;
}

}

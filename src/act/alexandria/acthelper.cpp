/*
 * This source file is part of the Alexandria Chemistry Toolkit.
 *
 * Copyright (C) 2014-2025
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
#include "acthelper.h"

#include <vector>
    
#include "acmfitnesscomputer.h"
#include "bayes.h"
#include "act/basics/dataset.h"

namespace alexandria
{
    
ACTHelper::ACTHelper(MsgHandler                *msghandler,
                     StaticIndividualInfo      *sii,
                     MolGen                    *mg,
                     double                     shellToler,
                     int                        shellMaxIter,
                     ChargeGenerationAlgorithm  algorithm)
{
    forceComp_ = new ForceComputer(shellToler, shellMaxIter);
    fitComp_   = new ACMFitnessComputer(msghandler, sii, mg, false, forceComp_, algorithm);
}

void ACTHelper::run(MsgHandler *msghandler)
{
    // H E L P E R   N O D E
    CalcDev cd = CalcDev::Compute;
    cd = fitComp_->distributeTasks(cd);
    while (CalcDev::Stop != cd)
    {
        if (debug)
        {
            fprintf(debug, "Received CalcDev %s\n", calcDevName(cd));
        }
        switch (cd)
        {
        case CalcDev::Parameters:
            {
                // All variables are set by the master or middlemen, but
                // we have to pass something.
                std::vector<double> dummy;
                std::set<int>       changed;
                // Get the new parameters to evaluate
                fitComp_->distributeParameters(msghandler, &dummy, changed);
            }
            break;
        case CalcDev::Compute:
            // Evaluate on my part of the dataset. The second flag
            // is not used inside the routine for helper, instead a
            // new dataset is distributed from the master or middlemen.
            (void) fitComp_->calcDeviation(msghandler, cd, iMolSelect::Train);
            break;
        case CalcDev::ComputeAll:
        case CalcDev::Stop:
            gmx_fatal(FARGS, "Incorrect calcdev %s in helper", calcDevName(cd));
            // This should never happen, but without the compiler will complain.
            break;
        }
        cd = fitComp_->distributeTasks(cd);
    }
}

} // namespace alexandria

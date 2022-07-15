/*
 * This source file is part of the Alexandria Chemistry Toolkit.
 *
 * Copyright (C) 2014-2022
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
    
ACTHelper::ACTHelper(StaticIndividualInfo *sii,
                     MolGen               *mg)
{
    forceComp_ = new ForceComputer(sii->poldata());
    fitComp_   = new ACMFitnessComputer(nullptr, false, sii, mg, false, forceComp_);
}

void ACTHelper::run()
{
    // H E L P E R   N O D E
    // All variables are set by the master or middlemen, but
    // we have to pass something.
    // If the result is less than zero (-1), we are done.
    std::vector<double> dummy;
    while (fitComp_->calcDeviation(&dummy, CalcDev::Parallel, iMolSelect::Train) >= 0)
    {
        ;
    }
}

} // namespace alexandria

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
/*! \internal \brief
 * Implements part of the alexandria program.
 * \author Julian Ramon Marrades Furquet <julian@marrad.es>
 */


#include "percentmutator.h"

#include "acmindividual.h"


namespace alexandria
{

void PercentMutator::mutate(ga::Genome *genome,
                            gmx_unused ga::Genome *bestGenome,
                            double      prMut)
{
    const std::vector<double> oldParam = genome->bases();
    const std::vector<double> lb = sii_->lowerBound();
    const std::vector<double> ub = sii_->upperBound();
    double newVal;

    for (size_t i = 0; i < genome->nBase(); i++)
    {
        if (randNum() <= prMut)
        {
            newVal = oldParam[i] + percent_*(2*randNum()-1)*(ub[i] - lb[i]);
            newVal = std::max(lb[i], std::min(ub[i], newVal));
            genome->setBase(i, newVal);
        }
    }
}


} //namespace alexandria

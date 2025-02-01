/*
 * This source file is part of the Alexandria Chemistry Toolkit.
 *
 * Copyright (C) 2021-2025
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


#include "terminator.h"

#include "gromacs/utility/basedefinitions.h"
#include "gromacs/utility/textwriter.h"

#include "gene_pool.h"

namespace ga
{


/* * * * * * * * * * * * * * * * * * * * * *
* BEGIN: GenerationTerminator              *
* * * * * * * * * * * * * * * * * * * * * */

bool GenerationTerminator::terminate(gmx::TextWriter *tw,
                                     const GenePool  *,
                                     const int        generationNumber)
{
    if (generationNumber >= maxGenerations_)
    {
        if (tw)
        {
            tw->writeStringFormatted("GenerationTerminator: evolution will be terminated as the maximum number of generations (%d) has been reached.\n",
                                     maxGenerations_);
        }
        return true;
    }
    else
    {
        return false;
    }
}

/* * * * * * * * * * * * * * * * * * * * * *
* END: GenerationTerminator                *
* * * * * * * * * * * * * * * * * * * * * */

/* * * * * * * * * * * * * * * * * * * * * *
* BEGIN: TestGenTerminator                 *
* * * * * * * * * * * * * * * * * * * * * */

bool TestGenTerminator::terminate(gmx::TextWriter *tw,
                                  const GenePool  *pool,
                                  const int        )
{
    remaining_ -= 1;
    
    // Get best genome in population
    const auto bestGenome = pool->getBest(iMolSelect::Train);
    // Check if its test fitness is better than the best one so far
    const double tmpTestFit = bestGenome.fitness(iMolSelect::Test);
    if (tmpTestFit < bestFitness_)
    {
        bestFitness_ = tmpTestFit;
        remaining_ = generations_;
    }

    if (remaining_ <= 0)
    {
        if (tw)
        {
            tw->writeStringFormatted("TestGenTerminator: evolution will be terminated as the best test fitness (%lf) has not improved in the last (%d) generation(s).\n",
                                     bestFitness_,
                                     generations_);
        }
        return true;
    }
    else
    {
        return false;
    }
}

/* * * * * * * * * * * * * * * * * * * * * *
* END: TestGenTerminator                   *
* * * * * * * * * * * * * * * * * * * * * */


} //namespace ga

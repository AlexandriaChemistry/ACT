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
 * \author Julian Ramon Marrades Furquet <julian.marrades@hotmail.es>
 */
#include <algorithm>

#include "npointcrossover.h"

namespace ga
{


void NPointCrossover::offspring(const ga::Genome  *parent1,
                                const ga::Genome  *parent2,
                                ga::Genome        *child1,
                                ga::Genome        *child2)
{

    // Iteration variable(s)
    size_t i;
    size_t j;

    // Fill in the crossover points
    // We start by shuffling the available indices
    std::shuffle(availableIndices_.begin(), availableIndices_.end(), gen);
    // Now we copy the first order_ elements to crossoverPoints_
    for (i = 0; i < order_; i++)
    {
        crossoverPoints_[i+1] = availableIndices_[i];
    }
    // Finally, sort the newly added crossover points
    std::sort(crossoverPoints_.begin()+1, crossoverPoints_.end()-1);
    // DONE! Sampled without replacement!

    // We now cross the genes of the genomes
    // Start by the regions that are not swapped
    for (i = 0; i < crossoverPoints_.size() - 1; i += 2)
    {
        for (j = crossoverPoints_[i]; j < crossoverPoints_[i+1]; j++)
        {
            child1->setBase(j, parent1->base(j));
            child2->setBase(j, parent2->base(j));
        }
    }
    // Now the regions that should be swapped
    for (i = 1; i < crossoverPoints_.size() - 1; i += 2)
    {
        for (j = crossoverPoints_[i]; j < crossoverPoints_[i+1]; j++)
        {
            child1->setBase(j, parent2->base(j));
            child2->setBase(j, parent1->base(j));
        }
    }
}


} //namespace ga

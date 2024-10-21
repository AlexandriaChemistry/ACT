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


#ifndef GA_NPOINTCROSSOVER_H
#define GA_NPOINTCROSSOVER_H


#include "crossover.h"


namespace ga
{


/*!
 * \brief Performs N-Point Crossover operation
 * See <a href="https://en.wikipedia.org/wiki/Crossover_(genetic_algorithm)#One-point_crossover">this</a> for details.
 *
 */
class NPointCrossover : public Crossover
{

private:

    // Random number stuff (for shuffling)
    std::random_device  rd;
    std::mt19937        gen;

    //! Amount of crossover cutting points
    size_t order_;
    //! List of all existing indices
    std::vector<size_t> availableIndices_;
    //! List of all crossover points. First element is 0, and last element is \p chromosomeLength_.
    std::vector<size_t> crossoverPoints_;

public:

    /*!
     * Property constructor
     * \param[in] chromosomeLength  the length of the chromosomes
     * \param[in] order             order of the crossover operator (amount of cutting points)
     * \param[in] seed              seed for the random number generator
     */
    NPointCrossover(const size_t chromosomeLength,
                    const size_t order,
                    const int seed)
    : ga::Crossover(chromosomeLength, seed), gen(rd()), 
      crossoverPoints_(order + 2)
    {
        gen.seed(seed);
        order_ = order;
        GMX_RELEASE_ASSERT(
            chromosomeLength >= 2,
            "There is 1 or less parameters to optimize, cannot crossover..."
        );
        availableIndices_.resize(chromosomeLength - 1);
        for (size_t i = 1; i < chromosomeLength; i++)
        {
            availableIndices_[i-1] = i;
        }
        crossoverPoints_[0] = 0;
        crossoverPoints_[crossoverPoints_.size() - 1] = chromosomeLength;
    }

    virtual void offspring(const ga::Genome *parent1,
                           const ga::Genome *parent2,
                           ga::Genome       *child1,
                           ga::Genome       *child2);

};


} //namespace ga


#endif //ALEXANDRIA_NPOINTCROSSOVER_H

/*
 * This source file is part of the Alexandria Chemistry Toolkit.
 *
 * Copyright (C) 2014-2022,2025
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


#ifndef GA_CROSSOVER_H
#define GA_CROSSOVER_H


#include <random>
#include <time.h>
#include <vector>

#include "genome.h"

namespace ga
{

/*! \brief Abstract class to perform crossover.
 * Given two \ref Individual (parents), it will mix their genomes to generate two new Individual (children)
 */
class Crossover
{

private:

    // Random number stuff
    std::random_device                     rd_base;
    std::mt19937                           gen_base;
    std::uniform_int_distribution<size_t>  dis_base;

protected:

    /*!
     * \brief Constructor
     * \param[in] chromosomeLength  length of the chromosome
     * \param[in] seed              seed for random number generation
     */
    Crossover(const size_t chromosomeLength,
              const int seed)
    : gen_base(rd_base()), dis_base(std::uniform_int_distribution<size_t>(1, chromosomeLength - 1))
    {
        gen_base.seed(seed);
    }

    //! Default destructor
    virtual ~Crossover() = default;

    /*! \brief Pick a random gene index
     * \return the selected index
     */
    size_t randIndex() { return dis_base(gen_base); }

public:

    /*!
     * \brief Perform crossover operation
     * \param[in] parent1   the first parent
     * \param[in] parent2   the second parent
     * \param[in] child1    the first child to write to
     * \param[in] child2    the second child to write to
     */
    virtual void offspring(const Genome *parent1,
                           const Genome *parent2,
                           Genome       *child1,
                           Genome       *child2) = 0;

};


} //namespace ga


#endif //GA_CROSSOVER_H

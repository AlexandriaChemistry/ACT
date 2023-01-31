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


#ifndef GA_SELECTOR_H
#define GA_SELECTOR_H

#include <random>
#include <time.h>
#include <vector>

#include "genome.h"

namespace ga
{


/*!
 * \brief Abstract class to select an Individual from the population based on its selection probability
 */
class Selector
{

public:

    /*!
     * Select an individual (by index) from the population
     * \param[in] pop   Pointer to the genes
     * \return          the index of the selected individual
     */
    virtual int select(const std::vector<Genome> *pop) = 0;

};


/*!
 * \brief Class for roulette-based selection. Uses cumulative probability to perform selection.
 */
class RouletteSelector : public Selector
{

private:

    // Random number stuff
    std::random_device                      rd;
    std::mt19937                            gen;
    std::uniform_real_distribution<double>  dis;

public:

    /*! \brief Constructor
     * \param[in] seed seed for random number generation
     */
    RouletteSelector(const int seed)
    : gen(rd()), dis(std::uniform_real_distribution<>(0.0, 1.0))
    {
        gen.seed(seed);
    }

    virtual int select(const std::vector<Genome> *pop);

};


} //namespace ga


#endif //GA_SELECTOR_H

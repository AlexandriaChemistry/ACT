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

/*! \internal \brief
 * Tests of the operators of the genetic algorithm: crossover, selection,
 * probability computation, termination and penalizing.
 * \author David van der Spoel <david.vanderspoel@icm.uu.se>
 */

#include <cmath>
#include <vector>

#include <gtest/gtest.h>

#include "act/ga/gene_pool.h"
#include "act/ga/genetic_algorithm.h"
#include "act/ga/initializer.h"
#include "act/ga/npointcrossover.h"
#include "act/ga/penalizer.h"
#include "act/ga/probability_computer.h"
#include "act/ga/selector.h"
#include "act/ga/terminator.h"

#include "gromacs/utility/exceptions.h"

namespace ga
{

namespace
{

//! Make a population with the training fitness given
std::vector<Genome> makePopulation(const std::vector<double> &fitness)
{
    std::vector<Genome> pop;
    for (size_t i = 0; i < fitness.size(); i++)
    {
        pop.push_back(Genome({ static_cast<double>(i) },
                             { { iMolSelect::Train, fitness[i] } }));
    }
    return pop;
}

//! \return the sum of the probabilities in a population
double sumOfProbabilities(const std::vector<Genome> &pop)
{
    double sum = 0;
    for (const auto &g : pop)
    {
        sum += g.probability();
    }
    return sum;
}

/*! \brief Check that the probabilities are normalized and decrease with the fitness
 * \param[in] pop    The population
 * \param[in] strict Whether a worse fitness gives a smaller probability,
 *                   otherwise it may be equal
 */
void checkProbabilities(const std::vector<Genome> &pop, bool strict = true)
{
    EXPECT_NEAR(1.0, sumOfProbabilities(pop), 1e-12);
    for (size_t i = 0; i < pop.size(); i++)
    {
        EXPECT_TRUE(std::isfinite(pop[i].probability())) << i;
        EXPECT_GE(pop[i].probability(), 0) << i;
        for (size_t j = 0; j < pop.size(); j++)
        {
            if (pop[i].fitness(iMolSelect::Train) < pop[j].fitness(iMolSelect::Train))
            {
                if (strict)
                {
                    EXPECT_GT(pop[i].probability(), pop[j].probability()) << i << " " << j;
                }
                else
                {
                    EXPECT_GE(pop[i].probability(), pop[j].probability()) << i << " " << j;
                }
            }
        }
    }
}

// Test: The Boltzmann probability at constant temperature
TEST(ProbabilityComputerTest, Boltzmann)
{
    std::vector<double> fitness = { 3, 1, 4, 2 };
    auto pop = makePopulation(fitness);
    double temp = 0.5;
    BoltzmannProbabilityComputer pc(temp, 10, 1, pop.size());
    pc.compute(&pop, 0);
    checkProbabilities(pop);
    double total = 0;
    for (auto f : fitness)
    {
        total += std::exp(1/(1e-4 + f)/temp);
    }
    for (size_t i = 0; i < pop.size(); i++)
    {
        EXPECT_NEAR(std::exp(1/(1e-4 + fitness[i])/temp)/total, pop[i].probability(), 1e-12) << i;
    }
}

// Test: When annealing to a low temperature the Boltzmann probabilities
// stay finite and the best genome gets all the probability
TEST(ProbabilityComputerTest, BoltzmannAnnealed)
{
    auto pop = makePopulation({ 3, 1, 4, 2 });
    int  maxGenerations = 100;
    BoltzmannProbabilityComputer pc(1, maxGenerations, 0.5, pop.size());
    for (int generation : { 0, 50, 75, 99, 100 })
    {
        pc.compute(&pop, generation);
        // Near zero temperature the probabilities of the worse genomes vanish
        checkProbabilities(pop, generation < 99);
    }
    EXPECT_EQ(1.0, pop[1].probability());
}

} // namespace

} // namespace ga

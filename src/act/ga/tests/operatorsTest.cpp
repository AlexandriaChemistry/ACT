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

/*! \brief Check that the children are a crossover of the parents
 * Each gene of the children comes from either parent, the children are
 * complementary, and the first child switches parent \p order times.
 */
void checkCrossover(const Genome &parent1, const Genome &parent2,
                    const Genome &child1,  const Genome &child2,
                    size_t        order)
{
    size_t switches = 0;
    bool   fromOne  = true;
    for (size_t j = 0; j < parent1.nBase(); j++)
    {
        bool one = child1.base(j) == parent1.base(j);
        if (one)
        {
            EXPECT_EQ(parent2.base(j), child2.base(j)) << j;
        }
        else
        {
            EXPECT_EQ(parent2.base(j), child1.base(j)) << j;
            EXPECT_EQ(parent1.base(j), child2.base(j)) << j;
        }
        if (one != fromOne)
        {
            switches += 1;
            fromOne   = one;
        }
    }
    // The first gene is never swapped
    EXPECT_EQ(child1.base(0), parent1.base(0));
    EXPECT_EQ(order, switches);
}

// Test: The children of an N-point crossover switch parent N times,
// for any number of crossover points
TEST(NPointCrossoverTest, SwitchesParentOrderTimes)
{
    const size_t length = 7;
    std::vector<double> bases1, bases2;
    for (size_t j = 0; j < length; j++)
    {
        bases1.push_back(j);
        bases2.push_back(-1.0 - j);
    }
    Genome parent1(bases1, {}), parent2(bases2, {});
    for (size_t order = 1; order < length; order++)
    {
        NPointCrossover crossover(length, order, 1993 + order);
        for (int k = 0; k < 20; k++)
        {
            Genome child1(std::vector<double>(length, 0), {});
            Genome child2(std::vector<double>(length, 0), {});
            crossover.offspring(&parent1, &parent2, &child1, &child2);
            checkCrossover(parent1, parent2, child1, child2, order);
        }
    }
}

// Test: The same seed gives the same children
TEST(NPointCrossoverTest, Reproducible)
{
    Genome parent1({ 1, 2, 3, 4, 5 }, {}), parent2({ 6, 7, 8, 9, 10 }, {});
    NPointCrossover a(5, 2, 42), b(5, 2, 42);
    for (int k = 0; k < 10; k++)
    {
        Genome a1(std::vector<double>(5, 0), {}), a2 = a1, b1 = a1, b2 = a1;
        a.offspring(&parent1, &parent2, &a1, &a2);
        b.offspring(&parent1, &parent2, &b1, &b2);
        EXPECT_EQ(a1.bases(), b1.bases());
        EXPECT_EQ(a2.bases(), b2.bases());
    }
}

// Test: The genomes are selected as often as their probability
TEST(RouletteSelectorTest, FrequencyMatchesProbability)
{
    std::vector<double> prob = { 0.1, 0.2, 0.3, 0.4 };
    auto pop = makePopulation({ 1, 2, 3, 4 });
    for (size_t i = 0; i < pop.size(); i++)
    {
        pop[i].setProbability(prob[i]);
    }
    RouletteSelector selector(1332);
    std::vector<int> count(pop.size(), 0);
    const int        ndraw = 100000;
    for (int k = 0; k < ndraw; k++)
    {
        int i = selector.select(&pop);
        ASSERT_GE(i, 0);
        ASSERT_LT(i, static_cast<int>(pop.size()));
        count[i] += 1;
    }
    for (size_t i = 0; i < pop.size(); i++)
    {
        EXPECT_NEAR(prob[i], static_cast<double>(count[i])/ndraw, 0.01) << i;
    }
}

// Test: The selected index is valid if the probabilities do not add up
// to one, and genomes with zero probability are not selected
TEST(RouletteSelectorTest, ValidIndex)
{
    auto pop = makePopulation({ 1, 2, 3 });
    pop[0].setProbability(0);
    pop[1].setProbability(0.2);
    pop[2].setProbability(0.3);
    RouletteSelector selector(506);
    for (int k = 0; k < 10000; k++)
    {
        int i = selector.select(&pop);
        EXPECT_TRUE(i == 1 || i == 2) << i;
    }
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

// Test: The probability is given by the rank, also if the population is not sorted
TEST(ProbabilityComputerTest, Rank)
{
    auto pop = makePopulation({ 3, 1, 4, 2 });
    RankProbabilityComputer pc(pop.size());
    pc.compute(&pop, 0);
    checkProbabilities(pop);
    std::vector<double> expected = { 0.2, 0.4, 0.1, 0.3 };
    for (size_t i = 0; i < pop.size(); i++)
    {
        EXPECT_NEAR(expected[i], pop[i].probability(), 1e-15) << i;
    }
    // A sorted population with a different size
    auto sorted = makePopulation({ 1, 2, 3, 4, 5 });
    pc.compute(&sorted, 0);
    for (size_t i = 0; i < sorted.size(); i++)
    {
        EXPECT_NEAR((5.0 - i)/15, sorted[i].probability(), 1e-15) << i;
    }
}

//! Make a pool with distinct genomes
GenePool makePool(size_t popSize)
{
    GenePool pool;
    for (size_t i = 0; i < popSize; i++)
    {
        pool.addGenome({ static_cast<double>(i), 2.0*i }, i);
    }
    return pool;
}

//! Terminator and penalizer that count how often they were asked
class CountingOperator : public Terminator, public Penalizer
{
public:
    //! Whether to answer yes
    bool answer_;
    //! The number of calls
    int  calls_ = 0;

    CountingOperator(bool answer) : answer_(answer) {}

    bool terminate(gmx::TextWriter *, const GenePool *, const int, bool) override
    {
        calls_ += 1;
        return answer_;
    }

    bool penalize(gmx::TextWriter *, GenePool *, const int) override
    {
        calls_ += 1;
        return answer_;
    }
};

//! A genetic algorithm that does not evolve
class StaticGA : public GeneticAlgorithm
{
public:
    StaticGA(std::vector<Terminator*> *terminators,
             std::vector<Penalizer*>  *penalizers)
    : GeneticAlgorithm(nullptr, nullptr, nullptr, nullptr, nullptr, nullptr,
                       terminators, penalizers, 4) {}

    bool evolve(alexandria::MsgHandler *, std::map<iMolSelect, Genome> *) override
    {
        return false;
    }
};

// Test: The genetic algorithm asks all terminators and penalizers,
// also after one of them said yes
TEST(GeneticAlgorithmTest, AskAll)
{
    CountingOperator yes(true), no(false);
    std::vector<Terminator*> terminators = { &yes, &no };
    std::vector<Penalizer*>  penalizers  = { &yes, &no };
    StaticGA ga(&terminators, &penalizers);
    auto pool = makePool(4);
    EXPECT_TRUE(ga.terminate(nullptr, &pool, 0, false));
    EXPECT_TRUE(ga.penalize(nullptr, &pool, 0));
    EXPECT_EQ(2, yes.calls_);
    EXPECT_EQ(2, no.calls_);
    std::vector<Terminator*> onlyNo = { &no };
    std::vector<Penalizer*>  noPen  = { &no };
    StaticGA ga2(&onlyNo, &noPen);
    EXPECT_FALSE(ga2.terminate(nullptr, &pool, 0, false));
    EXPECT_FALSE(ga2.penalize(nullptr, &pool, 0));
}

// Test: The gene pool can be updated with one of the same genome size only
TEST(GeneticAlgorithmTest, UpdateGenePool)
{
    CountingOperator         no(false);
    std::vector<Terminator*> terminators = { &no };
    std::vector<Penalizer*>  penalizers;
    StaticGA ga(&terminators, &penalizers);
    EXPECT_THROW(ga.updateGenePool(GenePool()), gmx::InvalidInputError);
    auto pool = makePool(4);
    ga.updateGenePool(pool);
    EXPECT_EQ(4u, ga.getLastPop().popSize());
    ga.updateGenePool(makePool(6));
    EXPECT_EQ(6u, ga.getLastPop().popSize());
    GenePool other;
    other.addGenome({ 1.0, 2.0, 3.0 }, 0);
    EXPECT_THROW(ga.updateGenePool(other), gmx::InvalidInputError);
}

} // namespace

} // namespace ga

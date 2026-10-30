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
 * Tests of the train_ff module with a population of four individuals on
 * four cores, using the HYBRID, GA and MCMC optimizers with input derived
 * from the benchmarks in ACTdata/Performance.
 * \author David van der Spoel <david.vanderspoel@icm.uu.se>
 */

#include "actpre.h"

#include <string>
#include <vector>

#include <gtest/gtest.h>

#include "trainffTestUtil.h"

#include "testutils/cmdlinetest.h"
#include "testutils/refdata.h"
#include "testutils/testasserts.h"

namespace alexandria
{

namespace
{

using trainfftest::TrainOutput;
using Options = std::vector<std::pair<std::string, std::string>>;

//! The population size, one individual per core
const int popSize = 4;

//! Options for the genetic algorithm, as in the benchmarks
Options gaOptions(const std::string &optimizer, int maxGenerations)
{
    return {
        { "-seed", "219463979" },
        { "-v", "3" },
        { "-optimizer", optimizer },
        { "-pop_size", std::to_string(popSize) },
        { "-max_generations", std::to_string(maxGenerations) },
        { "-max_test_generations", "10" },
        { "-random_init", "" },
        { "-n_crossovers", "1" },
        { "-ga_boltz_temp", "0.5" },
        { "-ga_boltz_anneal", "0.5" },
        { "-cp_gen_interval", "2" },
        { "-cp_pop_frac", "0.2" },
        { "-fc_bound", "10000" },
        { "-ener_boltz_temp", "3000" }
    };
}

//! Options for the MCMC part of HYBRID and MCMC
Options mcmcOptions(int maxIter)
{
    return {
        { "-max_iter", std::to_string(maxIter) },
        { "-temp", "0.025" },
        { "-anneal_begin", "0.1" },
        { "-anneal_end", "0.9" },
        { "-anneal_global_begin", "0.1" },
        { "-anneal_global_end", "0.9" }
    };
}

//! \return the options combined
Options combine(Options a, const Options &b)
{
    a.insert(a.end(), b.begin(), b.end());
    return a;
}

class TrainFFMpiTest : public gmx::test::CommandLineTestBase
{
protected:
    /*! \brief Run train_ff and check the output on the master
     * All ranks wait until the master has read the output, since the
     * temporary files are removed by each rank at the end of the test.
     * \param[in] ffName      Force field file
     * \param[in] molName     Base name of the molprop and selection file
     * \param[in] options     The options
     * \param[in] nFitness    Expected number of lines in the fitness file,
     *                        zero if there is none
     */
    void runAndCheck(const std::string &ffName,
                     const std::string &molName,
                     const Options     &options,
                     size_t             nFitness)
    {
        auto output = trainfftest::runTrainFF(&fileManager(), "train", ffName, molName, options);
        EXPECT_EQ(0, output.status);
        if (gmx_node_rank() == 0)
        {
            EXPECT_GT(output.finalChi2, 0);
            EXPECT_NE(std::string::npos, output.forceField.find("</alexandria_chemistry>"));
            EXPECT_EQ(nFitness, output.fitness.size());
            // Each line has the fitness of each individual
            for (const auto &line : output.fitness)
            {
                std::istringstream iss(line);
                double             f;
                int                n = 0;
                while (iss >> f)
                {
                    // The file has six decimals
                    EXPECT_GE(f, output.finalChi2*(1 - 1e-6) - 1e-6);
                    n += 1;
                }
                EXPECT_EQ(popSize, n) << line;
            }
            gmx::test::TestReferenceChecker checker(rootChecker());
            checker.setDefaultTolerance(gmx::test::relativeToleranceAsFloatingPoint(1.0, 1e-4));
            trainfftest::checkOutput(&checker, output);
        }
        trainfftest::barrier();
    }
};

// Test: HYBRID training of the hydrogen fluoride dimer, with virtual sites
// and a catastrophe in the second generation
TEST_F(TrainFFMpiTest, HybridHydrogenFluoride)
{
    auto options = combine(gaOptions("HYBRID", 3), mcmcOptions(2));
    options.push_back({ "-fit", "charge eta chi zeta vs2a sigma epsilon gamma delta" });
    options.push_back({ "-fc_allelec", "1" });
    options.push_back({ "-fc_exch", "1" });
    options.push_back({ "-fc_disp", "1" });
    options.push_back({ "-fc_mu", "1" });
    options.push_back({ "-fc_quad", "1" });
    // The initial population and one line per generation
    runAndCheck("4S-lj14-7-pg-hvdw.xml", "hf-hf", options, 4);
}

// Test: GA training of Lennard-Jones parameters on alcohol and water dimers
TEST_F(TrainFFMpiTest, GAInterBooze)
{
    auto options = gaOptions("GA", 5);
    options.push_back({ "-fit", "sigma epsilon" });
    options.push_back({ "-fc_allelec", "1" });
    options.push_back({ "-fc_exch", "1" });
    options.push_back({ "-fc_disp", "1" });
    // The initial population and one line per generation
    runAndCheck("lj-pc.xml", "booze-inter", options, 6);
}

// Test: MCMC training of bonded parameters of alcohols on four cores
TEST_F(TrainFFMpiTest, MCMCIntraBooze)
{
    auto options = combine(gaOptions("MCMC", 1), mcmcOptions(2));
    options.push_back({ "-fit", "De b c bondenergy kb bondlength angle kt c1 c2 c3 c4 c5" });
    options.push_back({ "-fc_epot", "1" });
    options.push_back({ "-fc_force", "0.1" });
    // MCMC does not write fitness files
    runAndCheck("halgren-PG.xml", "booze-mono", options, 0);
}

} // namespace

} // namespace alexandria

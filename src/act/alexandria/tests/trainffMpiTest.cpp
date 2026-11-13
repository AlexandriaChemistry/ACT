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

#include <algorithm>
#include <cstdio>
#include <cstdlib>
#include <fstream>
#include <sstream>
#include <string>
#include <vector>

#include <gtest/gtest.h>

#include "trainffTestUtil.h"

#include "act/gpu/gpu.h"

#include "testutils/cmdlinetest.h"
#include "testutils/gputests.h"
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

// Test: GA training continued from the gene pool of a previous run starts
// with the best individual of that run
TEST_F(TrainFFMpiTest, GAGenePool)
{
    auto options = gaOptions("GA", 2);
    options.push_back({ "-fit", "sigma epsilon" });
    options.push_back({ "-fc_allelec", "1" });
    options.push_back({ "-fc_exch", "1" });
    options.push_back({ "-fc_disp", "1" });
    auto first = trainfftest::runTrainFF(&fileManager(), "first", "lj-pc.xml", "booze-inter", options);
    EXPECT_EQ(0, first.status);
    // The pool is written by the master when there is a new best individual
    auto pool = fileManager().getTemporaryFilePath("first-pool.csv");
    auto continued = options;
    continued.push_back({ "-gpin", pool });
    auto second = trainfftest::runTrainFF(&fileManager(), "second", "lj-pc.xml", "booze-inter", continued);
    EXPECT_EQ(0, second.status);
    if (gmx_node_rank() == 0)
    {
        ASSERT_FALSE(second.fitness.empty());
        std::istringstream iss(second.fitness[0]);
        double             best = 1e30;
        int                n    = 0;
        for (double f; iss >> f; n++)
        {
            best = std::min(best, f);
        }
        EXPECT_EQ(popSize, n);
        // The fitness file has six decimals
        EXPECT_NEAR(first.finalChi2, best, 1e-6*first.finalChi2 + 1e-6);
        EXPECT_LE(second.finalChi2, first.finalChi2*(1 + 1e-6) + 1e-6);
    }
    trainfftest::barrier();
}

// Test: with fewer individuals than processes, the process left over is idle
// and leaves the communicator of ACT. The others must not use a collective
// of all processes then: the charge map was broadcast to all of them, which
// the barrier of the idle process matched, such that the processes hung at
// the end of the program, or here.
TEST_F(TrainFFMpiTest, IdleProcess)
{
    if (gmx_node_num() != popSize)
    {
        // Only with as many processes as the other tests have individuals
        return;
    }
    // The GA needs an even population
    auto options = combine(gaOptions("MCMC", 1), mcmcOptions(2));
    for (auto &o : options)
    {
        if (o.first == "-pop_size")
        {
            o.second = std::to_string(popSize - 1);
        }
    }
    options.push_back({ "-fit", "sigma epsilon" });
    options.push_back({ "-fc_allelec", "1" });
    options.push_back({ "-fc_exch", "1" });
    options.push_back({ "-fc_disp", "1" });
    auto output = trainfftest::runTrainFF(&fileManager(), "idle", "lj-pc.xml", "booze-inter", options);
    EXPECT_EQ(0, output.status);
    if (gmx_node_rank() == 0)
    {
        EXPECT_GT(output.finalChi2, 0);
        EXPECT_NE(std::string::npos, output.forceField.find("</alexandria_chemistry>"));
    }
    // A collective of all processes, as at the end of the program
    int rank = gmx_node_rank();
    int sum  = 0;
#if GMX_MPI
    MPI_Allreduce(&rank, &sum, 1, MPI_INT, MPI_SUM, MPI_COMM_WORLD);
#endif
    EXPECT_EQ(gmx_node_num()*(gmx_node_num() - 1)/2, sum);
    trainfftest::barrier();
}

//! \return the contents of a file
static std::string fileContents(const std::string &fileName)
{
    std::ifstream      in(fileName);
    std::ostringstream ss;
    ss << in.rdbuf();
    return ss.str();
}

// Test: The gene pool written every five generations and at the end is that
// written in each generation with a better individual. This training finds
// better individuals in generations 1, 2, 5, 6 and 7, so the last one is
// only written at the end.
TEST_F(TrainFFMpiTest, GenePoolInterval)
{
    auto options = combine(gaOptions("HYBRID", 7), mcmcOptions(2));
    options.push_back({ "-fit", "sigma epsilon" });
    options.push_back({ "-fc_allelec", "1" });
    options.push_back({ "-fc_exch", "1" });
    options.push_back({ "-fc_disp", "1" });
    auto each = options;
    each.push_back({ "-gpout_interval", "1" });
    auto five = options;
    five.push_back({ "-gpout_interval", "5" });
    auto first  = trainfftest::runTrainFF(&fileManager(), "each", "lj-pc.xml", "booze-inter", each);
    auto second = trainfftest::runTrainFF(&fileManager(), "five", "lj-pc.xml", "booze-inter", five);
    EXPECT_EQ(0, first.status);
    EXPECT_EQ(0, second.status);
    if (gmx_node_rank() == 0)
    {
        auto poolEach = fileContents(fileManager().getTemporaryFilePath("each-pool.csv"));
        auto poolFive = fileContents(fileManager().getTemporaryFilePath("five-pool.csv"));
        EXPECT_FALSE(poolEach.empty());
        EXPECT_EQ(poolEach, poolFive);
        EXPECT_EQ(first.finalChi2, second.finalChi2);
    }
    trainfftest::barrier();
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

// Test: MCMC with two individuals, on four cores each has a helper that
// computes part of the compounds. The reference data are from a run on
// two cores without helpers, such that the helpers give the same result.
TEST_F(TrainFFMpiTest, MCMCHelpers)
{
    auto options = combine(gaOptions("MCMC", 1), mcmcOptions(2));
    for (auto &o : options)
    {
        if (o.first == "-pop_size")
        {
            o.second = "2";
        }
    }
    options.push_back({ "-fit", "De b c bondenergy kb bondlength angle kt c1 c2 c3 c4 c5" });
    options.push_back({ "-fc_epot", "1" });
    options.push_back({ "-fc_force", "0.1" });
    runAndCheck("halgren-PG.xml", "booze-mono", options, 0);
}

// Test: The two-dimensional sensitivity analysis of a force field, with
// the chi-squared divided over the master and the middlemen. The reference
// data are a copy of those of TrainFFTest.Sensitivity2DJson, such that the
// distributed analysis has to give the same result as on a single core.
TEST_F(TrainFFMpiTest, Sensitivity2DMiddlemen)
{
    Options options = {
        { "-seed", "219463979" },
        { "-v", "3" },
        { "-optimizer", "MCMC" },
        { "-pop_size", std::to_string(popSize) },
        { "-max_iter", "1" },
        { "-nooptimize", "" },
        { "-norandom_init", "" },
        { "-sensitivity", "2D" },
        { "-fc_bound", "10000" },
        { "-ener_boltz_temp", "3000" },
        { "-fit", "sigma epsilon" },
        { "-fc_allelec", "1" },
        { "-fc_exch", "1" },
        { "-fc_disp", "1" }
    };
    auto output = trainfftest::runTrainFF(&fileManager(), "sens2d", "lj-pc-trained.xml",
                                          "booze-inter", options);
    EXPECT_EQ(0, output.status);
    if (gmx_node_rank() == 0)
    {
        gmx::test::TestReferenceChecker checker(rootChecker());
        trainfftest::checkSensitivity2D(&checker, output);
    }
    trainfftest::barrier();
}

// Test: HYBRID with the two-dimensional analysis as the test for a local
// minimum in each generation, the middlemen compute part of it
TEST_F(TrainFFMpiTest, HybridCurvature2D)
{
    auto options = combine(gaOptions("HYBRID", 3), mcmcOptions(2));
    options.push_back({ "-fit", "charge eta chi zeta vs2a sigma epsilon gamma delta" });
    options.push_back({ "-fc_allelec", "1" });
    options.push_back({ "-fc_exch", "1" });
    options.push_back({ "-fc_disp", "1" });
    options.push_back({ "-fc_mu", "1" });
    options.push_back({ "-fc_quad", "1" });
    options.push_back({ "-sensitivity", "2D" });
    options.push_back({ "-max_local_minimum_generations", "1" });
    auto output = trainfftest::runTrainFF(&fileManager(), "curv2d", "4S-lj14-7-pg-hvdw.xml",
                                          "hf-hf", options);
    EXPECT_EQ(0, output.status);
    if (gmx_node_rank() == 0)
    {
        // The test in each generation is two-dimensional
        int ngen = 0;
        for (const auto &line : output.log)
        {
            if (line.find("Curvature analysis (2D)") != std::string::npos)
            {
                EXPECT_NE(std::string::npos, line.find("pairs converged")) << line;
                ngen += 1;
            }
            EXPECT_EQ(std::string::npos, line.find("Curvature analysis (1D)")) << line;
        }
        EXPECT_GT(ngen, 0);
        // And so is the analysis after training
        ASSERT_TRUE(output.sensitivity.contains("sensitivity"));
        EXPECT_EQ("2D", output.sensitivity["sensitivity"]["mode"].get<std::string>());
        EXPECT_TRUE(output.sensitivity["sensitivity"].contains("pairs"));
    }
    trainfftest::barrier();
}

//! \return whether each rank trains polarizable models on the GPU
bool gpuTraining()
{
    int use = (!GMX_DOUBLE && gpu::useGpu()) ? 1 : 0;
    int all = 0;
    MPI_Allreduce(&use, &all, 1, MPI_INT, MPI_MIN, MPI_COMM_WORLD);
    return all == 1;
}

//! Options for short HYBRID runs of a polarizable model, as the alcohol benchmarks
Options polarizableOptions(int populationSize, int maxIter, const std::string &fit)
{
    return {
        { "-seed", "219463979" },
        { "-v", "3" },
        { "-optimizer", "HYBRID" },
        { "-pop_size", std::to_string(populationSize) },
        { "-max_generations", "1" },
        { "-random_init", "" },
        { "-n_crossovers", "1" },
        { "-max_iter", std::to_string(maxIter) },
        { "-temp", "0.025" },
        { "-fc_bound", "10000" },
        { "-fc_inter", "1" },
        { "-fit", fit }
    };
}

/*! \brief Train on the GPU twice in ways that have to give the same result
 * The shells are relaxed from the same positions each time, such that the
 * energies of the same parameters are the same, see -shell_keep.
 * \param[in] fileManager The test files
 * \param[in] label       Prefix of the output files
 * \param[in] options     The options of the first run
 * \param[in] other       The options of the second run
 * \param[in] otherEnv    The environment variable set during the second run, if any
 */
void compareGpuTrainings(gmx::test::TestFileManager *fileManager,
                         const std::string          &label,
                         const Options              &options,
                         const Options              &other,
                         const char                 *otherEnv = nullptr)
{
    auto first = trainfftest::runTrainFF(fileManager, label + "-1", "alc-pol-wbh.xml",
                                         "booze-inter", options);
    EXPECT_EQ(0, first.status);
    if (otherEnv)
    {
        setenv(otherEnv, "1", 1);
    }
    auto second = trainfftest::runTrainFF(fileManager, label + "-2", "alc-pol-wbh.xml",
                                          "booze-inter", other);
    if (otherEnv)
    {
        unsetenv(otherEnv);
    }
    EXPECT_EQ(0, second.status);
    if (gmx_node_rank() == 0)
    {
        for (const auto *output : { &first, &second })
        {
            EXPECT_TRUE(std::any_of(output->log.begin(), output->log.end(), [](const std::string &line)
                                    { return line.find("GPU resident relaxation") != std::string::npos; }))
                << "The GPU did not train";
        }
        ASSERT_GT(first.finalChi2, 0);
        // The log has six significant digits
        EXPECT_NEAR(first.finalChi2, second.finalChi2, 2e-5*first.finalChi2);
    }
    trainfftest::barrier();
}

// Test: With -lockstep, the GPU computes the individuals of a middleman
// together, each with the bounds of its own charges, as one at a time. The
// second charges of the polarizable model are trained with the charges, such
// that they differ between the individuals.
TEST_F(TrainFFMpiTest, LockstepChargesGpu)
{
    if (!gmx::test::gpuTestRuns(gpuTraining(), "no GPU or double precision"))
    {
        return;
    }
    // Four individuals on each of the three middlemen
    auto options = polarizableOptions(12, 2, "charge");
    options.push_back({ "-noshell_keep", "" });
    auto together = options, oneByOne = options;
    together.push_back({ "-lockstep", "4" });
    oneByOne.push_back({ "-lockstep", "1" });
    compareGpuTrainings(&fileManager(), "lockstep", together, oneByOne);
}

// Test: The molecules whose parameters did not change are not computed again
// on the GPU, their deviations, including the bounds of their charges, are
// those kept from before. ACT_NO_GPU_SKIP computes all of them each time.
TEST_F(TrainFFMpiTest, SkippedMoleculesGpu)
{
    if (!gmx::test::gpuTestRuns(gpuTraining(), "no GPU or double precision"))
    {
        return;
    }
    // Two individuals on each middleman, one after the other
    auto options = polarizableOptions(6, 6, "zeta charge alpha eta chi sigma epsilon gamma");
    options.push_back({ "-noshell_keep", "" });
    compareGpuTrainings(&fileManager(), "skipped", options, options, "ACT_NO_GPU_SKIP");
}

} // namespace

} // namespace alexandria

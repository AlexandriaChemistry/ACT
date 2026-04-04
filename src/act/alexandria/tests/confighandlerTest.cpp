/*
 * This source file is part of the Alexandria Chemistry Toolkit.
 *
 * Copyright (C) 2024-2026
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
 * Tests for configuration handler classes and enum conversions.
 * \author David van der Spoel <david.vanderspoel@icm.uu.se>
 */

#include "actpre.h"

#include <string>

#include <gtest/gtest.h>

#include "act/alexandria/confighandler.h"
#include "act/basics/msg_handler.h"
#include "gromacs/utility/exceptions.h"

#include "testutils/testasserts.h"

namespace alexandria
{

namespace
{

// ---- OptimizerAlg enum conversion tests ----

TEST(ConfigHandlerTest, OptimizerAlgRoundTrip)
{
    // All enum values should round-trip through string conversion
    auto check = [](OptimizerAlg alg)
    {
        std::string str = optimizerAlgToString(alg);
        EXPECT_EQ(stringToOptimizerAlg(str), alg);
    };
    check(OptimizerAlg::MCMC);
    check(OptimizerAlg::GA);
    check(OptimizerAlg::HYBRID);
}

TEST(ConfigHandlerTest, OptimizerAlgToStringValues)
{
    EXPECT_EQ(optimizerAlgToString(OptimizerAlg::MCMC), "MCMC");
    EXPECT_EQ(optimizerAlgToString(OptimizerAlg::GA), "GA");
    EXPECT_EQ(optimizerAlgToString(OptimizerAlg::HYBRID), "HYBRID");
}

TEST(ConfigHandlerTest, StringToOptimizerAlgInvalid)
{
    EXPECT_THROW(stringToOptimizerAlg("INVALID"), gmx::InvalidInputError);
}

// ---- ProbabilityComputerAlg enum conversion tests ----

TEST(ConfigHandlerTest, ProbabilityComputerAlgRoundTrip)
{
    auto check = [](ProbabilityComputerAlg alg)
    {
        std::string str = probabilityComputerAlgToString(alg);
        EXPECT_EQ(stringToProbabilityComputerAlg(str), alg);
    };
    check(ProbabilityComputerAlg::pcRANK);
    check(ProbabilityComputerAlg::pcFITNESS);
    check(ProbabilityComputerAlg::pcBOLTZMANN);
}

TEST(ConfigHandlerTest, ProbabilityComputerAlgToStringValues)
{
    EXPECT_EQ(probabilityComputerAlgToString(ProbabilityComputerAlg::pcRANK), "RANK");
    EXPECT_EQ(probabilityComputerAlgToString(ProbabilityComputerAlg::pcFITNESS), "FITNESS");
    EXPECT_EQ(probabilityComputerAlgToString(ProbabilityComputerAlg::pcBOLTZMANN), "BOLTZMANN");
}

TEST(ConfigHandlerTest, StringToProbabilityComputerAlgInvalid)
{
    EXPECT_THROW(stringToProbabilityComputerAlg("INVALID"), gmx::InvalidInputError);
}

// ---- eMinimizeAlgorithm enum conversion tests ----

TEST(ConfigHandlerTest, EMinimizeAlgorithmRoundTrip)
{
    auto check = [](eMinimizeAlgorithm alg)
    {
        std::string str = eMinimizeAlgorithmToString(alg);
        EXPECT_EQ(stringToEMinimizeAlgorithm(str), alg);
    };
    check(eMinimizeAlgorithm::LBFGS);
    check(eMinimizeAlgorithm::Steep);
    check(eMinimizeAlgorithm::Newton);
}

TEST(ConfigHandlerTest, EMinimizeAlgorithmToStringValues)
{
    EXPECT_EQ(eMinimizeAlgorithmToString(eMinimizeAlgorithm::LBFGS), "LBFGS");
    EXPECT_EQ(eMinimizeAlgorithmToString(eMinimizeAlgorithm::Steep), "Steep");
    EXPECT_EQ(eMinimizeAlgorithmToString(eMinimizeAlgorithm::Newton), "Newton");
}

// ---- GAConfigHandler getters/setters tests ----

TEST(ConfigHandlerTest, GAConfigHandlerDefaults)
{
    GAConfigHandler gaCfg;
    EXPECT_EQ(gaCfg.optimizer(), OptimizerAlg::GA);
    EXPECT_EQ(gaCfg.popSize(), 1);
    EXPECT_EQ(gaCfg.nElites(), 0);
    EXPECT_TRUE(gaCfg.randomInit());
    EXPECT_EQ(gaCfg.nCrossovers(), 1);
    EXPECT_TRUE(gaCfg.sort());
    EXPECT_EQ(gaCfg.probabilityComputerAlg(), ProbabilityComputerAlg::pcRANK);
    EXPECT_EQ(gaCfg.maxGenerations(), 10);
    EXPECT_FALSE(gaCfg.evaluateTestset());
}

TEST(ConfigHandlerTest, GAConfigHandlerSetters)
{
    GAConfigHandler gaCfg;
    gaCfg.setOptimizerAlg(OptimizerAlg::HYBRID);
    EXPECT_EQ(gaCfg.optimizer(), OptimizerAlg::HYBRID);

    gaCfg.setPopSize(50);
    EXPECT_EQ(gaCfg.popSize(), 50);

    gaCfg.setCrossovers(3);
    EXPECT_EQ(gaCfg.nCrossovers(), 3);

    gaCfg.setProbabilityComputerAlg(ProbabilityComputerAlg::pcBOLTZMANN);
    EXPECT_EQ(gaCfg.probabilityComputerAlg(), ProbabilityComputerAlg::pcBOLTZMANN);

    gaCfg.setPrMut(0.5);
    EXPECT_NEAR(gaCfg.prMut(), 0.5, 1e-10);

    gaCfg.setEvaluateTestset(true);
    EXPECT_TRUE(gaCfg.evaluateTestset());
}

// ---- BayesConfigHandler tests ----

TEST(ConfigHandlerTest, BayesConfigHandlerDefaults)
{
    BayesConfigHandler bayesCfg;
    EXPECT_EQ(bayesCfg.maxIter(), 100);
    EXPECT_EQ(bayesCfg.seed(), 0);
    EXPECT_NEAR(bayesCfg.step(), 0.02, 1e-10);
    EXPECT_NEAR(bayesCfg.temperature(), 5.0, 1e-10);
    EXPECT_FALSE(bayesCfg.temperatureWeighting());
    EXPECT_FALSE(bayesCfg.checkPoint());
}

TEST(ConfigHandlerTest, BayesConfigHandlerSetters)
{
    BayesConfigHandler bayesCfg;
    bayesCfg.setMaxIter(200);
    EXPECT_EQ(bayesCfg.maxIter(), 200);

    bayesCfg.setSeed(42);
    EXPECT_EQ(bayesCfg.seed(), 42);

    bayesCfg.setStep(0.05);
    EXPECT_NEAR(bayesCfg.step(), 0.05, 1e-10);

    bayesCfg.setTemperature(10.0);
    EXPECT_NEAR(bayesCfg.temperature(), 10.0, 1e-10);
}

TEST(ConfigHandlerTest, BayesTemperatureWithGlobalAnnealing)
{
    BayesConfigHandler bayesCfg;
    bayesCfg.setTemperature(10.0);

    // Without global annealing, temperature at generation 0 should be base temp
    real t0 = bayesCfg.temperature(0, 100);
    EXPECT_NEAR(t0, 10.0, 1e-10);
}

TEST(ConfigHandlerTest, BayesAnnealBehavior)
{
    BayesConfigHandler bayesCfg;
    // Default anneal is 1, so annealing should be off
    EXPECT_FALSE(bayesCfg.annealing());
    EXPECT_FALSE(bayesCfg.anneal(0, 0));

    bayesCfg.setAnneal(0.5);
    EXPECT_TRUE(bayesCfg.annealing());
    EXPECT_NEAR(bayesCfg.annealStart(), 0.5, 1e-10);
}

TEST(ConfigHandlerTest, BayesComputeBeta)
{
    BayesConfigHandler bayesCfg;
    bayesCfg.setTemperature(10.0);
    bayesCfg.setMaxIter(100);
    bayesCfg.setAnneal(0.5);

    // beta at iter 0 with no global annealing: should be 1/temp with annealing formula
    double beta = bayesCfg.computeBeta(0, 10, 0);
    EXPECT_GT(beta, 0.0);

    // When iter >= maxiter, temperature should be very small → beta very large
    double betaEnd = bayesCfg.computeBeta(0, 10, 100);
    EXPECT_GT(betaEnd, beta);
}

// ---- SimulationConfigHandler tests ----

TEST(ConfigHandlerTest, SimulationConfigHandlerDefaults)
{
    SimulationConfigHandler simCfg;
    EXPECT_EQ(simCfg.nsteps(), 0);
    EXPECT_FALSE(simCfg.minimize());
    EXPECT_FALSE(simCfg.singlePoint());
    EXPECT_EQ(simCfg.maxIter(), 100);
    EXPECT_GT(simCfg.deltat(), 0.0);
}

TEST(ConfigHandlerTest, SimulationConfigHandlerSetters)
{
    SimulationConfigHandler simCfg;
    simCfg.setMinimize(true);
    EXPECT_TRUE(simCfg.minimize());

    simCfg.setMinimizeAlgorithm(eMinimizeAlgorithm::Steep);
    EXPECT_EQ(simCfg.minAlg(), eMinimizeAlgorithm::Steep);

    simCfg.setForceTolerance(1e-8);
    EXPECT_NEAR(simCfg.forceTolerance(), 1e-8, 1e-15);

    simCfg.setMaxIter(500);
    EXPECT_EQ(simCfg.maxIter(), 500);
}

} // namespace

} // namespace alexandria

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
#include <vector>

#include <gtest/gtest.h>

#include "act/alexandria/confighandler.h"
#include "act/basics/msg_handler.h"
#include "gromacs/commandline/pargs.h"
#include "gromacs/fileio/oenv.h"
#include "gromacs/utility/exceptions.h"

#include "testutils/cmdlinetest.h"
#include "testutils/testasserts.h"

namespace alexandria
{

namespace
{

// ---- OptimizerAlg enum conversion tests ----

TEST(ConfigHandlerTest, OptimizerAlgRoundTrip)
{
    // All enum values should round-trip through string conversion
    const auto check = [](OptimizerAlg alg)
    {
        const std::string str = optimizerAlgToString(alg);
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

// ---- SensitivityMode enum conversion tests ----

TEST(ConfigHandlerTest, SensitivityModeRoundTrip)
{
    for (auto mode : { SensitivityMode::No, SensitivityMode::OneD, SensitivityMode::TwoD })
    {
        EXPECT_EQ(stringToSensitivityMode(sensitivityModeToString(mode)), mode);
    }
}

TEST(ConfigHandlerTest, SensitivityModeToStringValues)
{
    EXPECT_EQ(sensitivityModeToString(SensitivityMode::No), "No");
    EXPECT_EQ(sensitivityModeToString(SensitivityMode::OneD), "1D");
    EXPECT_EQ(sensitivityModeToString(SensitivityMode::TwoD), "2D");
}

TEST(ConfigHandlerTest, StringToSensitivityModeInvalid)
{
    EXPECT_THROW(stringToSensitivityMode("3D"), gmx::InvalidInputError);
    // The old boolean values are not accepted
    EXPECT_THROW(stringToSensitivityMode("yes"), gmx::InvalidInputError);
    EXPECT_THROW(stringToSensitivityMode(""), gmx::InvalidInputError);
}

// ---- ProbabilityComputerAlg enum conversion tests ----

TEST(ConfigHandlerTest, ProbabilityComputerAlgRoundTrip)
{
    const auto check = [](ProbabilityComputerAlg alg)
    {
        const std::string str = probabilityComputerAlgToString(alg);
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
    const auto check = [](eMinimizeAlgorithm alg)
    {
        const std::string str = eMinimizeAlgorithmToString(alg);
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

/*! \brief Parse the options of a GAConfigHandler from a command line
 * The MCMC optimizer is selected, since the default GA does not allow
 * the default population size of one.
 * \param[in] gach The config handler
 * \param[in] args The arguments after the MCMC optimizer
 */
static void parseGAOptions(GAConfigHandler                  *gach,
                           const std::vector<const char *> &args)
{
    std::vector<t_pargs>      pargs;
    std::vector<t_filenm>     fnms;
    gach->add_options(&pargs, &fnms);
    std::vector<const char *> cmdline = { "test", "-optimizer", "MCMC" };
    cmdline.insert(cmdline.end(), args.begin(), args.end());
    gmx::test::CommandLine    cmd(cmdline);
    gmx_output_env_t         *oenv = nullptr;
    ASSERT_TRUE(parse_common_args(&cmd.argc(), cmd.argv(), 0, 0, nullptr,
                                  pargs.size(), pargs.data(), 0, nullptr, 0, nullptr, &oenv));
    output_env_done(oenv);
    MsgHandler msghandler;
    gach->check_pargs(&msghandler);
}

TEST(ConfigHandlerTest, SensitivityOption)
{
    // The default is the one-dimensional analysis
    GAConfigHandler gaDefault;
    EXPECT_EQ(gaDefault.sensitivity(), SensitivityMode::OneD);
    parseGAOptions(&gaDefault, {});
    EXPECT_EQ(gaDefault.sensitivity(), SensitivityMode::OneD);
    for (auto mode : { SensitivityMode::No, SensitivityMode::OneD, SensitivityMode::TwoD })
    {
        GAConfigHandler gach;
        parseGAOptions(&gach, { "-sensitivity", sensitivityModeToString(mode).c_str() });
        EXPECT_EQ(gach.sensitivity(), mode);
    }
    // Parsing again without the option gives the default, not the last choice
    GAConfigHandler gaAgain;
    parseGAOptions(&gaAgain, {});
    EXPECT_EQ(gaAgain.sensitivity(), SensitivityMode::OneD);
}

TEST(ConfigHandlerTest, SensitivitySetter)
{
    GAConfigHandler gach;
    gach.setSensitivity(SensitivityMode::TwoD);
    EXPECT_EQ(gach.sensitivity(), SensitivityMode::TwoD);
    gach.setSensitivity(SensitivityMode::No);
    EXPECT_EQ(gach.sensitivity(), SensitivityMode::No);
}

TEST(ConfigHandlerTest, GAConfigHandlerDefaults)
{
    const GAConfigHandler gaCfg;
    EXPECT_EQ(gaCfg.optimizer(), OptimizerAlg::GA);
    EXPECT_EQ(gaCfg.popSize(), 1);
    EXPECT_EQ(gaCfg.nElites(), 0);
    EXPECT_TRUE(gaCfg.randomInit());
    EXPECT_EQ(gaCfg.nCrossovers(), 1);
    EXPECT_TRUE(gaCfg.sort());
    EXPECT_EQ(gaCfg.probabilityComputerAlg(), ProbabilityComputerAlg::pcRANK);
    EXPECT_EQ(gaCfg.maxGenerations(), 10);
    EXPECT_FALSE(gaCfg.evaluateTestset());
    EXPECT_EQ(gaCfg.nstFitness(), 1);
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
    EXPECT_NEAR(gaCfg.prMut(), 0.5, gmx::test::precisionTolerance(1e-10, 1e-7));

    gaCfg.setEvaluateTestset(true);
    EXPECT_TRUE(gaCfg.evaluateTestset());
}

// ---- BayesConfigHandler tests ----

TEST(ConfigHandlerTest, BayesConfigHandlerDefaults)
{
    const BayesConfigHandler bayesCfg;
    EXPECT_EQ(bayesCfg.maxIter(), 100);
    EXPECT_EQ(bayesCfg.seed(), 0);
    EXPECT_NEAR(bayesCfg.step(), 0.02, gmx::test::precisionTolerance(1e-10, 1e-7));
    EXPECT_NEAR(bayesCfg.temperature(), 5.0, gmx::test::precisionTolerance(1e-10, 1e-7));
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
    EXPECT_NEAR(bayesCfg.step(), 0.05, gmx::test::precisionTolerance(1e-10, 1e-7));

    bayesCfg.setTemperature(10.0);
    EXPECT_NEAR(bayesCfg.temperature(), 10.0, gmx::test::precisionTolerance(1e-10, 1e-7));
}

TEST(ConfigHandlerTest, BayesTemperatureWithGlobalAnnealing)
{
    BayesConfigHandler bayesCfg;
    double T0 = 10.0;
    bayesCfg.setTemperature(T0);

    // Without global annealing, temperature at generation 0 should be base temp
    real t0 = bayesCfg.temperature();
    EXPECT_NEAR(t0, T0, gmx::test::precisionTolerance(1e-10, 1e-7));
}

TEST(ConfigHandlerTest, BayesAnnealBehavior)
{
    BayesConfigHandler bayesCfg;
    // Default anneal is 1, so annealing should be off
    EXPECT_FALSE(bayesCfg.annealing());

    bayesCfg.setAnneal(0.5, 1.0);
    EXPECT_TRUE(bayesCfg.annealing());
}

TEST(ConfigHandlerTest, BayesComputeBeta)
{
    BayesConfigHandler bayesCfg;
    double T0 = 10;
    bayesCfg.setTemperature(T0);
    int maxiter = 100;
    int maxgeneration = 10;
    bayesCfg.setMaxIter(maxiter);
    bayesCfg.setAnneal(0.5, 1.0);

    // beta at iter 0 with no global annealing: should be 1/temp with annealing formula
    double beta = bayesCfg.computeBeta(0, maxgeneration, 0);
    EXPECT_GT(beta, 0.0);

    // When iter >= maxiter, temperature should be very small → beta very large
    double betaEnd = bayesCfg.computeBeta(0, maxgeneration, 100);
    EXPECT_GT(betaEnd, beta);

    // Max value that beta can adopt
    double betaMax = 1e6;
    // Test anneal ending
    bayesCfg.setAnneal(0.5, 0.8);
    beta = bayesCfg.computeBeta(0, maxgeneration, 90);
    EXPECT_NEAR(beta, betaMax, gmx::test::precisionTolerance(1e-8, 1e-6));
    beta = bayesCfg.computeBeta(0, maxgeneration, 65);
    EXPECT_NEAR(beta, 0.2, gmx::test::precisionTolerance(1e-8, 1e-6));

    // Test global annealing
    bayesCfg.setGlobalAnneal(0.0, 1.0);
    beta = bayesCfg.computeBeta(0, maxgeneration, 65);
    // At generation 0 annealing should be off, so starting T.
    EXPECT_NEAR(beta, 2/T0, gmx::test::precisionTolerance(1e-8, 1e-6));

    // Test global annealing more
    bayesCfg.setGlobalAnneal(0.2, 0.8);
    beta = bayesCfg.computeBeta(0, maxgeneration, 65);
    // At generation 0 annealing should be off, so starting T.
    EXPECT_NEAR(beta, 1/T0, gmx::test::precisionTolerance(1e-8, 1e-6));
    beta = bayesCfg.computeBeta(9, maxgeneration, 65);
    // At generation 9 annealing should be off, so final T
    EXPECT_NEAR(beta, betaMax, gmx::test::precisionTolerance(1e-8, 1e-6));

    // Now in between, before local annealing
    beta = bayesCfg.computeBeta(5, maxgeneration, 25);
    EXPECT_NEAR(beta, 0.2, gmx::test::precisionTolerance(1e-8, 1e-6));

    // Now in between, after local annealing
    beta = bayesCfg.computeBeta(5, maxgeneration, 85);
    EXPECT_NEAR(beta, betaMax, gmx::test::precisionTolerance(1e-8, 1e-6));

    // Now in between, during local annealing
    beta = bayesCfg.computeBeta(5, maxgeneration, 65);
    EXPECT_NEAR(beta, 0.4, gmx::test::precisionTolerance(1e-8, 1e-6));
}

// ---- SimulationConfigHandler tests ----

TEST(ConfigHandlerTest, SimulationConfigHandlerDefaults)
{
    const SimulationConfigHandler simCfg;
    EXPECT_EQ(simCfg.nsteps(), 0);
    EXPECT_FALSE(simCfg.minimize());
    EXPECT_FALSE(simCfg.singlePoint());
    EXPECT_EQ(simCfg.maxIter(), 0);
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

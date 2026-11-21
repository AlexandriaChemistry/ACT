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

#include <gtest/gtest.h>
#include <string>
#include <utility>
#include <vector>

#include "act/alexandria/confighandler.h"
#include "act/alexandria/sensitivity.h"
#include "act/basics/msg_handler.h"
#include "act/constraints/shake.h"
#include "act/properties/velocityhandler.h"
#include "gromacs/commandline/filenm.h"
#include "gromacs/commandline/pargs.h"
#include "gromacs/fileio/oenv.h"
#include "gromacs/utility/exceptions.h"
#include "gromacs/utility/real.h"
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

// ---- SimulationConfigHandler from the command line ----
//
// What the command line of minimize, nma and simulate gives, such that a
// change in how the options are turned into settings shows.

//! The option sets of the tools
enum class SimulationOptionSet
{
    //! minimize: the minimization options with -ws
    Minimize,
    //! nma: -sp, -minimize and the minimization options
    Nma,
    //! simulate: as nma, with the MD options
    Simulate
};

/*! \brief Parse a command line into simulation settings as a tool does
 * \param[in] sch  The settings
 * \param[in] set  The options that the tool has
 * \param[in] args The arguments
 */
static void parseSimulationOptions(SimulationConfigHandler         *sch,
                                   SimulationOptionSet              set,
                                   const std::vector<const char *> &args)
{
    std::vector<t_pargs>  pargs;
    std::vector<t_filenm> fnms;
    if (SimulationOptionSet::Minimize == set)
    {
        sch->add_minimize_options(&pargs, true);
    }
    else
    {
        sch->add_options(&pargs, &fnms);
        if (SimulationOptionSet::Simulate == set)
        {
            sch->add_MD_options(&pargs);
        }
    }
    std::vector<const char *> cmdline = { "test" };
    cmdline.insert(cmdline.end(), args.begin(), args.end());
    gmx::test::CommandLine    cmd(cmdline);
    gmx_output_env_t         *oenv = nullptr;
    ASSERT_TRUE(parse_common_args(&cmd.argc(), cmd.argv(), 0, 0, nullptr,
                                  pargs.size(), pargs.data(), 0, nullptr, 0, nullptr, &oenv));
    output_env_done(oenv);
    MsgHandler msghandler;
    sch->check_pargs(&msghandler);
}

/*! \brief Check that settings have the values of a tool without options
 * \param[in] sch The settings
 */
static void expectSimulationDefaults(const SimulationConfigHandler &sch)
{
    EXPECT_FALSE(sch.singlePoint());
    EXPECT_FALSE(sch.minimize());
    EXPECT_EQ(eMinimizeAlgorithm::LBFGS, sch.minAlg());
    EXPECT_EQ(0, sch.maxIter());
    EXPECT_EQ(100, sch.shellMaxIter());
    EXPECT_EQ(1, sch.minimizeRetries());
    EXPECT_DOUBLE_EQ(static_cast<real>(0.002), sch.minimizeDisplacement());
    EXPECT_FALSE(sch.forceReminimize());
    EXPECT_DOUBLE_EQ(static_cast<real>(1e-6), sch.forceTolerance());
    EXPECT_DOUBLE_EQ(1.0, sch.overRelax());
    EXPECT_EQ(0, sch.kFBPR());
    EXPECT_EQ(0, sch.r0FBPR());
    EXPECT_EQ(ConstraintType::No, sch.constraintType());
    EXPECT_FALSE(sch.constrain());
    EXPECT_EQ(1000, sch.constraintMaxIter());
    EXPECT_DOUBLE_EQ(static_cast<real>(1e-6), sch.constraintTolerance());
    EXPECT_FALSE(sch.writeShells());
    EXPECT_EQ(0, sch.nsteps());
    EXPECT_DOUBLE_EQ(static_cast<real>(0.0002), sch.deltat());
    EXPECT_EQ(0, sch.temperature());
    EXPECT_EQ(0, sch.seed());
    EXPECT_EQ(1, sch.nstxout());
    EXPECT_EQ(0, sch.nstvout());
    EXPECT_EQ(1, sch.nstener());
    EXPECT_EQ(COMMRemoval::Angular, sch.commRemoval());
    EXPECT_EQ(100, sch.nstcomm());
}

TEST(ConfigHandlerTest, SimulationOptionsDefaults)
{
    for (auto set : { SimulationOptionSet::Minimize, SimulationOptionSet::Nma, SimulationOptionSet::Simulate })
    {
        SCOPED_TRACE(static_cast<int>(set));
        SimulationConfigHandler sch;
        parseSimulationOptions(&sch, set, {});
        expectSimulationDefaults(sch);
    }
}

TEST(ConfigHandlerTest, SimulationOptionsAllSet)
{
    SimulationConfigHandler sch;
    parseSimulationOptions(&sch, SimulationOptionSet::Simulate,
                           { "-sp", "-minimize", "-minalg", "Steep", "-maxiter", "77",
                             "-shell_max_iter", "33", "-maxretries", "4", "-maxdisplacement", "0.01",
                             "-forcereminimize", "-toler", "0.0001", "-overrelax", "1.5",
                             "-kfbpr", "250", "-r0fbpr", "0.3", "-constrain", "angles",
                             "-constr_maxiter", "55", "-constr_toler", "0.0002",
                             "-nsteps", "1234", "-deltat", "0.0005", "-ws", "-temp", "310",
                             "-seed", "42", "-nstxout", "7", "-nstvout", "8", "-nstener", "9",
                             "-comm", "linear", "-nstcomm", "11" });
    EXPECT_TRUE(sch.singlePoint());
    EXPECT_TRUE(sch.minimize());
    EXPECT_EQ(eMinimizeAlgorithm::Steep, sch.minAlg());
    EXPECT_EQ(77, sch.maxIter());
    EXPECT_EQ(33, sch.shellMaxIter());
    EXPECT_EQ(4, sch.minimizeRetries());
    EXPECT_DOUBLE_EQ(static_cast<real>(0.01), sch.minimizeDisplacement());
    EXPECT_TRUE(sch.forceReminimize());
    EXPECT_DOUBLE_EQ(static_cast<real>(0.0001), sch.forceTolerance());
    EXPECT_DOUBLE_EQ(1.5, sch.overRelax());
    EXPECT_DOUBLE_EQ(250, sch.kFBPR());
    EXPECT_DOUBLE_EQ(static_cast<real>(0.3), sch.r0FBPR());
    EXPECT_EQ(ConstraintType::Angles, sch.constraintType());
    EXPECT_TRUE(sch.constrain());
    EXPECT_EQ(55, sch.constraintMaxIter());
    EXPECT_DOUBLE_EQ(static_cast<real>(0.0002), sch.constraintTolerance());
    EXPECT_EQ(1234, sch.nsteps());
    EXPECT_DOUBLE_EQ(static_cast<real>(0.0005), sch.deltat());
    EXPECT_TRUE(sch.writeShells());
    EXPECT_DOUBLE_EQ(310, sch.temperature());
    EXPECT_EQ(42, sch.seed());
    EXPECT_EQ(7, sch.nstxout());
    EXPECT_EQ(8, sch.nstvout());
    EXPECT_EQ(9, sch.nstener());
    EXPECT_EQ(COMMRemoval::Linear, sch.commRemoval());
    EXPECT_EQ(11, sch.nstcomm());
}

TEST(ConfigHandlerTest, SimulationOptionsMinimizationAlgorithms)
{
    const std::vector<std::pair<const char *, eMinimizeAlgorithm>> choices = {
        { "LBFGS", eMinimizeAlgorithm::LBFGS },
        { "Newton", eMinimizeAlgorithm::Newton },
        { "Steep", eMinimizeAlgorithm::Steep }
    };
    for (const auto &choice : choices)
    {
        for (auto set : { SimulationOptionSet::Minimize, SimulationOptionSet::Nma, SimulationOptionSet::Simulate })
        {
            SimulationConfigHandler sch;
            parseSimulationOptions(&sch, set, { "-minalg", choice.first });
            EXPECT_EQ(choice.second, sch.minAlg()) << choice.first;
        }
    }
}

TEST(ConfigHandlerTest, SimulationOptionsConstraints)
{
    const std::vector<std::pair<const char *, ConstraintType>> choices = {
        { "no", ConstraintType::No },
        { "bonds", ConstraintType::Bonds },
        { "angles", ConstraintType::Angles }
    };
    for (const auto &choice : choices)
    {
        for (auto set : { SimulationOptionSet::Minimize, SimulationOptionSet::Nma, SimulationOptionSet::Simulate })
        {
            SimulationConfigHandler sch;
            parseSimulationOptions(&sch, set, { "-constrain", choice.first });
            EXPECT_EQ(choice.second, sch.constraintType()) << choice.first;
        }
    }
}

TEST(ConfigHandlerTest, SimulationOptionsCenterOfMass)
{
    const std::vector<std::pair<const char *, COMMRemoval>> choices = {
        { "angular", COMMRemoval::Angular },
        { "linear", COMMRemoval::Linear },
        { "none", COMMRemoval::None }
    };
    for (const auto &choice : choices)
    {
        SimulationConfigHandler sch;
        parseSimulationOptions(&sch, SimulationOptionSet::Simulate, { "-comm", choice.first });
        EXPECT_EQ(choice.second, sch.commRemoval()) << choice.first;
    }
}

TEST(ConfigHandlerTest, SimulationOptionsNotKeptBetweenParses)
{
    // The choices of one command line do not carry over to the next one in
    // the same process, for any of the tools
    for (auto set : { SimulationOptionSet::Minimize, SimulationOptionSet::Nma, SimulationOptionSet::Simulate })
    {
        SCOPED_TRACE(static_cast<int>(set));
        SimulationConfigHandler first;
        parseSimulationOptions(&first, SimulationOptionSet::Simulate,
                               { "-minalg", "Newton", "-constrain", "bonds", "-comm", "linear" });
        SimulationConfigHandler second;
        parseSimulationOptions(&second, set, {});
        expectSimulationDefaults(second);
    }
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

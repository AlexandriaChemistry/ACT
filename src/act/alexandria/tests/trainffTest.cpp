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
 * Tests of the train_ff module with the MCMC optimizer on one core, using
 * input derived from the benchmarks in ACTdata/Performance: the
 * intermolecular and intramolecular training of alcohols and water, and
 * the hydrogen fluoride dimer with virtual sites.
 * \author David van der Spoel <david.vanderspoel@icm.uu.se>
 */

#include "actpre.h"

#include <cmath>
#include <limits>
#include <sstream>
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

//! Options for MCMC on a single core, as in the benchmarks
Options mcmcOptions(int maxIter, const std::string &seed = "219463979")
{
    return {
        { "-seed", seed },
        { "-v", "3" },
        { "-optimizer", "MCMC" },
        { "-pop_size", "1" },
        { "-max_iter", std::to_string(maxIter) },
        { "-temp", "0.025" },
        { "-anneal_begin", "0.1" },
        { "-anneal_end", "0.9" },
        { "-sensitivity", "" },
        { "-printSP", "" },
        { "-fc_bound", "10000" },
        { "-ener_boltz_temp", "3000" }
    };
}

//! Intermolecular terms of alcohol and water dimers
Options interOptions()
{
    return {
        { "-fit", "sigma epsilon" },
        { "-fc_allelec", "1" },
        { "-fc_exch", "1" },
        { "-fc_disp", "1" }
    };
}

//! Intramolecular terms of alcohols
Options intraOptions()
{
    return {
        { "-fit", "De b c bondenergy kb bondlength angle kt c1 c2 c3 c4 c5" },
        { "-fc_epot", "1" },
        { "-fc_force", "0.1" }
    };
}

//! Electrostatics, induction and van der Waals of hydrogen fluoride
Options hfOptions()
{
    return {
        { "-fit", "charge eta chi zeta vs2a sigma epsilon gamma delta" },
        { "-fc_allelec", "1" },
        { "-fc_exch", "1" },
        { "-fc_disp", "1" },
        { "-fc_mu", "1" },
        { "-fc_quad", "1" }
    };
}

//! \return the options combined
Options combine(Options a, const Options &b)
{
    a.insert(a.end(), b.begin(), b.end());
    return a;
}

class TrainFFTest : public gmx::test::CommandLineTestBase
{
protected:
    /*! \brief Run train_ff and check the output
     * \param[in] ffName  Force field file
     * \param[in] molName Base name of the molprop and selection file
     * \param[in] options The options
     */
    void runAndCheck(const std::string &ffName,
                     const std::string &molName,
                     const Options     &options)
    {
        auto output = trainfftest::runTrainFF(&fileManager(), "train", ffName, molName, options);
        ASSERT_EQ(0, output.status);
        // MCMC never accepts a worse best parameter vector
        ASSERT_GT(output.initialChi2, 0);
        ASSERT_GT(output.finalChi2, 0);
        EXPECT_LE(output.finalChi2, output.initialChi2);
        // The trained force field is written
        EXPECT_NE(std::string::npos, output.forceField.find("</alexandria_chemistry>"));
        gmx::test::TestReferenceChecker checker(rootChecker());
        checker.setDefaultTolerance(gmx::test::relativeToleranceAsFloatingPoint(1.0, 1e-4));
        trainfftest::checkOutput(&checker, output);
    }
};

// Test: Training Lennard-Jones parameters on alcohol and water dimers
TEST_F(TrainFFTest, MCMCInterBooze)
{
    runAndCheck("lj-pc.xml", "booze-inter", combine(mcmcOptions(10), interOptions()));
}

// Test: Training bonded parameters on the energies and forces of alcohols
TEST_F(TrainFFTest, MCMCIntraBooze)
{
    runAndCheck("halgren-PG.xml", "booze-mono", combine(mcmcOptions(3), intraOptions()));
}

// Test: Training the hydrogen fluoride dimer, with virtual sites, on the
// interaction energy components and the multipoles of the monomer
TEST_F(TrainFFTest, MCMCHydrogenFluoride)
{
    runAndCheck("4S-lj14-7-pg-hvdw.xml", "hf-hf", combine(mcmcOptions(3), hfOptions()));
}

// Test: The files with the correlation of the energy components get the
// prefix given with -compound_file
TEST_F(TrainFFTest, CompoundFilePrefix)
{
    auto prefix  = fileManager().getTemporaryFilePath("hf-");
    auto options = combine(mcmcOptions(1), hfOptions());
    options.push_back({ "-compound_file", prefix });
    auto output = trainfftest::runTrainFF(&fileManager(), "prefix", "4S-lj14-7-pg-hvdw.xml", "hf-hf", options);
    ASSERT_EQ(0, output.status);
    for (const std::string term : { "ALLELEC", "DISPERSION", "EXCHANGE" })
    {
        EXPECT_TRUE(std::filesystem::exists(prefix + term + ".xvg")) << prefix + term + ".xvg";
    }
}

// Test: The outliers are listed with their structure, from the largest
// to the smallest deviation
TEST_F(TrainFFTest, OutlierOverview)
{
    auto output = trainfftest::runTrainFF(&fileManager(), "outliers", "lj-pc.xml", "booze-inter",
                                          combine(mcmcOptions(1), interOptions()));
    ASSERT_EQ(0, output.status);
    int  noutlier = 0;
    bool inTable  = false;
    double previous = 0;
    for (const auto &line : output.log)
    {
        if (line.find("Name") == 0 && line.find("Structure") != std::string::npos)
        {
            inTable  = true;
            previous = std::numeric_limits<double>::max();
            continue;
        }
        if (inTable)
        {
            std::istringstream iss(line);
            std::string        name, structure;
            double             ref, act, diff;
            if (!(iss >> name >> structure >> ref >> act >> diff))
            {
                inTable = false;
                continue;
            }
            EXPECT_NEAR(act - ref, diff, 1e-3*std::abs(diff) + 1e-3) << line;
            EXPECT_LE(std::abs(diff), previous) << line;
            previous  = std::abs(diff);
            noutlier += 1;
        }
    }
    EXPECT_GT(noutlier, 1);
}

// Test: The same seed gives the same parameters, a different one does not
TEST_F(TrainFFTest, Reproducible)
{
    std::vector<nlohmann::json> params;
    int                         run = 0;
    for (const std::string seed : { "1993", "1993", "1997" })
    {
        auto options = combine(mcmcOptions(3, seed), hfOptions());
        auto output = trainfftest::runTrainFF(&fileManager(), "seed" + std::to_string(run++),
                                              "4S-lj14-7-pg-hvdw.xml", "hf-hf", options);
        ASSERT_EQ(0, output.status);
        params.push_back(output.json["train_ff"]["parameters"]);
    }
    EXPECT_EQ(params[0], params[1]);
    EXPECT_NE(params[0], params[2]);
}

} // namespace

} // namespace alexandria

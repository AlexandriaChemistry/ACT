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
/*! \internal \file
 * \brief
 * Tests of the settings of the generator of dimers, as the command lines of
 * alexandria b2 and simulate give them and as a program sets them.
 * \author David van der Spoel <david.vanderspoel@icm.uu.se>
 */
#include "actpre.h"

#include <gtest/gtest.h>
#include <string>
#include <utility>
#include <vector>

#include "act/properties/dimergenerator.h"
#include "gromacs/commandline/filenm.h"
#include "gromacs/commandline/pargs.h"
#include "gromacs/fileio/oenv.h"
#include "gromacs/utility/exceptions.h"
#include "testutils/cmdlinetest.h"
#include "testutils/testfilemanager.h"

namespace alexandria
{

namespace
{

/*! \brief Parse a command line into a generator of dimers as the tools do
 * \param[in] gendimers The generator
 * \param[in] args      The arguments
 * \param[in] cmd       The command line, which must exist as long as it is
 *                      parsed
 */
void parseDimerOptions(DimerGenerator                 *gendimers,
                       const std::vector<std::string> &args,
                       gmx::test::CommandLine         *cmd)
{
    std::vector<t_pargs>      pargs;
    std::vector<t_filenm>     fnms;
    std::vector<const char *> desc;
    gendimers->addOptions(&pargs, &fnms, &desc);
    cmd->append("test");
    for (const auto &a : args)
    {
        cmd->append(a);
    }
    gmx_output_env_t *oenv = nullptr;
    ASSERT_TRUE(parse_common_args(&cmd->argc(), cmd->argv(), 0, fnms.size(), fnms.data(),
                                  pargs.size(), pargs.data(), 0, nullptr, 0, nullptr, &oenv));
    output_env_done(oenv);
    gendimers->finishOptions(fnms);
}

TEST(DimerGeneratorOptionsTest, Defaults)
{
    DimerGenerator         gendimers;
    gmx::test::CommandLine cmd;
    parseDimerOptions(&gendimers, {}, &cmd);
    EXPECT_EQ(RandAlg::Sobol, gendimers.randAlg());
    EXPECT_EQ(1993, gendimers.seed());
    EXPECT_FALSE(gendimers.flexible());
    EXPECT_DOUBLE_EQ(static_cast<real>(0.1), gendimers.mindist());
    EXPECT_DOUBLE_EQ(static_cast<real>(4.0), gendimers.maxdist());
    // Bins of 2.5 pm from the shortest to the longest distance
    EXPECT_EQ(1 + std::lround((gendimers.maxdist() - gendimers.mindist())/0.0025), gendimers.ndist());
    EXPECT_FALSE(gendimers.hasTrajectory());
    EXPECT_EQ(nullptr, gendimers.outCoords());
}

TEST(DimerGeneratorOptionsTest, AllSet)
{
    const auto             traj = gmx::test::TestFileManager::getInputFilePath("molecules/water_dimer.pdb");
    DimerGenerator         gendimers;
    gmx::test::CommandLine cmd;
    parseDimerOptions(&gendimers, { "-ndist", "11", "-mindist", "0.2", "-maxdist", "1.2", "-dimerseed", "7",
                                    "-randalg", "Pseudo", "-flex", "-traj", traj, "-ox", "out.xyz" },
                      &cmd);
    EXPECT_EQ(11, gendimers.ndist());
    EXPECT_DOUBLE_EQ(static_cast<real>(0.2), gendimers.mindist());
    EXPECT_DOUBLE_EQ(static_cast<real>(1.2), gendimers.maxdist());
    EXPECT_NEAR(0.1, gendimers.binwidth(), 1e-6);
    EXPECT_EQ(7, gendimers.seed());
    EXPECT_EQ(RandAlg::Pseudo, gendimers.randAlg());
    EXPECT_TRUE(gendimers.flexible());
    EXPECT_TRUE(gendimers.hasTrajectory());
    EXPECT_EQ(traj, std::string(gendimers.trajname()));
    ASSERT_NE(nullptr, gendimers.outCoords());
    EXPECT_EQ("out.xyz", std::string(gendimers.outCoords()));
}

TEST(DimerGeneratorOptionsTest, RandomAlgorithms)
{
    const std::vector<std::pair<const char *, RandAlg>> choices = {
        { "Sobol", RandAlg::Sobol }, { "SobolOwen", RandAlg::SobolOwen },
        { "Pseudo", RandAlg::Pseudo }, { "Gauss", RandAlg::Gauss }
    };
    for (const auto &choice : choices)
    {
        DimerGenerator         gendimers;
        gmx::test::CommandLine cmd;
        parseDimerOptions(&gendimers, { "-randalg", choice.first }, &cmd);
        EXPECT_EQ(choice.second, gendimers.randAlg()) << choice.first;
    }
}

TEST(DimerGeneratorOptionsTest, NotKeptBetweenParses)
{
    // The choice of one command line does not carry over to the next one
    DimerGenerator first;
    {
        gmx::test::CommandLine cmd;
        parseDimerOptions(&first, { "-randalg", "Pseudo" }, &cmd);
    }
    DimerGenerator         second;
    gmx::test::CommandLine cmd;
    parseDimerOptions(&second, {}, &cmd);
    EXPECT_EQ(RandAlg::Sobol, second.randAlg());
}

TEST(DimerGeneratorOptionsTest, CopyOutlivesTheCommandLine)
{
    // The settings of a generator belong to it, not to the command line
    const auto     traj = gmx::test::TestFileManager::getInputFilePath("molecules/water_dimer.pdb");
    DimerGenerator copy;
    {
        DimerGenerator         gendimers;
        gmx::test::CommandLine cmd;
        parseDimerOptions(&gendimers, { "-traj", traj, "-ox", "out.xyz", "-randalg", "Pseudo" }, &cmd);
        copy = gendimers;
    }
    EXPECT_TRUE(copy.hasTrajectory());
    EXPECT_EQ(traj, std::string(copy.trajname()));
    EXPECT_EQ("out.xyz", std::string(copy.outCoords()));
    // A copy makes the same orientations as the original
    DimerGenerator original;
    original.setRandAlg(RandAlg::Pseudo);
    original.setSeed(11);
    original.finishSettings();
    DimerGenerator twin = original;
    original.generateRandomNumbers(16);
    twin.generateRandomNumbers(16);
    EXPECT_EQ(original.allRandom(), twin.allRandom());
}

TEST(DimerGeneratorOptionsTest, WithoutCommandLine)
{
    // As a program sets it up
    DimerGenerator gendimers;
    gendimers.setDistances(0.2, 1.2, 11);
    gendimers.setTrajectory("dimers.trr");
    gendimers.setOutputCoordinates("");
    gendimers.finishSettings();
    EXPECT_EQ(11, gendimers.ndist());
    EXPECT_NEAR(0.1, gendimers.binwidth(), 1e-6);
    EXPECT_TRUE(gendimers.hasTrajectory());
    EXPECT_EQ(nullptr, gendimers.outCoords());
}

TEST(DimerGeneratorOptionsTest, WrongSettings)
{
    // Each is refused, with the name of the option in the message
    const std::vector<std::pair<std::vector<std::string>, std::string>> wrong = {
        { { "-ndist", "-2" }, "-ndist" },
        { { "-mindist", "-0.1" }, "-mindist" },
        { { "-mindist", "1", "-maxdist", "0.5" }, "-maxdist" },
        { { "-randalg", "Gauss", "-nalpha", "0" }, "-nalpha" },
        { { "-randalg", "Gauss", "-nbeta", "-1" }, "-nbeta" }
    };
    for (const auto &w : wrong)
    {
        DimerGenerator         gendimers;
        gmx::test::CommandLine cmd;
        try
        {
            parseDimerOptions(&gendimers, w.first, &cmd);
            ADD_FAILURE() << "Accepted a wrong " << w.second;
        }
        catch (const gmx::InvalidInputError &ex)
        {
            EXPECT_NE(std::string::npos, std::string(ex.what()).find(w.second)) << ex.what();
        }
    }
    DimerGenerator set;
    set.setDistances(1, 0.5, 0);
    EXPECT_THROW(set.finishSettings(), gmx::InvalidInputError);
}

} // namespace

} // namespace alexandria

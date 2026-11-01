/*
 * This source file is part of the Alexandria program.
 *
 * Copyright (C) 2023,2026
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
 * Tests of the min_complex module.
 */

#include "actpre.h"

#include <filesystem>
#include <fstream>
#include <sstream>
#include <string>
#include <vector>

#include <gtest/gtest.h>

#include "act/alexandria/alex_modules.h"
#include "act/molprop/molprop.h"
#include "act/molprop/molprop_xml.h"

#include "trainffTestUtil.h"

#include "testutils/cmdlinetest.h"
#include "testutils/testfilemanager.h"

namespace alexandria
{

namespace
{

class MinComplexTest : public gmx::test::CommandLineTestBase
{
protected:
    /*! \brief Run min_complex in the temporary directory
     * \param[in]  ffName    Force field
     * \param[in]  molprops  The complexes
     * \param[out] log       The contents of the log file
     * \return the exit status
     */
    int run(const std::string &ffName, const std::string &molprops, std::string *log)
    {
        auto logName = fileManager().getTemporaryFilePath("min_complex.log");
        gmx::test::CommandLine cmd;
        cmd.append("min_complex");
        cmd.addOption("-ff", fileManager().getInputFilePath(ffName));
        cmd.addOption("-mp", molprops);
        cmd.addOption("-g", logName);
        int status;
        {
            // The structures are written to the working directory
            trainfftest::ScopedWorkingDirectory cwd(std::filesystem::path(logName).parent_path());
            status = min_complex(cmd.argc(), cmd.argv());
        }
        std::ifstream     in(logName);
        std::stringstream text;
        text << in.rdbuf();
        *log = text.str();
        return status;
    }

    //! \return a file with only the complexes of a file with the name
    std::string select(const std::string &fileName, const std::string &molname)
    {
        std::vector<MolProp> all, selected;
        MolPropRead(nullptr, fileManager().getInputFilePath(fileName).c_str(), &all);
        for (auto &mp : all)
        {
            if (mp.getMolname() == molname)
            {
                selected.push_back(std::move(mp));
            }
        }
        EXPECT_FALSE(selected.empty());
        auto out = fileManager().getTemporaryFilePath("selected.xml");
        MolPropWrite(out.c_str(), selected, false);
        return out;
    }
};

// Test: The water dimer is minimized near its structure, with a negative
// interaction energy, and written to a file named after it
TEST_F(MinComplexTest, WaterDimer)
{
    std::string log;
    auto input = select("trainff/booze-inter.xml", "water#water");
    ASSERT_EQ(0, run("trainff/halgren-PG.xml", input, &log));
    auto pos = log.find("water#water final energy:");
    ASSERT_NE(std::string::npos, pos) << log;
    double epot = 0, einter = 0, rmsd = 0;
    ASSERT_EQ(3, std::sscanf(log.c_str() + pos,
                             "water#water final energy: %lf. Interaction energy: %lf. RMSD wrt original structure %lf",
                             &epot, &einter, &rmsd));
    EXPECT_LT(einter, 0);
    EXPECT_LT(rmsd, 0.1);
    auto pdb = std::filesystem::path(fileManager().getTemporaryFilePath("x")).parent_path() / "water#water.pdb";
    EXPECT_TRUE(std::filesystem::exists(pdb));
    std::filesystem::remove(pdb);
}

// Test: Complexes without a topology are skipped, which is reported with
// the exit status
TEST_F(MinComplexTest, MissingAtomTypes)
{
    std::string log;
    // ACS-g has no atom type for the hydrogens of water
    auto input = select("trainff/booze-inter.xml", "water#water");
    EXPECT_EQ(1, run("ACS-g.xml", input, &log));
    EXPECT_NE(std::string::npos, log.find("Skipping 'water#water'"));
}

} // namespace

} // namespace alexandria

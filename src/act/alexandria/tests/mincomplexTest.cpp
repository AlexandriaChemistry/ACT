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
#include <map>
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
        const auto logName = fileManager().getTemporaryFilePath("min_complex.log");
        gmx::test::CommandLine cmd;
        cmd.append("min_complex");
        cmd.addOption("-ff", fileManager().getInputFilePath(ffName));
        cmd.addOption("-mp", molprops);
        cmd.addOption("-g", logName);
        int status;
        {
            // The structures are written to the working directory
            const trainfftest::ScopedWorkingDirectory cwd(std::filesystem::path(logName).parent_path());
            status = min_complex(cmd.argc(), cmd.argv());
        }
        const std::ifstream in(logName);
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
    const auto input = select("molprop/booze-inter.xml", "water#water");
    ASSERT_EQ(0, run("forcefields/halgren-PG.xml", input, &log));
    auto pos = log.find("water#water final energy:");
    ASSERT_NE(std::string::npos, pos) << log;
    double epot = 0, einter = 0, rmsd = 0;
    ASSERT_EQ(3, std::sscanf(log.c_str() + pos,
                             "water#water final energy: %lf. Interaction energy: %lf. RMSD wrt original structure %lf",
                             &epot, &einter, &rmsd));
    EXPECT_LT(einter, 0);
    EXPECT_LT(rmsd, 0.1);
    const auto pdb = std::filesystem::path(fileManager().getTemporaryFilePath("x")).parent_path() / "water#water.pdb";
    EXPECT_TRUE(std::filesystem::exists(pdb));
    std::filesystem::remove(pdb);
}

/*! \brief The results in a log file of min_complex
 * \param[in] log The contents of the log file
 * \return for each complex the interaction energy and the RMSD
 */
static std::map<std::string, std::pair<double, double>> results(const std::string &log)
{
    std::map<std::string, std::pair<double, double>> res;
    std::istringstream in(log);
    for (std::string line; std::getline(in, line); )
    {
        auto pos = line.find(" final energy:");
        if (pos == std::string::npos)
        {
            continue;
        }
        double epot, einter, rmsd;
        if (3 == std::sscanf(line.c_str() + pos,
                             " final energy: %lf. Interaction energy: %lf. RMSD wrt original structure %lf",
                             &epot, &einter, &rmsd))
        {
            res[line.substr(0, pos)] = { einter, rmsd };
        }
    }
    return res;
}

// Test: With a trained force field the complexes of alcohols and water are
// bound near their structures
TEST_F(MinComplexTest, TrainedForceField)
{
    std::string log;
    ASSERT_EQ(0, run("forcefields/lj-pc-trained.xml",
                     fileManager().getInputFilePath("molprop/booze-inter.xml"), &log));
    const auto res = results(log);
    ASSERT_EQ(3u, res.size()) << log;
    for (const auto &r : res)
    {
        EXPECT_LT(r.second.first, 0) << r.first;
        EXPECT_LT(r.second.second, 0.2) << r.first;
        const auto pdb = std::filesystem::path(fileManager().getTemporaryFilePath("x")).parent_path() / (r.first + ".pdb");
        EXPECT_TRUE(std::filesystem::exists(pdb)) << pdb;
        std::filesystem::remove(pdb);
    }
}

// Test: Newton steps and steepest descent bind the complexes as well, at the
// same minima, steepest descent needs thousands of steps. In single precision
// steepest descent does not converge within 10000 steps, since its step
// lengths follow from differences of the forces, which near the minimum are
// dominated by rounding, so it is only tested in double precision.
TEST_F(MinComplexTest, NewtonAndSteepestDescent)
{
    std::map<std::string, std::map<std::string, std::pair<double, double>>> all;
    std::vector<const char *> algorithms = { "Newton" };
    if (GMX_DOUBLE)
    {
        algorithms.push_back("Steep");
    }
    for (const char *alg : algorithms)
    {
        SCOPED_TRACE(alg);
        const auto logName = fileManager().getTemporaryFilePath(std::string(alg) + ".log");
        gmx::test::CommandLine cmd;
        cmd.append("min_complex");
        cmd.addOption("-ff", fileManager().getInputFilePath("forcefields/lj-pc-trained.xml"));
        cmd.addOption("-mp", fileManager().getInputFilePath("molprop/booze-inter.xml"));
        cmd.addOption("-g", logName);
        cmd.addOption("-minalg", alg);
        cmd.addOption("-maxiter", "10000");
        int status;
        {
            const trainfftest::ScopedWorkingDirectory cwd(std::filesystem::path(logName).parent_path());
            status = min_complex(cmd.argc(), cmd.argv());
        }
        EXPECT_EQ(0, status);
        const std::ifstream in(logName);
        std::stringstream text;
        text << in.rdbuf();
        all[alg] = results(text.str());
        EXPECT_EQ(3u, all[alg].size());
        for (const auto &r : all[alg])
        {
            EXPECT_LT(r.second.first, 0) << r.first;
            EXPECT_LT(r.second.second, 0.2) << r.first;
            std::filesystem::remove(std::filesystem::path(logName).parent_path() / (r.first + ".pdb"));
        }
    }
    if (GMX_DOUBLE)
    {
        for (const auto &r : all["Newton"])
        {
            EXPECT_NEAR(r.second.first, all["Steep"][r.first].first, 1e-2) << r.first;
        }
    }
}

// Test: A minimizer does not succeed far from the structure, where there are
// no forces, with the untrained force field. Steps that increase the energy
// moved atoms there before.
TEST_F(MinComplexTest, NewtonDoesNotDiverge)
{
    const auto logName = fileManager().getTemporaryFilePath("min_complex.log");
    gmx::test::CommandLine cmd;
    cmd.append("min_complex");
    cmd.addOption("-ff", fileManager().getInputFilePath("forcefields/halgren-PG.xml"));
    cmd.addOption("-mp", fileManager().getInputFilePath("molprop/booze-inter.xml"));
    cmd.addOption("-g", logName);
    cmd.addOption("-minalg", "Newton");
    // The molecules separate slowly where there are almost no forces,
    // until convergence that takes thousands of steps
    cmd.addOption("-maxiter", "100");
    {
        const trainfftest::ScopedWorkingDirectory cwd(std::filesystem::path(logName).parent_path());
        (void) min_complex(cmd.argc(), cmd.argv());
    }
    const std::ifstream in(logName);
    std::stringstream text;
    text << in.rdbuf();
    for (const auto &r : results(text.str()))
    {
        EXPECT_LT(r.second.second, 1.0) << r.first;
        const auto pdb = std::filesystem::path(logName).parent_path() / (r.first + ".pdb");
        std::filesystem::remove(pdb);
    }
}

// Test: Complexes without a topology are skipped, which is reported with
// the exit status
TEST_F(MinComplexTest, MissingAtomTypes)
{
    std::string log;
    // ACS-g has no atom type for the hydrogens of water
    const auto input = select("molprop/booze-inter.xml", "water#water");
    EXPECT_EQ(1, run("forcefields/ACS-g1.xml", input, &log));
    EXPECT_NE(std::string::npos, log.find("Skipping 'water#water'"));
}

} // namespace

} // namespace alexandria

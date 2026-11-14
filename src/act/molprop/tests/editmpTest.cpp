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
 * Tests of the edit_mp module.
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
#include "act/basics/msg_handler.h"
#include "act/molprop/molprop.h"
#include "act/molprop/molprop_xml.h"

#include "testutils/cmdlinetest.h"
#include "testutils/testfilemanager.h"

namespace alexandria
{

namespace
{

class EditMPTest : public gmx::test::CommandLineTestBase
{
protected:
    //! The log file of the last run
    std::string log_;

    /*! \brief Run edit_mp
     * \param[in] inputs  The molprop files
     * \param[in] options The other options
     * \return the molecules in the output file
     */
    std::vector<MolProp> run(const std::vector<std::string>                        &inputs,
                             const std::vector<std::pair<std::string, std::string>> &options = {})
    {
        const auto output = fileManager().getTemporaryFilePath("out.xml");
        log_        = fileManager().getTemporaryFilePath("edit_mp.log");
        gmx::test::CommandLine cmd;
        cmd.append("edit_mp");
        cmd.addOption("-mp");
        for (const auto &in : inputs)
        {
            cmd.append(in);
        }
        cmd.addOption("-o", output);
        cmd.addOption("-g", log_);
        for (const auto &o : options)
        {
            if (o.second.empty())
            {
                cmd.append(o.first);
            }
            else
            {
                cmd.addOption(o.first.c_str(), o.second);
            }
        }
        {
            // Some output is written in the working directory
            const auto old = std::filesystem::current_path();
            std::filesystem::current_path(std::filesystem::path(output).parent_path());
            EXPECT_EQ(0, edit_mp(cmd.argc(), cmd.argv()));
            std::filesystem::current_path(old);
        }
        std::vector<MolProp> mps;
        EXPECT_TRUE(std::filesystem::exists(output));
        MolPropRead(nullptr, output.c_str(), &mps);
        return mps;
    }

    //! \return the molecules in a file
    static std::vector<MolProp> read(const std::string &fileName)
    {
        std::vector<MolProp> mps;
        MolPropRead(nullptr, fileName.c_str(), &mps);
        return mps;
    }

    //! \return the number of experiments of each molecule
    static std::map<std::string, size_t> experiments(const std::vector<MolProp> &mps)
    {
        std::map<std::string, size_t> n;
        for (const auto &mp : mps)
        {
            n[mp.getMolname()] += mp.experimentConst().size();
        }
        return n;
    }
};

// Test: Writing the molecules of a file gives the same molecules
TEST_F(EditMPTest, RoundTrip)
{
    const auto input = fileManager().getInputFilePath("molprop/molprop.xml");
    auto ref   = read(input);
    auto mps   = run({ input });
    ASSERT_EQ(ref.size(), mps.size());
    EXPECT_EQ(experiments(ref), experiments(mps));
    for (size_t i = 0; i < ref.size(); i++)
    {
        EXPECT_EQ(ref[i].getMolname(), mps[i].getMolname());
        EXPECT_EQ(ref[i].NAtom(), mps[i].NAtom());
        EXPECT_EQ(ref[i].formula(), mps[i].formula());
    }
}

// Test: Molecules of different files are combined, sorted by name
TEST_F(EditMPTest, MergeFiles)
{
    const auto first  = fileManager().getInputFilePath("molprop/molprop.xml");
    const auto second = fileManager().getInputFilePath("molprop/booze-mono.xml");
    const auto mps    = run({ first, second });
    std::vector<std::string> names;
    for (const auto &mp : mps)
    {
        names.push_back(mp.getMolname());
    }
    const std::vector<std::string> expected = { "1-butanol", "acetone", "ethanol", "methanol", "water" };
    EXPECT_EQ(expected, names);
}

// Test: The same molecule in two files is one molecule with the
// experiments of both
TEST_F(EditMPTest, MergeDuplicates)
{
    const auto input = fileManager().getInputFilePath("molprop/molprop.xml");
    const auto ref   = experiments(read(input));
    const auto mps   = run({ input, input });
    auto n     = experiments(mps);
    ASSERT_EQ(ref.size(), mps.size());
    for (const auto &r : ref)
    {
        EXPECT_EQ(2*r.second, n[r.first]) << r.first;
    }
}

// Test: With a force field the compounds are checked, the results are in
// the log file, and the output has all their data
TEST_F(EditMPTest, CheckWithForceField)
{
    const auto input = fileManager().getInputFilePath("molprop/molprop.xml");
    const auto mps   = run({ input }, { { "-ff", fileManager().getInputFilePath("forcefields/ACS-g1.xml") } });
    // Checking does not change the compounds
    EXPECT_EQ(read(input).size(), mps.size());
    EXPECT_EQ(experiments(read(input)), experiments(mps));
    const std::ifstream in(log_);
    std::stringstream text;
    text << in.rdbuf();
    EXPECT_NE(std::string::npos, text.str().find("Force field file")) << text.str();
}

// Test: The histograms of the interaction energies are written for each
// compound that has them, the range is reported in the log file
TEST_F(EditMPTest, EnergyHistogram)
{
    const auto input = fileManager().getInputFilePath("molprop/hf-hf.xml");
    run({ input }, { { "-ehisto", "" } });
    const auto histo = std::filesystem::path(fileManager().getTemporaryFilePath("out.xml")).parent_path() /
        "hydrogen-fluoride#hydrogen-fluoride-InteractionEnergy.xvg";
    std::ifstream in(histo);
    ASSERT_TRUE(in.good()) << histo;
    int    npoints = 0;
    double total   = 0;
    for (std::string line; std::getline(in, line); )
    {
        if (line.empty() || line[0] == '#' || line[0] == '@')
        {
            continue;
        }
        std::istringstream iss(line);
        double x, y;
        ASSERT_TRUE(static_cast<bool>(iss >> x >> y)) << line;
        npoints++;
        total += y;
    }
    // Each energy is counted in one bin
    EXPECT_GT(npoints, 0);
    EXPECT_GT(total, 0);
    const std::ifstream logIn(log_);
    std::stringstream text;
    text << logIn.rdbuf();
    EXPECT_NE(std::string::npos, text.str().find("Range of energies for hydrogen-fluoride#hydrogen-fluoride")) << text.str();
    std::filesystem::remove(histo);
}

} // namespace

} // namespace alexandria

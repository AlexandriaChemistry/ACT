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
 * Tests of the geometry_ff module.
 */

#include "actpre.h"

#include <algorithm>
#include <filesystem>
#include <fstream>
#include <gtest/gtest.h>
#include <nlohmann/json.hpp>
#include <set>
#include <string>

#include "act/alexandria/alex_modules.h"
#include "act/basics/interactiontype.h"
#include "act/forcefield/forcefield.h"
#include "act/forcefield/forcefield_parameter.h"
#include "act/forcefield/forcefield_xml.h"
#include "testutils/cmdlinetest.h"
#include "testutils/testfilemanager.h"

namespace alexandria
{

namespace
{

class GeometryFFTest : public gmx::test::CommandLineTestBase
{
};

//! \return the parameters of a bond, given as text
static const ForceFieldParameterMap *bond(const ForceField &pd, const std::string &id)
{
    for (const auto &p : pd.findForcesConst(InteractionType::BONDS).parametersConst())
    {
        if (p.first.id() == id)
        {
            return &p.second;
        }
    }
    return nullptr;
}

// Test: The bond lengths of the compounds, with their standard deviation and
// number, are the new parameters, other bonds are removed, and the energies
// of harmonic bonds are negative
TEST_F(GeometryFFTest, AlcoholAndAcetone)
{
    const auto output = fileManager().getTemporaryFilePath("geometry.xml");
    gmx::test::CommandLine cmd;
    cmd.append("geometry_ff");
    cmd.addOption("-mp", fileManager().getInputFilePath("molprop/molprop.xml"));
    cmd.addOption("-ff", fileManager().getInputFilePath("forcefields/ACS-g1.xml"));
    cmd.addOption("-o", output);
    cmd.addOption("-g", fileManager().getTemporaryFilePath("geometry.log"));
    ASSERT_EQ(0, geometry_ff(cmd.argc(), cmd.argv()));
    ForceField pd;
    readForceField(output, &pd);
    auto ch = bond(pd, "c3_b~h_b");
    ASSERT_NE(nullptr, ch);
    // Nine C-H bonds in each of two structures of butanol, six in acetone
    auto &length = ch->find("bondlength")->second;
    EXPECT_EQ(24u, length.ntrain());
    EXPECT_NEAR(109.4, length.value(), 0.05);
    EXPECT_NEAR(0.5, length.uncertainty(), 0.05);
    EXPECT_LE(length.minimum(), length.value());
    EXPECT_GE(length.maximum(), length.value());
    auto &energy = ch->find("bondenergy")->second;
    EXPECT_LT(energy.value(), 0);
    EXPECT_FALSE(energy.nonNegative());
    EXPECT_LE(energy.minimum(), energy.value());
    EXPECT_GE(energy.maximum(), energy.value());
    // The bonds are made from the compounds only
    EXPECT_EQ(nullptr, bond(pd, "c3_b~cl_b"));
}


// Test: With -json the number of compounds, what was extracted and the force
// field that was written are described in JSON, and not without that option
TEST_F(GeometryFFTest, Json)
{
    const auto output   = fileManager().getTemporaryFilePath("geometry.xml");
    const auto jsonFile = fileManager().getTemporaryFilePath("geometry.json");
    std::filesystem::remove(jsonFile);
    std::filesystem::remove("geometry_ff.json");
    gmx::test::CommandLine cmd;
    cmd.append("geometry_ff");
    cmd.addOption("-mp", fileManager().getInputFilePath("molprop/molprop.xml"));
    cmd.addOption("-ff", fileManager().getInputFilePath("forcefields/ACS-g1.xml"));
    cmd.addOption("-o", output);
    cmd.addOption("-g", fileManager().getTemporaryFilePath("geometry.log"));
    {
        gmx::test::CommandLine without(cmd);
        ASSERT_EQ(0, geometry_ff(without.argc(), without.argv()));
        EXPECT_FALSE(std::filesystem::exists("geometry_ff.json"));
    }
    cmd.addOption("-json", jsonFile);
    ASSERT_EQ(0, geometry_ff(cmd.argc(), cmd.argv()));
    ForceField pd;
    readForceField(output, &pd);
    std::ifstream in(jsonFile);
    ASSERT_TRUE(in.good());
    const auto j = nlohmann::json::parse(in);
    ASSERT_TRUE(j.contains("geometry_ff"));
    const auto &root = j.at("geometry_ff");
    // Two structures of butanol and one of acetone
    EXPECT_EQ(3, root.at("compounds").get<int>());
    const auto &extracted = root.at("extracted");
    ASSERT_TRUE(extracted.is_object());
    const auto bondName = interactionTypeToString(InteractionType::BONDS);
    ASSERT_TRUE(extracted.contains(bondName));
    const auto nbonds = pd.findForcesConst(InteractionType::BONDS).parametersConst().size();
    EXPECT_EQ(nbonds, extracted.at(bondName).get<size_t>());
    // The force field is the one of the file
    const auto &ff = root.at("forcefield");
    EXPECT_EQ(pd.checkSum(), ff.at("checksum").get<std::string>());
    bool found = false;
    for (const auto &jfs : ff.at("interactions"))
    {
        if (jfs.at("type").get<std::string>() != bondName)
        {
            continue;
        }
        EXPECT_EQ(nbonds, jfs.at("parameterlists").size());
        for (const auto &jlist : jfs.at("parameterlists"))
        {
            if (jlist.at("identifier").get<std::string>() != "c3_b~h_b")
            {
                continue;
            }
            for (const auto &jp : jlist.at("parameters"))
            {
                if (jp.at("type").get<std::string>() == "bondlength")
                {
                    found = true;
                    EXPECT_EQ(24u, jp.at("ntrain").get<unsigned int>());
                    EXPECT_NEAR(109.4, jp.at("value").get<double>(), 0.05);
                }
            }
        }
    }
    EXPECT_TRUE(found);
    EXPECT_TRUE(root.at("validation").at("issues").is_array());
}

// Test: With -json the histograms of the bond lengths and angles are there,
// and with -samples each bin says, for some of its values, the compound,
// the calculation and the atoms they come from
TEST_F(GeometryFFTest, Histograms)
{
    const auto jsonFile = fileManager().getTemporaryFilePath("histograms.json");
    gmx::test::CommandLine cmd;
    cmd.append("geometry_ff");
    cmd.addOption("-mp", fileManager().getInputFilePath("molprop/molprop.xml"));
    cmd.addOption("-ff", fileManager().getInputFilePath("forcefields/ACS-g1.xml"));
    cmd.addOption("-o", fileManager().getTemporaryFilePath("histograms.xml"));
    cmd.addOption("-g", fileManager().getTemporaryFilePath("histograms.log"));
    cmd.addOption("-json", jsonFile);
    {
        // Without -samples: the numbers only
        gmx::test::CommandLine counts(cmd);
        ASSERT_EQ(0, geometry_ff(counts.argc(), counts.argv()));
        std::ifstream in(jsonFile);
        const auto j = nlohmann::json::parse(in);
        const auto &bonds = j.at("geometry_ff").at("histograms").at("BONDS");
        ASSERT_FALSE(bonds.empty());
        for (const auto &bin : bonds[0].at("bins"))
        {
            EXPECT_FALSE(bin.contains("samples"));
        }
    }
    cmd.addOption("-samples", 3);
    ASSERT_EQ(0, geometry_ff(cmd.argc(), cmd.argv()));
    std::ifstream in(jsonFile);
    const auto j    = nlohmann::json::parse(in);
    const auto &all = j.at("geometry_ff").at("histograms");
    ASSERT_TRUE(all.contains("BONDS"));
    ASSERT_TRUE(all.contains("ANGLES"));
    const nlohmann::json *ch = nullptr;
    for (const auto &h : all.at("BONDS"))
    {
        if (h.at("identifier").get<std::string>() == "c3_b~h_b")
        {
            ch = &h;
        }
    }
    ASSERT_NE(nullptr, ch);
    // Nine C-H bonds in each of two structures of butanol, six in acetone
    EXPECT_EQ(24, ch->at("count").get<int>());
    EXPECT_EQ("pm", ch->at("unit").get<std::string>());
    EXPECT_DOUBLE_EQ(1.0, ch->at("binwidth").get<double>());
    int    total = 0;
    double last  = -1;
    std::set<std::string> molecules;
    for (const auto &bin : ch->at("bins"))
    {
        const auto count = bin.at("count").get<int>();
        const auto value = bin.at("value").get<double>();
        EXPECT_GT(value, last);
        last   = value;
        total += count;
        const auto &samples = bin.at("samples");
        EXPECT_EQ(std::min(count, 3), static_cast<int>(samples.size()));
        EXPECT_EQ(count - static_cast<int>(samples.size()), bin.at("more").get<int>());
        for (const auto &sample : samples)
        {
            molecules.insert(sample.at("molecule").get<std::string>());
            EXPECT_FALSE(sample.at("calculation").get<std::string>().empty());
            const auto atoms = sample.at("atoms").get<std::vector<int>>();
            ASSERT_EQ(2u, atoms.size());
            EXPECT_GE(atoms[0], 1);
            // The names of the atoms with their numbers, a C and an H
            const auto names = sample.at("names").get<std::string>();
            const bool carbonFirst = names.rfind("C", 0) == 0;
            EXPECT_NE(std::string::npos, names.find((carbonFirst ? "C" : "H") + std::to_string(atoms[0]))) << names;
            EXPECT_NE(std::string::npos, names.find((carbonFirst ? "H" : "C") + std::to_string(atoms[1]))) << names;
            EXPECT_NEAR(value, sample.at("value").get<double>(), 0.5 + 1e-9);
        }
    }
    EXPECT_EQ(24, total);
    EXPECT_FALSE(molecules.empty());
    // The angles are in degrees
    EXPECT_EQ("degree", all.at("ANGLES")[0].at("unit").get<std::string>());
}

} // namespace

} // namespace alexandria

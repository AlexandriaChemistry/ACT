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
#include <vector>

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

// Test: The atoms that a value comes from are numbered as in the molprop
// file, also with a force field that puts a virtual site on bromine: the
// topology has it after the bromine, the molprop file does not have it
TEST_F(GeometryFFTest, AtomNumbersWithVirtualSites)
{
    const auto jsonFile = fileManager().getTemporaryFilePath("vsite.json");
    gmx::test::CommandLine cmd;
    cmd.append("geometry_ff");
    cmd.addOption("-mp", fileManager().getInputFilePath("molprop/hbr-vsite.xml"));
    cmd.addOption("-ff", fileManager().getInputFilePath("forcefields/ACS-pg-vs2.xml"));
    cmd.addOption("-o", fileManager().getTemporaryFilePath("vsite.xml"));
    cmd.addOption("-g", fileManager().getTemporaryFilePath("vsite.log"));
    cmd.addOption("-json", jsonFile);
    cmd.addOption("-samples", 3);
    ASSERT_EQ(0, geometry_ff(cmd.argc(), cmd.argv()));
    std::ifstream in(jsonFile);
    const auto j     = nlohmann::json::parse(in);
    const auto &root = j.at("geometry_ff");
    EXPECT_EQ(1, root.at("compounds").get<int>());
    const auto &bonds = root.at("histograms").at("BONDS");
    ASSERT_EQ(1u, bonds.size());
    EXPECT_EQ("br_b~h_b", bonds[0].at("identifier").get<std::string>());
    const auto &sample = bonds[0].at("bins")[0].at("samples")[0];
    EXPECT_EQ("hydrogen-bromide", sample.at("molecule").get<std::string>());
    EXPECT_EQ("0000", sample.at("calculation").get<std::string>());
    // Bromine is the first atom of the file, hydrogen the second
    auto atoms = sample.at("atoms").get<std::vector<int>>();
    std::sort(atoms.begin(), atoms.end());
    EXPECT_EQ((std::vector<int>{ 1, 2 }), atoms);
    const auto names = sample.at("names").get<std::string>();
    EXPECT_TRUE(names == "H2-Br1" || names == "Br1-H2") << names;
}

//! \return the lines of a file
static std::vector<std::string> linesOf(const std::string &filename)
{
    std::vector<std::string> lines;
    std::ifstream            in(filename);
    std::string              line;
    while (std::getline(in, line))
    {
        lines.push_back(line);
    }
    return lines;
}

// Test: With -analyse the bond lengths and the angles of each type are
// tested for normality: the log file warns for the types that fail and
// ends with their number, and the json file has the result for each type.
// Without -analyse the log file and the json file have none of this.
TEST_F(GeometryFFTest, Analyse)
{
    const auto jsonFile = fileManager().getTemporaryFilePath("analyse.json");
    const auto logFile  = fileManager().getTemporaryFilePath("analyse.log");
    gmx::test::CommandLine cmd;
    cmd.append("geometry_ff");
    cmd.addOption("-mp", fileManager().getInputFilePath("molprop/molprop.xml"));
    cmd.addOption("-ff", fileManager().getInputFilePath("forcefields/ACS-g1.xml"));
    cmd.addOption("-o", fileManager().getTemporaryFilePath("analyse.xml"));
    cmd.addOption("-g", logFile);
    cmd.addOption("-json", jsonFile);
    const std::string notNormal("WARNING not normal");
    {
        gmx::test::CommandLine without(cmd);
        ASSERT_EQ(0, geometry_ff(without.argc(), without.argv()));
        std::ifstream in(jsonFile);
        EXPECT_FALSE(nlohmann::json::parse(in).at("geometry_ff").contains("analysis"));
        for (const auto &line : linesOf(logFile))
        {
            EXPECT_EQ(std::string::npos, line.find(notNormal)) << line;
            EXPECT_EQ(std::string::npos, line.find("Normality")) << line;
        }
    }
    // The number of types that are not normal in the json file and the
    // number of warnings in the log file, for a significance level
    const auto run = [&](const char *alpha, nlohmann::json *analysis)
    {
        gmx::test::CommandLine with(cmd);
        with.append("-analyse");
        if (alpha)
        {
            with.addOption("-alpha", alpha);
        }
        EXPECT_EQ(0, geometry_ff(with.argc(), with.argv()));
        std::ifstream in(jsonFile);
        *analysis = nlohmann::json::parse(in).at("geometry_ff").at("analysis");
        int warnings = 0, summaries = 0;
        for (const auto &line : linesOf(logFile))
        {
            if (line.find(notNormal) != std::string::npos)
            {
                // A line of a bond or an angle, with what tells how
                EXPECT_TRUE(line.rfind("bond-", 0) == 0 || line.rfind("angle-", 0) == 0) << line;
                EXPECT_NE(std::string::npos, line.find("skewness")) << line;
                EXPECT_NE(std::string::npos, line.find("kurtosis")) << line;
                warnings += 1;
            }
            if (line.rfind("Normality of the distributions", 0) == 0)
            {
                summaries += 1;
            }
        }
        EXPECT_EQ(1, summaries);
        return warnings;
    };
    nlohmann::json analysis;
    const int      warnings = run(nullptr, &analysis);
    EXPECT_EQ("Anderson-Darling", analysis.at("test").get<std::string>());
    EXPECT_DOUBLE_EQ(0.01, analysis.at("alpha").get<double>());
    EXPECT_EQ(8, analysis.at("minimum").get<int>());
    ASSERT_TRUE(analysis.contains("BONDS"));
    ASSERT_TRUE(analysis.contains("ANGLES"));
    EXPECT_EQ("pm", analysis.at("BONDS").at("unit").get<std::string>());
    EXPECT_EQ("degree", analysis.at("ANGLES").at("unit").get<std::string>());
    EXPECT_DOUBLE_EQ(8, analysis.at("BONDS").at("tolerance").get<double>());
    int  total = 0;
    bool foundCH = false;
    for (const auto &itype : { "BONDS", "ANGLES" })
    {
        const auto &jtype = analysis.at(itype);
        int tested = 0, failed = 0, wide = 0;
        for (const auto &entry : jtype.at("entries"))
        {
            const auto count = entry.at("count").get<int>();
            EXPECT_GE(entry.at("sigma").get<double>(), 0);
            wide += entry.at("wide").get<bool>() ? 1 : 0;
            if (!entry.at("tested").get<bool>())
            {
                // Too few values, or all the same: no verdict
                EXPECT_FALSE(entry.contains("normal")) << entry.dump();
                EXPECT_FALSE(entry.contains("pvalue")) << entry.dump();
                continue;
            }
            tested += 1;
            EXPECT_GE(count, 8) << entry.dump();
            const auto p = entry.at("pvalue").get<double>();
            EXPECT_GE(p, 0);
            EXPECT_LE(p, 1);
            EXPECT_GT(entry.at("statistic").get<double>(), 0);
            EXPECT_TRUE(entry.contains("skewness"));
            EXPECT_TRUE(entry.contains("kurtosis"));
            EXPECT_EQ(p >= 0.01, entry.at("normal").get<bool>()) << entry.dump();
            failed += entry.at("normal").get<bool>() ? 0 : 1;
            if (entry.at("identifier").get<std::string>() == "c3_b~h_b")
            {
                // Nine C-H bonds in each of two structures of butanol, six in acetone
                EXPECT_EQ(24, count);
                EXPECT_NEAR(109.4, entry.at("average").get<double>(), 0.05);
                EXPECT_FALSE(entry.at("wide").get<bool>());
                foundCH = true;
            }
        }
        EXPECT_EQ(static_cast<int>(jtype.at("entries").size()), jtype.at("distributions").get<int>());
        EXPECT_EQ(tested, jtype.at("tested").get<int>());
        EXPECT_EQ(failed, jtype.at("not_normal").get<int>());
        EXPECT_EQ(wide, jtype.at("wide").get<int>());
        EXPECT_GT(tested, 0) << itype;
        total += failed;
    }
    EXPECT_TRUE(foundCH);
    // Each type that is not normal has its warning in the log file
    EXPECT_EQ(total, warnings);
    // With a significance level close to one almost nothing passes, with
    // one close to zero everything does
    nlohmann::json strict, loose;
    const int      many = run("0.999", &strict);
    EXPECT_DOUBLE_EQ(0.999, strict.at("alpha").get<double>());
    EXPECT_GT(many, warnings);
    EXPECT_EQ(many, strict.at("BONDS").at("not_normal").get<int>() + strict.at("ANGLES").at("not_normal").get<int>());
    EXPECT_EQ(0, run("1e-35", &loose));
    EXPECT_EQ(0, loose.at("BONDS").at("not_normal").get<int>() + loose.at("ANGLES").at("not_normal").get<int>());
    EXPECT_EQ(strict.at("BONDS").at("tested").get<int>(), loose.at("BONDS").at("tested").get<int>());
}

} // namespace

} // namespace alexandria

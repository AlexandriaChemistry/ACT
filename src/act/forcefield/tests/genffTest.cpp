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
 * Tests of the gen_ff module, with the atom types and virtual sites of
 * the hydrogen halides.
 */

#include "actpre.h"

#include <cstddef>
#include <filesystem>
#include <fstream>
#include <gtest/gtest.h>
#include <map>
#include <nlohmann/json.hpp>
#include <sstream>
#include <string>
#include <utility>
#include <vector>

#include "act/alexandria/alex_modules.h"
#include "act/basics/chargemodel.h"
#include "act/basics/interactiontype.h"
#include "act/basics/mutability.h"
#include "act/forcefield/forcefield.h"
#include "act/forcefield/forcefield_parameter.h"
#include "act/forcefield/forcefield_validation.h"
#include "act/forcefield/forcefield_xml.h"
#include "act/forcefield/gen_ff.h"
#include "act/forcefield/potential.h"
#include "act/testutils/output_checks.h"
#include "gromacs/utility/exceptions.h"
#include "testutils/cmdlinetest.h"
#include "testutils/refdata.h"
#include "testutils/testfilemanager.h"

namespace alexandria
{

namespace
{

class GenFFTest : public gmx::test::CommandLineTestBase
{
protected:
    /*! \brief Run gen_ff
     * \param[in] options Command line options, an empty value for a flag
     * \param[in] vsites  File with virtual sites in this directory or with an
     *                    absolute path, none if empty
     * \param[in] atomtypes The file with atom types in the test data
     * \return the command line for gen_ff
     */
    gmx::test::CommandLine commandLine(const std::vector<std::pair<std::string, std::string>> &options,
                                       const std::string &vsites,
                                       const std::string &atomtypes = "forcefields/genff-atomtypes.csv")
    {
        gmx::test::CommandLine cmd;
        cmd.append("gen_ff");
        cmd.addOption("-f", fileManager().getInputFilePath(atomtypes));
        if (!vsites.empty())
        {
            cmd.addOption("-vs", vsites[0] == '/' ? vsites : fileManager().getInputFilePath(vsites));
        }
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
        return cmd;
    }

    //! \return the force field generated with the options and virtual sites, see commandLine
    ForceField generate(const std::vector<std::pair<std::string, std::string>> &options,
                        const std::string &vsites = "forcefields/genff-vsites.json")
    {
        const auto output = fileManager().getTemporaryFilePath("generated.xml");
        auto cmd    = commandLine(options, vsites);
        cmd.addOption("-o", output);
        EXPECT_EQ(0, gen_ff(cmd.argc(), cmd.argv()));
        ForceField pd;
        EXPECT_TRUE(std::filesystem::exists(output));
        readForceField(output, &pd);
        return pd;
    }

    //! \return the value of an option of an interaction type
    static std::string option(const ForceField &pd, InteractionType itype, const std::string &key)
    {
        auto &fs = pd.findForcesConst(itype);
        return fs.optionExists(key) ? fs.optionValue(key) : "";
    }
};

// Test: The particle types, their charges and options, the default
// functions and the virtual sites
TEST_F(GenFFTest, Defaults)
{
    auto pd = generate({});
    EXPECT_EQ(55u, pd.particleTypesConst().size());
    ASSERT_TRUE(pd.hasParticleType("f"));
    auto f = pd.findParticleType("f");
    EXPECT_EQ("f_z", f->optionValue("acmtype"));
    EXPECT_EQ("f_b", f->optionValue("bondtype"));
    EXPECT_EQ("9", f->optionValue("atomnumber"));
    EXPECT_EQ("f", f->optionValue("vdwtype"));
    EXPECT_EQ(Mutability::ACM, f->parameterConst("charge").mutability());
    // The charge starts in the middle of its range
    const auto qhf = pd.findParticleType("hf")->parameterConst("charge");
    EXPECT_DOUBLE_EQ(0, qhf.minimum());
    EXPECT_DOUBLE_EQ(0.25, qhf.value());
    EXPECT_DOUBLE_EQ(0.5, qhf.maximum());
    std::vector<std::pair<InteractionType, Potential>> functions = {
        { InteractionType::BONDS,            Potential::HARMONIC_BONDS },
        { InteractionType::ANGLES,           Potential::HARMONIC_ANGLES },
        { InteractionType::PROPER_DIHEDRALS, Potential::FOURIER_DIHEDRALS },
        { InteractionType::VDW,              Potential::LJ12_6 },
        { InteractionType::VSITE2,           Potential::VSITE2 }
    };
    for (const auto &fn : functions)
    {
        ASSERT_TRUE(pd.interactionPresent(fn.first)) << interactionTypeToString(fn.first);
        EXPECT_EQ(fn.second, pd.findForcesConst(fn.first).potential()) << interactionTypeToString(fn.first);
    }
    // The charge distribution, the relative dielectric constant and the
    // number of exclusions are options of the particle types
    EXPECT_FALSE(pd.interactionPresent(InteractionType::ELECTROSTATICS));
    EXPECT_EQ(ChargeDistributionType::Gaussian, pd.chargeDistribution());
    EXPECT_EQ(Potential::COULOMB_GAUSSIAN, pd.potential(InteractionType::ELECTROSTATICS));
    EXPECT_DOUBLE_EQ(1, pd.epsilonr());
    EXPECT_EQ(2, pd.nexcl());
    EXPECT_FALSE(f->hasOption("zetatype"));
    EXPECT_TRUE(f->hasParameter("zeta"));
    // The parameter of a virtual site starts in the middle of its range
    const ForceFieldParameterMap *vsite = nullptr;
    for (const auto &p : pd.findForcesConst(InteractionType::VSITE2).parametersConst())
    {
        if (p.first.id() == "br_b~hbr_b!v2br1")
        {
            vsite = &p.second;
        }
    }
    ASSERT_NE(nullptr, vsite);
    auto &vs2a = vsite->find("vs2a")->second;
    EXPECT_DOUBLE_EQ(-1, vs2a.minimum());
    EXPECT_DOUBLE_EQ(-0.25, vs2a.value());
    EXPECT_DOUBLE_EQ(0.5, vs2a.maximum());
}

// Test: Other functions, charge distribution and options
TEST_F(GenFFTest, Options)
{
    const auto pd = generate({ { "-qdist", "POINT" }, { "-bondfn", "MORSE_BONDS" },
                               { "-anglefn", "UREY_BRADLEY_ANGLES" }, { "-vdwfn", "WANG_BUCKINGHAM" },
                               { "-epsilonr", "1.7" }, { "-nexcl", "3" } });
    EXPECT_EQ(ChargeDistributionType::Point, pd.chargeDistribution());
    EXPECT_EQ(Potential::MORSE_BONDS, pd.findForcesConst(InteractionType::BONDS).potential());
    EXPECT_EQ(Potential::UREY_BRADLEY_ANGLES, pd.findForcesConst(InteractionType::ANGLES).potential());
    EXPECT_EQ(Potential::WANG_BUCKINGHAM, pd.findForcesConst(InteractionType::VDW).potential());
    EXPECT_DOUBLE_EQ(1.7, pd.epsilonr());
    EXPECT_EQ(3, pd.nexcl());
    // Point charges have no width
    for (const auto &pt : pd.particleTypesConst())
    {
        EXPECT_EQ(0, pt.second.zeta()) << pt.first.id();
    }
}

// Test: The bond-angle cross term can be chosen for the angles
TEST_F(GenFFTest, CrossBondAngles)
{
    const auto pd = generate({ { "-anglefn", "CROSS_BOND_ANGLES" } });
    EXPECT_EQ(Potential::CROSS_BOND_ANGLES, pd.findForcesConst(InteractionType::ANGLES).potential());
}

// Test: The quartic potentials can be chosen for the bonds and the angles
TEST_F(GenFFTest, QuarticBondsAndAngles)
{
    const auto pd = generate({ { "-bondfn", "QUARTIC_BONDS" }, { "-anglefn", "QUARTIC_ANGLES" } });
    EXPECT_EQ(Potential::QUARTIC_BONDS, pd.findForcesConst(InteractionType::BONDS).potential());
    EXPECT_EQ(Potential::QUARTIC_ANGLES, pd.findForcesConst(InteractionType::ANGLES).potential());
}

// The force fields of the tool for several sets of options, as reference
// data, such that a change in how the options reach the work shows

class GenFFResultsTest : public GenFFTest
{
protected:
    /*! \brief Run gen_ff on the atom types with shells and compare the force field with the reference data
     * \param[in] options Command line options, an empty value for a flag
     */
    void checkGenerated(const std::vector<std::pair<std::string, std::string>> &options)
    {
        const auto output = fileManager().getTemporaryFilePath("generated.xml");
        auto       cmd    = commandLine(options, "", "forcefields/genff-atomtypes-shells.csv");
        cmd.addOption("-o", output);
        ASSERT_EQ(0, gen_ff(cmd.argc(), cmd.argv()));
        gmx::test::TestReferenceChecker checker(this->rootChecker());
        checker.setDefaultTolerance(buildTolerance(1e-6, 1e-8, 1e-6));
        // The version has the time when the file was written
        checkTextFile(&checker, output, { "<version" }, "forcefield");
    }
};

TEST_F(GenFFResultsTest, Defaults)
{
    checkGenerated({});
}

TEST_F(GenFFResultsTest, Options)
{
    checkGenerated({ { "-nexcl", "3" }, { "-epsilonr", "1.7" }, { "-shellmaxdistance", "0.03" },
                     { "-noqsymm", "" }, { "-qdist", "SLATER" }, { "-bondfn", "MORSE_BONDS" },
                     { "-anglefn", "UREY_BRADLEY_ANGLES" }, { "-dihfn", "PROPER_DIHEDRALS" },
                     { "-vdwfn", "WANG_BUCKINGHAM" }, { "-vdwcorr", "" }, { "-icfn", "MORSE_BONDS" },
                     { "-cr", "VANDERWAALS:sigma:Arithmetic" } });
}

TEST_F(GenFFResultsTest, VirtualSites)
{
    // What the virtual sites add to a force field
    const auto plain = generate({}, "");
    const auto sites = generate({}, "forcefields/genff-vsites.json");
    gmx::test::TestReferenceChecker checker(this->rootChecker());
    checker.setDefaultTolerance(buildTolerance(1e-6, 1e-8, 1e-6));
    checkJson(&checker, forceFieldChanges(plain, sites), "changes");
}

//! \return the contents of a force field file without its version, which has the time
static std::string forceFieldText(const std::string &fileName)
{
    std::ifstream in(fileName);
    std::string   text;
    for (std::string line; std::getline(in, line); )
    {
        if (line.find("<version") == std::string::npos)
        {
            text += line + "\n";
        }
    }
    return text;
}

// The force field without a command line, as a program that uses ACT makes it

TEST_F(GenFFTest, WithoutCommandLineSameAsTheProgram)
{
    // The defaults with virtual sites, and other options
    {
        const auto tool = fileManager().getTemporaryFilePath("tool-defaults.xml");
        auto       cmd  = commandLine({}, "forcefields/genff-vsites.json");
        cmd.addOption("-o", tool);
        ASSERT_EQ(0, gen_ff(cmd.argc(), cmd.argv()));
        GenFFOptions options;
        options.atomTypesFile = fileManager().getInputFilePath("forcefields/genff-atomtypes.csv");
        options.virtualSiteFiles.push_back(fileManager().getInputFilePath("forcefields/genff-vsites.json"));
        options.outputFile = fileManager().getTemporaryFilePath("direct-defaults.xml");
        runGenFF(options);
        EXPECT_EQ(forceFieldText(tool), forceFieldText(options.outputFile));
    }
    {
        const auto tool = fileManager().getTemporaryFilePath("tool-options.xml");
        auto       cmd  = commandLine({ { "-nexcl", "3" }, { "-epsilonr", "1.7" }, { "-shellmaxdistance", "0.03" },
                                        { "-noqsymm", "" }, { "-qdist", "SLATER" }, { "-bondfn", "MORSE_BONDS" },
                                        { "-anglefn", "UREY_BRADLEY_ANGLES" }, { "-dihfn", "PROPER_DIHEDRALS" },
                                        { "-vdwfn", "WANG_BUCKINGHAM" }, { "-vdwcorr", "" }, { "-icfn", "MORSE_BONDS" },
                                        { "-cr", "VANDERWAALS:sigma:Arithmetic" } },
                                      "", "forcefields/genff-atomtypes-shells.csv");
        cmd.addOption("-o", tool);
        ASSERT_EQ(0, gen_ff(cmd.argc(), cmd.argv()));
        GenFFOptions options;
        options.atomTypesFile       = fileManager().getInputFilePath("forcefields/genff-atomtypes-shells.csv");
        options.outputFile          = fileManager().getTemporaryFilePath("direct-options.xml");
        options.nexcl               = 3;
        options.epsilonr            = static_cast<real>(1.7);
        options.shellMaxDistance    = static_cast<real>(0.03);
        options.symmetricCharges    = false;
        options.chargeDistribution  = ChargeDistributionType::Slater;
        options.bondPotential       = Potential::MORSE_BONDS;
        options.anglePotential      = Potential::UREY_BRADLEY_ANGLES;
        options.dihedralPotential   = Potential::PROPER_DIHEDRALS;
        options.vdwPotential        = Potential::WANG_BUCKINGHAM;
        options.vdwCorrection       = true;
        options.inductionCorrection = Potential::MORSE_BONDS;
        options.combinationRules    = "VANDERWAALS:sigma:Arithmetic";
        runGenFF(options);
        EXPECT_EQ(forceFieldText(tool), forceFieldText(options.outputFile));
    }
}

TEST_F(GenFFTest, WithoutCommandLineWrongOptions)
{
    const auto valid = [this]()
    {
        GenFFOptions options;
        options.atomTypesFile = fileManager().getInputFilePath("forcefields/genff-atomtypes-shells.csv");
        options.outputFile    = fileManager().getTemporaryFilePath("wrong.xml");
        return options;
    };
    std::vector<std::pair<GenFFOptions, std::string>> wrong;
    {
        auto o = valid();
        o.atomTypesFile.clear();
        wrong.emplace_back(o, "-f");
    }
    {
        auto o = valid();
        o.outputFile.clear();
        wrong.emplace_back(o, "-o");
    }
    {
        auto o = valid();
        o.nexcl = -1;
        wrong.emplace_back(o, "-nexcl");
    }
    {
        auto o = valid();
        o.epsilonr = 0;
        wrong.emplace_back(o, "-epsilonr");
    }
    {
        auto o = valid();
        o.bondPotential = Potential::LJ12_6;
        wrong.emplace_back(o, "-bondfn");
    }
    {
        auto o = valid();
        o.vdwPotential = Potential::HARMONIC_BONDS;
        wrong.emplace_back(o, "-vdwfn");
    }
    for (const auto &w : wrong)
    {
        try
        {
            runGenFF(w.first);
            ADD_FAILURE() << "Accepted options with a wrong " << w.second;
        }
        catch (const gmx::InvalidInputError &ex)
        {
            EXPECT_NE(std::string::npos, std::string(ex.what()).find(w.second)) << ex.what();
        }
        EXPECT_FALSE(std::filesystem::exists(w.first.outputFile)) << w.second;
    }
}

// Test: The obsolete options for the exclusions are accepted if they are
// the same, otherwise the force field is not generated
TEST_F(GenFFTest, ObsoleteExclusions)
{
    const auto pd = generate({ { "-nexclqq", "3" }, { "-nexclvdw", "3" } });
    EXPECT_EQ(3, pd.nexcl());
    EXPECT_THROW(generate({ { "-nexclqq", "3" }, { "-nexclvdw", "1" } }), gmx::InvalidInputError);
}

// Test: Without the file of virtual sites there are none
TEST_F(GenFFTest, WithoutVirtualSites)
{
    const auto pd = generate({}, "");
    EXPECT_FALSE(pd.interactionPresent(InteractionType::VSITE2));
    EXPECT_EQ(55u, pd.particleTypesConst().size());
}


//! \return the JSON of a virtual site of type v2f1
static std::string vsite(const std::string &itype, const std::string &atoms, const std::string &params)
{
    return "{ \"type\": \"v2f1\", \"interaction\": \"" + itype + "\", \"atoms\": [" + atoms +
        "], \"parameters\": { " + params + " } }";
}

//! \return a file of virtual sites with the JSON given
static std::string vsiteFile(gmx::test::TestFileManager *fileManager, const std::string &name,
                             const std::string &sites)
{
    auto vsfile = fileManager->getTemporaryFilePath(name);
    std::ofstream out(vsfile);
    out << "{ \"description\": \"Test\", \"virtual_sites\": [ " << sites << " ] }\n";
    return vsfile;
}

// Test: All types of virtual sites with parameters, which start in the
// middle of their range, from two files
TEST_F(GenFFTest, VirtualSiteTypes)
{
    const std::string f3  = "\"f_b\", \"c3_b\", \"hf_b\"";
    const std::string f4  = "\"f_b\", \"c3_b\", \"hf_b\", \"h1_b\"";
    const std::string a   = "{ \"min\": 0, \"max\": 0.2 }";
    const std::string b   = "{ \"min\": 0.2, \"max\": 0.6 }";
    const std::string c   = "{ \"min\": -1, \"max\": 1 }";
    const auto first = vsiteFile(&fileManager(), "vsites1.json",
                                 vsite("VSITE2FD", "\"f_b\", \"hf_b\"", "\"vs2fd_a\": { \"min\": 0.1, \"max\": 0.3 }") + ", " +
                                 vsite("VSITE3", f3, "\"vs3a\": " + a + ", \"vs3b\": " + b) + ", " +
                                 vsite("VSITE3S", f3, "\"vs3sa\": " + a) + ", " +
                                 vsite("VSITE3FD", f3, "\"vs3fd_a\": " + a + ", \"vs3fd_b\": " + b) + ", " +
                                 vsite("VSITE3FAD", f3, "\"vs3fad_a\": " + a + ", \"vs3fad_b\": " + b));
    const auto second = vsiteFile(&fileManager(), "vsites2.json",
                                  vsite("VSITE3OUT", f3, "\"vs3out_a\": " + a + ", \"vs3out_b\": " + b + ", \"vs3out_c\": " + c) + ", " +
                                  vsite("VSITE3OUTS", f3, "\"vs3outs_a\": " + a + ", \"vs3outs_c\": " + c) + ", " +
                                  vsite("VSITE4", f4, "\"vs4a\": " + a + ", \"vs4b\": " + b + ", \"vs4c\": " + c) + ", " +
                                  vsite("VSITE4S", f4, "\"vs4sa\": " + a + ", \"vs4sb\": " + b) + ", " +
                                  vsite("VSITE4S3", f4, "\"vs4s3a\": " + a));
    auto cmd    = commandLine({}, first);
    const auto output = fileManager().getTemporaryFilePath("generated.xml");
    cmd.append(second);
    cmd.addOption("-o", output);
    ASSERT_EQ(0, gen_ff(cmd.argc(), cmd.argv()));
    ForceField pd;
    readForceField(output, &pd);
    // The middle of the ranges 0-0.2, 0.2-0.6 and -1-1
    std::vector<double> middle = { 0.1, 0.4, 0 };
    std::map<InteractionType, size_t> nparam = {
        { InteractionType::VSITE3,     2 },
        { InteractionType::VSITE3S,    1 },
        { InteractionType::VSITE3FD,   2 },
        { InteractionType::VSITE3FAD,  2 },
        { InteractionType::VSITE3OUT,  3 },
        { InteractionType::VSITE3OUTS, 2 },
        { InteractionType::VSITE4,     3 },
        { InteractionType::VSITE4S,    2 },
        { InteractionType::VSITE4S3,   1 }
    };
    for (const auto &np : nparam)
    {
        const auto name = interactionTypeToString(np.first);
        ASSERT_TRUE(pd.interactionPresent(np.first)) << name;
        auto &fs = pd.findForcesConst(np.first);
        ASSERT_EQ(1u, fs.parametersConst().size()) << name;
        auto &params = fs.parametersConst().begin()->second;
        auto  pnames = potentialToParameterName(fs.potential());
        ASSERT_EQ(np.second, pnames.size()) << name;
        for (size_t i = 0; i < pnames.size(); i++)
        {
            // The second parameter of VSITE3OUTS is the third of VSITE3OUT
            double expect = (np.first == InteractionType::VSITE3OUTS && i == 1) ? middle[2] : middle[i];
            auto &p = params.find(pnames[i])->second;
            EXPECT_NEAR(expect, p.value(), 1e-12) << name << " " << pnames[i];
            EXPECT_LE(p.minimum(), p.value()) << name << " " << pnames[i];
            EXPECT_GE(p.maximum(), p.value()) << name << " " << pnames[i];
        }
    }
    auto &vs2fd = pd.findForcesConst(InteractionType::VSITE2FD).parametersConst().begin()->second;
    EXPECT_NEAR(0.2, vs2fd.find("vs2fd_a")->second.value(), 1e-12);
}

// Test: Invalid virtual sites are an error: too few atoms, a parameter
// missing or unknown, a minimum larger than the maximum, an unknown
// interaction type, and a file that is not JSON
TEST_F(GenFFTest, VirtualSiteErrors)
{
    const std::string f2 = "\"f_b\", \"hf_b\"";
    const std::string ok = "{ \"min\": 0, \"max\": 0.2 }";
    for (const std::string &site : { vsite("VSITE3", f2, "\"vs3a\": " + ok + ", \"vs3b\": " + ok),
                                     vsite("VSITE2", f2, ""),
                                     vsite("VSITE2", f2, "\"vs2a\": " + ok + ", \"vs3a\": " + ok),
                                     vsite("VSITE2", f2, "\"vs2a\": { \"min\": 1, \"max\": 0 }"),
                                     vsite("VSITE2", f2, "\"vs2a\": { \"min\": 0 }"),
                                     vsite("BONDS", f2, "\"vs2a\": " + ok),
                                     std::string("{ \"interaction\": \"VSITE2\" }") })
    {
        const auto vsfile = vsiteFile(&fileManager(), "vsites.json", site);
        auto cmd = commandLine({ { "-o", fileManager().getTemporaryFilePath("generated.xml") } }, vsfile);
        EXPECT_THROW(gen_ff(cmd.argc(), cmd.argv()), gmx::InvalidInputError) << site;
    }
    const auto vsfile = fileManager().getTemporaryFilePath("notjson.json");
    {
        std::ofstream out(vsfile);
        out << "v2f1|VSITE2|f_b~hf_b|-1|0.5|||\n";
    }
    auto cmd = commandLine({ { "-o", fileManager().getTemporaryFilePath("generated.xml") } }, vsfile);
    EXPECT_THROW(gen_ff(cmd.argc(), cmd.argv()), gmx::InvalidInputError);
}

// Test: A force field with shells is described as polarizable in the JSON
// output, as it is when the file is read, and one without is not
TEST_F(GenFFTest, JsonPolarizable)
{
    for (const auto &input : std::vector<std::pair<std::string, bool>>{
             { "forcefields/genff-atomtypes-shells.csv", true },
             { "forcefields/genff-atomtypes.csv", false } })
    {
        SCOPED_TRACE(input.first);
        const auto output   = fileManager().getTemporaryFilePath("polarizable.xml");
        const auto jsonFile = fileManager().getTemporaryFilePath("polarizable.json");
        auto cmd = commandLine({ { "-o", output }, { "-json", jsonFile } }, "", input.first);
        ASSERT_EQ(0, gen_ff(cmd.argc(), cmd.argv()));
        ForceField pd;
        readForceField(output, &pd);
        EXPECT_EQ(input.second, pd.polarizable());
        std::ifstream in(jsonFile);
        ASSERT_TRUE(in.good());
        const auto ff = nlohmann::json::parse(in).at("gen_ff").at("forcefield");
        EXPECT_EQ(input.second, ff.at("polarizable").get<bool>());
        if (input.second)
        {
            // An atom with a shell is a point charge, the shell is not
            ASSERT_TRUE(pd.hasParticleType("c3_s"));
            EXPECT_EQ("c3_s", pd.findParticleType("c3")->optionValue("poltype"));
        }
    }
}

/*! \brief Copy the atom types of the tests with other bounds for chi of fluorine
 * \param[in] fileManager The test files
 * \param[in] name        Name of the copy
 * \param[in] min         The minimum of chi
 * \param[in] max         The maximum of chi
 * \param[in] nonneg      Whether chi must be non-negative
 * \return the name of the copy with its path
 */
static std::string chiFile(gmx::test::TestFileManager *fileManager, const std::string &name,
                           const std::string &min, const std::string &max, bool nonneg)
{
    std::ifstream in(fileManager->getInputFilePath("forcefields/genff-atomtypes.csv"));
    auto          copy = fileManager->getTemporaryFilePath(name);
    std::ofstream out(copy);
    const std::string prefix = "ELECTRONEGATIVITYEQUALIZATION:chi:";
    std::map<std::string, size_t> column;
    std::string line;
    while (std::getline(in, line))
    {
        std::vector<std::string> words;
        std::istringstream       stream(line);
        std::string              word;
        while (std::getline(stream, word, '|'))
        {
            words.push_back(word);
        }
        if (column.empty())
        {
            for (size_t i = 0; i < words.size(); i++)
            {
                column[words[i]] = i;
            }
        }
        else if (words.size() > column[prefix + "nonneg"] && words[column["acmtype"]] == "f_z")
        {
            words[column[prefix + "min"]]    = min;
            words[column[prefix + "max"]]    = max;
            words[column[prefix + "nonneg"]] = nonneg ? "True" : "False";
            line.clear();
            for (size_t i = 0; i < words.size(); i++)
            {
                line += (i > 0 ? "|" : "") + words[i];
            }
        }
        out << line << "\n";
    }
    return copy;
}

// Test: A parameter starts at its minimum, also one that may be negative,
// and a new force field has no parameters that were corrected. A minimum
// that is larger than the maximum is an error.
TEST_F(GenFFTest, ParametersStartAtMinimum)
{
    const auto output = fileManager().getTemporaryFilePath("generated.xml");
    {
        gmx::test::CommandLine cmd;
        cmd.append("gen_ff");
        cmd.addOption("-f", chiFile(&fileManager(), "negative.csv", "-3", "5", false));
        cmd.addOption("-o", output);
        ASSERT_EQ(0, gen_ff(cmd.argc(), cmd.argv()));
        ForceField pd;
        readForceField(output, &pd);
        bool found = false;
        for (const auto &p : pd.findForcesConst(InteractionType::ELECTRONEGATIVITYEQUALIZATION).parametersConst())
        {
            if (p.first.id() == "f_z")
            {
                const auto &chi = p.second.find("chi")->second;
                EXPECT_DOUBLE_EQ(-3, chi.value());
                EXPECT_DOUBLE_EQ(-3, chi.minimum());
                EXPECT_DOUBLE_EQ(5, chi.maximum());
                EXPECT_FALSE(chi.nonNegative());
                found = true;
            }
            else if (p.first.id() == "hf_z")
            {
                // As it is in the file of the tests
                EXPECT_DOUBLE_EQ(4, p.second.find("chi")->second.value());
            }
        }
        EXPECT_TRUE(found);
        EXPECT_TRUE(validateForceField(pd).empty());
    }
    {
        gmx::test::CommandLine cmd;
        cmd.append("gen_ff");
        cmd.addOption("-f", chiFile(&fileManager(), "wrong.csv", "6", "5", true));
        cmd.addOption("-o", output);
        EXPECT_THROW(gen_ff(cmd.argc(), cmd.argv()), gmx::InvalidInputError);
    }
}

// Test: With -json the force field that was written is described in JSON,
// with the result of the validation, and not without that option
TEST_F(GenFFTest, Json)
{
    const auto jsonFile = fileManager().getTemporaryFilePath("generated.json");
    std::filesystem::remove(jsonFile);
    std::filesystem::remove("gen_ff.json");
    const auto pd = generate({ { "-json", jsonFile } });
    EXPECT_FALSE(std::filesystem::exists("gen_ff.json"));
    std::ifstream in(jsonFile);
    ASSERT_TRUE(in.good());
    const auto j = nlohmann::json::parse(in);
    ASSERT_TRUE(j.contains("gen_ff"));
    const auto &ff = j.at("gen_ff").at("forcefield");
    // It is the force field of the file
    EXPECT_EQ(pd.checkSum(), ff.at("checksum").get<std::string>());
    EXPECT_EQ(pd.polarizable(), ff.at("polarizable").get<bool>());
    const auto &ptypes = ff.at("particletypes").at("types");
    ASSERT_EQ(pd.particleTypesConst().size(), ptypes.size());
    bool foundHF = false;
    for (const auto &jpt : ptypes)
    {
        if (jpt.at("identifier").get<std::string>() != "hf")
        {
            continue;
        }
        foundHF = true;
        for (const auto &jp : jpt.at("parameters"))
        {
            if (jp.at("type").get<std::string>() == "charge")
            {
                EXPECT_DOUBLE_EQ(0.25, jp.at("value").get<double>());
                EXPECT_DOUBLE_EQ(0.5, jp.at("maximum").get<double>());
            }
        }
    }
    EXPECT_TRUE(foundHF);
    // The interactions, with the virtual sites of the file
    std::map<std::string, std::string> functions;
    for (const auto &jfs : ff.at("interactions"))
    {
        functions[jfs.at("type").get<std::string>()] = jfs.at("function").get<std::string>();
    }
    for (const auto itype : { InteractionType::BONDS, InteractionType::VDW, InteractionType::VSITE2 })
    {
        const auto name = interactionTypeToString(itype);
        ASSERT_EQ(1u, functions.count(name)) << name;
        EXPECT_EQ(potentialToString(pd.findForcesConst(itype).potential()), functions[name]);
    }
    // A new force field has no problems
    const auto &validation = j.at("gen_ff").at("validation");
    EXPECT_TRUE(validation.at("ok").get<bool>());
    EXPECT_TRUE(validation.at("issues").empty());
    // No file without the option
    std::filesystem::remove(jsonFile);
    generate({});
    EXPECT_FALSE(std::filesystem::exists(jsonFile));
    EXPECT_FALSE(std::filesystem::exists("gen_ff.json"));
}

} // namespace

} // namespace alexandria

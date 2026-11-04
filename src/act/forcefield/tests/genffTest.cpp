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

#include <filesystem>
#include <fstream>
#include <string>
#include <vector>

#include <gtest/gtest.h>

#include "act/alexandria/alex_modules.h"
#include "act/basics/interactiontype.h"
#include "act/forcefield/forcefield.h"
#include "act/forcefield/forcefield_xml.h"
#include "act/forcefield/forcefield_parametername.h"
#include "act/forcefield/potential.h"
#include "gromacs/utility/exceptions.h"

#include "testutils/cmdlinetest.h"
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
     * \return the command line for gen_ff
     */
    gmx::test::CommandLine commandLine(const std::vector<std::pair<std::string, std::string>> &options,
                                       const std::string &vsites)
    {
        gmx::test::CommandLine cmd;
        cmd.append("gen_ff");
        cmd.addOption("-f", fileManager().getInputFilePath("genff-atomtypes.csv"));
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
                        const std::string &vsites = "genff-vsites.csv")
    {
        auto output = fileManager().getTemporaryFilePath("generated.xml");
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
    auto qhf = pd.findParticleType("hf")->parameterConst("charge");
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
    auto pd = generate({ { "-qdist", "POINT" }, { "-bondfn", "MORSE_BONDS" },
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

// Test: The obsolete options for the exclusions are accepted if they are
// the same, otherwise the force field is not generated
TEST_F(GenFFTest, ObsoleteExclusions)
{
    auto pd = generate({ { "-nexclqq", "3" }, { "-nexclvdw", "3" } });
    EXPECT_EQ(3, pd.nexcl());
    EXPECT_THROW(generate({ { "-nexclqq", "3" }, { "-nexclvdw", "1" } }), gmx::InvalidInputError);
}

// Test: Without the file of virtual sites there are none
TEST_F(GenFFTest, WithoutVirtualSites)
{
    auto pd = generate({}, "");
    EXPECT_FALSE(pd.interactionPresent(InteractionType::VSITE2));
    EXPECT_EQ(55u, pd.particleTypesConst().size());
}


// Test: All types of virtual sites with parameters, which start in the
// middle of their range, and comments and empty lines in the file
TEST_F(GenFFTest, VirtualSiteTypes)
{
    auto vsfile = fileManager().getTemporaryFilePath("vsites.csv");
    {
        std::ofstream out(vsfile);
        out << "# Virtual sites\n"
            << "v2f1|VSITE2FD|f_b~hf_b|0.1|0.3|||\n"
            << "\n"
            << "v2f1|VSITE3|f_b~c3_b~hf_b|0|0.2|0.2|0.6|\n"
            << "v2f1|VSITE3S|f_b~c3_b~hf_b|0|0.2|0|0|||\n"
            << "v2f1|VSITE3FD|f_b~c3_b~hf_b|0|0.2|0.2|0.6\n"
            << "v2f1|VSITE3FAD|f_b~c3_b~hf_b|0|0.2|0.2|0.6\n"
            << "v2f1|VSITE3OUT|f_b~c3_b~hf_b|0|0.2|0.2|0.6|-1|1\n"
            << "v2f1|VSITE3OUTS|f_b~c3_b~hf_b|0|0.2|-1|1|1\n"
            << "v2f1|VSITE4|f_b~c3_b~hf_b~h1_b|0|0.2|0.2|0.6|-1|1\n"
            << "v2f1|VSITE4S|f_b~c3_b~hf_b~h1_b|0|0.2|0.2|0.6|||\n"
            << "v2f1|VSITE4S3|f_b~c3_b~hf_b~h1_b|0|0.2|||\n";
    }
    auto pd = generate({}, vsfile);
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
        auto name = interactionTypeToString(np.first);
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

// Test: Lines in the file of virtual sites with too few items are an error
TEST_F(GenFFTest, VirtualSiteTooFewItems)
{
    for (const char *line : { "v2f1|VSITE3OUT|f_b~c3_b~hf_b|0|0.2|0.2|0.6", "v2f1|VSITE2" })
    {
        auto vsfile = fileManager().getTemporaryFilePath("vsites.csv");
        {
            std::ofstream out(vsfile);
            out << line << "\n";
        }
        auto cmd = commandLine({ { "-o", fileManager().getTemporaryFilePath("generated.xml") } }, vsfile);
        EXPECT_THROW(gen_ff(cmd.argc(), cmd.argv()), gmx::InvalidInputError) << line;
    }
}
} // namespace

} // namespace alexandria

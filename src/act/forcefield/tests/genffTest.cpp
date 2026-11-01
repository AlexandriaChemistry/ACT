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
#include <string>
#include <vector>

#include <gtest/gtest.h>

#include "act/alexandria/alex_modules.h"
#include "act/basics/interactiontype.h"
#include "act/forcefield/forcefield.h"
#include "act/forcefield/forcefield_xml.h"
#include "act/forcefield/potential.h"

#include "testutils/cmdlinetest.h"
#include "testutils/testfilemanager.h"

namespace alexandria
{

namespace
{

class GenFFTest : public gmx::test::CommandLineTestBase
{
protected:
    //! \return the force field generated with the options, an empty value for a flag
    ForceField generate(const std::vector<std::pair<std::string, std::string>> &options,
                        bool vsites = true)
    {
        auto output = fileManager().getTemporaryFilePath("generated.xml");
        gmx::test::CommandLine cmd;
        cmd.append("gen_ff");
        cmd.addOption("-f", fileManager().getInputFilePath("genff-atomtypes.csv"));
        if (vsites)
        {
            cmd.addOption("-vs", fileManager().getInputFilePath("genff-vsites.csv"));
        }
        cmd.addOption("-o", output);
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
    auto &qhf = pd.findParticleType("hf")->parameterConst("charge");
    EXPECT_DOUBLE_EQ(0, qhf.minimum());
    EXPECT_DOUBLE_EQ(0.25, qhf.value());
    EXPECT_DOUBLE_EQ(0.5, qhf.maximum());
    std::vector<std::pair<InteractionType, Potential>> functions = {
        { InteractionType::BONDS,            Potential::HARMONIC_BONDS },
        { InteractionType::ANGLES,           Potential::HARMONIC_ANGLES },
        { InteractionType::PROPER_DIHEDRALS, Potential::FOURIER_DIHEDRALS },
        { InteractionType::VDW,              Potential::LJ12_6 },
        { InteractionType::ELECTROSTATICS,   Potential::COULOMB_GAUSSIAN },
        { InteractionType::VSITE2,           Potential::VSITE2 }
    };
    for (const auto &fn : functions)
    {
        ASSERT_TRUE(pd.interactionPresent(fn.first)) << interactionTypeToString(fn.first);
        EXPECT_EQ(fn.second, pd.findForcesConst(fn.first).potential()) << interactionTypeToString(fn.first);
    }
    EXPECT_DOUBLE_EQ(1, std::stod(option(pd, InteractionType::ELECTROSTATICS, "epsilonr")));
    EXPECT_EQ("2", option(pd, InteractionType::ELECTROSTATICS, "nexcl"));
    EXPECT_EQ("2", option(pd, InteractionType::VDW, "nexcl"));
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
                         { "-epsilonr", "1.7" }, { "-nexclqq", "3" }, { "-nexclvdw", "1" } });
    EXPECT_EQ(Potential::COULOMB_POINT, pd.findForcesConst(InteractionType::ELECTROSTATICS).potential());
    EXPECT_EQ(Potential::MORSE_BONDS, pd.findForcesConst(InteractionType::BONDS).potential());
    EXPECT_EQ(Potential::UREY_BRADLEY_ANGLES, pd.findForcesConst(InteractionType::ANGLES).potential());
    EXPECT_EQ(Potential::WANG_BUCKINGHAM, pd.findForcesConst(InteractionType::VDW).potential());
    EXPECT_DOUBLE_EQ(1.7, std::stod(option(pd, InteractionType::ELECTROSTATICS, "epsilonr")));
    EXPECT_EQ("3", option(pd, InteractionType::ELECTROSTATICS, "nexcl"));
    EXPECT_EQ("1", option(pd, InteractionType::VDW, "nexcl"));
}

// Test: Without the file of virtual sites there are none
TEST_F(GenFFTest, WithoutVirtualSites)
{
    auto pd = generate({}, false);
    EXPECT_FALSE(pd.interactionPresent(InteractionType::VSITE2));
    EXPECT_EQ(55u, pd.particleTypesConst().size());
}

} // namespace

} // namespace alexandria

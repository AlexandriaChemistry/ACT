/*
 * This source file is part of the Alexandria Chemistry Toolkit.
 *
 * Copyright (C) 2026
 *
 * Developers:
 *             Mohammad Mehdi Ghahremanpour,
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
 * Tests of shells with a fixed force constant, of which the polarizability
 * or the charge follows from the force constant, see generateShellParameters.
 * \author David van der Spoel <david.vanderspoel@icm.uu.se>
 */
#include <algorithm>
#include <cmath>
#include <gtest/gtest.h>
#include <string>
#include <vector>

#include "act/alexandria/alex_modules.h"
#include "act/basics/identifier.h"
#include "act/basics/interactiontype.h"
#include "act/basics/msg_handler.h"
#include "act/basics/mutability.h"
#include "act/forcefield/edit_ff.h"
#include "act/forcefield/forcefield.h"
#include "act/forcefield/forcefield_parameter.h"
#include "act/forcefield/forcefield_xml.h"
#include "act/forcefield/gen_ff.h"
#include "act/forcefield/generate_dependent.h"
#include "act/forcefield/particletype.h"
#include "act/testutils/forcefield_utils.h"
#include "act/utility/units.h"
#include "gromacs/math/units.h"
#include "gromacs/utility/exceptions.h"
#include "testutils/cmdlinetest.h"
#include "testutils/testfilemanager.h"

namespace alexandria
{

namespace
{

//! The shell type of the tests, of the bromide ion in ACS-pg1
const char *shellType = "Br-_s";

class ShellForceConstantTest : public ::testing::Test
{
protected:
    //! A copy of the force field, such that the tests do not change the shared one
    ForceField pd_;

    void SetUp() override
    {
        pd_ = *getForceField("ACS-pg1");
    }

    //! \return the polarization parameter of the shell type
    ForceFieldParameter *polarization(ForceField *pd, const std::string &name)
    {
        auto &pmap = pd->findForces(InteractionType::POLARIZATION)->parametersNoRemoval()->find(Identifier(shellType))->second;
        return &pmap.find(name)->second;
    }

    //! \return the polarizability
    ForceFieldParameter *alpha(ForceField *pd = nullptr) { return polarization(pd ? pd : &pd_, "alpha"); }

    //! \return the force constant
    ForceFieldParameter *kshell(ForceField *pd = nullptr) { return polarization(pd ? pd : &pd_, kshellName); }

    //! \return the charge of the shell
    ForceFieldParameter *charge(ForceField *pd = nullptr)
    {
        return (pd ? pd : &pd_)->findParticleType(shellType)->parameter("charge");
    }

    //! Give a parameter a value and a mutability, with the value as bounds
    static void set(ForceFieldParameter *p, double value, Mutability m)
    {
        p->setMutability(Mutability::Bounded);
        p->setMaximum(value);
        (void) p->setMinimum(value);
        p->forceSetValue(value);
        p->setMutability(m);
    }
};

TEST_F(ShellForceConstantTest, DependentForceConstantChangesNothing)
{
    EXPECT_EQ(Mutability::Dependent, kshell()->mutability());
    double a      = alpha()->value();
    double q      = charge()->value();
    const auto changed  = generateShellParameters(&pd_);
    EXPECT_FALSE(changed.alpha);
    EXPECT_FALSE(changed.charge);
    EXPECT_EQ(a, alpha()->value());
    EXPECT_EQ(q, charge()->value());
    EXPECT_TRUE(fixedShellForceConstants(&pd_).empty());
}

TEST_F(ShellForceConstantTest, AlphaFollowsFromCharge)
{
    set(kshell(), 400000, Mutability::Fixed);
    alpha()->setMutability(Mutability::Dependent);
    double q     = charge()->value();
    auto changed = generateShellParameters(&pd_);
    EXPECT_TRUE(changed.alpha);
    EXPECT_FALSE(changed.charge);
    // In nm3, the units of the program
    double expect = ONE_4PI_EPS0*q*q/400000;
    EXPECT_NEAR(expect, alpha()->internalValue(), 1e-12*expect);
    EXPECT_EQ(alpha()->value(), alpha()->minimum());
    EXPECT_EQ(alpha()->value(), alpha()->maximum());
    EXPECT_EQ(Mutability::Dependent, alpha()->mutability());
    EXPECT_EQ(q, charge()->value());
    // The second time nothing changes
    changed = generateShellParameters(&pd_);
    EXPECT_FALSE(changed.alpha);
}

TEST_F(ShellForceConstantTest, ChargeFollowsFromAlpha)
{
    set(kshell(), 400000, Mutability::Fixed);
    charge()->setMutability(Mutability::Dependent);
    double a     = alpha()->internalValue();
    const auto changed = generateShellParameters(&pd_);
    EXPECT_FALSE(changed.alpha);
    EXPECT_TRUE(changed.charge);
    // The charge of the shell stays negative
    double expect = -std::sqrt(400000*a/ONE_4PI_EPS0);
    EXPECT_NEAR(expect, charge()->value(), 1e-12*std::abs(expect));
    EXPECT_FALSE(charge()->nonNegative());
    EXPECT_EQ(Mutability::Dependent, charge()->mutability());
    // As the force computation has it
    double q = charge()->value();
    EXPECT_NEAR(400000, ONE_4PI_EPS0*q*q/a, 1e-6);
}

TEST_F(ShellForceConstantTest, PositiveChargeStaysPositive)
{
    set(kshell(), 400000, Mutability::Fixed);
    set(charge(), 2, Mutability::Dependent);
    (void) generateShellParameters(&pd_);
    EXPECT_GT(charge()->value(), 0);
}

TEST_F(ShellForceConstantTest, TrainedParameterChangesDependentOne)
{
    set(kshell(), 400000, Mutability::Fixed);
    charge()->setMutability(Mutability::Dependent);
    (void) generateShellParameters(&pd_);
    double q = charge()->value();
    // As training changes alpha
    set(alpha(), 2*alpha()->value(), Mutability::Bounded);
    const auto changed = generateShellParameters(&pd_);
    EXPECT_TRUE(changed.charge);
    EXPECT_NEAR(std::sqrt(2.0)*q, charge()->value(), 1e-9*std::abs(q));
}

TEST_F(ShellForceConstantTest, NeitherDependentThrows)
{
    set(kshell(), 400000, Mutability::Fixed);
    EXPECT_THROW(generateShellParameters(&pd_), gmx::InvalidInputError);
}

TEST_F(ShellForceConstantTest, BothDependentThrows)
{
    set(kshell(), 400000, Mutability::Fixed);
    alpha()->setMutability(Mutability::Dependent);
    charge()->setMutability(Mutability::Dependent);
    EXPECT_THROW(generateShellParameters(&pd_), gmx::InvalidInputError);
}

TEST_F(ShellForceConstantTest, ZeroForceConstantThrows)
{
    set(kshell(), 0, Mutability::Fixed);
    alpha()->setMutability(Mutability::Dependent);
    EXPECT_THROW(generateShellParameters(&pd_), gmx::InvalidInputError);
}

TEST_F(ShellForceConstantTest, ListsFixedShells)
{
    set(kshell(), 400000, Mutability::Fixed);
    alpha()->setMutability(Mutability::Dependent);
    (void) generateShellParameters(&pd_);
    auto shells = fixedShellForceConstants(&pd_);
    ASSERT_EQ(1u, shells.size());
    EXPECT_EQ(shellType, shells[0].shellType);
    EXPECT_EQ(400000, shells[0].kshell);
    EXPECT_EQ(kshell()->unit(), shells[0].kshellUnit);
    EXPECT_EQ(alpha()->value(), shells[0].alpha);
    EXPECT_EQ(alpha()->unit(), shells[0].alphaUnit);
    EXPECT_EQ(charge()->value(), shells[0].charge);
    EXPECT_TRUE(shells[0].alphaDependent);
}

TEST_F(ShellForceConstantTest, WriteAndRead)
{
    set(kshell(), 400000, Mutability::Fixed);
    charge()->setMutability(Mutability::Dependent);
    (void) generateShellParameters(&pd_);
    gmx::test::TestFileManager files;
    const auto fn = files.getTemporaryFilePath("shell.xml");
    writeForceField(fn, &pd_, false);
    ForceField pd;
    readForceField(fn, &pd);
    EXPECT_EQ(Mutability::Fixed, kshell(&pd)->mutability());
    EXPECT_EQ(Mutability::Dependent, charge(&pd)->mutability());
    EXPECT_NEAR(charge()->value(), charge(&pd)->value(), 1e-5*std::abs(charge()->value()));
    EXPECT_EQ(1u, fixedShellForceConstants(&pd).size());
}

TEST_F(ShellForceConstantTest, ShellTypes)
{
    const auto types = shellTypes(&pd_);
    EXPECT_NE(types.end(), std::find(types.begin(), types.end(), Identifier(shellType)));
    // The atom has no shell type
    EXPECT_EQ(types.end(), std::find(types.begin(), types.end(), Identifier("Br-")));
}

TEST_F(ShellForceConstantTest, FixKeepsCurrentForceConstant)
{
    const double a = alpha()->internalValue();
    const double q = charge()->value();
    fixShellForceConstant(&pd_, Identifier(shellType), ShellFollows::Charge, 0);
    const double k = ONE_4PI_EPS0*q*q/a;
    EXPECT_NEAR(k, kshell()->value(), 1e-12*k);
    EXPECT_EQ(Mutability::Fixed, kshell()->mutability());
    EXPECT_EQ(Mutability::Dependent, charge()->mutability());
    EXPECT_EQ(Mutability::Fixed, alpha()->mutability());
    EXPECT_NEAR(q, charge()->value(), 1e-12*std::abs(q));
}

TEST_F(ShellForceConstantTest, FixWithValueAndSwitch)
{
    fixShellForceConstant(&pd_, Identifier(shellType), ShellFollows::Alpha, 400000);
    EXPECT_EQ(400000, kshell()->value());
    EXPECT_EQ(Mutability::Dependent, alpha()->mutability());
    const double q = charge()->value();
    EXPECT_NEAR(ONE_4PI_EPS0*q*q/400000, alpha()->internalValue(), 1e-12*alpha()->internalValue());
    // Now the charge follows, from the same force constant
    fixShellForceConstant(&pd_, Identifier(shellType), ShellFollows::Charge, 0);
    EXPECT_NEAR(400000, kshell()->value(), 1e-6);
    EXPECT_EQ(Mutability::Fixed, alpha()->mutability());
    EXPECT_EQ(Mutability::Dependent, charge()->mutability());
}

TEST_F(ShellForceConstantTest, Release)
{
    fixShellForceConstant(&pd_, Identifier(shellType), ShellFollows::Alpha, 400000);
    EXPECT_TRUE(releaseShellForceConstant(&pd_, Identifier(shellType)));
    EXPECT_EQ(Mutability::Dependent, kshell()->mutability());
    EXPECT_EQ(Mutability::Fixed, alpha()->mutability());
    EXPECT_NEAR(400000, kshell()->value(), 1e-6);
    EXPECT_TRUE(fixedShellForceConstants(&pd_).empty());
    EXPECT_FALSE(releaseShellForceConstant(&pd_, Identifier(shellType)));
}

TEST_F(ShellForceConstantTest, NotAShellThrows)
{
    EXPECT_THROW(fixShellForceConstant(&pd_, Identifier("Br-"), ShellFollows::Charge, 0), gmx::InvalidInputError);
    EXPECT_THROW(releaseShellForceConstant(&pd_, Identifier("Br-")), gmx::InvalidInputError);
}

TEST_F(ShellForceConstantTest, GenFFFixesShells)
{
    gmx::test::TestFileManager files;
    const auto out = files.getTemporaryFilePath("genff.xml");
    gmx::test::CommandLine cmd;
    cmd.append("gen_ff");
    cmd.addOption("-f", files.getInputFilePath("forcefields/genff-atomtypes-shells.csv"));
    cmd.addOption("-o", out);
    cmd.addOption("-fixkshell", "charge");
    cmd.addOption("-kshell", "300000");
    ASSERT_EQ(0, gen_ff(cmd.argc(), cmd.argv()));
    ForceField pd;
    readForceField(out, &pd);
    const auto types  = shellTypes(&pd);
    const auto shells = fixedShellForceConstants(&pd);
    ASSERT_FALSE(types.empty());
    ASSERT_EQ(types.size(), shells.size());
    for (const auto &s : shells)
    {
        EXPECT_NEAR(300000, s.kshell, 1e-6) << s.shellType;
        EXPECT_FALSE(s.alphaDependent) << s.shellType;
        EXPECT_EQ(Mutability::Dependent, pd.findParticleType(s.shellType)->parameterConst("charge").mutability());
    }
}

TEST_F(ShellForceConstantTest, GenFFForceConstantOfTheRanges)
{
    gmx::test::TestFileManager files;
    GenFFOptions options;
    options.atomTypesFile = files.getInputFilePath("forcefields/genff-atomtypes-shells.csv");
    options.outputFile    = files.getTemporaryFilePath("genff-ranges.xml");
    options.shellFollows  = ShellFollows::Alpha;
    runGenFF(options);
    ForceField pd;
    readForceField(options.outputFile, &pd);
    for (const auto &s : fixedShellForceConstants(&pd))
    {
        // The charge does not change, alpha starts where it gives the
        // force constant, in the middle of its range
        const auto  &a = pd.findForcesConst(InteractionType::POLARIZATION).parametersConst().find(Identifier(s.shellType))->second.find("alpha")->second;
        const double q = s.charge;
        EXPECT_TRUE(s.alphaDependent) << s.shellType;
        EXPECT_NEAR(s.kshell, ONE_4PI_EPS0*q*q/a.internalValue(), 1e-4*s.kshell) << s.shellType;
    }
    EXPECT_FALSE(fixedShellForceConstants(&pd).empty());
}

TEST_F(ShellForceConstantTest, GenFFKshellNeedsFix)
{
    GenFFOptions options;
    options.atomTypesFile = "x.csv";
    options.outputFile    = "x.xml";
    options.kshell        = 1000;
    EXPECT_THROW(options.check(), gmx::InvalidInputError);
}

TEST_F(ShellForceConstantTest, EditFFFixesAndReleases)
{
    gmx::test::TestFileManager files;
    const auto input = files.getInputFilePath("forcefields/ACS-pol.xml");
    const auto fixed = files.getTemporaryFilePath("editff-fixed.xml");
    {
        gmx::test::CommandLine cmd;
        cmd.append("edit_ff");
        cmd.addOption("-ff", input);
        cmd.addOption("-o", fixed);
        cmd.addOption("-fixkshell", "alpha");
        cmd.addOption("-a", "c3_s");
        cmd.addOption("-kshell", "1000000");
        ASSERT_EQ(0, edit_ff(cmd.argc(), cmd.argv()));
    }
    ForceField pd;
    readForceField(fixed, &pd);
    auto shells = fixedShellForceConstants(&pd);
    ASSERT_EQ(1u, shells.size());
    EXPECT_EQ("c3_s", shells[0].shellType);
    EXPECT_NEAR(1000000, shells[0].kshell, 1e-6);
    EXPECT_TRUE(shells[0].alphaDependent);
    // As an edit without the command line
    EditFFOptions options;
    options.forceFieldFile = fixed;
    options.outputFile     = files.getTemporaryFilePath("editff-released.xml");
    options.releaseShells  = true;
    MsgHandler msghandler;
    msghandler.setPrintLevel(ACTStatus::Fatal);
    EXPECT_EQ(0, runEditFF(&msghandler, options));
    ForceField released;
    readForceField(options.outputFile, &released);
    EXPECT_TRUE(fixedShellForceConstants(&released).empty());
    const auto &alpha = released.findForcesConst(InteractionType::POLARIZATION).parametersConst().find(Identifier("c3_s"))->second.find("alpha")->second;
    EXPECT_EQ(Mutability::Fixed, alpha.mutability());
}

TEST_F(ShellForceConstantTest, EditFFWrongOptions)
{
    EditFFOptions options;
    options.forceFieldFile = "x.xml";
    options.kshell         = 1000;
    EXPECT_THROW(options.check(), gmx::InvalidInputError);
    options.shellFollows   = ShellFollows::Charge;
    EXPECT_NO_THROW(options.check());
    options.parameter      = "chi";
    EXPECT_THROW(options.check(), gmx::InvalidInputError);
}

} // namespace

} // namespace alexandria

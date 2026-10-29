/*
 * This source file is part of the Alexandria program.
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

/*! \internal \brief
 * Tests for replacing virtual sites on top of atoms by second charges,
 * see ForceField::convertVsite1.
 * \author David van der Spoel <david.vanderspoel@icm.uu.se>
 */
#include <string>

#include <gtest/gtest.h>

#include "act/basics/msg_handler.h"
#include "act/forcefield/forcefield.h"
#include "act/forcefield/forcefield_xml.h"

#include "testutils/testfilemanager.h"

namespace alexandria
{

namespace
{

class Vsite1ConversionTest : public ::testing::Test
{
protected:
    //! The force field with virtual sites on the ions
    ForceField pd_;

    void SetUp() override
    {
        auto fn = gmx::test::TestFileManager::getInputFilePath("../../alexandria/tests/PG-pol-ions.xml");
        readForceField(fn, &pd_, nullptr, false);
        ASSERT_TRUE(pd_.hasParticleType("v1Na+"));
        ASSERT_TRUE(pd_.hasParticleType("v1K+"));
    }

    //! Convert and check that there were no errors
    void convert()
    {
        MsgHandler msghandler;
        pd_.convertVsite1(&msghandler);
        EXPECT_TRUE(msghandler.ok());
    }

    //! \return A parameter of a virtual site type
    ForceFieldParameter *parameter(InteractionType itype, const std::string &id, const std::string &type)
    {
        auto params = pd_.findForces(itype)->parametersNoRemoval();
        auto pp     = params->find(Identifier(id));
        EXPECT_TRUE(params->end() != pp) << id;
        return &pp->second.find(type)->second;
    }

    //! Check that the virtual site on an atom became a second charge
    void expectConverted(const std::string &atom, const std::string &vsite,
                         double charge, double zeta)
    {
        EXPECT_FALSE(pd_.hasParticleType(vsite)) << vsite;
        auto pa = pd_.findParticleType(atom);
        EXPECT_TRUE(pa->hasChargeVS()) << atom;
        EXPECT_DOUBLE_EQ(charge, pa->chargeVS()) << atom;
        EXPECT_DOUBLE_EQ(zeta, pa->zetaVS()) << atom;
    }

    //! Check that a virtual site is kept, and its atom has no second charge
    void expectKept(const std::string &atom, const std::string &vsite)
    {
        EXPECT_TRUE(pd_.hasParticleType(vsite)) << vsite;
        EXPECT_FALSE(pd_.findParticleType(atom)->hasChargeVS()) << atom;
    }
};

// Test: Virtual sites with only a charge become second charges with the
// charge and distribution width of the virtual site, which are removed
TEST_F(Vsite1ConversionTest, ChargeOnlyVsitesConverted)
{
    double q    = pd_.findParticleType("v1Na+")->charge();
    double zeta = parameter(InteractionType::ELECTROSTATICS, "v1Na+_z", "zeta")->value();
    convert();
    expectConverted("Na+", "v1Na+", q, zeta);
    EXPECT_TRUE(pd_.findParticleType("K+")->hasChargeVS());
    EXPECT_FALSE(pd_.findForces(InteractionType::ELECTROSTATICS)->parameterExists(Identifier("v1Na+_z")));
    EXPECT_FALSE(pd_.findForces(InteractionType::VDW)->parameterExists(Identifier("v1Na+")));
    EXPECT_FALSE(pd_.interactionPresent(InteractionType::VSITE1));
}

// Test: A virtual site with Van der Waals interactions is kept
TEST_F(Vsite1ConversionTest, VsiteWithVdwKept)
{
    parameter(InteractionType::VDW, "v1Na+", "epsilon")->forceSetValue(0.5);
    convert();
    expectKept("Na+", "v1Na+");
    EXPECT_TRUE(pd_.findForces(InteractionType::VSITE1)->parameterExists(Identifier(InteractionType::VSITE1, "Na+!v1Na+", CanSwap::Yes)));
    // The others are converted
    EXPECT_TRUE(pd_.findParticleType("K+")->hasChargeVS());
}

// Test: A virtual site with Van der Waals parameters that are zero but can
// be trained is kept
TEST_F(Vsite1ConversionTest, VsiteWithTrainableVdwKept)
{
    parameter(InteractionType::VDW, "v1Na+", "epsilon")->setMutability(Mutability::Bounded);
    convert();
    expectKept("Na+", "v1Na+");
    EXPECT_TRUE(pd_.findParticleType("K+")->hasChargeVS());
}

// Test: A virtual site with a shell is kept
TEST_F(Vsite1ConversionTest, VsiteWithShellKept)
{
    pd_.findParticleType("v1Na+")->setOption("poltype", "v1Na+");
    convert();
    expectKept("Na+", "v1Na+");
    EXPECT_TRUE(pd_.findParticleType("K+")->hasChargeVS());
}

// Test: A virtual site that is part of another interaction is kept
TEST_F(Vsite1ConversionTest, VsiteWithBondKept)
{
    ForceFieldParameter kb("kJ/mol/nm2", 1000, 0, 0, 0, 2000, Mutability::Fixed, true, true);
    pd_.findForces(InteractionType::BONDS)->addParameter(Identifier({ "Na+", "v1Na+" }, { 1 }, CanSwap::Yes),
                                                         "kb", kb);
    convert();
    expectKept("Na+", "v1Na+");
    EXPECT_TRUE(pd_.findParticleType("K+")->hasChargeVS());
}

// Test: A virtual site type on two atom types becomes a second charge of both
TEST_F(Vsite1ConversionTest, SharedVsiteConverted)
{
    ForceFieldParameter vs1a("", 1, 0, 0, 1, 1, Mutability::Fixed, true, false);
    pd_.findForces(InteractionType::VSITE1)->addParameter(Identifier(InteractionType::VSITE1, "Li+!v1Na+", CanSwap::Yes),
                                                          "vs1a", vs1a);
    // Li+ has its own virtual site, which cannot be on it as well
    pd_.findForces(InteractionType::VSITE1)->parametersNoRemoval()->erase(Identifier(InteractionType::VSITE1, "Li+!v1Li+", CanSwap::Yes));
    double q    = pd_.findParticleType("v1Na+")->charge();
    double zeta = parameter(InteractionType::ELECTROSTATICS, "v1Na+_z", "zeta")->value();
    convert();
    expectConverted("Na+", "v1Na+", q, zeta);
    expectConverted("Li+", "v1Na+", q, zeta);
}

// Test: Without combination rules for the Van der Waals interactions,
// the pair parameters of virtual sites cannot be checked, all are kept
TEST_F(Vsite1ConversionTest, NoCombinationRulesKept)
{
    auto                    vdw = pd_.findForces(InteractionType::VDW);
    ForceFieldParameterList noRules(potentialToString(vdw->potential()), vdw->canSwap());
    for(const auto &p : vdw->parametersConst())
    {
        for(const auto &pp : p.second)
        {
            noRules.addParameter(p.first, pp.first, pp.second);
        }
    }
    (*pd_.forces())[InteractionType::VDW] = noRules;
    MsgHandler msghandler;
    pd_.convertVsite1(&msghandler);
    expectKept("Na+", "v1Na+");
    expectKept("K+", "v1K+");
}

} // namespace

} // namespace alexandria

/*
 * This source file is part of the Alexandria program.
 *
 * Copyright (C) 2020-2025
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
 * Implements part of the alexandria program.
 * \author David van der Spoel <david.vanderspoel@icm.uu.se>
 */
#include <gtest/gtest.h>

#include "../interactiontype.h"

#include "testutils/cmdlinetest.h"
#include "testutils/refdata.h"
#include "testutils/testasserts.h"
#include "testutils/testfilemanager.h"

#include "gromacs/utility/exceptions.h"

namespace alexandria
{

namespace
{

TEST(InteractionTypeTest, stringToInteractionType)
{
    InteractionType iType;
    EXPECT_TRUE(stringToInteractionType("BONDCORRECTIONS", &iType));
    EXPECT_TRUE(iType == InteractionType::BONDCORRECTIONS);
    EXPECT_TRUE(stringToInteractionType("BONDS", &iType));
    EXPECT_FALSE(iType == InteractionType::ANGLES);
    EXPECT_FALSE(stringToInteractionType("FOO", &iType));
}

TEST(InteractionTypeTest, stringToInteractionTypeAllValid)
{
    // Verify that every string produced by interactionTypeToString is
    // accepted back by stringToInteractionType and yields the original enum.
    const std::vector<InteractionType> allTypes = {
        InteractionType::BONDS,
        InteractionType::ANGLES,
        InteractionType::LINEAR_ANGLES,
        InteractionType::PROPER_DIHEDRALS,
        InteractionType::IMPROPER_DIHEDRALS,
        InteractionType::VDW,
        InteractionType::DISPERSION,
        InteractionType::EXCHANGE,
        InteractionType::VDWCORRECTION,
        InteractionType::ELECTROSTATICS,
        InteractionType::POLARIZATION,
        InteractionType::POSITION_RESTRAINT,
        InteractionType::INDUCTION,
        InteractionType::INDUCTIONCORRECTION,
        InteractionType::CHARGETRANSFER,
        InteractionType::ALLELEC,
        InteractionType::EXCHIND,
        InteractionType::EPOT,
        InteractionType::CONSTR,
        InteractionType::VSITE1,
        InteractionType::VSITE2,
        InteractionType::VSITE2FD,
        InteractionType::VSITE3,
        InteractionType::VSITE3S,
        InteractionType::VSITE3FD,
        InteractionType::VSITE3FAD,
        InteractionType::VSITE3OUT,
        InteractionType::VSITE3OUTS,
        InteractionType::VSITE4,
        InteractionType::VSITE4S,
        InteractionType::VSITE4S3,
        InteractionType::BONDCORRECTIONS,
        InteractionType::ELECTRONEGATIVITYEQUALIZATION,
        InteractionType::CHARGE
    };
    for (auto it : allTypes)
    {
        InteractionType result;
        const auto &name = interactionTypeToString(it);
        EXPECT_TRUE(stringToInteractionType(name, &result))
            << "Failed to parse string '" << name << "' back to InteractionType";
        EXPECT_EQ(it, result)
            << "Round-trip mismatch for '" << name << "'";
    }
}

TEST(InteractionTypeTest, interactionTypeToStringKnownValues)
{
    EXPECT_EQ("BONDS",              interactionTypeToString(InteractionType::BONDS));
    EXPECT_EQ("ANGLES",             interactionTypeToString(InteractionType::ANGLES));
    EXPECT_EQ("PROPER_DIHEDRALS",   interactionTypeToString(InteractionType::PROPER_DIHEDRALS));
    EXPECT_EQ("IMPROPER_DIHEDRALS", interactionTypeToString(InteractionType::IMPROPER_DIHEDRALS));
    EXPECT_EQ("VANDERWAALS",        interactionTypeToString(InteractionType::VDW));
    EXPECT_EQ("COULOMB",            interactionTypeToString(InteractionType::ELECTROSTATICS));
    EXPECT_EQ("BONDCORRECTIONS",    interactionTypeToString(InteractionType::BONDCORRECTIONS));
    EXPECT_EQ("CHARGE",             interactionTypeToString(InteractionType::CHARGE));
    EXPECT_EQ("EPOT",               interactionTypeToString(InteractionType::EPOT));
}

TEST(InteractionTypeTest, interactionTypeToDescriptionNotEmpty)
{
    const std::vector<InteractionType> allTypes = {
        InteractionType::BONDS,
        InteractionType::ANGLES,
        InteractionType::LINEAR_ANGLES,
        InteractionType::PROPER_DIHEDRALS,
        InteractionType::IMPROPER_DIHEDRALS,
        InteractionType::VDW,
        InteractionType::DISPERSION,
        InteractionType::EXCHANGE,
        InteractionType::VDWCORRECTION,
        InteractionType::ELECTROSTATICS,
        InteractionType::POLARIZATION,
        InteractionType::POSITION_RESTRAINT,
        InteractionType::INDUCTION,
        InteractionType::INDUCTIONCORRECTION,
        InteractionType::CHARGETRANSFER,
        InteractionType::ALLELEC,
        InteractionType::EXCHIND,
        InteractionType::EPOT,
        InteractionType::CONSTR,
        InteractionType::VSITE1,
        InteractionType::VSITE2,
        InteractionType::VSITE2FD,
        InteractionType::VSITE3,
        InteractionType::VSITE3S,
        InteractionType::VSITE3FD,
        InteractionType::VSITE3FAD,
        InteractionType::VSITE3OUT,
        InteractionType::VSITE3OUTS,
        InteractionType::VSITE4,
        InteractionType::VSITE4S,
        InteractionType::VSITE4S3,
        InteractionType::BONDCORRECTIONS,
        InteractionType::ELECTRONEGATIVITYEQUALIZATION,
        InteractionType::CHARGE
    };
    for (auto it : allTypes)
    {
        EXPECT_FALSE(interactionTypeToDescription(it).empty())
            << "Description is empty for " << interactionTypeToString(it);
    }
}

TEST(InteractionTypeTest, interactionTypeToNatomsReturnsFive)
{
    EXPECT_EQ(5, interactionTypeToNatoms(InteractionType::VSITE4));
    EXPECT_EQ(5, interactionTypeToNatoms(InteractionType::VSITE4S));
    EXPECT_EQ(5, interactionTypeToNatoms(InteractionType::VSITE4S3));
}

TEST(InteractionTypeTest, interactionTypeToNatomsReturnsFour)
{
    EXPECT_EQ(4, interactionTypeToNatoms(InteractionType::PROPER_DIHEDRALS));
    EXPECT_EQ(4, interactionTypeToNatoms(InteractionType::IMPROPER_DIHEDRALS));
    EXPECT_EQ(4, interactionTypeToNatoms(InteractionType::VSITE3));
    EXPECT_EQ(4, interactionTypeToNatoms(InteractionType::VSITE3S));
    EXPECT_EQ(4, interactionTypeToNatoms(InteractionType::VSITE3FD));
    EXPECT_EQ(4, interactionTypeToNatoms(InteractionType::VSITE3FAD));
    EXPECT_EQ(4, interactionTypeToNatoms(InteractionType::VSITE3OUT));
    EXPECT_EQ(4, interactionTypeToNatoms(InteractionType::VSITE3OUTS));
}

TEST(InteractionTypeTest, interactionTypeToNatomsReturnsThree)
{
    EXPECT_EQ(3, interactionTypeToNatoms(InteractionType::ANGLES));
    EXPECT_EQ(3, interactionTypeToNatoms(InteractionType::LINEAR_ANGLES));
    EXPECT_EQ(3, interactionTypeToNatoms(InteractionType::VSITE2));
    EXPECT_EQ(3, interactionTypeToNatoms(InteractionType::VSITE2FD));
}

TEST(InteractionTypeTest, interactionTypeToNatomsReturnsTwo)
{
    EXPECT_EQ(2, interactionTypeToNatoms(InteractionType::BONDS));
    EXPECT_EQ(2, interactionTypeToNatoms(InteractionType::CONSTR));
    EXPECT_EQ(2, interactionTypeToNatoms(InteractionType::VSITE1));
    EXPECT_EQ(2, interactionTypeToNatoms(InteractionType::BONDCORRECTIONS));
}

TEST(InteractionTypeTest, interactionTypeToNatomsReturnsOne)
{
    EXPECT_EQ(1, interactionTypeToNatoms(InteractionType::VDW));
    EXPECT_EQ(1, interactionTypeToNatoms(InteractionType::DISPERSION));
    EXPECT_EQ(1, interactionTypeToNatoms(InteractionType::EXCHANGE));
    EXPECT_EQ(1, interactionTypeToNatoms(InteractionType::VDWCORRECTION));
    EXPECT_EQ(1, interactionTypeToNatoms(InteractionType::INDUCTIONCORRECTION));
    EXPECT_EQ(1, interactionTypeToNatoms(InteractionType::CHARGETRANSFER));
    EXPECT_EQ(1, interactionTypeToNatoms(InteractionType::POLARIZATION));
    EXPECT_EQ(1, interactionTypeToNatoms(InteractionType::INDUCTION));
    EXPECT_EQ(1, interactionTypeToNatoms(InteractionType::ELECTROSTATICS));
    EXPECT_EQ(1, interactionTypeToNatoms(InteractionType::ALLELEC));
    EXPECT_EQ(1, interactionTypeToNatoms(InteractionType::EXCHIND));
    EXPECT_EQ(1, interactionTypeToNatoms(InteractionType::ELECTRONEGATIVITYEQUALIZATION));
}

TEST(InteractionTypeTest, interactionTypeToNatomsReturnsZero)
{
    EXPECT_EQ(0, interactionTypeToNatoms(InteractionType::CHARGE));
    EXPECT_EQ(0, interactionTypeToNatoms(InteractionType::EPOT));
}

TEST(InteractionTypeTest, isVsiteReturnsTrueForVsiteTypes)
{
    EXPECT_TRUE(isVsite(InteractionType::VSITE1));
    EXPECT_TRUE(isVsite(InteractionType::VSITE2));
    EXPECT_TRUE(isVsite(InteractionType::VSITE2FD));
    EXPECT_TRUE(isVsite(InteractionType::VSITE3));
    EXPECT_TRUE(isVsite(InteractionType::VSITE3S));
    EXPECT_TRUE(isVsite(InteractionType::VSITE3FD));
    EXPECT_TRUE(isVsite(InteractionType::VSITE3FAD));
    EXPECT_TRUE(isVsite(InteractionType::VSITE3OUT));
    EXPECT_TRUE(isVsite(InteractionType::VSITE3OUTS));
    EXPECT_TRUE(isVsite(InteractionType::VSITE4));
    EXPECT_TRUE(isVsite(InteractionType::VSITE4S));
    EXPECT_TRUE(isVsite(InteractionType::VSITE4S3));
}

TEST(InteractionTypeTest, isVsiteReturnsFalseForNonVsiteTypes)
{
    EXPECT_FALSE(isVsite(InteractionType::BONDS));
    EXPECT_FALSE(isVsite(InteractionType::ANGLES));
    EXPECT_FALSE(isVsite(InteractionType::PROPER_DIHEDRALS));
    EXPECT_FALSE(isVsite(InteractionType::VDW));
    EXPECT_FALSE(isVsite(InteractionType::ELECTROSTATICS));
    EXPECT_FALSE(isVsite(InteractionType::CHARGE));
    EXPECT_FALSE(isVsite(InteractionType::EPOT));
    EXPECT_FALSE(isVsite(InteractionType::BONDCORRECTIONS));
}

}

}

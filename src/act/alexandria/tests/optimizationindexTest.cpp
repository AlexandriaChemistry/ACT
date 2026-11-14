/*
 * This source file is part of the Alexandria Chemistry Toolkit.
 *
 * Copyright (C) 2024-2026
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
 * Tests for the OptimizationIndex class.
 * \author David van der Spoel <david.vanderspoel@icm.uu.se>
 */

#include "actpre.h"

#include <string>

#include <gtest/gtest.h>

#include "act/alexandria/optimizationindex.h"
#include "act/basics/identifier.h"
#include "act/basics/interactiontype.h"

#include "testutils/testasserts.h"

namespace alexandria
{

namespace
{

TEST(OptimizationIndexTest, DefaultConstructor)
{
    const OptimizationIndex oi;
    // Default interaction type is CHARGE
    EXPECT_EQ(oi.iType(), InteractionType::CHARGE);
    EXPECT_TRUE(oi.particleType().empty());
    EXPECT_TRUE(oi.parameterType().empty());
}

TEST(OptimizationIndexTest, InteractionTypeConstructor)
{
    const Identifier id("h");
    const OptimizationIndex oi(InteractionType::BONDS, id, "bondlength");
    EXPECT_EQ(oi.iType(), InteractionType::BONDS);
    EXPECT_EQ(oi.id().id(), "h");
    EXPECT_EQ(oi.parameterType(), "bondlength");
    EXPECT_TRUE(oi.particleType().empty());
}

TEST(OptimizationIndexTest, ParticleTypeConstructor)
{
    const OptimizationIndex oi("opls_116", "charge");
    EXPECT_EQ(oi.particleType(), "opls_116");
    EXPECT_EQ(oi.parameterType(), "charge");
    // Default iType should be CHARGE since not set by this constructor
    EXPECT_EQ(oi.iType(), InteractionType::CHARGE);
}

TEST(OptimizationIndexTest, NameForChargeType)
{
    const OptimizationIndex oi("opls_116", "charge");
    const std::string name = oi.name();
    // For CHARGE type: "particleType parameterType"
    EXPECT_EQ(name, "opls_116 charge");
}

TEST(OptimizationIndexTest, NameForInteractionType)
{
    const Identifier id("h");
    const OptimizationIndex oi(InteractionType::BONDS, id, "bondlength");
    const std::string name = oi.name();
    // For non-CHARGE: "parameterId parameterType"
    EXPECT_EQ(name, "h bondlength");
}

TEST(OptimizationIndexTest, ParameterName)
{
    const Identifier id("c-o");
    const OptimizationIndex oi(InteractionType::BONDS, id, "force_constant");
    EXPECT_EQ(oi.parameterName(), "c-o");
}

TEST(OptimizationIndexTest, ForceFieldParameterInitiallyNull)
{
    OptimizationIndex oi;
    EXPECT_EQ(oi.forceFieldParameter(), nullptr);
}

TEST(OptimizationIndexTest, MultiAtomIdentifier)
{
    const std::vector<std::string> atoms = {"C", "O"};
    const std::vector<double> bondorders = {2.0};
    const Identifier id(atoms, bondorders, CanSwap::Yes);
    const OptimizationIndex oi(InteractionType::BONDS, id, "bondlength");

    EXPECT_EQ(oi.iType(), InteractionType::BONDS);
    EXPECT_EQ(oi.parameterType(), "bondlength");
    // The name should include the identifier
    const std::string name = oi.name();
    EXPECT_FALSE(name.empty());
}

TEST(OptimizationIndexTest, DifferentInteractionTypes)
{
    const Identifier id("c");
    const OptimizationIndex oiBonds(InteractionType::BONDS, id, "bondlength");
    EXPECT_EQ(oiBonds.iType(), InteractionType::BONDS);

    const OptimizationIndex oiAngles(InteractionType::ANGLES, id, "angle");
    EXPECT_EQ(oiAngles.iType(), InteractionType::ANGLES);

    const OptimizationIndex oiVdw(InteractionType::VDW, id, "sigma");
    EXPECT_EQ(oiVdw.iType(), InteractionType::VDW);
}

} // namespace

} // namespace alexandria

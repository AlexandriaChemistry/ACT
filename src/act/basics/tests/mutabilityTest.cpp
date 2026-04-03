/*
 * This source file is part of the Alexandria Chemistry Toolkit.
 *
 * Copyright (C) 2020,2025
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
 * Tests for mutability.h/cpp (mutabilityName and nameToMutability functions).
 * \author David van der Spoel <david.vanderspoel@icm.uu.se>
 */
#include <gtest/gtest.h>

#include "../mutability.h"

#include "testutils/testasserts.h"

namespace alexandria
{

namespace
{

TEST(MutabilityTest, mutabilityName)
{
    EXPECT_EQ("Fixed",     mutabilityName(Mutability::Fixed));
    EXPECT_EQ("Dependent", mutabilityName(Mutability::Dependent));
    EXPECT_EQ("ACM",       mutabilityName(Mutability::ACM));
    EXPECT_EQ("Bounded",   mutabilityName(Mutability::Bounded));
    EXPECT_EQ("Free",      mutabilityName(Mutability::Free));
}

TEST(MutabilityTest, nameToMutabilityValid)
{
    Mutability mut;

    EXPECT_TRUE(nameToMutability("Fixed", &mut));
    EXPECT_EQ(Mutability::Fixed, mut);

    EXPECT_TRUE(nameToMutability("Dependent", &mut));
    EXPECT_EQ(Mutability::Dependent, mut);

    EXPECT_TRUE(nameToMutability("ACM", &mut));
    EXPECT_EQ(Mutability::ACM, mut);

    EXPECT_TRUE(nameToMutability("Bounded", &mut));
    EXPECT_EQ(Mutability::Bounded, mut);

    EXPECT_TRUE(nameToMutability("Free", &mut));
    EXPECT_EQ(Mutability::Free, mut);
}

TEST(MutabilityTest, nameToMutabilityInvalid)
{
    Mutability mut;
    EXPECT_FALSE(nameToMutability("InvalidMutability", &mut));
    // Input is case-sensitive
    EXPECT_FALSE(nameToMutability("fixed",             &mut));
    EXPECT_FALSE(nameToMutability("acm",               &mut));
    EXPECT_FALSE(nameToMutability("",                  &mut));
}

TEST(MutabilityTest, roundTrip)
{
    for (auto m : { Mutability::Fixed,
                    Mutability::Dependent,
                    Mutability::ACM,
                    Mutability::Bounded,
                    Mutability::Free })
    {
        Mutability out;
        EXPECT_TRUE(nameToMutability(mutabilityName(m), &out));
        EXPECT_EQ(m, out);
    }
}

} // namespace

} // namespace alexandria

/*
 * This source file is part of the Alexandria Chemistry Toolkit.
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
 * Tests that the random numbers are the same on all platforms. The
 * expected values follow from the Mersenne twister, of which the output
 * is specified by the C++ standard, and were checked with an independent
 * implementation of it.
 * \author David van der Spoel <david.vanderspoel@icm.uu.se>
 */

#include "actpre.h"

#include <algorithm>
#include <vector>

#include <gtest/gtest.h>

#include "act/utility/portable_random.h"

namespace alexandria
{

namespace
{

// Test: The seeds are the upper 31 bits of the output of the engine
TEST(PortableRandomTest, Seed)
{
    // The first output of the default engine is 3499211612
    std::mt19937 gen(5489u);
    EXPECT_EQ(1749605806, portableSeed(&gen));
    EXPECT_EQ(290934651, portableSeed(&gen));
    EXPECT_EQ(1945173367, portableSeed(&gen));
}

// Test: The indices are in range and the same on all platforms
TEST(PortableRandomTest, Index)
{
    std::mt19937        gen(1993);
    std::vector<size_t> expected = { 1, 2, 2, 9, 9 };
    for (auto e : expected)
    {
        EXPECT_EQ(e, portableIndex(&gen, 10));
    }
    for (int i = 0; i < 1000; i++)
    {
        EXPECT_LT(portableIndex(&gen, 7), 7u);
    }
    EXPECT_EQ(0u, portableIndex(&gen, 1));
}

// Test: A shuffle is a permutation that is the same on all platforms
TEST(PortableRandomTest, Shuffle)
{
    std::mt19937     gen(1993);
    std::vector<int> v = { 0, 1, 2, 3, 4, 5, 6, 7, 8, 9 };
    portableShuffle(&v, &gen);
    EXPECT_EQ(std::vector<int>({ 3, 8, 4, 7, 0, 5, 6, 9, 2, 1 }), v);
    std::vector<int> empty;
    portableShuffle(&empty, &gen);
    EXPECT_TRUE(empty.empty());
}

} // namespace

} // namespace alexandria

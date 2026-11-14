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
 * Tests of how evenly the middlemen of train_ff share the individuals of a
 * generation, see queuebalance.h.
 * \author David van der Spoel <david.vanderspoel@icm.uu.se>
 */
#include "act/alexandria/queuebalance.h"

#include <array>

#include <gtest/gtest.h>

namespace alexandria
{

namespace
{

// Test: 256 individuals over 63 middlemen, one at a time: some compute 5,
// most 4, and 32 or 64 middlemen would share them evenly
TEST(QueueBalanceTest, OneAtATime)
{
    const auto qb = queueBalance(256, 1, 63);
    EXPECT_EQ(256u, qb.groups);
    EXPECT_EQ(5u, qb.rounds);
    EXPECT_NEAR(1 - 256.0/315, qb.idle, 1e-12);
    EXPECT_EQ(32u, qb.fewer);
    EXPECT_EQ(64u, qb.more);
}

// Test: With -lockstep 4 there are 64 groups, so one of 63 middlemen computes
// two of them while the others wait half the time
TEST(QueueBalanceTest, LockstepOneMiddlemanShort)
{
    const auto qb = queueBalance(256, 4, 63);
    EXPECT_EQ(64u, qb.groups);
    EXPECT_EQ(2u, qb.rounds);
    EXPECT_NEAR(1 - 64.0/126, qb.idle, 1e-12);
    EXPECT_EQ(32u, qb.fewer);
    EXPECT_EQ(64u, qb.more);
}

// Test: As many middlemen as groups, or a divisor of them, do not wait
TEST(QueueBalanceTest, Even)
{
    for (auto c : { std::array<size_t, 3>{ 256, 4, 64 }, { 256, 2, 128 }, { 256, 1, 64 }, { 480, 4, 40 } })
    {
        const auto qb = queueBalance(c[0], c[1], c[2]);
        EXPECT_EQ(0.0, qb.idle) << c[0] << " " << c[1] << " " << c[2];
    }
}

// Test: More middlemen than groups leave some without work
TEST(QueueBalanceTest, MoreMiddlemenThanGroups)
{
    const auto qb = queueBalance(256, 4, 128);
    EXPECT_EQ(64u, qb.groups);
    EXPECT_EQ(1u, qb.rounds);
    EXPECT_NEAR(0.5, qb.idle, 1e-12);
    EXPECT_EQ(64u, qb.fewer);
}

// Test: A population that is not a multiple of -lockstep has a smaller last group
TEST(QueueBalanceTest, LastGroupSmaller)
{
    const auto qb = queueBalance(10, 4, 3);
    EXPECT_EQ(3u, qb.groups);
    EXPECT_EQ(1u, qb.rounds);
    EXPECT_EQ(0.0, qb.idle);
}

// Test: Nothing to share
TEST(QueueBalanceTest, Empty)
{
    const auto qb = queueBalance(0, 4, 63);
    EXPECT_EQ(0u, qb.groups);
    EXPECT_EQ(0.0, qb.idle);
}

} // namespace

} // namespace alexandria

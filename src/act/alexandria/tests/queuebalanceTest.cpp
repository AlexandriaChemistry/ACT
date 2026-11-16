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

#include <algorithm>
#include <array>
#include <vector>

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

// Test: The settings of the benchmarks in the manual, a population of 128
TEST(LockstepPlanTest, Benchmarks)
{
    // Four GH200 (arrhenius): groups of 16 on 8 middlemen
    auto plan = planLockstep(128, 129, 4);
    EXPECT_EQ(16u, plan.lockstep);
    EXPECT_EQ(8u, plan.middlemen);
    EXPECT_EQ(9u, plan.processes);
    EXPECT_EQ(120u, plan.idle);
    EXPECT_EQ(1u, plan.rounds);
    EXPECT_EQ(0u, plan.better);
    // Eight MI250X GCDs (dardel): groups of 8 on 16 middlemen
    plan = planLockstep(128, 129, 8);
    EXPECT_EQ(8u, plan.lockstep);
    EXPECT_EQ(16u, plan.middlemen);
    EXPECT_EQ(17u, plan.processes);
}

// Test: A large population on one GPU takes several rounds if the
// processes are too few, and says how many would do it in one
TEST(LockstepPlanTest, TooFewProcesses)
{
    auto plan = planLockstep(512, 5, 1);
    EXPECT_EQ(16u, plan.lockstep);
    EXPECT_EQ(4u, plan.middlemen);
    EXPECT_EQ(8u, plan.rounds);
    EXPECT_EQ(0u, plan.idle);
    EXPECT_EQ(33u, plan.better);
    plan = planLockstep(512, 64, 1);
    EXPECT_EQ(16u, plan.lockstep);
    EXPECT_EQ(32u, plan.middlemen);
    EXPECT_EQ(1u, plan.rounds);
    EXPECT_EQ(31u, plan.idle);
}

// Test: A small population still uses at least two middlemen per GPU
TEST(LockstepPlanTest, SmallPopulation)
{
    auto plan = planLockstep(12, 4, 1);
    EXPECT_EQ(6u, plan.lockstep);
    EXPECT_EQ(2u, plan.middlemen);
    EXPECT_EQ(1u, plan.idle);
    // Groups of equal size: 20 in 2 groups of 10, not 16 and 4
    plan = planLockstep(20, 9, 1);
    EXPECT_EQ(10u, plan.lockstep);
    EXPECT_EQ(2u, plan.middlemen);
}

// Test: The number of middlemen per GPU stays below the limit of the
// multi-process service, even if the groups get several rounds
TEST(LockstepPlanTest, ManyIndividualsPerGpu)
{
    const auto plan = planLockstep(4096, 257, 1);
    EXPECT_EQ(16u, plan.lockstep);
    // 256 groups take 7 rounds on 40 middlemen, as on 37
    EXPECT_LE(plan.middlemen, maxMiddlemenPerGpu);
    EXPECT_EQ(37u, plan.middlemen);
    EXPECT_EQ(7u, plan.rounds);
}

// Test: Without GPU or with one process nothing is planned
TEST(LockstepPlanTest, NoGpu)
{
    auto plan = planLockstep(128, 33, 0);
    EXPECT_EQ(1u, plan.lockstep);
    EXPECT_EQ(32u, plan.middlemen);
    EXPECT_EQ(0u, plan.idle);
    plan = planLockstep(128, 1, 1);
    EXPECT_EQ(1u, plan.lockstep);
}

// Test: The working processes are spread over the GPUs, also when the
// processes were bound to the GPUs in blocks, as with srun
// --gpu-bind=single:4, and when the local ranks alternate between them
TEST(LockstepPlanTest, SpreadOverGpus)
{
    // 16 processes on 4 GPUs in blocks of 4, 9 of them work
    std::vector<int> blocks;
    for (int r = 0; r < 16; r++)
    {
        blocks.push_back(r/4);
    }
    auto working = spreadOverGpus(blocks, 9);
    std::array<int, 4> count = { 0, 0, 0, 0 };
    for (size_t r = 0; r < blocks.size(); r++)
    {
        count[blocks[r]] += working[r];
    }
    EXPECT_EQ(1, working[0]);
    EXPECT_EQ((std::array<int, 4>{ 3, 2, 2, 2 }), count);
    // Ranks on the GPUs in turn, the lowest ones work, as without spreading
    std::vector<int> alternating;
    for (int r = 0; r < 16; r++)
    {
        alternating.push_back(r % 4);
    }
    working = spreadOverGpus(alternating, 9);
    for (size_t r = 0; r < alternating.size(); r++)
    {
        EXPECT_EQ(r < 9 ? 1 : 0, working[r]) << "rank " << r;
    }
    // All of them if as many work as there are
    working = spreadOverGpus(blocks, 16);
    EXPECT_EQ(16, std::count(working.begin(), working.end(), 1));
}

} // namespace alexandria

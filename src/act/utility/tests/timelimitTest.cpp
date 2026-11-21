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
 * Tests of the limit on the wall time of a program.
 * \author David van der Spoel <david.vanderspoel@icm.uu.se>
 */
#include "actpre.h"

#include <chrono>
#include <gtest/gtest.h>

#include "act/utility/time_limit.h"

namespace alexandria
{

namespace
{

// Test: Without a limit the time is never up
TEST(TimeLimitTest, NoLimit)
{
    for (double hours : { 0.0, -1.0 })
    {
        const TimeLimit limit(hours);
        EXPECT_FALSE(limit.limited());
        EXPECT_FALSE(limit.expired());
        EXPECT_EQ(0, limit.hours());
    }
    const TimeLimit byDefault;
    EXPECT_FALSE(byDefault.limited());
    EXPECT_FALSE(byDefault.expired());
    // Nor is it up when the limit counts from long ago
    const TimeLimit old(0, TimeLimit::Clock::now() - std::chrono::hours(100));
    EXPECT_FALSE(old.expired());
}

// Test: The time is up once the limit has passed since the moment it counts from
TEST(TimeLimitTest, CountsFromAMoment)
{
    const auto now = TimeLimit::Clock::now();
    const TimeLimit ample(2, now);
    EXPECT_TRUE(ample.limited());
    EXPECT_DOUBLE_EQ(2, ample.hours());
    EXPECT_FALSE(ample.expired());
    EXPECT_GE(ample.elapsedHours(), 0);
    EXPECT_LT(ample.elapsedHours(), 1);

    // Two hours, of which three have passed, and of which one has
    const TimeLimit over(2, now - std::chrono::hours(3));
    EXPECT_TRUE(over.expired());
    EXPECT_GE(over.elapsedHours(), 3);
    const TimeLimit under(2, now - std::chrono::hours(1));
    EXPECT_FALSE(under.expired());
    EXPECT_GE(under.elapsedHours(), 1);
    EXPECT_LT(under.elapsedHours(), 2);
}

// Test: A limit without a moment counts from the start of the program,
// which is the same for every limit and lies in the past
TEST(TimeLimitTest, CountsFromTheStartOfTheProgram)
{
    const auto start = TimeLimit::programStart();
    EXPECT_EQ(start, TimeLimit::programStart());
    EXPECT_LE(start, TimeLimit::Clock::now());
    // The program has run for a nanosecond, not for a year
    EXPECT_TRUE(TimeLimit(1e-9/3600).expired());
    EXPECT_FALSE(TimeLimit(24*365).expired());
    const TimeLimit first(1), second(1);
    EXPECT_NEAR(first.elapsedHours(), second.elapsedHours(), 1.0/3600);
    // A copy has the same limit, as the processes of a parallel run have
    const TimeLimit copy = TimeLimit(1e-9/3600);
    EXPECT_TRUE(copy.expired());
}

} // namespace

} // namespace alexandria

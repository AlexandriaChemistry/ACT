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
 * Tests of the integer coordinates.
 * \author David van der Spoel <david.vanderspoel@icm.uu.se>
 */
#include "actpre.h"

#include "act/basics/coordinate.h"

#include <cmath>
#include <limits>

#include <gtest/gtest.h>

namespace alexandria
{

namespace
{

// Test: The unit is a power of two, such that scaling is exact
TEST(CoordinateTest, Unit)
{
    EXPECT_EQ(std::ldexp(1.0, -56), coordinateUnit);
    EXPECT_EQ(1.0, coordinateUnit*coordinateScale);
    EXPECT_EQ(std::ldexp(1.0, 63), coordinateMax*coordinateScale);
}

// Test: Coordinates are rounded to the nearest unit, also negative ones
TEST(CoordinateTest, Rounding)
{
    EXPECT_EQ(0, toCoordinate(0.0));
    EXPECT_EQ(Coordinate(1) << 56, toCoordinate(1.0));
    EXPECT_EQ(-(Coordinate(1) << 56), toCoordinate(-1.0));
    EXPECT_EQ(3, toCoordinate(3*coordinateUnit));
    EXPECT_EQ(-3, toCoordinate(-3*coordinateUnit));
    EXPECT_EQ(3, toCoordinate(3.4*coordinateUnit));
    EXPECT_EQ(-4, toCoordinate(-3.6*coordinateUnit));
    // Halfway cases to even
    EXPECT_EQ(2, toCoordinate(2.5*coordinateUnit));
    EXPECT_EQ(-2, toCoordinate(-2.5*coordinateUnit));
}

// Test: Converting back gives the coordinate within half a unit, and
// coordinates that are a multiple of the unit exactly
TEST(CoordinateTest, RoundTrip)
{
    for (double x : { 0.0, 1e-9, 0.1, -0.3, 1.234567891, -12.3456789, 100.0, -127.0 })
    {
        EXPECT_NEAR(x, fromCoordinate(toCoordinate(x)), 0.5*coordinateUnit*(1 + 1e-6)) << x;
    }
    for (double x : { 0.5, -0.125, 3.0, 64.0 })
    {
        EXPECT_EQ(x, fromCoordinate(toCoordinate(x))) << x;
    }
}

// Test: Coordinates too far from the origin or not finite are an error
TEST(CoordinateTest, Range)
{
    EXPECT_NO_THROW(toCoordinate(0.999*coordinateMax));
    EXPECT_NO_THROW(toCoordinate(-0.999*coordinateMax));
    EXPECT_THROW(toCoordinate(coordinateMax), gmx::RangeError);
    EXPECT_THROW(toCoordinate(-1.001*coordinateMax), gmx::RangeError);
    EXPECT_THROW(toCoordinate(std::nan("")), gmx::RangeError);
    EXPECT_THROW(toCoordinate(HUGE_VAL), gmx::RangeError);
}

// Test: Differences are exact wherever the coordinates are, also across the
// ends of the range, as long as they are less than the range apart
TEST(CoordinateTest, Difference)
{
    EXPECT_EQ(3, coordinateDifference(5, 2));
    EXPECT_EQ(-3, coordinateDifference(2, 5));
    Coordinate big   = toCoordinate(120.0);
    Coordinate small = toCoordinate(0.1);
    EXPECT_EQ(small, coordinateDifference(big, big - small));
    // A difference of 0.1 nm far from the origin is that near the origin
    EXPECT_EQ(fromCoordinate(small), fromCoordinate(coordinateDifference(coordinateSum(big, small), big)));
    EXPECT_EQ(big, coordinateSum(coordinateDifference(big, small), small));
    EXPECT_EQ(toCoordinate(-120.0), coordinateDifference(toCoordinate(-60.0), toCoordinate(60.0)));
}

} // namespace

} // namespace alexandria

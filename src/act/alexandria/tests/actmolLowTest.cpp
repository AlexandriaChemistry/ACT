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
 * Tests for low-level molecule utility functions.
 * \author David van der Spoel <david.vanderspoel@icm.uu.se>
 */

#include "actpre.h"

#include <cmath>
#include <string>
#include <vector>

#include <gtest/gtest.h>

#include "act/alexandria/actmol_low.h"
#include "act/basics/msg_handler.h"
#include "gromacs/math/vec.h"

#include "testutils/testasserts.h"

namespace alexandria
{

namespace
{

class ActmolLowTest : public ::testing::Test
{
protected:
    MsgHandler msghandler_;

    ActmolLowTest()
    {
        msghandler_.setPrintLevel(ACTStatus::Warning);
    }
};

TEST_F(ActmolLowTest, IsPlanarTrue)
{
    // Four atoms in the XY plane
    rvec xi = { 0.0, 0.0, 0.0 };
    rvec xj = { 1.0, 0.0, 0.0 };
    rvec xk = { 1.0, 1.0, 0.0 };
    rvec xl = { 0.0, 1.0, 0.0 };
    EXPECT_TRUE(is_planar(xi, xj, xk, xl, 10.0));
}

TEST_F(ActmolLowTest, IsPlanarFalse)
{
    // Fourth atom out of plane
    rvec xi = { 0.0, 0.0, 0.0 };
    rvec xj = { 1.0, 0.0, 0.0 };
    rvec xk = { 1.0, 1.0, 0.0 };
    rvec xl = { 0.0, 1.0, 1.0 };
    EXPECT_FALSE(is_planar(xi, xj, xk, xl, 10.0));
}

TEST_F(ActmolLowTest, IsLinearTrue)
{
    // Three atoms in a line along X axis
    rvec xi = { 0.0, 0.0, 0.0 };
    rvec xj = { 1.0, 0.0, 0.0 };
    rvec xk = { 2.0, 0.0, 0.0 };
    // is_linear returns true when angle > th_toler OR angle < 180-th_toler
    // For collinear atoms, angle is 180. With th_toler=175, 180>175 is true
    EXPECT_TRUE(is_linear(&msghandler_, xi, xj, xk, 175.0));
}

TEST_F(ActmolLowTest, IsLinearFalse)
{
    // Three atoms at a 90 degree angle - not linear
    rvec xi = { 0.0, 0.0, 0.0 };
    rvec xj = { 1.0, 0.0, 0.0 };
    rvec xk = { 1.0, 1.0, 0.0 };
    // Angle is 90 degrees. 90 < 180-175=5 is false, 90 > 175 is false
    EXPECT_FALSE(is_linear(&msghandler_, xi, xj, xk, 175.0));
}

TEST_F(ActmolLowTest, GetDoublesSimple)
{
    std::string s = "1.5 2.3 3.7";
    auto result = getDoubles(s);
    ASSERT_EQ(result.size(), 3u);
    EXPECT_NEAR(result[0], 1.5, 1e-10);
    EXPECT_NEAR(result[1], 2.3, 1e-10);
    EXPECT_NEAR(result[2], 3.7, 1e-10);
}

TEST_F(ActmolLowTest, GetDoublesSingle)
{
    std::string s = "42.0";
    auto result = getDoubles(s);
    ASSERT_EQ(result.size(), 1u);
    EXPECT_NEAR(result[0], 42.0, 1e-10);
}

TEST_F(ActmolLowTest, GetDoublesNegative)
{
    std::string s = "-1.5 -2.3";
    auto result = getDoubles(s);
    ASSERT_EQ(result.size(), 2u);
    EXPECT_NEAR(result[0], -1.5, 1e-10);
    EXPECT_NEAR(result[1], -2.3, 1e-10);
}

TEST_F(ActmolLowTest, GetDoublesEmpty)
{
    std::string s = "";
    auto result = getDoubles(s);
    EXPECT_TRUE(result.empty());
}

TEST_F(ActmolLowTest, CalcRotmatrixAligned)
{
    // When target and reference are the same direction, should give an outer product
    rvec target = { 1.0, 0.0, 0.0 };
    rvec ref    = { 1.0, 0.0, 0.0 };
    matrix rotmat;
    calc_rotmatrix(target, ref, rotmat);
    // For aligned unit vectors: bu = ref/|ref| = (1,0,0), au = target/|target| = (1,0,0)
    // rotmat[i][j] = bu[i]*au[j]
    EXPECT_NEAR(rotmat[0][0], 1.0, 1e-10);
    EXPECT_NEAR(rotmat[0][1], 0.0, 1e-10);
    EXPECT_NEAR(rotmat[0][2], 0.0, 1e-10);
    EXPECT_NEAR(rotmat[1][0], 0.0, 1e-10);
    EXPECT_NEAR(rotmat[1][1], 0.0, 1e-10);
    EXPECT_NEAR(rotmat[2][2], 0.0, 1e-10);
}

TEST_F(ActmolLowTest, CalcRotmatrixPerpendicular)
{
    // Perpendicular vectors
    rvec target = { 1.0, 0.0, 0.0 };
    rvec ref    = { 0.0, 1.0, 0.0 };
    matrix rotmat;
    calc_rotmatrix(target, ref, rotmat);
    // bu = (0,1,0), au = (1,0,0)
    // rotmat[i][j] = bu[i]*au[j]
    EXPECT_NEAR(rotmat[0][0], 0.0, 1e-10);
    EXPECT_NEAR(rotmat[0][1], 0.0, 1e-10);
    EXPECT_NEAR(rotmat[1][0], 1.0, 1e-10);
    EXPECT_NEAR(rotmat[1][1], 0.0, 1e-10);
}

TEST_F(ActmolLowTest, PutInBoxSetsBoxDimensions)
{
    rvec x[3] = { {0.0, 0.0, 0.0}, {1.0, 2.0, 3.0}, {0.5, 1.0, 1.5} };
    matrix box;
    clear_mat(box);
    real dbox = 0.5;
    put_in_box(3, box, x, dbox);
    // box[m][m] = dbox + xmax[m] - xmin[m]
    // X: 0.5 + 1.0 - 0.0 = 1.5
    // Y: 0.5 + 2.0 - 0.0 = 2.5
    // Z: 0.5 + 3.0 - 0.0 = 3.5
    EXPECT_NEAR(box[XX][XX], 1.5, 1e-10);
    EXPECT_NEAR(box[YY][YY], 2.5, 1e-10);
    EXPECT_NEAR(box[ZZ][ZZ], 3.5, 1e-10);
}

} // namespace

} // namespace alexandria

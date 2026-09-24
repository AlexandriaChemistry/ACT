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
 * Tests for the Sensitivity class and the SensitivityAnalysis() routine.
 * \author David van der Spoel <david.vanderspoel@icm.uu.se>
 */

#include "actpre.h"

#include <cmath>
#include <string>

#include <gtest/gtest.h>

#include "act/alexandria/loss_function.h"

namespace alexandria
{

namespace
{

// ============================================================================
// LossFunction -- test all variants.
// ============================================================================

TEST(LossFunctionTest, MSE)
{
    auto lf = alexandria::LossFunction::MSE;
    EXPECT_TRUE(loss(lf, 1,  2) == 4.0);
    EXPECT_TRUE(loss(lf, 1, -2) == 4.0);
    EXPECT_TRUE(loss(lf, 2,  4) == 16.0);
    EXPECT_TRUE(loss(lf, 0, -4) == 16.0);
}

TEST(LossFunctionTest, MAE)
{
    auto lf = alexandria::LossFunction::MAE;
    EXPECT_TRUE(loss(lf, 1,  2) == 2.0);
    EXPECT_TRUE(loss(lf, 1, -2) == 2.0);
    EXPECT_TRUE(loss(lf, 2,  4) == 4.0);
    EXPECT_TRUE(loss(lf, 0, -4) == 4.0);
}

TEST(LossFunctionTest, Huber)
{
    auto lf = alexandria::LossFunction::Huber;
    EXPECT_TRUE(loss(lf, 1,  2) == std::sqrt(5.0)-1.0);
    EXPECT_TRUE(loss(lf, 1, -2) == std::sqrt(5.0)-1.0);
    EXPECT_TRUE(loss(lf, 2,  4) == 4.0*(std::sqrt(5.0)-1.0));
    EXPECT_TRUE(loss(lf, 4, -4) == 16.0*(std::sqrt(2.0)-1.0));
}

TEST(LossFunctionTest, Asinh)
{
    double toler = 1e-4;
    auto lf = alexandria::LossFunction::Asinh;
    EXPECT_TRUE(std::abs(loss(lf, 1,  2) - 1.651203) < toler);
    EXPECT_TRUE(std::abs(loss(lf, 1, -2) - 1.651203) < toler);
    EXPECT_TRUE(std::abs(loss(lf, 2,  4) - 6.604812) < toler);
    EXPECT_TRUE(std::abs(loss(lf, 4, -4) - 7.474560) < toler);
}

} // namespace

} // namespace alexandria

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
 * Tests for charge symmetrization functions.
 * \author David van der Spoel <david.vanderspoel@icm.uu.se>
 */

#include "actpre.h"

#include <vector>

#include <gtest/gtest.h>

#include "act/alexandria/symmetrize_charges.h"
#include "gromacs/utility/exceptions.h"

#include "testutils/testasserts.h"

namespace alexandria
{

namespace
{

TEST(SymmetrizeChargesTest, ApplySymmetrizedChargesNoSymmetry)
{
    // Each atom in its own group (no symmetry)
    std::vector<double> q = { 0.1, 0.2, 0.3 };
    std::vector<int> sym  = { 0, 1, 2 };
    apply_symmetrized_charges(&q, sym);
    EXPECT_NEAR(q[0], 0.1, 1e-10);
    EXPECT_NEAR(q[1], 0.2, 1e-10);
    EXPECT_NEAR(q[2], 0.3, 1e-10);
}

TEST(SymmetrizeChargesTest, ApplySymmetrizedChargesAllSame)
{
    // All atoms in same group
    std::vector<double> q = { 0.1, 0.2, 0.3 };
    std::vector<int> sym  = { 0, 0, 0 };
    apply_symmetrized_charges(&q, sym);
    double expected = (0.1 + 0.2 + 0.3) / 3.0;
    EXPECT_NEAR(q[0], expected, 1e-10);
    EXPECT_NEAR(q[1], expected, 1e-10);
    EXPECT_NEAR(q[2], expected, 1e-10);
}

TEST(SymmetrizeChargesTest, ApplySymmetrizedChargesPartialGroups)
{
    // Water-like: O by itself, two H's symmetrized
    std::vector<double> q = { -0.8, 0.35, 0.45 };
    std::vector<int> sym  = { 0, 1, 1 };
    apply_symmetrized_charges(&q, sym);
    EXPECT_NEAR(q[0], -0.8, 1e-10);
    double hAvg = (0.35 + 0.45) / 2.0;
    EXPECT_NEAR(q[1], hAvg, 1e-10);
    EXPECT_NEAR(q[2], hAvg, 1e-10);
}

TEST(SymmetrizeChargesTest, ApplySymmetrizedChargesMethane)
{
    // Methane-like: C and 4 H's
    std::vector<double> q = { -0.4, 0.08, 0.12, 0.10, 0.10 };
    std::vector<int> sym  = { 0, 1, 1, 1, 1 };
    apply_symmetrized_charges(&q, sym);
    EXPECT_NEAR(q[0], -0.4, 1e-10);
    double hAvg = (0.08 + 0.12 + 0.10 + 0.10) / 4.0;
    EXPECT_NEAR(q[1], hAvg, 1e-10);
    EXPECT_NEAR(q[2], hAvg, 1e-10);
    EXPECT_NEAR(q[3], hAvg, 1e-10);
    EXPECT_NEAR(q[4], hAvg, 1e-10);
}

TEST(SymmetrizeChargesTest, ApplySymmetrizedChargesMultipleGroups)
{
    // Two separate symmetric groups
    std::vector<double> q = { 0.1, 0.2, 0.5, 0.3, 0.4 };
    std::vector<int> sym  = { 0, 0, 2, 3, 3 };
    apply_symmetrized_charges(&q, sym);
    double avg01 = (0.1 + 0.2) / 2.0;
    EXPECT_NEAR(q[0], avg01, 1e-10);
    EXPECT_NEAR(q[1], avg01, 1e-10);
    EXPECT_NEAR(q[2], 0.5, 1e-10);
    double avg34 = (0.3 + 0.4) / 2.0;
    EXPECT_NEAR(q[3], avg34, 1e-10);
    EXPECT_NEAR(q[4], avg34, 1e-10);
}

TEST(SymmetrizeChargesTest, ApplySymmetrizedChargesSizeMismatch)
{
    // Size mismatch should throw
    std::vector<double> q = { 0.1, 0.2, 0.3 };
    std::vector<int> sym  = { 0, 1 };
    EXPECT_THROW(apply_symmetrized_charges(&q, sym), gmx::InternalError);
}

TEST(SymmetrizeChargesTest, ApplySymmetrizedChargesPreservesTotal)
{
    // Total charge should be preserved after symmetrization
    std::vector<double> q = { -0.82, 0.39, 0.43 };
    std::vector<int> sym  = { 0, 1, 1 };
    double totalBefore = q[0] + q[1] + q[2];
    apply_symmetrized_charges(&q, sym);
    double totalAfter = q[0] + q[1] + q[2];
    EXPECT_NEAR(totalAfter, totalBefore, 1e-10);
}

} // namespace

} // namespace alexandria

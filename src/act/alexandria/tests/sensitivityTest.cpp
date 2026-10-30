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

#include <map>
#include <string>

#include <gtest/gtest.h>

#include "act/alexandria/sensitivity.h"
#include "act/alexandria/acmfitnesscomputer.h"
#include "act/alexandria/staticindividualinfo.h"
#include "act/forcefield/forcefield_utils.h"
#include "act/ga/genome.h"
#include "act/basics/msg_handler.h"
#include "act/utility/communicationrecord.h"
#include "act/utility/jsontree.h"

#include "testutils/testfilemanager.h"

namespace alexandria
{

namespace
{

// ============================================================================
// Sensitivity -- the parabola fit that decides "is it a minimum?"
// ============================================================================

/*! \brief chi2 = p^2 - 2p + 1 must fit to a >= 0 (upward parabola => minimum).
 */
TEST(SensitivityTest, FitsUpwardParabola)
{
    Sensitivity s;
    //add 3 points; 
    s.add(0, 1);
    s.add(1, 0);
    s.add(2, 1);
    s.computeForceConstants(nullptr);
    // Expected result is 1.0
    EXPECT_TRUE(std::abs(s.a() - 1.0) < 1e-7);
}

/*! \brief chi2 = -p^2 + 2p - 1 must fit to a < 0 (downward parabola => not min).
 */
TEST(SensitivityTest, FitsDownwardParabola)
{
    Sensitivity s;
    //add 3 points; 
    s.add(0, 1);
    s.add(1, 2);
    s.add(2, 1);
    s.computeForceConstants(nullptr);
    // Expected results is -1.0
    EXPECT_TRUE(std::abs(s.a() + 1.0) < 1e-7);
}

/*! \brief With fewer than 3 points the coefficients stay at their zero initial
 * values.
 */
TEST(SensitivityTest, InsufficientPointsLeavesConstantsZero)
{
    Sensitivity s;
    //add 2 points; 
    s.add(0, 1);
    s.add(1, 2);
    // TODO: Sensitivity s; add 2 points; 
    s.computeForceConstants(nullptr);
    // TODO: 
    EXPECT_TRUE(s.a()==0);
    EXPECT_TRUE(s.b()==0);
    EXPECT_TRUE(s.c()==0);
}

// Test: The minimum of the parabola is in the output, also without a text writer
TEST(SensitivityTest, PrintsMinimumWithoutTextWriter)
{
    // chi2 = 3 (p - 2)^2 + 5 = 3 p^2 - 12 p + 17
    Sensitivity s;
    for(double p : { 1.5, 2.25, 3.0 })
    {
        s.add(p, 3*p*p - 12*p + 17);
    }
    s.computeForceConstants(nullptr);
    JsonTree tree("sensitivity");
    s.print(nullptr, &tree, "0", "c3_z chi");
    ASSERT_EQ(1u, tree.objects().size());
    std::map<std::string, std::string> values;
    for(const auto &obj : tree.objects()[0].objects())
    {
        values[obj.key()] = obj.value();
    }
    EXPECT_EQ("c3_z", values["particle"]);
    EXPECT_EQ("chi", values["parameter"]);
    EXPECT_NEAR(2.25, std::stod(values["p_orig"]), 1e-6);
    EXPECT_NEAR(2, std::stod(values["p_min"]), 1e-6);
    EXPECT_NEAR(5, std::stod(values["chi2_min"]), 1e-6);
}

} // namespace

} // namespace alexandria

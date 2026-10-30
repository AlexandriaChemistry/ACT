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

/*! \brief A parabola through three points around a parameter
 * \param[in] a     Curvature
 * \param[in] pmin  Minimum, or maximum if a < 0
 * \param[in] p     The parameter, the second point
 * \param[in] delta Distance between the points
 * \return the fitted parabola
 */
static Sensitivity parabola(double a, double pmin, double p, double delta)
{
    Sensitivity s;
    for(double x : { p - delta, p, p + delta })
    {
        s.add(x, a*(x - pmin)*(x - pmin) + 1, x == p);
    }
    s.computeForceConstants(nullptr);
    return s;
}

// Test: A parameter at the minimum of an upward parabola is converged
TEST(SensitivityTest, ConvergedAtMinimum)
{
    auto s = parabola(2, 5, 5, 0.05);
    EXPECT_TRUE(s.checkConvergence(0, 10, 0.5));
}

// Test: Within 5% of the range of the minimum is converged, beyond it not
TEST(SensitivityTest, ConvergedWithinTolerance)
{
    EXPECT_TRUE(parabola(2, 5, 5.4, 0.05).checkConvergence(0, 10, 0.5));
    EXPECT_TRUE(parabola(2, 5, 4.6, 0.05).checkConvergence(0, 10, 0.5));
    EXPECT_FALSE(parabola(2, 5, 5.6, 0.05).checkConvergence(0, 10, 0.5));
    EXPECT_FALSE(parabola(2, 5, 4.4, 0.05).checkConvergence(0, 10, 0.5));
}

// Test: At the maximum of a downward parabola or on a straight line the
// parameter is not converged
TEST(SensitivityTest, NotConvergedWithoutUpwardCurvature)
{
    EXPECT_FALSE(parabola(-2, 5, 5, 0.05).checkConvergence(0, 10, 0.5));
    EXPECT_FALSE(parabola(0, 5, 5, 0.05).checkConvergence(0, 10, 0.5));
}

// Test: A parameter at a bound is converged if the minimum is beyond the bound
TEST(SensitivityTest, ConvergedAtBound)
{
    // The minimum is far below the lower bound of 4
    Sensitivity s;
    for(double x : { 4.0, 4.05, 4.1 })
    {
        s.add(x, 2*(x - 1)*(x - 1), x == 4.0);
    }
    s.computeForceConstants(nullptr);
    EXPECT_TRUE(s.checkConvergence(4, 14, 0.5));
    // Not so with a lower bound that leaves room for the minimum
    EXPECT_FALSE(s.checkConvergence(0, 10, 0.5));
}

// Test: Whether the parameter is converged is in the output
TEST(SensitivityTest, PrintsConvergence)
{
    for(bool conv : { true, false })
    {
        auto     s = parabola(2, 5, conv ? 5.1 : 6, 0.05);
        EXPECT_EQ(conv, s.checkConvergence(0, 10, 0.5));
        JsonTree tree("sensitivity");
        s.print(nullptr, &tree, "0", "c3_z chi");
        ASSERT_EQ(1u, tree.objects().size());
        std::string value;
        for(const auto &obj : tree.objects()[0].objects())
        {
            if (obj.key() == "converged")
            {
                value = obj.value();
            }
        }
        EXPECT_EQ(conv ? "yes" : "no", value);
    }
}

} // namespace

} // namespace alexandria

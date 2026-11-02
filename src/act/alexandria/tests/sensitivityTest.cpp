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
#include <nlohmann/json.hpp>

#include "act/alexandria/sensitivity.h"
#include "act/alexandria/acmfitnesscomputer.h"
#include "act/alexandria/staticindividualinfo.h"
#include "act/forcefield/forcefield_utils.h"
#include "act/ga/genome.h"
#include "act/basics/msg_handler.h"
#include "act/utility/communicationrecord.h"

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

// Test: The parabola is accurate also when the points are close together
// compared to the parameter value, as for a parameter with narrow bounds.
// A fit in the parameter itself loses the curvature completely for the
// larger ratios.
TEST(SensitivityTest, FitIsAccurateForNarrowBounds)
{
    // Parameter value and range, the points are half a percent of the
    // range apart, as in the sensitivity analysis
    struct Case { double p0, range; };
    for (const auto &cs : { Case{ 109, 20 }, Case{ 0.1, 0.02 }, Case{ 3000, 10 },
                            Case{ -3000, 10 }, Case{ 1.5e5, 1e4 }, Case{ 1e6, 200 } })
    {
        double delta = cs.range/200;
        double pmin  = cs.p0 + 0.3*delta;
        // The chi2 changes by 1e-4 over one step
        double a     = 1e-4/(delta*delta);
        Sensitivity s;
        for (double x : { cs.p0 - delta, cs.p0, cs.p0 + delta })
        {
            s.add(x, a*(x - pmin)*(x - pmin) + 1, x == cs.p0);
        }
        s.computeForceConstants(nullptr);
        EXPECT_NEAR(a, s.a(), 1e-8*a) << "p0 " << cs.p0;
        EXPECT_NEAR(pmin, s.pMin(), 1e-6*delta) << "p0 " << cs.p0;
        EXPECT_NEAR(1, s.chi2Min(), 1e-10) << "p0 " << cs.p0;
        // The coefficients describe the same parabola
        double chi2 = s.a()*cs.p0*cs.p0 + s.b()*cs.p0 + s.c();
        EXPECT_NEAR(a*(cs.p0 - pmin)*(cs.p0 - pmin) + 1, chi2, 1e-6) << "p0 " << cs.p0;
        // Within 5% of the range of the minimum
        double lower = cs.p0 - cs.range/2;
        EXPECT_TRUE(s.checkConvergence(lower, lower + cs.range, 0.05*cs.range)) << "p0 " << cs.p0;
    }
}

// Test: Points that all have the same parameter value cannot be fitted
TEST(SensitivityTest, IdenticalPointsLeaveConstantsZero)
{
    Sensitivity s;
    for (double chi2 : { 1.0, 2.0, 3.0 })
    {
        s.add(5, chi2);
    }
    s.computeForceConstants(nullptr);
    EXPECT_EQ(0, s.a());
    EXPECT_EQ(0, s.b());
    EXPECT_EQ(0, s.c());
    EXPECT_EQ(0, s.pMin());
    EXPECT_FALSE(s.checkConvergence(0, 10, 0.5));
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
    auto js = s.toJson(0, "c3_z chi");
    EXPECT_EQ(0, js["index"].get<int>());
    EXPECT_EQ("c3_z chi", js["name"].get<std::string>());
    EXPECT_EQ("c3_z", js["particle"].get<std::string>());
    EXPECT_EQ("chi", js["parameter"].get<std::string>());
    EXPECT_NEAR(2.25, js["p_orig"].get<double>(), 1e-12);
    EXPECT_NEAR(2, js["p_min"].get<double>(), 1e-12);
    EXPECT_NEAR(5, js["chi2_min"].get<double>(), 1e-12);
    EXPECT_NEAR(3, js["a"].get<double>(), 1e-12);
    EXPECT_NEAR(-12, js["b"].get<double>(), 1e-12);
    EXPECT_NEAR(17, js["c"].get<double>(), 1e-12);
    EXPECT_TRUE(js["fit_ok"].get<bool>());
    // The points are in the output
    ASSERT_EQ(3u, js["points"].size());
    EXPECT_EQ(1.5, js["points"][0]["p"].get<double>());
    EXPECT_EQ(3*1.5*1.5 - 12*1.5 + 17, js["points"][0]["chi2"].get<double>());
    // Without checkConvergence there are no bounds and no convergence
    EXPECT_FALSE(js.contains("lower"));
    EXPECT_FALSE(js.contains("converged"));
}

// Test: Small numbers are written with full precision, not as zero
TEST(SensitivityTest, JsonKeepsSmallNumbers)
{
    // The curvature of a C12 like parameter
    Sensitivity s;
    for(double p : { 4.95e4, 5e4, 5.05e4 })
    {
        s.add(p, 4e-12*(p - 5.01e4)*(p - 5.01e4) + 0.482, p == 5e4);
    }
    s.computeForceConstants(nullptr);
    auto js = s.toJson(0, "c12 epsilon");
    EXPECT_NEAR(4e-12, js["a"].get<double>(), 1e-20);
    // Through a string, as in the file
    auto again = nlohmann::ordered_json::parse(js.dump());
    EXPECT_EQ(js["a"].get<double>(), again["a"].get<double>());
}

// Test: The minimum of a straight line is null, not a number
TEST(SensitivityTest, JsonStraightLineHasNoMinimum)
{
    Sensitivity s;
    for(double p : { 1.0, 2.0, 3.0 })
    {
        s.add(p, 2*p + 1, p == 2.0);
    }
    s.computeForceConstants(nullptr);
    EXPECT_FALSE(s.checkConvergence(0, 10, 0.5));
    auto js = s.toJson(3, "h chi");
    EXPECT_TRUE(js["fit_ok"].get<bool>());
    EXPECT_TRUE(js["p_min"].is_null());
    EXPECT_TRUE(js["chi2_min"].is_null());
    EXPECT_TRUE(js["p_min_bounded"].is_null());
    EXPECT_FALSE(js["converged"].get<bool>());
    EXPECT_EQ(std::string::npos, js.dump().find("nan"));
}

// Test: Without enough points there is no fit
TEST(SensitivityTest, JsonWithoutFit)
{
    Sensitivity s;
    s.add(1, 1);
    s.add(2, 2, true);
    s.computeForceConstants(nullptr);
    auto js = s.toJson(0, "h chi");
    EXPECT_FALSE(js["fit_ok"].get<bool>());
    EXPECT_FALSE(js.contains("a"));
    EXPECT_TRUE(js["p_min"].is_null());
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
        auto js = s.toJson(0, "c3_z chi");
        EXPECT_EQ(conv, js["converged"].get<bool>());
        EXPECT_EQ(0, js["lower"].get<double>());
        EXPECT_EQ(10, js["upper"].get<double>());
        // The minimum is within the bounds
        EXPECT_NEAR(5, js["p_min_bounded"].get<double>(), 1e-10);
    }
}

} // namespace

} // namespace alexandria

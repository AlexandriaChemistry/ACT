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

#include <array>
#include <cmath>
#include <functional>
#include <limits>
#include <map>
#include <string>
#include <vector>

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

// Test: A parameter that does not change the chi2 is flat, with neither
// curvature nor slope, also with differences at the level of rounding
// errors. A straight line and a parabola are not flat.
TEST(SensitivityTest, FlatParameter)
{
    double chi2 = 564174.9941317956;
    double eps  = std::numeric_limits<double>::epsilon();
    for (const auto &values : { std::vector<double>{ chi2, chi2, chi2 },
                                std::vector<double>{ chi2*(1 + eps), chi2, chi2*(1 - 2*eps) } })
    {
        Sensitivity s;
        std::vector<double> p = { 0, 0.01, 0.02 };
        for (size_t k = 0; k < 3; k++)
        {
            s.add(p[k], values[k], k == 0);
        }
        s.computeForceConstants(nullptr);
        EXPECT_TRUE(s.flat());
        EXPECT_EQ(0, s.a());
        EXPECT_EQ(0, s.b());
        EXPECT_FALSE(s.checkConvergence(0, 2, 0.1));
        auto js = s.toJson(0, "f sigma");
        EXPECT_TRUE(js["flat"].get<bool>());
        EXPECT_TRUE(js["p_min"].is_null());
    }
    Sensitivity line, parabola;
    for (double p : { 1.0, 2.0, 3.0 })
    {
        line.add(p, 2*p + 1, p == 2.0);
        parabola.add(p, (p - 2)*(p - 2) + 1, p == 2.0);
    }
    line.computeForceConstants(nullptr);
    parabola.computeForceConstants(nullptr);
    EXPECT_FALSE(line.flat());
    EXPECT_NE(0, line.b());
    EXPECT_FALSE(parabola.flat());
    EXPECT_FALSE(parabola.toJson(0, "h chi")["flat"].get<bool>());
}

// Test: A pair with a flat parameter is degenerate
TEST(SensitivityTest, FlatParameterMakesDegeneratePair)
{
    // chi2 = 5 + y^2 does not depend on x, which is at its lower bound
    const auto  g = [](double y) { return 5 + y*y; };
    Sensitivity si, sj;
    for (double x : { 0.0, 0.01, 0.02 })
    {
        si.add(x, g(0), x == 0);
    }
    for (double y : { -0.1, 0.0, 0.1 })
    {
        sj.add(y, g(y), y == 0);
    }
    si.computeForceConstants(nullptr);
    sj.computeForceConstants(nullptr);
    ASSERT_TRUE(si.flat());
    Sensitivity2D s2(si, sj);
    for (double x : { 0.01, 0.02 })
    {
        for (double y : { -0.1, 0.1 })
        {
            s2.addCorner(x, y, g(y));
        }
    }
    s2.computeCrossTerm();
    EXPECT_EQ(PairShape::Degenerate, s2.shape());
    EXPECT_FALSE(s2.checkConvergence({ 0, 2, -1, 1 }, { 0.1, 0.1 }));
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

// ============================================================================
// Sensitivity2D -- the chi2 surface of a pair of parameters
// ============================================================================

//! A chi-squared surface of two parameters
using Surface = std::function<double(double, double)>;

/*! \brief The three grid values of a parameter, as in the sensitivity analysis
 * \param[in] p     The parameter
 * \param[in] delta The step
 * \param[in] side  0 for p - delta, p, p + delta, 1 for p, p + delta, p + 2 delta
 *                  (at a lower bound) and -1 for p - 2 delta, p - delta, p
 *                  (at an upper bound)
 * \return the values in increasing order
 */
static std::vector<double> gridValues(double p, double delta, int side)
{
    if (side > 0)
    {
        return { p, p + delta, p + 2*delta };
    }
    else if (side < 0)
    {
        return { p - 2*delta, p - delta, p };
    }
    return { p - delta, p, p + delta };
}

/*! \brief The analysis of a pair of parameters on a surface
 * The one-dimensional parabolas go through the grid values of each
 * parameter at the original value of the other one, the corners are
 * the grid points where both differ from their original value.
 * \param[in] f      The surface
 * \param[in] pi     First parameter
 * \param[in] pj     Second parameter
 * \param[in] deltai Step of the first parameter
 * \param[in] deltaj Step of the second parameter
 * \param[in] sidei  Side of the grid of the first parameter, see gridValues
 * \param[in] sidej  Same for the second parameter
 * \return the analysis with the coupling computed
 */
static Sensitivity2D pairOnSurface(const Surface &f, double pi, double pj,
                                   double deltai, double deltaj,
                                   int sidei = 0, int sidej = 0)
{
    const auto  xs = gridValues(pi, deltai, sidei);
    const auto  ys = gridValues(pj, deltaj, sidej);
    Sensitivity si, sj;
    for (double x : xs)
    {
        si.add(x, f(x, pj), x == pi);
    }
    for (double y : ys)
    {
        sj.add(y, f(pi, y), y == pj);
    }
    si.computeForceConstants(nullptr);
    sj.computeForceConstants(nullptr);
    Sensitivity2D s2(si, sj);
    for (double x : xs)
    {
        for (double y : ys)
        {
            if (x != pi && y != pj)
            {
                s2.addCorner(x, y, f(x, y));
            }
        }
    }
    s2.computeCrossTerm();
    return s2;
}

//! A quadratic surface with a coupling term
static double quadratic(double x, double y)
{
    return 3*x*x + 2*y*y - 1.5*x*y + 0.7*x - 0.2*y + 4;
}

// Test: The coupling of a quadratic surface is exact and the corners are on it
TEST(Sensitivity2DTest, CrossTermExactForQuadratic)
{
    const auto s2 = pairOnSurface(quadratic, 0.3, -0.4, 0.05, 0.02);
    EXPECT_NEAR(-1.5, s2.aij(), 1e-8);
    EXPECT_NEAR(0, s2.cornerRms(), 1e-12);
    auto H = s2.hessian();
    EXPECT_NEAR(6, H[0][0], 1e-8);
    EXPECT_NEAR(4, H[1][1], 1e-8);
    EXPECT_EQ(H[0][1], H[1][0]);
    // The gradient at (0.3, -0.4)
    auto g = s2.gradient();
    EXPECT_NEAR(6*0.3 - 1.5*(-0.4) + 0.7, g[0], 1e-8);
    EXPECT_NEAR(4*(-0.4) - 1.5*0.3 - 0.2, g[1], 1e-8);
    // The surface is the quadratic one, also away from the grid
    EXPECT_NEAR(quadratic(1, 2), s2.value(1, 2), 1e-8);
    EXPECT_EQ(PairShape::Minimum, s2.shape());
}

// Test: The same at bounds, where the grid is on one side of the parameters
TEST(Sensitivity2DTest, CrossTermOneSidedGrid)
{
    for (auto sides : { std::array<int, 2>{ 1, -1 }, std::array<int, 2>{ -1, 1 },
                        std::array<int, 2>{ 1, 0 }, std::array<int, 2>{ 0, -1 } })
    {
        const auto s2 = pairOnSurface(quadratic, 0.3, -0.4, 0.05, 0.02, sides[0], sides[1]);
        EXPECT_NEAR(-1.5, s2.aij(), 1e-8) << sides[0] << " " << sides[1];
        EXPECT_NEAR(0, s2.cornerRms(), 1e-12);
    }
}

// Test: On a symmetric grid the coupling is the mixed central difference,
// also for a surface that is not quadratic, which shows in the corners
TEST(Sensitivity2DTest, CrossTermIsCentralDifference)
{
    const Surface f = [](double x, double y) { return std::exp(x + 2*y) + x*x*y; };
    double  pi = 0.1, pj = 0.2, di = 0.05, dj = 0.03;
    const auto s2 = pairOnSurface(f, pi, pj, di, dj);
    double  cd = (f(pi + di, pj + dj) - f(pi + di, pj - dj) -
                  f(pi - di, pj + dj) + f(pi - di, pj - dj))/(4*di*dj);
    EXPECT_NEAR(cd, s2.aij(), 1e-9*std::abs(cd));
    EXPECT_GT(s2.cornerRms(), 0);
}

// Test: Both one-dimensional parabolas have their minimum at the
// parameters, but the pair is at a saddle point
TEST(Sensitivity2DTest, SaddleNotConverged)
{
    const Surface f  = [](double x, double y) { return x*x + y*y + 3*x*y + 1; };
    auto    s2 = pairOnSurface(f, 0, 0, 0.05, 0.05);
    EXPECT_EQ(PairShape::Saddle, s2.shape());
    EXPECT_FALSE(s2.checkConvergence({ -5, 5, -5, 5 }, { 0.5, 0.5 }));
    EXPECT_FALSE(s2.boxMinimum({ -5, 5, -5, 5 }).has_value());
    // The rho is beyond one
    ASSERT_TRUE(s2.rho().has_value());
    EXPECT_NEAR(1.5, *s2.rho(), 1e-8);
}

/*! \brief A valley along x = y with its minimum at (1, 1)
 * Along each parameter on its own the minimum is close to zero.
 */
static double valley(double x, double y)
{
    return 100*(x - y)*(x - y) + (x + y - 2)*(x + y - 2) + 1;
}

// Test: In a correlated valley each parameter on its own is within the
// tolerance of its minimum, but the pair is not within the tolerance of
// the joint minimum
TEST(Sensitivity2DTest, CorrelatedValleyNotConverged)
{
    // Range 10, tolerance 5% of the range
    std::array<double, 4> bounds = { -5, 5, -5, 5 };
    double                tol    = 0.5;
    auto s2 = pairOnSurface(valley, 0, 0, 0.05, 0.05);
    // The one-dimensional minima are at 2/101
    Sensitivity si;
    for (double x : { -0.05, 0.0, 0.05 })
    {
        si.add(x, valley(x, 0), x == 0);
    }
    si.computeForceConstants(nullptr);
    EXPECT_NEAR(2.0/101, si.pMin(), 1e-8);
    EXPECT_TRUE(si.checkConvergence(bounds[0], bounds[1], tol));
    // The pair is a minimum at (1, 1) with a strong negative coupling
    EXPECT_EQ(PairShape::Minimum, s2.shape());
    ASSERT_TRUE(s2.rho().has_value());
    EXPECT_NEAR(-198.0/202, *s2.rho(), 1e-8);
    auto sp = s2.stationaryPoint();
    ASSERT_TRUE(sp.has_value());
    EXPECT_NEAR(1, (*sp)[0], 1e-6);
    EXPECT_NEAR(1, (*sp)[1], 1e-6);
    EXPECT_FALSE(s2.checkConvergence(bounds, { tol, tol }));
}

// Test: At the joint minimum, and within the tolerance of it on each
// axis, the pair is converged; just outside it is not
TEST(Sensitivity2DTest, ConvergedWithinTolerance)
{
    std::array<double, 4> bounds = { -5, 5, -5, 5 };
    EXPECT_TRUE(pairOnSurface(valley, 1, 1, 0.05, 0.05).checkConvergence(bounds, { 0.5, 0.5 }));
    // Along the valley
    EXPECT_TRUE(pairOnSurface(valley, 0.55, 0.55, 0.05, 0.05).checkConvergence(bounds, { 0.5, 0.5 }));
    EXPECT_FALSE(pairOnSurface(valley, 0.45, 0.45, 0.05, 0.05).checkConvergence(bounds, { 0.5, 0.5 }));
    // One parameter close, the other not
    EXPECT_FALSE(pairOnSurface(valley, 1, 0.45, 0.05, 0.05).checkConvergence(bounds, { 0.5, 0.5 }));
    // A tolerance per parameter
    EXPECT_TRUE(pairOnSurface(valley, 1, 0.45, 0.05, 0.05).checkConvergence(bounds, { 0.5, 0.6 }));
}

// Test: With the minimum beyond one bound, the minimum within the bounds
// is on that edge, where the surface is a parabola in the other parameter
TEST(Sensitivity2DTest, BoxMinimumOnEdge)
{
    // Minimum at (3, 0.5), the upper bound of x is 2. Along x = 2 the
    // minimum is at y = 0.5 + 0.25.
    const Surface f = [](double x, double y)
    { return (x - 3)*(x - 3) + (y - 0.5)*(y - 0.5) + 0.5*(x - 3)*(y - 0.5) + 1; };
    std::array<double, 4> bounds = { 0, 2, 0, 1 };
    // The parameter x is at its upper bound, the grid is below it
    auto s2 = pairOnSurface(f, 2, 0.75, 0.01, 0.005, -1, 0);
    auto bm = s2.boxMinimum(bounds);
    ASSERT_TRUE(bm.has_value());
    EXPECT_NEAR(2, (*bm)[0], 1e-8);
    EXPECT_NEAR(0.75, (*bm)[1], 1e-8);
    EXPECT_TRUE(s2.checkConvergence(bounds, { 0.1, 0.05 }));
    // Away from the edge minimum along the edge
    auto s2b = pairOnSurface(f, 2, 0.6, 0.01, 0.005, -1, 0);
    EXPECT_FALSE(s2b.checkConvergence(bounds, { 0.1, 0.05 }));
    // Without the bound the minimum is far away
    EXPECT_FALSE(s2.checkConvergence({ 0, 5, 0, 1 }, { 0.1, 0.05 }));
}

// Test: With the minimum beyond two bounds it is in the corner
TEST(Sensitivity2DTest, BoxMinimumAtCorner)
{
    const Surface f = [](double x, double y)
    { return (x - 3)*(x - 3) + (y - 2)*(y - 2) + 0.3*(x - 3)*(y - 2) + 1; };
    std::array<double, 4> bounds = { 0, 2, 0, 1 };
    auto s2 = pairOnSurface(f, 2, 1, 0.01, 0.005, -1, -1);
    auto bm = s2.boxMinimum(bounds);
    ASSERT_TRUE(bm.has_value());
    EXPECT_NEAR(2, (*bm)[0], 1e-8);
    EXPECT_NEAR(1, (*bm)[1], 1e-8);
    EXPECT_TRUE(s2.checkConvergence(bounds, { 0.1, 0.05 }));
}

// Test: A maximum, and surfaces without curvature in a direction
TEST(Sensitivity2DTest, MaximumAndDegenerate)
{
    std::array<double, 4> bounds = { -5, 5, -5, 5 };
    auto maximum = pairOnSurface([](double x, double y) { return 10 - x*x - y*y - 0.5*x*y; },
                                 0, 0, 0.05, 0.05);
    EXPECT_EQ(PairShape::Maximum, maximum.shape());
    EXPECT_FALSE(maximum.checkConvergence(bounds, { 0.5, 0.5 }));
    EXPECT_FALSE(maximum.rho().has_value());
    // No dependence on y
    auto flat = pairOnSurface([](double x, double) { return x*x + 1; }, 0, 0, 0.05, 0.05);
    EXPECT_EQ(PairShape::Degenerate, flat.shape());
    EXPECT_FALSE(flat.checkConvergence(bounds, { 0.5, 0.5 }));
    // A valley without a minimum along x = y
    auto line = pairOnSurface([](double x, double y) { return (x - y)*(x - y) + 1; }, 0.3, 0.1, 0.05, 0.05);
    EXPECT_EQ(PairShape::Degenerate, line.shape());
    EXPECT_FALSE(line.stationaryPoint().has_value());
    EXPECT_FALSE(line.checkConvergence(bounds, { 0.5, 0.5 }));
    for (const auto &s2 : { maximum, flat, line })
    {
        auto js = s2.toJson(0, 1, "a chi", "b chi");
        EXPECT_TRUE(js["p_min"].is_null());
        EXPECT_TRUE(js["chi2_min"].is_null());
        EXPECT_TRUE(js["p_min_bounded"].is_null());
        EXPECT_FALSE(js["converged"].get<bool>());
        EXPECT_EQ(std::string::npos, js.dump().find("nan"));
    }
}

// Test: Without corners there is no coupling and no shape
TEST(Sensitivity2DTest, WithoutCorners)
{
    Sensitivity si, sj;
    for (double x : { -0.1, 0.0, 0.1 })
    {
        si.add(x, x*x, x == 0);
        sj.add(x, x*x, x == 0);
    }
    si.computeForceConstants(nullptr);
    sj.computeForceConstants(nullptr);
    Sensitivity2D s2(si, sj);
    s2.computeCrossTerm();
    EXPECT_EQ(PairShape::Degenerate, s2.shape());
    auto js = s2.toJson(0, 1, "a chi", "b chi");
    EXPECT_TRUE(js["a_ij"].is_null());
    EXPECT_TRUE(js["grid"]["chi2"][0][0].is_null());
}

// Test: The output has all fields, with the grid as it was computed
TEST(Sensitivity2DTest, Json)
{
    double pi = 0.3, pj = -0.4, di = 0.05, dj = 0.02;
    auto   s2 = pairOnSurface(quadratic, pi, pj, di, dj);
    (void) s2.checkConvergence({ -5, 5, -5, 5 }, { 0.5, 0.5 });
    auto js = s2.toJson(2, 4, "c3 chi", "h chi");
    EXPECT_EQ(2u, js["i"].get<size_t>());
    EXPECT_EQ(4u, js["j"].get<size_t>());
    EXPECT_EQ("c3 chi", js["names"][0].get<std::string>());
    EXPECT_EQ("h chi", js["names"][1].get<std::string>());
    auto xs = gridValues(pi, di, 0);
    auto ys = gridValues(pj, dj, 0);
    ASSERT_EQ(3u, js["grid"]["p_i"].size());
    ASSERT_EQ(3u, js["grid"]["p_j"].size());
    ASSERT_EQ(3u, js["grid"]["chi2"].size());
    for (size_t k = 0; k < 3; k++)
    {
        EXPECT_EQ(xs[k], js["grid"]["p_i"][k].get<double>());
        EXPECT_EQ(ys[k], js["grid"]["p_j"][k].get<double>());
        ASSERT_EQ(3u, js["grid"]["chi2"][k].size());
        for (size_t l = 0; l < 3; l++)
        {
            EXPECT_EQ(quadratic(xs[k], ys[l]), js["grid"]["chi2"][k][l].get<double>()) << k << " " << l;
        }
    }
    EXPECT_NEAR(-1.5, js["a_ij"].get<double>(), 1e-8);
    EXPECT_NEAR(-1.5/(2*std::sqrt(3.0*2.0)), js["rho"].get<double>(), 1e-8);
    EXPECT_TRUE(js["corner_rms"].is_number());
    EXPECT_NEAR(6, js["hessian"][0][0].get<double>(), 1e-8);
    EXPECT_NEAR(-1.5, js["hessian"][0][1].get<double>(), 1e-8);
    EXPECT_EQ(2u, js["gradient"].size());
    EXPECT_EQ("minimum", js["shape"].get<std::string>());
    // The minimum of the quadratic: 6x - 1.5y = -0.7, 4y - 1.5x = 0.2
    double det = 6*4 - 1.5*1.5;
    double xm  = (-0.7*4 + 1.5*0.2)/det;
    double ym  = (6*0.2 - 1.5*0.7)/det;
    EXPECT_NEAR(xm, js["p_min"][0].get<double>(), 1e-8);
    EXPECT_NEAR(ym, js["p_min"][1].get<double>(), 1e-8);
    EXPECT_NEAR(quadratic(xm, ym), js["chi2_min"].get<double>(), 1e-8);
    EXPECT_NEAR(xm, js["p_min_bounded"][0].get<double>(), 1e-8);
    EXPECT_TRUE(js["converged"].is_boolean());
    EXPECT_EQ("minimum", pairShapeName(PairShape::Minimum));
    EXPECT_EQ("saddle", pairShapeName(PairShape::Saddle));
    EXPECT_EQ("maximum", pairShapeName(PairShape::Maximum));
    EXPECT_EQ("degenerate", pairShapeName(PairShape::Degenerate));
}

// ============================================================================
// HessianDiagnostic -- the full Hessian of all parameters
// ============================================================================

// Test: A positive definite Hessian, scaled by the ranges, and the Newton
// step to its minimum relative to the ranges
TEST(HessianDiagnosticTest, PositiveDefinite)
{
    // chi2 = (x - 1)^2 + 4 (y + 1)^2 around (0, 0)
    auto hd = diagnoseHessian({ { 2, 0 }, { 0, 8 } }, { -2, 8 }, { 1, 0.5 });
    EXPECT_TRUE(hd.positiveDefinite);
    // Scaled by the ranges the Hessian is twice the unit matrix
    EXPECT_NEAR(2, hd.eigenvalueMin, 1e-12);
    EXPECT_NEAR(2, hd.eigenvalueMax, 1e-12);
    // The step is (1, -1), twice the range of y
    ASSERT_TRUE(hd.maxNewtonStepFraction.has_value());
    EXPECT_NEAR(2, *hd.maxNewtonStepFraction, 1e-12);
    auto js = hd.toJson();
    EXPECT_TRUE(js["positive_definite"].get<bool>());
    EXPECT_NEAR(2, js["max_newton_step_fraction"].get<double>(), 1e-12);
}

// Test: With a coupling the Newton step goes to the joint minimum
TEST(HessianDiagnosticTest, NewtonStepWithCoupling)
{
    // chi2 = 3x^2 + 2y^2 - 1.5xy + 0.7x - 0.2y around (0, 0)
    auto   hd  = diagnoseHessian({ { 6, -1.5 }, { -1.5, 4 } }, { 0.7, -0.2 }, { 10, 10 });
    double det = 6*4 - 1.5*1.5;
    double xm  = (-0.7*4 + 1.5*0.2)/det;
    double ym  = (6*0.2 - 1.5*0.7)/det;
    ASSERT_TRUE(hd.maxNewtonStepFraction.has_value());
    EXPECT_NEAR(std::max(std::abs(xm), std::abs(ym))/10, *hd.maxNewtonStepFraction, 1e-12);
}

// Test: An indefinite Hessian has no Newton step
TEST(HessianDiagnosticTest, Indefinite)
{
    const auto hd = diagnoseHessian({ { 2, 3 }, { 3, 2 } }, { 0, 0 }, { 1, 1 });
    EXPECT_FALSE(hd.positiveDefinite);
    EXPECT_NEAR(-1, hd.eigenvalueMin, 1e-12);
    EXPECT_NEAR(5, hd.eigenvalueMax, 1e-12);
    EXPECT_FALSE(hd.maxNewtonStepFraction.has_value());
    EXPECT_TRUE(hd.toJson()["max_newton_step_fraction"].is_null());
}

// Test: Each pair of three parameters is positive definite, all three are not
TEST(HessianDiagnosticTest, PairsPositiveDefiniteButNotAll)
{
    const std::vector<std::vector<double>> H = { { 1, -0.6, -0.6 }, { -0.6, 1, -0.6 }, { -0.6, -0.6, 1 } };
    // Each 2x2 block has determinant 1 - 0.36 > 0
    const auto hd = diagnoseHessian(H, { 0, 0, 0 }, { 1, 1, 1 });
    EXPECT_FALSE(hd.positiveDefinite);
    EXPECT_NEAR(1 - 2*0.6, hd.eigenvalueMin, 1e-12);
    EXPECT_NEAR(1 + 0.6, hd.eigenvalueMax, 1e-12);
}

// Test: A singular Hessian is not positive definite
TEST(HessianDiagnosticTest, Singular)
{
    const auto hd = diagnoseHessian({ { 2, -2 }, { -2, 2 } }, { 0, 0 }, { 1, 1 });
    EXPECT_FALSE(hd.positiveDefinite);
}

} // namespace

} // namespace alexandria

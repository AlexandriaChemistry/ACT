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
 * Tests of the explicit SIMD math functions.
 * \author David van der Spoel <david.vanderspoel@icm.uu.se>
 */
#include "actpre.h"

#include "../simd_math.h"

#include <cmath>

#include <gtest/gtest.h>

#ifdef ACT_HAVE_SIMD

namespace alexandria
{

namespace
{

using simd::SimdDouble;

/*! \brief Largest error of a SIMD function compared to the library
 * \param[in] func     The SIMD function
 * \param[in] ref      The library function
 * \param[in] xmin     Start of the range
 * \param[in] xmax     End of the range
 * \param[in] relative Whether to return the error in units in the last place
 * \return The largest error
 */
template <typename Func, typename Ref>
double maxError(Func func, Ref ref, double xmin, double xmax, bool relative)
{
    const int n      = 200000;
    double    maxErr = 0;
    for (int i = 0; i < n; i += simd::simdWidth)
    {
        double x[simd::simdWidth], y[simd::simdWidth];
        for (size_t k = 0; k < simd::simdWidth; k++)
        {
            x[k] = xmin + (xmax-xmin)*(i+k)/n;
        }
        simd::store(y, func(simd::load(x)));
        for (size_t k = 0; k < simd::simdWidth; k++)
        {
            double yref = ref(x[k]);
            double err  = std::abs(y[k] - yref);
            if (relative)
            {
                err /= std::nextafter(std::abs(yref), INFINITY) - std::abs(yref);
            }
            maxErr = std::max(maxErr, err);
        }
    }
    return maxErr;
}

TEST(SimdMathTest, Exp)
{
    const auto func = [](SimdDouble x) { return simd::exp(x); };
    const auto ref  = [](double x) { return std::exp(x); };
    EXPECT_LE(maxError(func, ref, -708, 709.7, true), 2);
    EXPECT_LE(maxError(func, ref, -50, 5, true), 2);
}

TEST(SimdMathTest, ExpLimits)
{
    double x[simd::simdWidth], y[simd::simdWidth];
    for (size_t k = 0; k < simd::simdWidth; k++)
    {
        x[k] = (k % 2 == 0) ? -800 : 800;
    }
    simd::store(y, simd::exp(simd::load(x)));
    for (size_t k = 0; k < simd::simdWidth; k++)
    {
        if (k % 2 == 0)
        {
            EXPECT_EQ(0, y[k]);
        }
        else
        {
            EXPECT_TRUE(std::isinf(y[k]));
        }
    }
}

TEST(SimdMathTest, Erf)
{
    const auto func = [](SimdDouble x) { return simd::erf(x); };
    const auto ref  = [](double x) { return std::erf(x); };
    EXPECT_LE(maxError(func, ref, -10, 10, false), 5e-16);
    EXPECT_LE(maxError(func, ref, -1, 1, true), 8);
}

TEST(SimdMathTest, ErfAndExp)
{
    // The error function is the same as without the exponential
    const auto ferf = [](SimdDouble x) { SimdDouble e; return simd::erfAndExp(x, &e); };
    const auto rerf = [](double x) { return std::erf(x); };
    EXPECT_LE(maxError(ferf, rerf, -10, 10, false), 5e-16);
    // The exponential is that of simd::exp(-x^2) up to |x| = 8
    const auto fexp = [](SimdDouble x) { SimdDouble e; simd::erfAndExp(x, &e); return e; };
    const auto rexp = [](double x) { return std::exp(-x*x); };
    EXPECT_LE(maxError(fexp, rexp, -8, 8, true), 2);
}

TEST(SimdMathTest, Asin)
{
    const auto func = [](SimdDouble x) { return simd::asin(x); };
    const auto ref  = [](double x) { return std::asin(x); };
    EXPECT_LE(maxError(func, ref, -1, 1, true), 2);
    EXPECT_LE(maxError(func, ref, -1e-6, 1e-6, true), 2);
}

TEST(SimdMathTest, Acos)
{
    const auto func = [](SimdDouble x) { return simd::acos(x); };
    const auto ref  = [](double x) { return std::acos(x); };
    EXPECT_LE(maxError(func, ref, -1, 1, true), 2);
    EXPECT_LE(maxError(func, ref, 1-1e-6, 1, false), 1e-16);
}

TEST(SimdMathTest, Atan2)
{
    // Points on circles of different radius around the origin, in all
    // quadrants, and the axes, which are on the grid of angles as well
    for (double radius : { 1e-5, 0.3, 1.0, 70.0 })
    {
        const auto func = [radius](SimdDouble phi)
        {
            double p[simd::simdWidth], x[simd::simdWidth], y[simd::simdWidth];
            simd::store(p, phi);
            for (size_t k = 0; k < simd::simdWidth; k++)
            {
                x[k] = radius*std::cos(p[k]);
                y[k] = radius*std::sin(p[k]);
            }
            return simd::atan2(simd::load(y), simd::load(x));
        };
        const auto ref  = [radius](double phi)
        {
            return std::atan2(radius*std::sin(phi), radius*std::cos(phi));
        };
        EXPECT_LE(maxError(func, ref, -M_PI, M_PI, true), 2) << "radius " << radius;
    }
    const double x[] = { 0, 1, 0, -1, 0, 3e-200, -2, 5 };
    const double y[] = { 0, 0, 1, 0, -1, 1, 1e-300, -1e-300 };
    for (size_t i = 0; i < sizeof(x)/sizeof(x[0]); i += simd::simdWidth)
    {
        double xs[simd::simdWidth], ys[simd::simdWidth], r[simd::simdWidth];
        for (size_t k = 0; k < simd::simdWidth; k++)
        {
            xs[k] = x[(i+k) % 8];
            ys[k] = y[(i+k) % 8];
        }
        simd::store(r, simd::atan2(simd::load(ys), simd::load(xs)));
        for (size_t k = 0; k < simd::simdWidth; k++)
        {
            EXPECT_EQ(std::atan2(ys[k], xs[k]), r[k]) << "x " << xs[k] << " y " << ys[k];
        }
    }
}

using simd::SimdFloat;

/*! \brief The largest error of a single precision SIMD function
 * The reference is computed in double precision, the relative errors
 * are in units of the last place of the reference rounded to float.
 * \param[in] func     The SIMD function
 * \param[in] ref      The reference function
 * \param[in] xmin     Smallest argument
 * \param[in] xmax     Largest argument
 * \param[in] relative Whether the error is relative, in units in the
 *                     last place, or absolute
 * \return the largest error
 */
template <typename Func, typename Ref>
double maxErrorFloat(Func func, Ref ref, double xmin, double xmax, bool relative)
{
    constexpr size_t width  = simd::simdFloatWidth;
    const int        n      = 200000;
    double           maxErr = 0;
    for (int i = 0; i < n; i += width)
    {
        float x[width], y[width];
        for (size_t k = 0; k < width; k++)
        {
            x[k] = xmin + (xmax-xmin)*(i+k)/n;
        }
        simd::store(y, func(simd::load(x)));
        for (size_t k = 0; k < width; k++)
        {
            double yref = ref(static_cast<double>(x[k]));
            double err  = std::abs(y[k] - yref);
            if (relative)
            {
                float fref = std::abs(static_cast<float>(yref));
                err /= std::nextafter(fref, INFINITY) - fref;
            }
            maxErr = std::max(maxErr, err);
        }
    }
    return maxErr;
}

TEST(SimdMathFloatTest, Exp)
{
    const auto func = [](SimdFloat x) { return simd::exp(x); };
    const auto ref  = [](double x) { return std::exp(x); };
    EXPECT_LE(maxErrorFloat(func, ref, -87, 88.7, true), 2);
    EXPECT_LE(maxErrorFloat(func, ref, -50, 5, true), 2);
}

TEST(SimdMathFloatTest, ExpLimits)
{
    constexpr size_t width = simd::simdFloatWidth;
    float            x[width], y[width];
    for (size_t k = 0; k < width; k++)
    {
        x[k] = (k % 2 == 0) ? -100 : 100;
    }
    simd::store(y, simd::exp(simd::load(x)));
    for (size_t k = 0; k < width; k++)
    {
        if (k % 2 == 0)
        {
            EXPECT_EQ(0, y[k]);
        }
        else
        {
            EXPECT_TRUE(std::isinf(y[k]));
        }
    }
}

TEST(SimdMathFloatTest, Erf)
{
    const auto func = [](SimdFloat x) { return simd::erf(x); };
    const auto ref  = [](double x) { return std::erf(x); };
    EXPECT_LE(maxErrorFloat(func, ref, -10, 10, false), 2e-7);
    EXPECT_LE(maxErrorFloat(func, ref, -1, 1, true), 3);
}

TEST(SimdMathFloatTest, ErfAndExp)
{
    // The exponential of -x^2 with x^2 in single precision, as computed
    const auto fexp = [](SimdFloat x) { SimdFloat e; simd::erfAndExp(x, &e); return e; };
    const auto rexp = [](double x) { return std::exp(-static_cast<double>(static_cast<float>(x*x))); };
    EXPECT_LE(maxErrorFloat(fexp, rexp, -8, 8, true), 2);
}

TEST(SimdMathFloatTest, Asin)
{
    const auto func = [](SimdFloat x) { return simd::asin(x); };
    const auto ref  = [](double x) { return std::asin(x); };
    EXPECT_LE(maxErrorFloat(func, ref, -1, 1, true), 2);
    EXPECT_LE(maxErrorFloat(func, ref, -1e-3, 1e-3, true), 2);
}

TEST(SimdMathFloatTest, Acos)
{
    const auto func = [](SimdFloat x) { return simd::acos(x); };
    const auto ref  = [](double x) { return std::acos(x); };
    EXPECT_LE(maxErrorFloat(func, ref, -1, 1, true), 2);
    EXPECT_LE(maxErrorFloat(func, ref, 1-1e-3, 1, false), 1e-8);
}

TEST(SimdMathFloatTest, Atan2)
{
    // Points on circles of different radius around the origin, the
    // reference from the coordinates rounded to float
    for (double radius : { 1e-5, 0.3, 1.0, 70.0 })
    {
        const auto func = [radius](SimdFloat phi)
        {
            constexpr size_t width = simd::simdFloatWidth;
            float            p[width], x[width], y[width];
            simd::store(p, phi);
            for (size_t k = 0; k < width; k++)
            {
                x[k] = radius*std::cos(p[k]);
                y[k] = radius*std::sin(p[k]);
            }
            return simd::atan2(simd::load(y), simd::load(x));
        };
        const auto ref = [radius](double phi)
        {
            float x = radius*std::cos(static_cast<float>(phi));
            float y = radius*std::sin(static_cast<float>(phi));
            return std::atan2(static_cast<double>(y), static_cast<double>(x));
        };
        EXPECT_LE(maxErrorFloat(func, ref, -M_PI, M_PI, true), 2) << "radius " << radius;
    }
}

} // namespace

} // namespace alexandria

#endif

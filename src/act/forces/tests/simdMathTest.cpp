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
    auto func = [](SimdDouble x) { return simd::exp(x); };
    auto ref  = [](double x) { return std::exp(x); };
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
    auto func = [](SimdDouble x) { return simd::erf(x); };
    auto ref  = [](double x) { return std::erf(x); };
    EXPECT_LE(maxError(func, ref, -10, 10, false), 5e-16);
    EXPECT_LE(maxError(func, ref, -1, 1, true), 8);
}

} // namespace

} // namespace alexandria

#endif

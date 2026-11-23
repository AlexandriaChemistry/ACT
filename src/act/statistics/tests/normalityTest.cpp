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
 * Tests of the test for normality.
 * \author David van der Spoel <david.vanderspoel@icm.uu.se>
 */

#include "actpre.h"

#include <cmath>
#include <gtest/gtest.h>
#include <vector>

#include "act/statistics/normality.h"

namespace alexandria
{

namespace
{

/*! \brief Values that follow a normal distribution as well as n values can
 * \param[in] n       The number of values
 * \param[in] average The average
 * \param[in] sigma   The standard deviation
 * \return the values at the probabilities (i+0.5)/n
 */
std::vector<double> normalValues(size_t n, double average, double sigma)
{
    std::vector<double> values;
    for (size_t i = 0; i < n; i++)
    {
        const double p = (i + 0.5)/n;
        // The inverse of the cumulative distribution by bisection
        double lo = -10, hi = 10;
        for (int k = 0; k < 100; k++)
        {
            const double mid = 0.5*(lo + hi);
            if (0.5*std::erfc(-mid/std::sqrt(2.0)) < p)
            {
                lo = mid;
            }
            else
            {
                hi = mid;
            }
        }
        values.push_back(average + sigma*0.5*(lo + hi));
    }
    return values;
}

// Test: The statistic, the p-value and the moments are those that scipy
// (stats.anderson, with the correction for the number of values, skew and
// kurtosis) gives for the same values
TEST(NormalityTest, ReferenceValues)
{
    const auto y = testNormality({ 13, 51, 18, 21, 16, 27, 26, 40, 36, 19, 42 });
    EXPECT_EQ(11u, y.n);
    EXPECT_TRUE(y.tested);
    EXPECT_NEAR(28.0909090909, y.average, 1e-9);
    EXPECT_NEAR(12.3972137632, y.sigma, 1e-9);
    EXPECT_NEAR(0.5255680614, y.skewness, 1e-9);
    EXPECT_NEAR(-0.9980905756, y.kurtosis, 1e-9);
    EXPECT_NEAR(0.4040674461, y.statistic, 1e-9);
    EXPECT_NEAR(0.3546521, y.pvalue, 1e-6);
    EXPECT_TRUE(y.normal(0.01));
    EXPECT_TRUE(y.normal(0.05));
}

// Test: Two groups of bond lengths ten pm apart, as when one bond type
// covers two kinds of bonds, are not normal: no skewness, negative kurtosis
TEST(NormalityTest, TwoPeaks)
{
    const auto b = testNormality({ 100.1, 100.4, 99.8, 100.0, 100.3, 99.9, 100.2, 99.7,
                                   110.2, 110.0, 109.8, 110.3, 109.9, 110.1, 110.4, 109.7 });
    EXPECT_TRUE(b.tested);
    EXPECT_NEAR(2.4899656729, b.statistic, 1e-9);
    EXPECT_NEAR(2.743722e-06, b.pvalue, 1e-11);
    EXPECT_NEAR(0, b.skewness, 1e-9);
    EXPECT_NEAR(-1.9916318232, b.kurtosis, 1e-9);
    EXPECT_FALSE(b.normal(0.01));
}

// Test: A tail to large values is not normal: positive skewness and kurtosis.
// The order of the values does not matter.
TEST(NormalityTest, Tail)
{
    const std::vector<double> values = { 12.0, 1, 1.1, 7.0, 1.2, 1.25, 1.3, 1.4, 1.5, 1.7, 2.0, 2.4, 3.1, 4.5 };
    const auto                s      = testNormality(values);
    EXPECT_TRUE(s.tested);
    EXPECT_NEAR(2.0163135020, s.statistic, 1e-9);
    EXPECT_NEAR(3.939723e-05, s.pvalue, 1e-10);
    EXPECT_NEAR(2.0866676790, s.skewness, 1e-9);
    EXPECT_NEAR(3.3971220093, s.kurtosis, 1e-9);
    EXPECT_FALSE(s.normal(0.01));
}

// Test: Values from a normal distribution pass, with few and with many
// values, whatever their average and width
TEST(NormalityTest, NormalValues)
{
    for (const size_t n : { 8u, 25u, 1000u, 100000u })
    {
        const auto r = testNormality(normalValues(n, 109.5, 0.7));
        EXPECT_TRUE(r.tested) << n;
        EXPECT_GT(r.pvalue, 0.9) << n;
        EXPECT_TRUE(r.normal(0.05)) << n;
        EXPECT_NEAR(109.5, r.average, 1e-6) << n;
        EXPECT_NEAR(0, r.skewness, 1e-6) << n;
        if (n >= 1000)
        {
            EXPECT_NEAR(0.7, r.sigma, 0.01) << n;
            EXPECT_NEAR(0, r.kurtosis, 0.05) << n;
        }
    }
}

// Test: A far outlier among many normal values is found, and the p-value
// stays a probability for a statistic of any size
TEST(NormalityTest, Outlier)
{
    auto values = normalValues(500, 120, 2);
    values.push_back(190);
    const auto r = testNormality(values);
    EXPECT_TRUE(r.tested);
    EXPECT_FALSE(r.normal(0.01));
    EXPECT_GT(r.skewness, 1);
    EXPECT_GT(r.kurtosis, 10);
    // Half of the values at each of two places far apart: a very large statistic
    std::vector<double> two(2000, 100.0);
    for (size_t i = 0; i < 1000; i++)
    {
        two[i] = 150 + 1e-3*i;
    }
    const auto t = testNormality(two);
    EXPECT_TRUE(t.tested);
    EXPECT_GT(t.statistic, 100);
    EXPECT_EQ(0, t.pvalue);
    EXPECT_FALSE(t.normal(1e-6));
}

// Test: Too few values, or values that are all the same, are not tested and
// do not count as not normal
TEST(NormalityTest, NotTested)
{
    const auto none = testNormality({});
    EXPECT_EQ(0u, none.n);
    EXPECT_FALSE(none.tested);
    const auto one = testNormality({ 109.4 });
    EXPECT_FALSE(one.tested);
    EXPECT_NEAR(109.4, one.average, 1e-12);
    EXPECT_EQ(0, one.sigma);
    const auto few = testNormality({ 1, 2, 3, 4, 5, 6, 100 });
    EXPECT_EQ(normalityMinimum - 1, few.n);
    EXPECT_FALSE(few.tested);
    EXPECT_TRUE(few.normal(0.01));
    EXPECT_GT(few.skewness, 1);
    EXPECT_EQ(1, few.pvalue);
    const auto same = testNormality(std::vector<double>(50, 109.4));
    EXPECT_EQ(50u, same.n);
    EXPECT_FALSE(same.tested);
    EXPECT_TRUE(same.normal(0.01));
    EXPECT_EQ(0, same.skewness);
}

} // namespace

} // namespace alexandria

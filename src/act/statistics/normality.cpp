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
 * Implements part of the alexandria program.
 * \author David van der Spoel <david.vanderspoel@icm.uu.se>
 */

#include "actpre.h"

#include <algorithm>
#include <cmath>
#include <cstddef>
#include <limits>
#include <vector>

#include "normality.h"

namespace alexandria
{

/*! \brief The logarithm of the cumulative normal distribution
 * \param[in] z The standardized value
 * \return ln(Phi(z)), finite far out in the tails as well
 */
static double logNormalCdf(double z)
{
    // One over the square root of two
    const double invSqrt2 = 0.70710678118654752440;
    const double cdf      = 0.5*std::erfc(-z*invSqrt2);
    return std::log(std::max(cdf, std::numeric_limits<double>::min()));
}

/*! \brief The p-value of the Anderson-Darling statistic for a normal
 * distribution of which the average and the standard deviation are
 * estimated
 * \param[in] astar The statistic, corrected for the number of values
 * \return the p-value
 */
static double andersonDarlingPvalue(double astar)
{
    const double a = astar;
    double       p = 0;
    if (a >= 13)
    {
        // Beyond this value the first formula does not hold: below 1e-30
        p = 0;
    }
    else if (a >= 0.6)
    {
        p = std::exp(1.2937 - 5.709*a + 0.0186*a*a);
    }
    else if (a >= 0.34)
    {
        p = std::exp(0.9177 - 4.279*a - 1.38*a*a);
    }
    else if (a >= 0.2)
    {
        p = 1 - std::exp(-8.318 + 42.796*a - 59.938*a*a);
    }
    else
    {
        p = 1 - std::exp(-13.436 + 101.14*a - 223.73*a*a);
    }
    return std::min(1.0, std::max(0.0, p));
}

NormalityResult testNormality(std::vector<double> values)
{
    NormalityResult result;
    const size_t    n = values.size();
    result.n          = n;
    if (n == 0)
    {
        return result;
    }
    double sum = 0;
    for (const auto v : values)
    {
        sum += v;
    }
    result.average = sum/n;
    // The central moments
    double m2 = 0, m3 = 0, m4 = 0;
    for (const auto v : values)
    {
        const double d  = v - result.average;
        const double d2 = d*d;
        m2 += d2;
        m3 += d2*d;
        m4 += d2*d2;
    }
    if (n < 2)
    {
        return result;
    }
    result.sigma = std::sqrt(m2/(n-1));
    m2 /= n;
    m3 /= n;
    m4 /= n;
    // Values that are the same but for rounding have no shape
    const double scale = std::max(std::abs(result.average), std::numeric_limits<double>::min());
    if (result.sigma <= 1e-12*scale)
    {
        return result;
    }
    result.skewness = m3/(m2*std::sqrt(m2));
    result.kurtosis = m4/(m2*m2) - 3;
    if (n < normalityMinimum)
    {
        return result;
    }
    std::sort(values.begin(), values.end());
    double a2 = 0;
    for (size_t i = 0; i < n; i++)
    {
        const double zlow  = (values[i] - result.average)/result.sigma;
        const double zhigh = (values[n-1-i] - result.average)/result.sigma;
        // ln(1-Phi(z)) is ln(Phi(-z))
        a2 += (2.0*i + 1)*(logNormalCdf(zlow) + logNormalCdf(-zhigh));
    }
    a2               = -static_cast<double>(n) - a2/n;
    result.statistic = a2*(1 + 0.75/n + 2.25/(static_cast<double>(n)*n));
    result.pvalue    = andersonDarlingPvalue(result.statistic);
    result.tested    = true;
    return result;
}

} // namespace alexandria

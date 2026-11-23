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

#ifndef ACT_STATISTICS_NORMALITY_H
#define ACT_STATISTICS_NORMALITY_H

#include <cstddef>
#include <vector>

namespace alexandria
{

//! The fewest values for which the test for normality is done
constexpr size_t normalityMinimum = 8;

/*! \brief How a set of values is distributed, and whether that can be a
 * normal distribution
 */
struct NormalityResult
{
    //! The number of values
    size_t n         = 0;
    //! Their average
    double average   = 0;
    //! Their standard deviation, with n-1 in the denominator
    double sigma     = 0;
    //! The skewness: zero for a symmetric distribution, positive with a
    //! tail to large values
    double skewness  = 0;
    //! The excess kurtosis: zero for a normal distribution, positive with
    //! heavy tails or outliers, negative for a flat distribution or one
    //! with two peaks
    double kurtosis  = 0;
    //! Whether the test was done: not with fewer than normalityMinimum
    //! values or with values that are all the same
    bool   tested    = false;
    //! The statistic A2 of the Anderson-Darling test, corrected for the
    //! number of values
    double statistic = 0;
    //! The probability of a statistic at least this large if the values
    //! come from a normal distribution
    double pvalue    = 1;

    /*! \brief Whether the values can come from a normal distribution
     * \param[in] alpha The significance level
     * \return false if the test was done and its p-value is below alpha
     */
    bool normal(double alpha) const { return !tested || pvalue >= alpha; }
};

/*! \brief Test whether values come from a normal distribution
 *
 * The Anderson-Darling test with the average and the standard deviation
 * estimated from the values. The statistic is corrected for the number of
 * values and the p-value follows from the formulae of D'Agostino and
 * Stephens (Goodness-of-fit techniques, Marcel Dekker, 1986, Table 4.9),
 * which hold for eight values or more. A p-value below 1e-30 is given as
 * zero.
 * With many values the test finds small deviations from a normal
 * distribution as well: the skewness and the kurtosis tell how large the
 * deviation is.
 * \param[in] values The values, in any order
 * \return the moments of the values and the result of the test
 */
NormalityResult testNormality(std::vector<double> values);

} // namespace alexandria

#endif

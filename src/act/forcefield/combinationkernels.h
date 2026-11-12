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
/*! \internal \file
 * \brief
 * The combination rules as functions for the processor and the GPU.
 *
 * The formulas of the combination rules exist only here. The processor
 * uses them through combinationrules.h, with the functions of the standard
 * library (StdMath), and GPU code can instantiate them with those of SYCL.
 * Therefore these functions do not throw and use no strings or containers:
 * the functions of combinationrules.h check the input and report errors.
 * The expressions are those that ACT used before, in the same order, such
 * that results do not change.
 *
 * \author David van der Spoel <david.vanderspoel@icm.uu.se>
 * \ingroup group_forcefield
 */
#ifndef ACT_FORCEFIELD_COMBINATIONKERNELS_H
#define ACT_FORCEFIELD_COMBINATIONKERNELS_H

#include <cmath>

namespace alexandria
{

//! The combination rules, see the manual for their equations
enum class CombRule {
    Geometric, Arithmetic, Volumetric, InverseSquare,
    HogervorstEpsilon, HogervorstSigma, Yang,
    WaldmanSigma, WaldmanEpsilon, HalgrenEpsilon,
    QiSigma, QiEpsilon, MasonGamma, Kronecker,
    GeneralizedMean
};

namespace combination
{

//! The mathematical functions of the standard library, for the processor
struct StdMath
{
    //! \return the square root of x
    static double sqrt(double x) { return std::sqrt(x); }
    //! \return x to the power y
    static double pow(double x, double y) { return std::pow(x, y); }
    //! \return the absolute value of x
    static double abs(double x) { return std::abs(x); }
};

/*! \brief Combine two values with a rule that only needs these values
 * \tparam M     The mathematical functions, such as StdMath
 * \tparam T     The type of the values
 * \param[in]  rule The combination rule
 * \param[in]  x1   First value
 * \param[in]  x2   Second value
 * \param[out] ok   Whether the rule is one of these
 * \return The combined value, zero if not ok
 */
template <class M, typename T>
inline T two(CombRule rule, T x1, T x2, bool *ok)
{
    *ok = true;
    switch (rule)
    {
    case CombRule::Geometric:
        return M::sqrt(x1*x2);
    case CombRule::Arithmetic:
        return 0.5 * (x1 + x2);
    case CombRule::HogervorstEpsilon:
        return (2.0 * x1 * x2)/( x1 + x2 );
    case CombRule::Yang:
        if (x1*x2 == 0)
        {
            return 0;
        }
        return ( x1 * x2 ) * ( x1 + x2 ) / (x1*x1 + x2*x2);
    case CombRule::WaldmanSigma:
        return M::pow(0.5*(M::pow(x1, 6) + M::pow(x2, 6)), (1.0/6.0));
    case CombRule::Volumetric:
        // Kriz2024a Eqn. 20
        return M::pow(0.5*(x1*x1*x1 + x2*x2*x2), (1.0/3.0));
    case CombRule::InverseSquare:
        // Kriz2024a Eqn. 21
        if (x1 == 0 || x2 == 0)
        {
            return 0;
        }
        else
        {
            T denom = 1/(x1*x1) + 1/(x2*x2);
            return M::sqrt(2/denom);
        }
    case CombRule::QiSigma:
        if (x1 == 0 && x2 == 0)
        {
            return 0;
        }
        // Qi2016a Eqn. 2
        return (x1*x1*x1+x2*x2*x2)/(x1*x1+x2*x2);
    case CombRule::HalgrenEpsilon:
    {
        // Halgren1992a, Eqn. 14
        T s = M::sqrt(x1) + M::sqrt(x2);
        return 4*x1*x2/(s*s);
    }
    default:
        *ok = false;
        return 0;
    }
}

/*! \brief The sigma of Hogervorst1971a, for gamma1 and gamma2 larger than 6, which the caller checks
 * \tparam M The mathematical functions
 * \tparam T The type of the values
 * \return The combined sigma
 */
template <class M, typename T>
inline T hogervorstSigma(T e1, T e2, T g1, T g2, T s1, T s2)
{
    T    tempi = M::abs(e1 * g1 * (M::pow(s1, 6)) /(g1 - 6 ));
    T    tempj = M::abs(e2 * g2 * (M::pow(s2, 6)) /(g2 - 6  ));
    bool ok;
    T    gam12 = two<M>(CombRule::Arithmetic, g1, g2, &ok);
    T    eps12 = two<M>(CombRule::HogervorstEpsilon, e1, e2, &ok);
    if (0 == eps12)
    {
        return 0;
    }
    return M::pow((M::sqrt( tempi * tempj ) )* M::abs(gam12 - 6) / (gam12 * eps12), 1.0/6.0);
}

/*! \brief The epsilon of Waldman and Hagler, which Qi et al. use too (Qi2016a Eqn. 3)
 * \tparam M The mathematical functions
 * \tparam T The type of the values
 * \return The combined epsilon
 */
template <class M, typename T>
inline T waldmanEpsilon(T e1, T e2, T s1, T s2)
{
    T s13 = s1*s1*s1;
    T s23 = s2*s2*s2;
    if (s13 == 0 && s23 == 0)
    {
        return 0;
    }
    return M::sqrt( e1 * e2 ) * ( 2 * s13 * s23 / (s13*s13 + s23*s23));
}

/*! \brief The generalized mean, see https://en.wikipedia.org/wiki/Generalized_mean
 * The paper by Hohm points out that for very small exponents, double
 * precision is not good enough, however the analytical result converges to
 * the geometric mean, which is therefore used for small exponents.
 * \tparam M The mathematical functions
 * \tparam T The type of the values
 * \return The combined value
 */
template <class M, typename T>
inline T generalizedMean(T x1, T x2, T exponent)
{
    const T tolerance = 0.1;
    if (M::abs(exponent) < tolerance)
    {
        return M::sqrt(x1*x2);
    }
    T sum = 0;
    if (x1 > 0)
    {
        sum += M::pow(x1, exponent);
    }
    if (x2 > 0)
    {
        sum += M::pow(x2, exponent);
    }
    if (sum > 0)
    {
        return M::pow(0.5*sum, 1/exponent);
    }
    return 0;
}

/*! \brief The gamma of Mason1955a
 * \tparam M The mathematical functions
 * \tparam T The type of the values
 * \return The combined gamma
 */
template <class M, typename T>
inline T masonGamma(T g1, T g2, T s1, T s2)
{
    T sigmaIJ = M::sqrt(s1*s2);
    return sigmaIJ * (0.5*((g1/s1)+(g2/s2)));
}

} // namespace combination

} // namespace alexandria

#endif

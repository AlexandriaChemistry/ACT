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
 * Explicit SIMD in double precision.
 *
 * Selects the SIMD primitives chosen with the CMake option ACT_SIMD,
 * AVX-512 (simd_avx512.h), AVX2 with FMA (simd_avx2.h) or ARM NEON
 * (simd_neon.h), if the compiler targets that instruction set, and defines
 * ACT_HAVE_SIMD if there are any. Both provide the type SimdDouble with
 * simdWidth doubles and the same operations, such that the code using them
 * does not depend on the processor.
 *
 * The math functions needed by the force routines are implemented here
 * from these primitives, such that no library functions are called inside
 * loops. Vectorized library functions force the compiler to save and
 * restore all registers around each call, which costs more than the math,
 * and some of them compute one element at a time.
 * \author David van der Spoel <david.vanderspoel@icm.uu.se>
 */
#ifndef ACT_FORCES_SIMD_MATH_H
#define ACT_FORCES_SIMD_MATH_H

#include "config.h"

#if ACT_SIMD_AVX_512 && defined(__AVX512F__)
#    include "act/forces/simd_avx512.h"
#    define ACT_HAVE_SIMD 1
#elif ACT_SIMD_AVX2 && defined(__AVX2__) && defined(__FMA__)
#    include "act/forces/simd_avx2.h"
#    define ACT_HAVE_SIMD 1
#elif ACT_SIMD_ARM_NEON && defined(__ARM_NEON) && defined(__aarch64__)
#    include "act/forces/simd_neon.h"
#    define ACT_HAVE_SIMD 1
#endif

#ifdef ACT_HAVE_SIMD

#include <cmath>

namespace alexandria
{

namespace simd
{

//! \return 1/sqrt(a), computed exactly as in the scalar code
static inline SimdDouble invsqrt(SimdDouble a) { return 1.0/sqrt(a); }

/*! \brief Evaluate a polynomial with the Horner scheme
 * \param[in] x The argument
 * \param[in] c The coefficients, highest order first
 * \return c[0]*x^(N-1) + c[1]*x^(N-2) + ... + c[N-1]
 */
template <size_t N>
static inline SimdDouble polynomial(SimdDouble x, const double (&c)[N])
{
    SimdDouble p(c[0]);
    for (size_t i = 1; i < N; i++)
    {
        p = fma(p, x, c[i]);
    }
    return p;
}

/*! \brief Exponential function
 * The algorithm of the Cephes library: exp(x) = 2^n exp(r) with
 * |r| <= ln(2)/2 and exp(r) from a Pade approximation. The relative
 * error is at most two units in the last place. Results below the
 * smallest normal number are flushed to zero.
 */
static inline SimdDouble exp(SimdDouble x)
{
    const double P[] = { 1.26177193074810590878E-4, 3.02994407707441961300E-2, 9.99999999999999999910E-1 };
    const double Q[] = { 3.00198505138664455042E-6, 2.52448340349684104192E-3, 2.27265548208155028766E-1, 2.00000000000000000009E0 };
    // ln(2) split in two parts, such that n*ln2Hi is exact
    const double ln2Hi  = 6.93145751953125E-1;
    const double ln2Lo  = 1.42860682030941723212E-6;
    const double log2e  = 1.4426950408889634073599;
    const double minArg = -708.39;
    const double maxArg = 709.78;
    SimdDouble   xc     = min(max(x, minArg), maxArg);
    SimdDouble   n      = floor(fma(xc, log2e, 0.5));
    SimdDouble   r      = fma(-n, ln2Lo, fma(-n, ln2Hi, xc));
    SimdDouble   r2     = r*r;
    SimdDouble   p      = r*polynomial(r2, P);
    SimdDouble   er     = 1.0 + 2.0*p/(polynomial(r2, Q) - p);
    // Multiply by 2^n as 2^n1 2^n2 with n1 + n2 = n, since 2^1024 does not exist
    SimdDouble   n1     = floor(0.5*n);
    SimdDouble   result = er*pow2(n1)*pow2(n - n1);
    return select(x < minArg, 0.0, select(x > maxArg, HUGE_VAL, result));
}

/*! \brief Error function
 * The rational approximations of the Cephes library: for |x| < 1
 * erf(x) = x T(x^2)/U(x^2), otherwise erf(x) = 1 - erfc(x) with
 * erfc(x) = exp(-x^2) P(x)/Q(x). Both are computed and the right one
 * selected, using a single division. The absolute error is at most a
 * few times 1e-16.
 */
static inline SimdDouble erf(SimdDouble x)
{
    const double T[] = { 9.60497373987051638749E0, 9.00260197203842689217E1, 2.23200534594684319226E3,
                         7.00332514112805075473E3, 5.55923013010394962768E4 };
    const double U[] = { 1.0, 3.35617141647503099647E1, 5.21357949780152679795E2, 4.59432382970980127987E3,
                         2.26290000613890934246E4, 4.92673942608635921086E4 };
    const double P[] = { 2.46196981473530512524E-10, 5.64189564831068821977E-1, 7.46321056442269912687E0,
                         4.86371970985681366614E1, 1.96520832956077098242E2, 5.26445194995477358631E2,
                         9.34528527171957607540E2, 1.02755188689515710272E3, 5.57535335369399327526E2 };
    const double Q[] = { 1.0, 1.32281951154744992508E1, 8.67072140885989742329E1, 3.54937778887819891062E2,
                         9.75708501743205489753E2, 1.82390916687909736289E3, 2.24633760818710981792E3,
                         1.65666309194161350182E3, 5.57535340817727675546E2 };
    SimdDouble ax    = abs(x);
    SimdBool   small = ax < 1.0;
    SimdDouble x2    = x*x;
    // Beyond 8 erfc is zero in double precision
    SimdDouble a     = min(ax, 8.0);
    SimdDouble num   = select(small, x*polynomial(x2, T), exp(-a*a)*polynomial(a, P));
    SimdDouble den   = select(small, polynomial(x2, U), polynomial(a, Q));
    SimdDouble ratio = num/den;
    SimdDouble large = 1.0 - ratio;
    return select(small, ratio, select(x < 0.0, -large, large));
}

} // namespace simd

} // namespace alexandria

#endif

#endif

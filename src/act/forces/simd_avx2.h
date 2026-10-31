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
 * SIMD primitives in double precision for x86 AVX2 with FMA.
 *
 * A SimdDouble holds four doubles. This is the best choice for processors
 * that have AVX2 but not AVX-512, like AMD Zen 2 and Zen 3 (for instance
 * EPYC 7xx2 and 7xx3) and Intel processors before Skylake-X. The
 * interface is the same as in simd_neon.h. Do not include this file
 * directly, but simd_math.h.
 * \author David van der Spoel <david.vanderspoel@icm.uu.se>
 */
#ifndef ACT_FORCES_SIMD_AVX2_H
#define ACT_FORCES_SIMD_AVX2_H

#include <immintrin.h>

#include <cstddef>
#include <cstdint>

namespace alexandria
{

namespace simd
{

//! Number of doubles in a SimdDouble
constexpr size_t simdWidth = 4;

/*! \brief Default unroll factor of the loops over conformations of pair interactions
 * See computePairConformations, it can be changed with the CMake option
 * ACT_PAIR_UNROLL for benchmarking.
 */
constexpr size_t pairUnrollDefault = 4;

//! Four doubles
class SimdDouble
{
public:
    //! The data
    __m256d v;
    //! Uninitialized
    SimdDouble() = default;
    /*! \brief Set all elements to a
     * \param[in] a The value
     */
    SimdDouble(double a) : v(_mm256_set1_pd(a)) {}
    /*! \brief From AVX type
     * \param[in] a The AVX vector
     */
    explicit SimdDouble(__m256d a) : v(a) {}
};

//! Four booleans, the result of comparing SimdDoubles
class SimdBool
{
public:
    //! The data, all bits set for true
    __m256d v;
    /*! \brief From AVX type
     * \param[in] a The AVX vector
     */
    explicit SimdBool(__m256d a) : v(a) {}
};

//! Load four doubles from memory
static inline SimdDouble load(const double *p) { return SimdDouble(_mm256_loadu_pd(p)); }

//! Store four doubles to memory
static inline void store(double *p, SimdDouble a) { _mm256_storeu_pd(p, a.v); }

static inline SimdDouble operator+(SimdDouble a, SimdDouble b) { return SimdDouble(_mm256_add_pd(a.v, b.v)); }
static inline SimdDouble operator-(SimdDouble a, SimdDouble b) { return SimdDouble(_mm256_sub_pd(a.v, b.v)); }
static inline SimdDouble operator*(SimdDouble a, SimdDouble b) { return SimdDouble(_mm256_mul_pd(a.v, b.v)); }
static inline SimdDouble operator/(SimdDouble a, SimdDouble b) { return SimdDouble(_mm256_div_pd(a.v, b.v)); }
//! Flip the sign bit, such that -0 is correct as well
static inline SimdDouble operator-(SimdDouble a) { return SimdDouble(_mm256_xor_pd(a.v, _mm256_set1_pd(-0.0))); }
static inline SimdDouble &operator+=(SimdDouble &a, SimdDouble b) { a = a + b; return a; }
static inline SimdDouble &operator-=(SimdDouble &a, SimdDouble b) { a = a - b; return a; }
static inline SimdDouble &operator*=(SimdDouble &a, SimdDouble b) { a = a * b; return a; }

static inline SimdBool operator<(SimdDouble a, SimdDouble b) { return SimdBool(_mm256_cmp_pd(a.v, b.v, _CMP_LT_OQ)); }
static inline SimdBool operator>(SimdDouble a, SimdDouble b) { return SimdBool(_mm256_cmp_pd(a.v, b.v, _CMP_GT_OQ)); }
static inline SimdBool operator&&(SimdBool a, SimdBool b) { return SimdBool(_mm256_and_pd(a.v, b.v)); }
static inline SimdBool operator!(SimdBool a)
{
    return SimdBool(_mm256_xor_pd(a.v, _mm256_castsi256_pd(_mm256_set1_epi64x(-1))));
}

//! \return a*b + c
static inline SimdDouble fma(SimdDouble a, SimdDouble b, SimdDouble c) { return SimdDouble(_mm256_fmadd_pd(a.v, b.v, c.v)); }

//! \return a where m is true, b otherwise
static inline SimdDouble select(SimdBool m, SimdDouble a, SimdDouble b) { return SimdDouble(_mm256_blendv_pd(b.v, a.v, m.v)); }

static inline SimdDouble min(SimdDouble a, SimdDouble b) { return SimdDouble(_mm256_min_pd(a.v, b.v)); }
static inline SimdDouble max(SimdDouble a, SimdDouble b) { return SimdDouble(_mm256_max_pd(a.v, b.v)); }
//! Clear the sign bit
static inline SimdDouble abs(SimdDouble a) { return SimdDouble(_mm256_andnot_pd(_mm256_set1_pd(-0.0), a.v)); }
static inline SimdDouble floor(SimdDouble a) { return SimdDouble(_mm256_floor_pd(a.v)); }
static inline SimdDouble sqrt(SimdDouble a) { return SimdDouble(_mm256_sqrt_pd(a.v)); }

/*! \brief Power of two
 * \param[in] n Integer values from -1022 to 1023
 * \return 2^n, made by putting n in the exponent bits
 */
static inline SimdDouble pow2(SimdDouble n)
{
    // Converting to 64-bit integers directly needs AVX-512, so go via 32 bits
    __m256i n64  = _mm256_cvtepi32_epi64(_mm256_cvtpd_epi32(n.v));
    __m256i bits = _mm256_slli_epi64(_mm256_add_epi64(n64, _mm256_set1_epi64x(1023)), 52);
    return SimdDouble(_mm256_castsi256_pd(bits));
}

} // namespace simd

} // namespace alexandria

#endif

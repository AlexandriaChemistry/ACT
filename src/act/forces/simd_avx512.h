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
 * SIMD primitives in double precision for x86 AVX-512.
 *
 * A SimdDouble holds eight doubles. Only the AVX512F instructions are
 * used, which all processors with AVX-512 have. The interface is the same
 * as in simd_neon.h. Do not include this file directly, but simd_math.h.
 * \author David van der Spoel <david.vanderspoel@icm.uu.se>
 */
#ifndef ACT_FORCES_SIMD_AVX512_H
#define ACT_FORCES_SIMD_AVX512_H

#include <immintrin.h>

#include <cstddef>
#include <cstdint>

namespace alexandria
{

namespace simd
{

//! Number of doubles in a SimdDouble
constexpr size_t simdWidth = 8;

//! Eight doubles
class SimdDouble
{
public:
    //! The data
    __m512d v;
    //! Uninitialized
    SimdDouble() = default;
    //! Set all elements to a
    SimdDouble(double a) : v(_mm512_set1_pd(a)) {}
    //! From AVX-512 type
    explicit SimdDouble(__m512d a) : v(a) {}
};

//! Eight booleans, the result of comparing SimdDoubles
class SimdBool
{
public:
    //! The data, one bit per element
    __mmask8 v;
    //! From AVX-512 type
    explicit SimdBool(__mmask8 a) : v(a) {}
};

//! Load eight doubles from memory
static inline SimdDouble load(const double *p) { return SimdDouble(_mm512_loadu_pd(p)); }

//! Store eight doubles to memory
static inline void store(double *p, SimdDouble a) { _mm512_storeu_pd(p, a.v); }

static inline SimdDouble operator+(SimdDouble a, SimdDouble b) { return SimdDouble(_mm512_add_pd(a.v, b.v)); }
static inline SimdDouble operator-(SimdDouble a, SimdDouble b) { return SimdDouble(_mm512_sub_pd(a.v, b.v)); }
static inline SimdDouble operator*(SimdDouble a, SimdDouble b) { return SimdDouble(_mm512_mul_pd(a.v, b.v)); }
static inline SimdDouble operator/(SimdDouble a, SimdDouble b) { return SimdDouble(_mm512_div_pd(a.v, b.v)); }
//! Flip the sign bit, such that -0 is correct as well
static inline SimdDouble operator-(SimdDouble a)
{
    return SimdDouble(_mm512_castsi512_pd(_mm512_xor_si512(_mm512_castpd_si512(a.v), _mm512_set1_epi64(INT64_MIN))));
}
static inline SimdDouble &operator+=(SimdDouble &a, SimdDouble b) { a = a + b; return a; }
static inline SimdDouble &operator-=(SimdDouble &a, SimdDouble b) { a = a - b; return a; }
static inline SimdDouble &operator*=(SimdDouble &a, SimdDouble b) { a = a * b; return a; }

static inline SimdBool operator<(SimdDouble a, SimdDouble b) { return SimdBool(_mm512_cmp_pd_mask(a.v, b.v, _CMP_LT_OQ)); }
static inline SimdBool operator>(SimdDouble a, SimdDouble b) { return SimdBool(_mm512_cmp_pd_mask(a.v, b.v, _CMP_GT_OQ)); }
static inline SimdBool operator&&(SimdBool a, SimdBool b) { return SimdBool(static_cast<__mmask8>(a.v & b.v)); }
static inline SimdBool operator!(SimdBool a) { return SimdBool(static_cast<__mmask8>(~a.v)); }

//! \return a*b + c
static inline SimdDouble fma(SimdDouble a, SimdDouble b, SimdDouble c) { return SimdDouble(_mm512_fmadd_pd(a.v, b.v, c.v)); }

//! \return a where m is true, b otherwise
static inline SimdDouble select(SimdBool m, SimdDouble a, SimdDouble b) { return SimdDouble(_mm512_mask_blend_pd(m.v, b.v, a.v)); }

static inline SimdDouble min(SimdDouble a, SimdDouble b) { return SimdDouble(_mm512_min_pd(a.v, b.v)); }
static inline SimdDouble max(SimdDouble a, SimdDouble b) { return SimdDouble(_mm512_max_pd(a.v, b.v)); }
static inline SimdDouble abs(SimdDouble a) { return SimdDouble(_mm512_abs_pd(a.v)); }
static inline SimdDouble floor(SimdDouble a) { return SimdDouble(_mm512_roundscale_pd(a.v, _MM_FROUND_TO_NEG_INF | _MM_FROUND_NO_EXC)); }
static inline SimdDouble sqrt(SimdDouble a) { return SimdDouble(_mm512_sqrt_pd(a.v)); }

/*! \brief Power of two
 * \param[in] n Integer values from -1022 to 1023
 * \return 2^n, made by putting n in the exponent bits
 */
static inline SimdDouble pow2(SimdDouble n)
{
    // Converting to 64-bit integers directly needs AVX512DQ, so go via 32 bits
    __m512i n64  = _mm512_cvtepi32_epi64(_mm512_cvtpd_epi32(n.v));
    __m512i bits = _mm512_slli_epi64(_mm512_add_epi64(n64, _mm512_set1_epi64(1023)), 52);
    return SimdDouble(_mm512_castsi512_pd(bits));
}

} // namespace simd

} // namespace alexandria

#endif

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
 * SIMD primitives in double and single precision for x86 AVX-512.
 *
 * A SimdDouble holds eight doubles, a SimdFloat sixteen floats. Only the AVX512F instructions are
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

/*! \brief Default unroll factor of the loops over conformations of pair interactions
 * See computePairConformations, it can be changed with the CMake option
 * ACT_PAIR_UNROLL for benchmarking.
 */
constexpr size_t pairUnrollDefault = 4;

//! Eight doubles
class SimdDouble
{
public:
    //! The data
    __m512d v;
    //! Uninitialized
    SimdDouble() = default;
    /*! \brief Set all elements to a
     * \param[in] a The value
     */
    SimdDouble(double a) : v(_mm512_set1_pd(a)) {}
    /*! \brief From AVX-512 type
     * \param[in] a The AVX-512 vector
     */
    explicit SimdDouble(__m512d a) : v(a) {}
};

//! Eight booleans, the result of comparing SimdDoubles
class SimdBool
{
public:
    //! The data, one bit per element
    __mmask8 v;
    /*! \brief From AVX-512 type
     * \param[in] a The AVX-512 mask
     */
    explicit SimdBool(__mmask8 a) : v(a) {}
};

//! Load eight doubles from memory
static inline SimdDouble load(const double *p) { return SimdDouble(_mm512_loadu_pd(p)); }

//! Store eight doubles to memory
static inline void store(double *p, SimdDouble a) { _mm512_storeu_pd(p, a.v); }

//! \return a + b, element-wise
static inline SimdDouble operator+(SimdDouble a, SimdDouble b) { return SimdDouble(_mm512_add_pd(a.v, b.v)); }
//! \return a - b, element-wise
static inline SimdDouble operator-(SimdDouble a, SimdDouble b) { return SimdDouble(_mm512_sub_pd(a.v, b.v)); }
//! \return a * b, element-wise
static inline SimdDouble operator*(SimdDouble a, SimdDouble b) { return SimdDouble(_mm512_mul_pd(a.v, b.v)); }
//! \return a / b, element-wise
static inline SimdDouble operator/(SimdDouble a, SimdDouble b) { return SimdDouble(_mm512_div_pd(a.v, b.v)); }
//! Flip the sign bit, such that -0 is correct as well
static inline SimdDouble operator-(SimdDouble a)
{
    return SimdDouble(_mm512_castsi512_pd(_mm512_xor_si512(_mm512_castpd_si512(a.v), _mm512_set1_epi64(INT64_MIN))));
}
//! Add b to a element-wise, \return a
static inline SimdDouble &operator+=(SimdDouble &a, SimdDouble b) { a = a + b; return a; }
//! Subtract b from a element-wise, \return a
static inline SimdDouble &operator-=(SimdDouble &a, SimdDouble b) { a = a - b; return a; }
//! Multiply a by b element-wise, \return a
static inline SimdDouble &operator*=(SimdDouble &a, SimdDouble b) { a = a * b; return a; }

//! \return a < b, element-wise
static inline SimdBool operator<(SimdDouble a, SimdDouble b) { return SimdBool(_mm512_cmp_pd_mask(a.v, b.v, _CMP_LT_OQ)); }
//! \return a > b, element-wise
static inline SimdBool operator>(SimdDouble a, SimdDouble b) { return SimdBool(_mm512_cmp_pd_mask(a.v, b.v, _CMP_GT_OQ)); }
//! \return Logical and of a and b, element-wise
static inline SimdBool operator&&(SimdBool a, SimdBool b) { return SimdBool(static_cast<__mmask8>(a.v & b.v)); }
//! \return Logical not of a, element-wise
static inline SimdBool operator!(SimdBool a) { return SimdBool(static_cast<__mmask8>(~a.v)); }

//! \return a*b + c
static inline SimdDouble fma(SimdDouble a, SimdDouble b, SimdDouble c) { return SimdDouble(_mm512_fmadd_pd(a.v, b.v, c.v)); }

//! \return a where m is true, b otherwise
static inline SimdDouble select(SimdBool m, SimdDouble a, SimdDouble b) { return SimdDouble(_mm512_mask_blend_pd(m.v, b.v, a.v)); }

//! \return Element-wise minimum of a and b
static inline SimdDouble min(SimdDouble a, SimdDouble b) { return SimdDouble(_mm512_min_pd(a.v, b.v)); }
//! \return Element-wise maximum of a and b
static inline SimdDouble max(SimdDouble a, SimdDouble b) { return SimdDouble(_mm512_max_pd(a.v, b.v)); }
//! \return Element-wise absolute value of a
static inline SimdDouble abs(SimdDouble a) { return SimdDouble(_mm512_abs_pd(a.v)); }
//! \return Element-wise a rounded down to an integer value
static inline SimdDouble floor(SimdDouble a) { return SimdDouble(_mm512_roundscale_pd(a.v, _MM_FROUND_TO_NEG_INF | _MM_FROUND_NO_EXC)); }
//! \return Element-wise square root of a
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

//! Number of floats in a SimdFloat
constexpr size_t simdFloatWidth = 16;

//! Sixteen floats
class SimdFloat
{
public:
    //! The data
    __m512 v;
    //! Uninitialized
    SimdFloat() = default;
    /*! \brief Set all elements to a
     * \param[in] a The value
     */
    SimdFloat(float a) : v(_mm512_set1_ps(a)) {}
    /*! \brief From AVX-512 type
     * \param[in] a The AVX-512 vector
     */
    explicit SimdFloat(__m512 a) : v(a) {}
};

//! Sixteen booleans, the result of comparing SimdFloats
class SimdFloatBool
{
public:
    //! The data, one bit per element
    __mmask16 v;
    /*! \brief From AVX-512 type
     * \param[in] a The mask
     */
    explicit SimdFloatBool(__mmask16 a) : v(a) {}
};

//! Load sixteen floats from memory
static inline SimdFloat load(const float *p) { return SimdFloat(_mm512_loadu_ps(p)); }

//! Store sixteen floats to memory
static inline void store(float *p, SimdFloat a) { _mm512_storeu_ps(p, a.v); }

/*! \brief Convert eight 64-bit integers to doubles, rounded once
 * The conversion of 64-bit integers needs AVX512DQ, the parts are
 * converted as in simd_avx2.h instead.
 * \param[in] v The integers
 * \return the doubles
 */
static inline __m512d int64ToDouble(__m512i v)
{
    const __m512i magicLo  = _mm512_set1_epi64(0x4330000000000000LL); // 2^52
    const __m512i magicHi  = _mm512_set1_epi64(0x4530000080000000LL); // 2^84 + 2^63
    const __m512d magicAll = _mm512_set1_pd(19342822341709703277445120.0); // 2^84 + 2^63 + 2^52
    __m512i       lo       = _mm512_mask_blend_epi32(0x5555, magicLo, v);
    __m512i       hi       = _mm512_xor_si512(_mm512_srli_epi64(v, 32), magicHi);
    __m512d       hiD      = _mm512_sub_pd(_mm512_castsi512_pd(hi), magicAll);
    return _mm512_add_pd(hiD, _mm512_castsi512_pd(lo));
}

/*! \brief The differences of two times eight integers, as doubles
 * \param[in] a The first integers
 * \param[in] b The second integers
 * \return a - b, computed with wrap-around in 64 bits
 */
static inline SimdDouble integerDifferenceDouble(const std::int64_t *a, const std::int64_t *b)
{
    __m512i d = _mm512_sub_epi64(_mm512_loadu_si512(a), _mm512_loadu_si512(b));
    return SimdDouble(int64ToDouble(d));
}

/*! \brief The differences of two times sixteen integers, as floats
 * Converted to double first, then rounded once to float.
 * \param[in] a The first integers
 * \param[in] b The second integers
 * \return a - b, computed with wrap-around in 64 bits
 */
static inline SimdFloat integerDifferenceFloat(const std::int64_t *a, const std::int64_t *b)
{
    __m256 lo = _mm512_cvtpd_ps(integerDifferenceDouble(a, b).v);
    __m256 hi = _mm512_cvtpd_ps(integerDifferenceDouble(a + 8, b + 8).v);
    return SimdFloat(_mm512_castpd_ps(_mm512_insertf64x4(_mm512_castps_pd(_mm512_castps256_ps512(lo)),
                                                         _mm256_castps_pd(hi), 1)));
}

/*! \brief Add sixteen floats to sixteen doubles in memory
 * \param[inout] p The doubles, two SimdDouble
 * \param[in]    a The floats
 */
static inline void addToDouble(double *p, SimdFloat a)
{
    // The upper half through doubles, since AVX512F has no extract of eight floats
    __m256 hi = _mm256_castpd_ps(_mm512_extractf64x4_pd(_mm512_castps_pd(a.v), 1));
    _mm512_storeu_pd(p, _mm512_add_pd(_mm512_loadu_pd(p), _mm512_cvtps_pd(_mm512_castps512_ps256(a.v))));
    _mm512_storeu_pd(p + 8, _mm512_add_pd(_mm512_loadu_pd(p + 8), _mm512_cvtps_pd(hi)));
}

//! \return a + b, element-wise
static inline SimdFloat operator+(SimdFloat a, SimdFloat b) { return SimdFloat(_mm512_add_ps(a.v, b.v)); }
//! \return a - b, element-wise
static inline SimdFloat operator-(SimdFloat a, SimdFloat b) { return SimdFloat(_mm512_sub_ps(a.v, b.v)); }
//! \return a * b, element-wise
static inline SimdFloat operator*(SimdFloat a, SimdFloat b) { return SimdFloat(_mm512_mul_ps(a.v, b.v)); }
//! \return a / b, element-wise
static inline SimdFloat operator/(SimdFloat a, SimdFloat b) { return SimdFloat(_mm512_div_ps(a.v, b.v)); }
//! \return -a, element-wise
static inline SimdFloat operator-(SimdFloat a)
{
    return SimdFloat(_mm512_castsi512_ps(_mm512_xor_si512(_mm512_castps_si512(a.v), _mm512_set1_epi32(INT32_MIN))));
}
//! Add b to a element-wise, \return a
static inline SimdFloat &operator+=(SimdFloat &a, SimdFloat b) { a = a + b; return a; }
//! Subtract b from a element-wise, \return a
static inline SimdFloat &operator-=(SimdFloat &a, SimdFloat b) { a = a - b; return a; }
//! Multiply a by b element-wise, \return a
static inline SimdFloat &operator*=(SimdFloat &a, SimdFloat b) { a = a * b; return a; }

//! \return a < b, element-wise
static inline SimdFloatBool operator<(SimdFloat a, SimdFloat b) { return SimdFloatBool(_mm512_cmp_ps_mask(a.v, b.v, _CMP_LT_OQ)); }
//! \return a > b, element-wise
static inline SimdFloatBool operator>(SimdFloat a, SimdFloat b) { return SimdFloatBool(_mm512_cmp_ps_mask(a.v, b.v, _CMP_GT_OQ)); }
//! \return Logical and of a and b, element-wise
static inline SimdFloatBool operator&&(SimdFloatBool a, SimdFloatBool b) { return SimdFloatBool(static_cast<__mmask16>(a.v & b.v)); }
//! \return Logical not of a, element-wise
static inline SimdFloatBool operator!(SimdFloatBool a) { return SimdFloatBool(static_cast<__mmask16>(~a.v)); }

//! \return a*b + c
static inline SimdFloat fma(SimdFloat a, SimdFloat b, SimdFloat c) { return SimdFloat(_mm512_fmadd_ps(a.v, b.v, c.v)); }

//! \return a where m is true, b otherwise
static inline SimdFloat select(SimdFloatBool m, SimdFloat a, SimdFloat b) { return SimdFloat(_mm512_mask_blend_ps(m.v, b.v, a.v)); }

//! \return Element-wise minimum of a and b
static inline SimdFloat min(SimdFloat a, SimdFloat b) { return SimdFloat(_mm512_min_ps(a.v, b.v)); }
//! \return Element-wise maximum of a and b
static inline SimdFloat max(SimdFloat a, SimdFloat b) { return SimdFloat(_mm512_max_ps(a.v, b.v)); }
//! \return Element-wise absolute value of a
static inline SimdFloat abs(SimdFloat a) { return SimdFloat(_mm512_abs_ps(a.v)); }
//! \return Element-wise a rounded down to an integer value
static inline SimdFloat floor(SimdFloat a) { return SimdFloat(_mm512_roundscale_ps(a.v, _MM_FROUND_TO_NEG_INF | _MM_FROUND_NO_EXC)); }
//! \return Element-wise square root of a
static inline SimdFloat sqrt(SimdFloat a) { return SimdFloat(_mm512_sqrt_ps(a.v)); }

/*! \brief Power of two
 * \param[in] n Integer values from -126 to 127
 * \return 2^n, made by putting n in the exponent bits
 */
static inline SimdFloat pow2(SimdFloat n)
{
    __m512i bits = _mm512_slli_epi32(_mm512_add_epi32(_mm512_cvtps_epi32(n.v), _mm512_set1_epi32(127)), 23);
    return SimdFloat(_mm512_castsi512_ps(bits));
}

} // namespace simd

} // namespace alexandria

#endif

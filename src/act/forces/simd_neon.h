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
 * SIMD primitives in double and single precision for ARM NEON (aarch64).
 *
 * A SimdDouble holds two doubles, a SimdFloat four floats. Only the operations that map directly
 * onto instructions are here, the math functions built from them are in
 * simd_math.h. Do not include this file directly, but simd_math.h.
 * \author David van der Spoel <david.vanderspoel@icm.uu.se>
 */
#ifndef ACT_FORCES_SIMD_NEON_H
#define ACT_FORCES_SIMD_NEON_H

#include <arm_neon.h>

#include <cstddef>
#include <cstdint>

namespace alexandria
{

namespace simd
{

//! Number of doubles in a SimdDouble
constexpr size_t simdWidth = 2;

/*! \brief Default unroll factor of the loops over conformations of pair interactions
 * See computePairConformations, it can be changed with the CMake option
 * ACT_PAIR_UNROLL for benchmarking.
 */
constexpr size_t pairUnrollDefault = 4;

//! Two doubles
class SimdDouble
{
public:
    //! The data
    float64x2_t v;
    //! Uninitialized
    SimdDouble() = default;
    /*! \brief Set all elements to a
     * \param[in] a The value
     */
    SimdDouble(double a) : v(vdupq_n_f64(a)) {}
    /*! \brief From NEON type
     * \param[in] a The NEON vector
     */
    explicit SimdDouble(float64x2_t a) : v(a) {}
};

//! Two booleans, the result of comparing SimdDoubles
class SimdBool
{
public:
    //! The data, all bits set for true
    uint64x2_t v;
    /*! \brief From NEON type
     * \param[in] a The NEON vector
     */
    explicit SimdBool(uint64x2_t a) : v(a) {}
};

//! Load two doubles from memory
static inline SimdDouble load(const double *p) { return SimdDouble(vld1q_f64(p)); }

//! Store two doubles to memory
static inline void store(double *p, SimdDouble a) { vst1q_f64(p, a.v); }

//! \return a + b, element-wise
static inline SimdDouble operator+(SimdDouble a, SimdDouble b) { return SimdDouble(vaddq_f64(a.v, b.v)); }
//! \return a - b, element-wise
static inline SimdDouble operator-(SimdDouble a, SimdDouble b) { return SimdDouble(vsubq_f64(a.v, b.v)); }
//! \return a * b, element-wise
static inline SimdDouble operator*(SimdDouble a, SimdDouble b) { return SimdDouble(vmulq_f64(a.v, b.v)); }
//! \return a / b, element-wise
static inline SimdDouble operator/(SimdDouble a, SimdDouble b) { return SimdDouble(vdivq_f64(a.v, b.v)); }
//! \return -a, element-wise
static inline SimdDouble operator-(SimdDouble a) { return SimdDouble(vnegq_f64(a.v)); }
//! Add b to a element-wise, \return a
static inline SimdDouble &operator+=(SimdDouble &a, SimdDouble b) { a = a + b; return a; }
//! Subtract b from a element-wise, \return a
static inline SimdDouble &operator-=(SimdDouble &a, SimdDouble b) { a = a - b; return a; }
//! Multiply a by b element-wise, \return a
static inline SimdDouble &operator*=(SimdDouble &a, SimdDouble b) { a = a * b; return a; }

//! \return a < b, element-wise
static inline SimdBool operator<(SimdDouble a, SimdDouble b) { return SimdBool(vcltq_f64(a.v, b.v)); }
//! \return a > b, element-wise
static inline SimdBool operator>(SimdDouble a, SimdDouble b) { return SimdBool(vcgtq_f64(a.v, b.v)); }
//! \return Logical and of a and b, element-wise
static inline SimdBool operator&&(SimdBool a, SimdBool b) { return SimdBool(vandq_u64(a.v, b.v)); }
//! \return Logical not of a, element-wise
static inline SimdBool operator!(SimdBool a) { return SimdBool(vreinterpretq_u64_u32(vmvnq_u32(vreinterpretq_u32_u64(a.v)))); }

//! \return a*b + c
static inline SimdDouble fma(SimdDouble a, SimdDouble b, SimdDouble c) { return SimdDouble(vfmaq_f64(c.v, a.v, b.v)); }

//! \return a where m is true, b otherwise
static inline SimdDouble select(SimdBool m, SimdDouble a, SimdDouble b) { return SimdDouble(vbslq_f64(m.v, a.v, b.v)); }

//! \return Element-wise minimum of a and b
static inline SimdDouble min(SimdDouble a, SimdDouble b) { return SimdDouble(vminq_f64(a.v, b.v)); }
//! \return Element-wise maximum of a and b
static inline SimdDouble max(SimdDouble a, SimdDouble b) { return SimdDouble(vmaxq_f64(a.v, b.v)); }
//! \return Element-wise absolute value of a
static inline SimdDouble abs(SimdDouble a) { return SimdDouble(vabsq_f64(a.v)); }
//! \return Element-wise a rounded down to an integer value
static inline SimdDouble floor(SimdDouble a) { return SimdDouble(vrndmq_f64(a.v)); }
//! \return Element-wise square root of a
static inline SimdDouble sqrt(SimdDouble a) { return SimdDouble(vsqrtq_f64(a.v)); }

/*! \brief Power of two
 * \param[in] n Integer values from -1022 to 1023
 * \return 2^n, made by putting n in the exponent bits
 */
static inline SimdDouble pow2(SimdDouble n)
{
    const int64x2_t bits = vshlq_n_s64(vaddq_s64(vcvtq_s64_f64(n.v), vdupq_n_s64(1023)), 52);
    return SimdDouble(vreinterpretq_f64_s64(bits));
}

//! Number of floats in a SimdFloat
constexpr size_t simdFloatWidth = 4;

//! Four floats
class SimdFloat
{
public:
    //! The data
    float32x4_t v;
    //! Uninitialized
    SimdFloat() = default;
    /*! \brief Set all elements to a
     * \param[in] a The value
     */
    SimdFloat(float a) : v(vdupq_n_f32(a)) {}
    /*! \brief From NEON type
     * \param[in] a The NEON vector
     */
    explicit SimdFloat(float32x4_t a) : v(a) {}
};

//! Four booleans, the result of comparing SimdFloats
class SimdFloatBool
{
public:
    //! The data, all bits set for true
    uint32x4_t v;
    /*! \brief From NEON type
     * \param[in] a The NEON vector
     */
    explicit SimdFloatBool(uint32x4_t a) : v(a) {}
};

//! Load four floats from memory
static inline SimdFloat load(const float *p) { return SimdFloat(vld1q_f32(p)); }

//! Store four floats to memory
static inline void store(float *p, SimdFloat a) { vst1q_f32(p, a.v); }

/*! \brief The differences of two times two integers, as doubles
 * \param[in] a The first integers
 * \param[in] b The second integers
 * \return a - b, computed with wrap-around in 64 bits
 */
static inline SimdDouble integerDifferenceDouble(const std::int64_t *a, const std::int64_t *b)
{
    return SimdDouble(vcvtq_f64_s64(vsubq_s64(vld1q_s64(a), vld1q_s64(b))));
}

/*! \brief The differences of two times four integers, as floats
 * Converted to double first, then rounded once to float.
 * \param[in] a The first integers
 * \param[in] b The second integers
 * \return a - b, computed with wrap-around in 64 bits
 */
static inline SimdFloat integerDifferenceFloat(const std::int64_t *a, const std::int64_t *b)
{
    const float64x2_t lo = vcvtq_f64_s64(vsubq_s64(vld1q_s64(a), vld1q_s64(b)));
    const float64x2_t hi = vcvtq_f64_s64(vsubq_s64(vld1q_s64(a + 2), vld1q_s64(b + 2)));
    return SimdFloat(vcvt_high_f32_f64(vcvt_f32_f64(lo), hi));
}

/*! \brief Add four floats to four doubles in memory
 * \param[inout] p The doubles, two SimdDouble
 * \param[in]    a The floats
 */
static inline void addToDouble(double *p, SimdFloat a)
{
    vst1q_f64(p, vaddq_f64(vld1q_f64(p), vcvt_f64_f32(vget_low_f32(a.v))));
    vst1q_f64(p + 2, vaddq_f64(vld1q_f64(p + 2), vcvt_high_f64_f32(a.v)));
}

//! \return a + b, element-wise
static inline SimdFloat operator+(SimdFloat a, SimdFloat b) { return SimdFloat(vaddq_f32(a.v, b.v)); }
//! \return a - b, element-wise
static inline SimdFloat operator-(SimdFloat a, SimdFloat b) { return SimdFloat(vsubq_f32(a.v, b.v)); }
//! \return a * b, element-wise
static inline SimdFloat operator*(SimdFloat a, SimdFloat b) { return SimdFloat(vmulq_f32(a.v, b.v)); }
//! \return a / b, element-wise
static inline SimdFloat operator/(SimdFloat a, SimdFloat b) { return SimdFloat(vdivq_f32(a.v, b.v)); }
//! \return -a, element-wise
static inline SimdFloat operator-(SimdFloat a) { return SimdFloat(vnegq_f32(a.v)); }
//! Add b to a element-wise, \return a
static inline SimdFloat &operator+=(SimdFloat &a, SimdFloat b) { a = a + b; return a; }
//! Subtract b from a element-wise, \return a
static inline SimdFloat &operator-=(SimdFloat &a, SimdFloat b) { a = a - b; return a; }
//! Multiply a by b element-wise, \return a
static inline SimdFloat &operator*=(SimdFloat &a, SimdFloat b) { a = a * b; return a; }

//! \return a < b, element-wise
static inline SimdFloatBool operator<(SimdFloat a, SimdFloat b) { return SimdFloatBool(vcltq_f32(a.v, b.v)); }
//! \return a > b, element-wise
static inline SimdFloatBool operator>(SimdFloat a, SimdFloat b) { return SimdFloatBool(vcgtq_f32(a.v, b.v)); }
//! \return Logical and of a and b, element-wise
static inline SimdFloatBool operator&&(SimdFloatBool a, SimdFloatBool b) { return SimdFloatBool(vandq_u32(a.v, b.v)); }
//! \return Logical not of a, element-wise
static inline SimdFloatBool operator!(SimdFloatBool a) { return SimdFloatBool(vmvnq_u32(a.v)); }

//! \return a*b + c
static inline SimdFloat fma(SimdFloat a, SimdFloat b, SimdFloat c) { return SimdFloat(vfmaq_f32(c.v, a.v, b.v)); }

//! \return a where m is true, b otherwise
static inline SimdFloat select(SimdFloatBool m, SimdFloat a, SimdFloat b) { return SimdFloat(vbslq_f32(m.v, a.v, b.v)); }

//! \return Element-wise minimum of a and b
static inline SimdFloat min(SimdFloat a, SimdFloat b) { return SimdFloat(vminq_f32(a.v, b.v)); }
//! \return Element-wise maximum of a and b
static inline SimdFloat max(SimdFloat a, SimdFloat b) { return SimdFloat(vmaxq_f32(a.v, b.v)); }
//! \return Element-wise absolute value of a
static inline SimdFloat abs(SimdFloat a) { return SimdFloat(vabsq_f32(a.v)); }
//! \return Element-wise a rounded down to an integer value
static inline SimdFloat floor(SimdFloat a) { return SimdFloat(vrndmq_f32(a.v)); }
//! \return Element-wise square root of a
static inline SimdFloat sqrt(SimdFloat a) { return SimdFloat(vsqrtq_f32(a.v)); }

/*! \brief Power of two
 * \param[in] n Integer values from -126 to 127
 * \return 2^n, made by putting n in the exponent bits
 */
static inline SimdFloat pow2(SimdFloat n)
{
    const int32x4_t bits = vshlq_n_s32(vaddq_s32(vcvtq_s32_f32(n.v), vdupq_n_s32(127)), 23);
    return SimdFloat(vreinterpretq_f32_s32(bits));
}

} // namespace simd

} // namespace alexandria

#endif

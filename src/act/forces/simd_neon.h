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
 * SIMD primitives in double precision for ARM NEON (aarch64).
 *
 * A SimdDouble holds two doubles. Only the operations that map directly
 * onto instructions are here, the math functions built from them are in
 * simd_math.h. Do not include this file directly, but simd_math.h.
 * \author David van der Spoel <david.vanderspoel@icm.uu.se>
 */
#ifndef ACT_FORCES_SIMD_NEON_H
#define ACT_FORCES_SIMD_NEON_H

#include <arm_neon.h>

#include <cstddef>

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

static inline SimdDouble operator+(SimdDouble a, SimdDouble b) { return SimdDouble(vaddq_f64(a.v, b.v)); }
static inline SimdDouble operator-(SimdDouble a, SimdDouble b) { return SimdDouble(vsubq_f64(a.v, b.v)); }
static inline SimdDouble operator*(SimdDouble a, SimdDouble b) { return SimdDouble(vmulq_f64(a.v, b.v)); }
static inline SimdDouble operator/(SimdDouble a, SimdDouble b) { return SimdDouble(vdivq_f64(a.v, b.v)); }
static inline SimdDouble operator-(SimdDouble a) { return SimdDouble(vnegq_f64(a.v)); }
static inline SimdDouble &operator+=(SimdDouble &a, SimdDouble b) { a = a + b; return a; }
static inline SimdDouble &operator-=(SimdDouble &a, SimdDouble b) { a = a - b; return a; }
static inline SimdDouble &operator*=(SimdDouble &a, SimdDouble b) { a = a * b; return a; }

static inline SimdBool operator<(SimdDouble a, SimdDouble b) { return SimdBool(vcltq_f64(a.v, b.v)); }
static inline SimdBool operator>(SimdDouble a, SimdDouble b) { return SimdBool(vcgtq_f64(a.v, b.v)); }
static inline SimdBool operator&&(SimdBool a, SimdBool b) { return SimdBool(vandq_u64(a.v, b.v)); }
static inline SimdBool operator!(SimdBool a) { return SimdBool(vreinterpretq_u64_u32(vmvnq_u32(vreinterpretq_u32_u64(a.v)))); }

//! \return a*b + c
static inline SimdDouble fma(SimdDouble a, SimdDouble b, SimdDouble c) { return SimdDouble(vfmaq_f64(c.v, a.v, b.v)); }

//! \return a where m is true, b otherwise
static inline SimdDouble select(SimdBool m, SimdDouble a, SimdDouble b) { return SimdDouble(vbslq_f64(m.v, a.v, b.v)); }

static inline SimdDouble min(SimdDouble a, SimdDouble b) { return SimdDouble(vminq_f64(a.v, b.v)); }
static inline SimdDouble max(SimdDouble a, SimdDouble b) { return SimdDouble(vmaxq_f64(a.v, b.v)); }
static inline SimdDouble abs(SimdDouble a) { return SimdDouble(vabsq_f64(a.v)); }
static inline SimdDouble floor(SimdDouble a) { return SimdDouble(vrndmq_f64(a.v)); }
static inline SimdDouble sqrt(SimdDouble a) { return SimdDouble(vsqrtq_f64(a.v)); }

/*! \brief Power of two
 * \param[in] n Integer values from -1022 to 1023
 * \return 2^n, made by putting n in the exponent bits
 */
static inline SimdDouble pow2(SimdDouble n)
{
    int64x2_t bits = vshlq_n_s64(vaddq_s64(vcvtq_s64_f64(n.v), vdupq_n_s64(1023)), 52);
    return SimdDouble(vreinterpretq_f64_s64(bits));
}

} // namespace simd

} // namespace alexandria

#endif

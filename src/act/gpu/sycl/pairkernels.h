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
 * Pair potentials for the GPU kernels, as the explicit SIMD routines in
 * forces_explicit_simd.cpp compute them, in single precision. Each
 * computes the energy terms and the scalar force divided by the distance,
 * which is added to the force on the first particle times the distance
 * vector from the second to the first, and subtracted from the other.
 * \author David van der Spoel <david.vanderspoel@icm.uu.se>
 */
#ifndef ACT_GPU_SYCL_PAIRKERNELS_H
#define ACT_GPU_SYCL_PAIRKERNELS_H

#include <sycl/sycl.hpp>

#include "gpu.h"

namespace alexandria
{

namespace gpu
{

//! The unit of the coordinates in nm, 2^-56, see coordinate.h
constexpr double coordinateUnitNm = 1.0/72057594037927936.0;

//! The number of units of the coordinates in a nm, 2^56
constexpr double coordinateScaleNm = 72057594037927936.0;

//! 1/sqrt(pi)
constexpr float invSqrtPi = 0.564189583547756287f;

/*! \brief Scalar force divided by r of two Gaussian charges, without charges
 * As gaussianForceOverR in forces_simd.cpp in single precision: below
 * x^2 = 0.16 the series is more accurate than the difference of the terms.
 * \param[in] zeta Width of the combined distribution
 * \param[in] rinv Inverse distance
 * \param[in] x    zeta times the distance
 * \param[in] erfx The error function of x
 * \return the scalar force divided by r
 */
inline float gaussianForceOverR(float zeta, float rinv, float x, float erfx)
{
    float x2 = x*x;
    if (x2 < 0.16f)
    {
        return zeta*zeta*zeta*2*invSqrtPi*
            (2.0f/3.0f + x2*(-2.0f/5.0f + x2*(1.0f/7.0f + x2*(-1.0f/27.0f + x2*(1.0f/132.0f)))));
    }
    return (erfx*rinv - 2*invSqrtPi*zeta*sycl::exp(-x2))*rinv*rinv;
}

/*! \brief Energy and scalar force divided by r of one term of a Coulomb interaction
 * As Coulomb::compute in forces_explicit_simd.cpp: particles on top of each
 * other only have an energy with Gaussian charges, and no force.
 * \param[in]  qq    The charge product
 * \param[in]  zeta  The combined width, zero for point charges
 * \param[in]  rsq   The distance squared
 * \param[in]  rinv  The inverse distance, not used if rsq is zero
 * \param[out] fscal The scalar force divided by r, added to
 * \return the energy
 */
inline float coulombTerm(float qq, float zeta, float rsq, float rinv, float *fscal)
{
    if (zeta == 0)
    {
        if (rsq > 0)
        {
            *fscal += qq*rinv*rinv*rinv;
            return qq*rinv;
        }
        return 0;
    }
    if (rsq > 0)
    {
        float rz   = zeta*rsq*rinv;
        float erfz = sycl::erf(rz);
        *fscal    += qq*gaussianForceOverR(zeta, rinv, rz, erfz);
        return qq*erfz*rinv;
    }
    return qq*2*zeta*invSqrtPi;
}

/*! \brief Compute one pair interaction
 * The parameters are those of PairEntry for the kernel.
 * \param[in]  kernel The potential
 * \param[in]  p      The parameters
 * \param[in]  rsq    The distance squared
 * \param[out] e      The energy terms, see PairKernel
 * \return the scalar force divided by the distance, zero for particles on
 *         top of each other
 */
inline float pairInteraction(PairKernel kernel, const float *p, float rsq, float *e)
{
    e[0] = 0;
    e[1] = 0;
    e[2] = 0;
    float fscal = 0;
    switch (kernel)
    {
    case PairKernel::Coulomb:
    {
        // The charge products in p[0..3], the widths in p[4..7]
        float rinv = rsq > 0 ? sycl::rsqrt(rsq) : 0.0f;
        e[0] = coulombTerm(p[0], p[4], rsq, rinv, &fscal);
        for (int k = 1; k < 4; k++)
        {
            if (p[k] != 0)
            {
                e[0] += coulombTerm(p[k], p[4+k], rsq, rinv, &fscal);
            }
        }
        return rsq > 0 ? fscal : 0.0f;
    }
    case PairKernel::Polarization:
    {
        // Force constant, distance and force constant of the quartic term
        if (!(rsq > 0))
        {
            // A shell on top of its core has no energy and no force
            return 0;
        }
        float ksh = p[0];
        e[0]  = 0.5f*ksh*rsq;
        fscal = -ksh;
        if (p[2] > 0)
        {
            float dr   = sycl::sqrt(rsq);
            float ddr  = sycl::fmax(dr - p[1], 0.0f);
            float ddr2 = ddr*ddr;
            e[0]  += p[2]*ddr2*ddr2;
            fscal -= 4*p[2]*ddr2*ddr/dr;
        }
        return fscal;
    }
    case PairKernel::LJ12_6:
    {
        // 4 epsilon sigma^6 and 4 epsilon sigma^12
        float rinv2 = 1.0f/rsq;
        float rinv6 = rinv2*rinv2*rinv2;
        float disp  = -p[0]*rinv6;
        float rep   = p[1]*rinv6*rinv6;
        e[0] = rep;
        e[1] = disp;
        return (12*rep + 6*disp)*rinv2;
    }
    case PairKernel::LJ14_7:
    {
        // 1/sigma, epsilon, gamma, delta
        float gamma1  = p[2] + 1;
        float rinv    = sycl::rsqrt(rsq);
        float rstar   = rsq*rinv*p[0];
        float rstar2  = rstar*rstar;
        float rstar6  = rstar2*rstar2*rstar2;
        float deltars = p[3] + rstar;
        float invdrs  = 1.0f/deltars;
        float t       = (p[3] + 1)*invdrs;
        float t2      = t*t;
        float repfac  = p[1]*t2*t2*t2*t;
        float invgr7  = 1.0f/(p[2] + rstar6*rstar);
        e[0] = repfac*gamma1*invgr7;
        e[1] = -2*repfac;
        float f147 = 7*repfac*(gamma1*rstar6*deltars*invgr7*invgr7 + gamma1*invgr7 - 2)*p[0]*invdrs;
        return f147*rinv;
    }
    case PairKernel::WangBuckingham:
    {
        // sigma^6, 3/(gamma+3), dispersion prefactor, gamma, 1/sigma
        float rinv     = sycl::rsqrt(rsq);
        float r        = rsq*rinv;
        float r5       = rsq*rsq*r;
        float inv6     = 1.0f/(p[0] + r5*r);
        float erep_exp = p[1]*sycl::exp(p[3]*(1 - r*p[4]));
        float vdisp    = -p[2]*p[0]*inv6;
        float fdisp    = vdisp*6*r5*inv6;
        e[0] = -vdisp*erep_exp;
        e[1] = vdisp;
        return (-fdisp*erep_exp - vdisp*erep_exp*p[3]*p[4] + fdisp)*rinv;
    }
    case PairKernel::DoubleExponential:
    {
        // A2 - A1 and the exponent, as bornMayer
        float rinv = sycl::rsqrt(rsq);
        float eexp = p[0]*sycl::exp(-p[1]*rsq*rinv);
        e[0] = eexp;
        return p[1]*eexp*rinv;
    }
    }
    return 0;
}

} // namespace gpu

} // namespace alexandria

#endif

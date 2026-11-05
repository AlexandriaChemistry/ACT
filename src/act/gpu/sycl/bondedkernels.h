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
 * Bonded potentials and virtual sites for the GPU kernels, as the explicit
 * SIMD routines in forces_explicit_simd.cpp and the routines of
 * vsitehandler.cpp compute them, in single precision.
 * \author David van der Spoel <david.vanderspoel@icm.uu.se>
 */
#ifndef ACT_GPU_SYCL_BONDEDKERNELS_H
#define ACT_GPU_SYCL_BONDEDKERNELS_H

#include <cfloat>

#include <sycl/sycl.hpp>

#include "gpu.h"

namespace alexandria
{

namespace gpu
{

//! Cross product
inline void cross(const float *a, const float *b, float *axb)
{
    axb[0] = a[1]*b[2] - a[2]*b[1];
    axb[1] = a[2]*b[0] - a[0]*b[2];
    axb[2] = a[0]*b[1] - a[1]*b[0];
}

//! Inner product
inline float dot(const float *a, const float *b)
{
    return a[0]*b[0] + a[1]*b[1] + a[2]*b[2];
}

/*! \brief The forces of a dihedral, as dihedralForces in forces_explicit_simd.cpp
 * \param[in]  r_ij  Vector from j to i
 * \param[in]  r_kj  Vector from j to k
 * \param[in]  r_kl  Vector from l to k
 * \param[in]  m     r_ij x r_kj
 * \param[in]  n     r_kj x r_kl
 * \param[in]  iprm  |m|^2
 * \param[in]  iprn  |n|^2
 * \param[in]  nrkj2 |r_kj|^2
 * \param[in]  nrkj  |r_kj|
 * \param[in]  s     1/(|m| |n|)
 * \param[in]  ddphi The derivative of the energy to the angle
 * \param[out] f     The forces on i, j, k and l
 */
inline void dihedralForces(const float *r_ij, const float *r_kj, const float *r_kl,
                           const float *m, const float *n, float iprm, float iprn,
                           float nrkj2, float nrkj, float s, float ddphi, float f[4][3])
{
    float toler  = nrkj2*FLT_EPSILON;
    bool  ok     = iprm > toler && iprn > toler;
    float dd     = (ok ? ddphi : 0.0f)*nrkj;
    float s2     = s*s;
    float a      = -dd*(s2*iprn);
    float b      = dd*(s2*iprm);
    float nrkj_2 = 1.0f/(ok ? nrkj2 : 1.0f);
    float p      = dot(r_ij, r_kj)*nrkj_2;
    float q      = dot(r_kl, r_kj)*nrkj_2;
    for (int d = 0; d < 3; d++)
    {
        float f_i  = a*m[d];
        float f_l  = b*n[d];
        float svec = p*f_i - q*f_l;
        f[0][d] = f_i;
        f[1][d] = -(f_i - svec);
        f[2][d] = -(f_l + svec);
        f[3][d] = f_l;
    }
}

/*! \brief Compute one bonded interaction
 * \param[in]  kernel The potential
 * \param[in]  p      The parameters, see BondedKernel
 * \param[in]  r      Positions of the second to fourth particle relative
 *                    to the first, r[0] of particle 1 and so on
 * \param[out] f      Forces on the particles, zero for those not used
 * \return the energy
 */
inline float bondedInteraction(BondedKernel kernel, const float *p, const float r[3][3], float f[4][3])
{
    for (int k = 0; k < 4; k++)
    {
        for (int d = 0; d < 3; d++)
        {
            f[k][d] = 0;
        }
    }
    switch (kernel)
    {
    case BondedKernel::HarmonicAngle:
    {
        // Atoms i j k, r_ij = x_i - x_j, r_kj = x_k - x_j
        float rij[3], rkj[3];
        for (int d = 0; d < 3; d++)
        {
            rij[d] = -r[0][d];
            rkj[d] = r[1][d] - r[0][d];
        }
        float ip    = dot(rij, rkj);
        float ipa   = dot(rij, rij);
        float ipb   = dot(rkj, rkj);
        float ipab  = ipa*ipb;
        float costh = ipab > 0 ? ip*sycl::rsqrt(ipab) : 1.0f;
        costh       = sycl::fmin(sycl::fmax(costh, -1.0f), 1.0f);
        float dth   = sycl::acos(costh) - p[1];
        float costh2 = costh*costh;
        float st     = costh2 < 1 ? (-p[0]*dth)*sycl::rsqrt(1 - costh2) : 0.0f;
        float sth    = st*costh;
        float nrij_1 = sycl::rsqrt(ipa);
        float nrkj_1 = sycl::rsqrt(ipb);
        float cik    = st*nrij_1*nrkj_1;
        float cii    = sth*nrij_1*nrij_1;
        float ckk    = sth*nrkj_1*nrkj_1;
        for (int d = 0; d < 3; d++)
        {
            float f_i = -(cik*rkj[d] - cii*rij[d]);
            float f_k = -(cik*rij[d] - ckk*rkj[d]);
            f[0][d] = f_i;
            f[1][d] = -f_i - f_k;
            f[2][d] = f_k;
        }
        return 0.5f*p[0]*dth*dth;
    }
    case BondedKernel::LinearAngle:
    {
        // As linearAngleOne in forces_simd.cpp
        float a = p[0], klin = p[1], b = 1 - a;
        float dr2 = 0;
        for (int d = 0; d < 3; d++)
        {
            float rij = -r[0][d];
            float rkj = r[1][d] - r[0][d];
            float dr  = -a*rij - b*rkj;
            dr2      += dr*dr;
            float f_i = a*klin*dr;
            float f_k = b*klin*dr;
            f[0][d] = f_i;
            f[1][d] = -(f_i + f_k);
            f[2][d] = f_k;
        }
        return 0.5f*klin*dr2;
    }
    case BondedKernel::FourierDihedral:
    case BondedKernel::HarmonicDihedral:
    {
        // Atoms i j k l, as dihedralVectors
        float r_ij[3], r_kj[3], r_kl[3], m[3], n[3];
        for (int d = 0; d < 3; d++)
        {
            r_ij[d] = -r[0][d];
            r_kj[d] = r[1][d] - r[0][d];
            r_kl[d] = r[1][d] - r[2][d];
        }
        cross(r_ij, r_kj, m);
        cross(r_kj, r_kl, n);
        float iprm  = dot(m, m);
        float iprn  = dot(n, n);
        float nrkj2 = dot(r_kj, r_kj);
        float mn    = dot(m, n);
        float ipr   = dot(r_ij, n);
        float mn2   = iprm*iprn;
        float s     = sycl::rsqrt(mn2 > 0 ? mn2 : 1.0f);
        float nrkj  = sycl::sqrt(nrkj2);
        if (kernel == BondedKernel::HarmonicDihedral)
        {
            float phi = sycl::atan2(nrkj*ipr, mn);
            dihedralForces(r_ij, r_kj, r_kl, m, n, iprm, iprn, nrkj2, nrkj, s, p[0]*phi, f);
            return 0.5f*p[0]*phi*phi;
        }
        // Fourier: c0..c5 in p[0..5], the number of them in p[6]
        int   nparam  = static_cast<int>(p[6]);
        float cos_phi = -(mn2 > 0 ? mn*s : 1.0f);
        float cosfac  = 1;
        float ener    = p[0];
        float dvdphi  = 0;
        for (int j = 1; j < nparam; j++)
        {
            dvdphi += j*p[j]*cosfac;
            cosfac *= cos_phi;
            ener   += p[j]*cosfac;
        }
        dvdphi *= nrkj*ipr*s;
        dihedralForces(r_ij, r_kj, r_kl, m, n, iprm, iprn, nrkj2, nrkj, s, dvdphi, f);
        return ener;
    }
    }
    return 0;
}

/*! \brief The position of a virtual site relative to its first parent
 * As constructSlots and the construction routines in vsitehandler.cpp,
 * with the parents ordered as those routines get them, see VsiteEntry.
 * \param[in]  kind The construction, see VsiteKind
 * \param[in]  p    The parameters
 * \param[in]  r    The positions of the other parents relative to the first,
 *                  r[0] of the second parent and so on
 * \param[out] x    The position of the site relative to the first parent
 */
inline void vsitePosition(VsiteKind kind, const float *p, const float r[3][3], float *x)
{
    switch (kind)
    {
    case VsiteKind::Copy:
        x[0] = x[1] = x[2] = 0;
        return;
    case VsiteKind::Linear2:
        for (int d = 0; d < 3; d++)
        {
            x[d] = p[0]*r[0][d];
        }
        return;
    case VsiteKind::Fd2:
    {
        float factor = p[0]*sycl::rsqrt(dot(r[0], r[0]));
        for (int d = 0; d < 3; d++)
        {
            x[d] = factor*r[0][d];
        }
        return;
    }
    case VsiteKind::Linear3:
        // c*xj + a*xi + b*xk with xi the first parent: a*0 + c*rj + b*rk
        for (int d = 0; d < 3; d++)
        {
            x[d] = (1 - p[0] - p[1])*r[0][d] + p[1]*r[1][d];
        }
        return;
    case VsiteKind::Fd3:
    {
        // xij = xi - xj, xjk = xk - xj, temp = xij + a xjk, x = xj + c temp
        float temp[3];
        for (int d = 0; d < 3; d++)
        {
            temp[d] = -r[0][d] + p[0]*(r[1][d] - r[0][d]);
        }
        float c = p[1]*sycl::rsqrt(dot(temp, temp));
        for (int d = 0; d < 3; d++)
        {
            x[d] = r[0][d] + c*temp[d];
        }
        return;
    }
    case VsiteKind::Out3:
    {
        // x = xi + a xij + b xik + c (xij x xik), xij = xj - xi, xik = xk - xi
        float temp[3];
        cross(r[0], r[1], temp);
        for (int d = 0; d < 3; d++)
        {
            x[d] = p[0]*r[0][d] + p[1]*r[1][d] + p[2]*temp[d];
        }
        return;
    }
    case VsiteKind::Linear4:
        // d*xj + a*xi + b*xk + c*xl with xi first
        for (int d = 0; d < 3; d++)
        {
            x[d] = (1 - p[0] - p[1] - p[2])*r[0][d] + p[1]*r[1][d] + p[2]*r[2][d];
        }
        return;
    }
}

} // namespace gpu

} // namespace alexandria

#endif

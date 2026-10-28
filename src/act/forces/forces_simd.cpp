/*
 * This source file is part of the Alexandria Chemistry Toolkit.
 *
 * Copyright (C) 2014-2026
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
#include "forces_simd.h"

#include <algorithm>
#include <array>
#include <type_traits>

#include "act/basics/msg_handler.h"
#include "act/coulombintegrals/gaussian_integrals.h"
#include "act/coulombintegrals/slater_integrals.h"
#include "act/forces/conformationdata.h"
#include "act/forcefield/forcefield_parametername.h"
#include "act/forces/forces_explicit_simd.h"
#include "act/forces/simd_information.h"
#include "act/forces/forcecomputerutils.h"
#include "act/molprop/topologyentry.h"
#include "act/topology/actmol.h"
#include "act/topology/coords_energy.h"
#include "act/topology/topology.h"
#include "gromacs/math/vec.h"
#include "gromacs/math/units.h"

namespace alexandria
{

#define MAX_ENER 3
static_assert(MAX_ENER == explicitSimdEnergyTerms);

/*! \brief Force inlining of the functions computing the interactions.
 * The loops over conformations can only be vectorized when everything
 * called inside them is inlined.
 */
#if defined(__GNUC__) || defined(__clang__)
#define ACT_ALWAYS_INLINE inline __attribute__((always_inline))
#else
#define ACT_ALWAYS_INLINE inline
#endif

/*! \brief Prevent inlining of a function.
 * Used to keep loops that are almost the same in different functions,
 * since the compiler may otherwise merge them into a loop that it
 * cannot vectorize.
 */
#if defined(__GNUC__) || defined(__clang__)
#define ACT_NOINLINE __attribute__((noinline))
#else
#define ACT_NOINLINE
#endif

/* The energy and forces of each potential are computed for one
 * interaction at a time by a function xxxOne. These functions are used
 * both for one conformation at a time (see calcSome) and for all
 * conformations at once (see calcConformations). They are templates in
 * the storage of the coordinates X and forces F, and the atom indices I
 * used to access them. In the latter case, the coordinates and forces of
 * the atoms in one interaction are stored in small local arrays, the
 * atom indices into those are 0, 1, ...  The global atom indices, needed
 * for looking up atom properties, are always taken from the topology entry.
 * The energies are added to e[0], e[1] and e[2].
 */

/*! \brief Spread scalar force over two atoms
 * \param[in]    fscalar The scalar force
 * \param[in]    dx      The distance vector
 * \param[in]    indices The atom indices
 * \param[inout] forces  The forces to be updated
 */
template <typename I, typename F>
static ACT_ALWAYS_INLINE void pairforces(double           fscalar,
                              const gmx::RVec &dx,
                              const I         &indices,
                              F               *forces)
{
    auto  &f = *forces;

    auto fijX          = fscalar*dx[XX];
    f[indices[0]][XX] += fijX;
    f[indices[1]][XX] -= fijX;
    auto fijY          = fscalar*dx[YY];
    f[indices[0]][YY] += fijY;
    f[indices[1]][YY] -= fijY;
    auto fijZ          = fscalar*dx[ZZ];
    f[indices[0]][ZZ] += fijZ;
    f[indices[1]][ZZ] -= fijZ;
}

//! Flat-bottom position restraint for one atom
template <typename I, typename X, typename F>
static ACT_ALWAYS_INLINE void fbposreOne(gmx_unused MsgHandler *msghandler,
                              const TopologyEntry   &b,
                              const I               &idx,
                              const X               &x,
                              F                     *forces,
                              double                *e)
{
    // Get the parameters. We have to know their names to do this.
    auto &params = b.params();
    auto K       = params[fbprK];
    auto R0      = params[fbprR0];
    // Get the atom indices
    auto ai      = idx[0];
    auto dr2     = iprod(x[ai], x[ai]);
    if (dr2 > R0*R0)
    {
        auto r1   = std::sqrt(dr2);
        auto dr   = r1-R0;
        e[0]     += 0.5*K*dr*dr;
        // F = -dV/dr x/r
        auto fpr  = -K*dr/r1;
        for (int m = 0; m < DIM; m++)
        {
            (*forces)[ai][m] += fpr*x[ai][m];
        }
#ifdef DEBUG
        if (msghandler && msghandler->debug())
        {
            msghandler->writeDebug(gmx::formatString("ACT FBPOSRES ai %d epr: %10g K: %g R0 %g dr2 %g\n",
                                                     ai, e[0], K, R0, dr2));
        }
#endif
    }
}

//! Lennard-Jones 12-6 for one pair
template <typename I, typename X, typename F>
static ACT_ALWAYS_INLINE void lj12_6One(gmx_unused MsgHandler *msghandler,
                             const TopologyEntry   &b,
                             const I               &indices,
                             const X               &x,
                             F                     *forces,
                             double                *e)
{
    // Get the parameters. We have to know their names to do this.
    auto &params    = b.params();
    auto sig_ij     = params[lj12_6SIGMA];
    auto eps_ij     = params[lj12_6EPSILON];
    auto sig6       = gmx::square(sig_ij*sig_ij*sig_ij);
    auto c6         = 4*eps_ij*sig6;
    auto c12        = c6*sig6;
    // Get the atom indices
    auto ai         = indices[0];
    auto aj         = indices[1];
    rvec dx;
    rvec_sub(x[ai], x[aj], dx);
    auto dr2        = iprod(dx, dx);
    auto rinv       = gmx::invsqrt(dr2);
    auto rinv2      = rinv*rinv;
    auto rinv6      = rinv2*rinv2*rinv2;
    auto vvdw_disp  = -c6*rinv6;
    auto vvdw_rep   = c12*rinv6*rinv6;
    auto flj        = (12*vvdw_rep + 6*vvdw_disp)*rinv2;
#ifdef DEBUG
    if (msghandler && msghandler->debug())
    {
        msghandler->writeDebug(gmx::formatString("ACT ai %d aj %d vvdw_rep: %10g vvdw_disp: %10g c6: %10g c12: %10g",
                                                 ai, aj, vvdw_rep, vvdw_disp, c6, c12));
    }
#endif
    pairforces(flj, dx, indices, forces);
    e[0] += vvdw_rep;
    e[1] += vvdw_disp;
}

//! Lennard-Jones 12-6-4 for one pair
template <typename I, typename X, typename F>
static ACT_ALWAYS_INLINE void lj12_6_4One(gmx_unused MsgHandler *msghandler,
                               const TopologyEntry   &b,
                               const I               &indices,
                               const X               &x,
                               F                     *forces,
                               double                *e)
{
    // Get the parameters. We have to know their names to do this.
    auto &params    = b.params();
    auto sig_ij     = params[lj12_6_4SIGMA];
    auto eps_ij     = params[lj12_6_4EPSILON];
    auto gamma      = params[lj12_6_4GAMMA];
    auto sig6       = gmx::square(sig_ij*sig_ij*sig_ij);
    auto c6         = 4*eps_ij*sig6;
    auto c12        = c6*sig6;
    // Get the atom indices
    auto ai         = indices[0];
    auto aj         = indices[1];
    rvec dx;
    rvec_sub(x[ai], x[aj], dx);
    auto dr2        = iprod(dx, dx);
    auto rinv       = gmx::invsqrt(dr2);
    auto rinv2      = rinv*rinv;
    auto rinv4      = rinv2*rinv2;
    auto rinv6      = rinv2*rinv2*rinv2;
    // Calculate energy terms
    auto vvdw_mid   = -gamma*rinv4;
    auto vvdw_disp  = -c6*rinv6;
    auto vvdw_rep   = c12*rinv6*rinv6;
    // Calculate force
    auto flj        = (12*vvdw_rep + 6*vvdw_disp + 4*vvdw_mid)*rinv2;
#ifdef DEBUG
    if (msghandler && msghandler->debug())
    {
        msghandler->writeDebug(gmx::formatString("ACT ai %d aj %d vvdw_rep: %10g vvdw_disp: %10g vvdw_mid: 10%g c6: %10g c12: %10g gamma: %10g",
                                                 ai, aj, vvdw_rep, vvdw_disp, vvdw_mid, c6, c12, gamma));
    }
#endif
    pairforces(flj, dx, indices, forces);
    // Pauli repulsion (1/r^12), dipole-dipole (1/r^6) and monopole-induced dipole (1/r^4)
    e[0] += vvdw_rep;
    e[1] += vvdw_disp;
    e[2] += vvdw_mid;
}

//! Lennard-Jones 8-6 for one pair
template <typename I, typename X, typename F>
static ACT_ALWAYS_INLINE void lj8_6One(gmx_unused MsgHandler *msghandler,
                            const TopologyEntry   &b,
                            const I               &indices,
                            const X               &x,
                            F                     *forces,
                            double                *e)
{
    // Get the parameters. We have to know their names to do this.
    auto &params    = b.params();
    auto sig_ij     = params[lj8_6SIGMA];
    auto eps_ij     = params[lj8_6EPSILON];
    auto sig6       = gmx::square(sig_ij*sig_ij*sig_ij);
    auto c6         = 4*eps_ij*sig6;
    auto c8         = 0.75*c6*sig_ij*sig_ij;
    // Get the atom indices
    auto ai         = indices[0];
    auto aj         = indices[1];
    rvec dx;
    rvec_sub(x[ai], x[aj], dx);
    auto dr2        = iprod(dx, dx);
    auto rinv       = gmx::invsqrt(dr2);
    auto rinv2      = rinv*rinv;
    auto rinv6      = rinv2*rinv2*rinv2;
    auto vvdw_disp  = -c6*rinv6;
    auto vvdw_rep   = c8*rinv6*rinv2;
    auto flj        = (8*vvdw_rep + 6*vvdw_disp)*rinv2;
#ifdef DEBUG
    if (msghandler && msghandler->debug())
    {
        msghandler->writeDebug(gmx::formatString("ACT ai %d aj %d vrep: %10g vdisp: %10g c6: %10g c8: %10g", ai, aj, vvdw_rep, vvdw_disp, c6, c8));
    }
#endif
    e[0] += vvdw_rep;
    e[1] += vvdw_disp;
    pairforces(flj, dx, indices, forces);
}

//! Return x^7
static ACT_ALWAYS_INLINE double power7(double x)
{
    return x*gmx::power6(x);
}

//! Buffered 14-7 for one pair
template <typename I, typename X, typename F>
static ACT_ALWAYS_INLINE void lj14_7One(gmx_unused MsgHandler *msghandler,
                             const TopologyEntry   &b,
                             const I               &indices,
                             const X               &x,
                             F                     *forces,
                             double                *e)
{
    // Get the parameters. We have to know their names to do this.
    auto &params    = b.params();
    auto sigma      = params[lj14_7SIGMA];
    auto epsilon    = params[lj14_7EPSILON];
    real f147       = 0;
    // Get the atom indices
    auto ai         = indices[0];
    auto aj         = indices[1];
    rvec dx;
    rvec_sub(x[ai], x[aj], dx);
    auto dr2        = iprod(dx, dx);
    auto rinv       = gmx::invsqrt(dr2);
    real eerep = 0, eedisp = 0;
    if (epsilon > 0 && sigma > 0)
    {
        auto gamma      = params[lj14_7GAMMA];
        auto delta      = params[lj14_7DELTA];
        // Integer powers are written out, std::pow is much slower
        real rstar      = dr2*rinv/sigma;
        real rstar6     = gmx::power6(rstar);
        real delta1     = delta + 1;
        real gamma1     = gamma + 1;
        real deltars    = delta + rstar;
        real repfac     = epsilon * power7(delta1/deltars);
        real gamrstar7  = gamma + rstar6*rstar;
        eerep           = repfac * gamma1/gamrstar7;
        eedisp          = -2 * repfac;
        // Since epsilon*delta1^7/deltars^8 = repfac/deltars
        f147            = 7*repfac*(gamma1*rstar6*deltars/gmx::square(gamrstar7) + gamma1/gamrstar7 - 2)/(sigma*deltars);
#ifdef DEBUG
        if (msghandler && msghandler->debug())
        {
            msghandler->writeDebug(gmx::formatString("ACT ai %d aj %d vvdw: %10g epsilon: %10g gamma: %10g sigma: %10g delta: %10g",
                                                     ai, aj, eerep + eedisp, epsilon, gamma, sigma, delta));
        }
#endif
        e[0]    += eerep;
        e[1]    += eedisp;
    }
    real fbond  = f147*rinv;
    pairforces(fbond, dx, indices, forces);
}

/*! \brief Low-level utility to compute Born-Mayer energy and forces
 * \param[in]    msghandler  For warnings and errors
 * \param[in]    indices The atom identifiers and parameters
 * \param[in]    x       The coordinates of all particles
 * \param[in]    aexp    The prefactor in the energy expression
 * \param[in]    bexp    The exponent in the energy expression
 * \param[inout] forces      The forces on all particles
 * \return total energy
 */
template <typename I, typename X, typename F>
static ACT_ALWAYS_INLINE double lowBornMayer(gmx_unused MsgHandler *msghandler,
                                  const I               &indices,
                                  const X               &x,
                                  double                 aexp,
                                  double                 bexp,
                                  F                     *forces)
{
    // Get the atom indices
    auto ai         = indices[0];
    auto aj         = indices[1];
    rvec dx;
    rvec_sub(x[ai], x[aj], dx);
    auto dr2        = iprod(dx, dx);
    auto rinv       = gmx::invsqrt(dr2);

    auto eeexp  = aexp*std::exp(-bexp*dr2*rinv);
#ifdef DEBUG
    if (msghandler && msghandler->debug())
    {
        msghandler->writeDebug(gmx::formatString("lowBornMayer ai %d aj %d r %g eeexp %g aexp %g bexp %g",
                                                 ai, aj, dr2*rinv, eeexp, aexp, bexp));
    }
#endif
    real fexp  = bexp*eeexp*rinv;
    pairforces(fexp, dx, indices, forces);

    return eeexp;
}

//! Born-Mayer for one pair
template <typename I, typename X, typename F>
static ACT_ALWAYS_INLINE void bornMayerOne(MsgHandler          *msghandler,
                                const TopologyEntry &b,
                                const I             &indices,
                                const X             &x,
                                F                   *forces,
                                double              *e)
{
    // Get the parameters. We have to know their names to do this.
    auto &params    = b.params();
    // Charge transfer correction, according to Eqn. 20, Walker et al.
    // https://doi.org/10.1002/jcc.26954
    e[0] += lowBornMayer(msghandler, indices, x, -params[expA],
                         params[expB], forces);
}

/*! \brief Low-level utility to compute Slater-ISA energy and forces
 * \param[in]    msghandler  For warnings and errors
 * \param[in]    indices The atom identifiers and parameters
 * \param[in]    x       The coordinates of all particles
 * \param[in]    aexp    The prefactor in the energy expression
 * \param[in]    bexp    The exponent in the energy expression
 * \param[inout] forces      The forces on all particles
 * \return total energy
 */
template <typename I, typename X, typename F>
static ACT_ALWAYS_INLINE double lowSlaterISA(gmx_unused MsgHandler *msghandler,
                                  const I               &indices,
                                  const X               &x,
                                  double                 aexp,
                                  double                 bexp,
                                  F                     *forces)
{
    const double third = 1.0/3.0;
    // Get the atom indices
    auto ai         = indices[0];
    auto aj         = indices[1];
    rvec dx;
    rvec_sub(x[ai], x[aj], dx);
    auto dr2        = iprod(dx, dx);
    auto rinv       = gmx::invsqrt(dr2);
    auto dr         = rinv*dr2;
    // Beyond Born-Mayer: Improved Models for Short-Range Repulsion in ab Initio Force Fields
    // Mary J. Van Vleet, Alston J. Misquitta, Anthony J. Stone, and J. R. Schmidt
    // J. Chem. Theory Comput. 12 (2016) 3851-3870
    real br    = bexp*dr;
    real Pbr   = (br*br*third + br + 1);
    real aterm = aexp*std::exp(-br);
    auto eeexp = aterm*Pbr;
#ifdef DEBUG
    if (msghandler && msghandler->debug())
    {
        msghandler->writeDebug(gmx::formatString("SlaterISA r  %g  eeexp %g aexp %g bexp %g aterm %g P(br) %g",
                                                 dr, eeexp, aexp, bexp, aterm, Pbr));
    }
#endif
    real fexp  = bexp*eeexp - bexp*aterm*(2*br*third + 1);
    real fbond = fexp*rinv;
    pairforces(fbond, dx, indices, forces);

    return eeexp;
}

//! Slater-ISA for one pair
template <typename I, typename X, typename F>
static ACT_ALWAYS_INLINE void slaterISAOne(MsgHandler          *msghandler,
                                const TopologyEntry &b,
                                const I             &indices,
                                const X             &x,
                                F                   *forces,
                                double              *e)
{
    // Get the parameters. We have to know their names to do this.
    auto &params    = b.params();
    e[0] += lowSlaterISA(msghandler, indices, x, params[expA], params[expB], forces);
}

//! Double exponential for one pair
template <typename I, typename X, typename F>
static ACT_ALWAYS_INLINE void doubleExponentialOne(gmx_unused MsgHandler *msghandler,
                                        const TopologyEntry   &b,
                                        const I               &indices,
                                        const X               &x,
                                        F                     *forces,
                                        double                *e)
{
    // Get the parameters. We have to know their names to do this.
    auto &params    = b.params();
    real aexp       = params[dexpA1] - params[dexpA2];
    real bexp       = params[dexpB];
    // Get the atom indices
    auto ai         = indices[0];
    auto aj         = indices[1];
    rvec dx;
    rvec_sub(x[ai], x[aj], dx);
    auto dr2        = iprod(dx, dx);
    auto rinv       = gmx::invsqrt(dr2);
    auto eeexp      = -aexp*std::exp(-bexp*dr2*rinv);
#ifdef DEBUG
    if (msghandler && msghandler->debug())
    {
        msghandler->writeDebug(gmx::formatString("r  %g  dexp %g aexp %g bexp %g",
                                                 dr2*rinv, eeexp, aexp, bexp));
    }
#endif
    real fexp  = bexp*eeexp;
    e[0]      += eeexp;

    real fbond  = fexp*rinv;
    pairforces(fbond, dx, indices, forces);
}

/*! \brief Compute one wang_buckingham interaction
 * \param[in]  sigma       Parameter
 * \param[in]  epsilon     Parameter
 * \param[in]  gamma       Parameter
 * \param[in]  rsq         Square distance
 * \param[in]  rinv        Inverse distance
 * \param[out] vrepulsion  Repulsion energy
 * \param[out] vdispersion Dispersion energy
 * \param[out] fvdw        Scalar force
 */
static ACT_ALWAYS_INLINE void wang_buckingham(real sigma, real epsilon, real gamma,
                                   real rsq, real rinv,
                                   real *vrepulsion,
                                   real *vdispersion,
                                   real *fvdw)
{
    /* Modified Buckingham: JCTC  Volume: 9  Page: 452  Year: 2012 */
    real r           = rsq*rinv;
    real r5          = rsq*rsq*r;
    real r6          = r5*r;
    real sigma2      = sigma*sigma;
    real sigma6      = sigma2*sigma2*sigma2;
    real sigma6_r6   = sigma6 + r6;
    real gamma_3     = 3.0/(gamma + 3.0);
    real gamma3_inv  = 1.0/(1.0 - gamma_3);
    real disp_pre    = 2 * epsilon * gamma3_inv;
    real erep_exp    = gamma_3*std::exp(gamma*(1-(r/sigma)));

    real vvdw_disp   = -disp_pre * (sigma6 / sigma6_r6);
    *vdispersion     = vvdw_disp;
    *vrepulsion      = -vvdw_disp*erep_exp;

    real fvdw_disp   = - disp_pre * sigma6 * 6 * r5 / (sigma6_r6 * sigma6_r6);
    real fvdw_rep    = - fvdw_disp*erep_exp - vvdw_disp*erep_exp*(gamma/sigma);
    *fvdw            = fvdw_rep + fvdw_disp;
}

//! Wang-Buckingham for one pair
template <typename I, typename X, typename F>
static ACT_ALWAYS_INLINE void wbhOne(gmx_unused MsgHandler *msghandler,
                          const TopologyEntry   &b,
                          const I               &indices,
                          const X               &x,
                          F                     *forces,
                          double                *e)
{
    // Get the parameters. We have to know their names to do this.
    auto &params    = b.params();
    auto sigma      = params[wbhSIGMA];
    auto epsilon    = params[wbhEPSILON];
    auto gamma      = params[wbhGAMMA];
    if (epsilon > 0 && gamma > 0 && sigma > 0)
    {
        // Get the atom indices
        rvec dx;
        rvec_sub(x[indices[0]], x[indices[1]], dx);
        auto dr2        = iprod(dx, dx);
        auto rinv       = gmx::invsqrt(dr2);
        real eerep = 0, eedisp = 0, fwbh = 0;
        wang_buckingham(sigma, epsilon, gamma, dr2, rinv, &eerep, &eedisp, &fwbh);
        e[0]     += eerep;
        e[1]     += eedisp;
#ifdef DEBUG
        if (msghandler && msghandler->debug())
        {
            msghandler->writeDebug(gmx::formatString("WBHAM ai %d aj %d sigma %g epsilon %g gamma %g erep: %g edisp: %g",
                                                     indices[0], indices[1], sigma, epsilon, gamma, eerep, eedisp));
        }
#endif
        real fbond  = fwbh*rinv;
        pairforces(fbond, dx, indices, forces);
    }
}

//! Buckingham for one pair
template <typename I, typename X, typename F>
static ACT_ALWAYS_INLINE void buckinghamOne(gmx_unused MsgHandler *msghandler,
                                 const TopologyEntry   &b,
                                 const I               &indices,
                                 const X               &x,
                                 F                     *forces,
                                 double                *e)
{
    // Get the parameters. We have to know their names to do this.
    auto &params    = b.params();
    auto Abh  = params[bhA];
    auto bbh  = params[bhB];
    auto c6bh = params[bhC6];
    if (Abh > 0 && bbh > 0 && c6bh > 0)
    {
        // Get the atom indices
        rvec dx;
        rvec_sub(x[indices[0]], x[indices[1]], dx);
        auto dr2    = iprod(dx, dx);
        auto rinv   = gmx::invsqrt(dr2);
        auto rinv2  = rinv*rinv;
        real eerep  = Abh*std::exp(-bbh*dr2*rinv);
        real eedisp = -c6bh*rinv2*rinv2*rinv2;
        e[0]     += eerep;
        e[1]     += eedisp;
#ifdef DEBUG
        if (msghandler && msghandler->debug())
        {
            msghandler->writeDebug(gmx::formatString("BHAM ai %d aj %d A %g b %g c6 %g erep: %g edisp: %g",
                                                     indices[0], indices[1], Abh, bbh, c6bh, eerep, eedisp));
        }
#endif
        real fbond  = (bbh*eerep + 6*eedisp*rinv)*rinv;
        pairforces(fbond, dx, indices, forces);
    }
}

/*! \brief Low-level utility to compute Tang-Toennies dispersion energy and forces
 * \param[in]    msghandler  For warnings and errors
 * \param[in]    indices The atom identifiers and parameters
 * \param[in]    x       The coordinates of all particles
 * \param[in]    bDisp   The exponent in the energy expression
 * \param[in]    ctt     The dispersion coefficient C6, C8, C10
 * \param[inout] forces  The forces on all particles
 * \return total energy
 */
template <typename I, typename X, typename F>
static ACT_ALWAYS_INLINE double TT_Dispersion(gmx_unused MsgHandler *msghandler,
                                   const I               &indices,
                                   const X               &x,
                                   double                 bDisp,
                                   const double           ctt[3],
                                   F                     *forces)
{
    rvec dx;
    rvec_sub(x[indices[0]], x[indices[1]], dx);
    auto dr2     = iprod(dx, dx);
    auto rinv1   = gmx::invsqrt(dr2);
    auto rinv2   = rinv1*rinv1;
    real bDispr  = bDisp*dr2*rinv1;
    real ebrDisp = std::exp(-bDispr);
    real eedisp  = 0;
    real fdisp   = 0;
    real rinvn   = rinv2*rinv2*rinv2;
    // The damping function for C_2m is f_2m = 1 - exp(-br) sum_k=0^2m (br)^k/k!
    // The sum is accumulated in fk, with term holding (br)^k/k!
    real fk      = 1;
    real term    = 1;
    int  k       = 0;
#pragma unroll
    for (int m = 3; m <= 5; m++)
    {
        int  nn  = 2*m;
#pragma unroll
        for (; k < nn; k++)
        {
            term *= bDispr/(k+1);
            fk   += term;
        }
        real ed = -(1-ebrDisp*fk)*ctt[m-3]*rinvn;
        eedisp += ed;
        fdisp  += ed*nn*rinv1;
        // The derivative of the sum is b times the sum up to 2m-1,
        // the terms up to 2m-1 cancel against the derivative of exp(-br).
        fdisp  += ctt[m-3]*rinvn*bDisp*ebrDisp*term;
        rinvn  *= rinv2;
    }
#ifdef DEBUG
    if (msghandler && msghandler->debug())
    {
        msghandler->writeDebug(gmx::formatString("TT_Dispersion ai %d aj %d dr %g bDisp %g c6 %g c8 %g c10 %g edisp: %g fdisp: %g",
                                                 indices[0], indices[1], 1/rinv1, bDisp, ctt[0], ctt[1], ctt[2], eedisp, fdisp));
    }
#endif
    real fbond = fdisp*rinv1;
    pairforces(fbond, dx, indices, forces);

    return eedisp;
}

//! Tang-Toennies for one pair
template <typename I, typename X, typename F>
static ACT_ALWAYS_INLINE void tangToenniesOne(MsgHandler          *msghandler,
                                   const TopologyEntry &b,
                                   const I             &indices,
                                   const X             &x,
                                   F                   *forces,
                                   double              *e)
{
    // Get the parameters. We have to know their names to do this.
    auto &params    = b.params();
    // Call low level routines
    e[0] += lowBornMayer(msghandler, indices, x, params[ttA], params[ttB], forces);
    const double ctt[3] = { params[ttC6], params[ttC8], params[ttC10] };
    e[1] += TT_Dispersion(msghandler, indices, x, params[ttB], ctt, forces);
}

//! Tang-Toennies with separate exponents for one pair
template <typename I, typename X, typename F>
static ACT_ALWAYS_INLINE void tt2bOne(MsgHandler          *msghandler,
                           const TopologyEntry &b,
                           const I             &indices,
                           const X             &x,
                           F                   *forces,
                           double              *e)
{
    // Get the parameters. We have to know their names to do this.
    auto &params    = b.params();
    // Call low level routines
    e[0] += lowBornMayer(msghandler, indices, x, params[tt2bA], params[tt2bBexch], forces);
    const double ctt[3] = { params[tt2bC6], params[tt2bC8], params[tt2bC10] };
    e[1] += TT_Dispersion(msghandler, indices, x, params[tt2bBdisp], ctt, forces);
}

//! Slater-ISA exchange and Tang-Toennies dispersion for one pair
template <typename I, typename X, typename F>
static ACT_ALWAYS_INLINE void slaterISATTOne(MsgHandler          *msghandler,
                                  const TopologyEntry &b,
                                  const I             &indices,
                                  const X             &x,
                                  F                   *forces,
                                  double              *e)
{
    // Get the parameters. We have to know their names to do this.
    auto &params    = b.params();
    // Call low level routines
    e[0] += lowSlaterISA(msghandler, indices, x, params[SIttA], params[SIttBexch], forces);
    const double ctt[3] = { params[SIttC6], params[SIttC8], params[SIttC10] };
    e[1] += TT_Dispersion(msghandler, indices, x, params[SIttBdisp], ctt, forces);
}

/*! \brief Generalized Buckingham for one pair
 * Source: Werhahn et al. Chem. Phys. Lett. 619 (2015) 133-138, Eqn. 6.
 * Code was double checked here https://github.com/AlexandriaChemistry/ACT/issues/1299
 */
template <typename I, typename X, typename F>
static ACT_ALWAYS_INLINE void generalizedBhamOne(gmx_unused MsgHandler *msghandler,
                                      const TopologyEntry   &b,
                                      const I               &indices,
                                      const X               &x,
                                      F                     *forces,
                                      double                *e)
{
    // Get the parameters. We have to know their names to do this.
    auto &params    = b.params();
    auto rmin       = params[gbhRMIN];
    auto epsilon    = params[gbhEPSILON];
    auto gamma      = params[gbhGAMMA];
    auto delta      = params[gbhDELTA];
    // \todo Move this check to creation of pairlist? Will that work during training?
    if (epsilon > 0 && gamma > 0 && rmin > 0 && delta > 0)
    {
        // Get the atom indices
        rvec dx;
        rvec_sub(x[indices[0]], x[indices[1]], dx);
        auto dr2        = iprod(dx, dx);
        auto rinv       = gmx::invsqrt(dr2);

        real delta6     = 6+delta;
        real delta6gam2 = delta6 + 2*gamma;
        real rstar      = dr2*rinv/rmin;
        // Integer powers are written out, std::pow is much slower
        real rstar5     = gmx::power5(rstar);
        real rstardelta = std::pow(rstar, delta);
        real sixterm    = 1 + rstar5*rstar;
        real delterm    = 1 + rstardelta;
        real expterm    = std::exp(gamma*(1 - rstar));
        real sixdenom   = 1/(2*gamma*sixterm);
        real eerep      = epsilon*delta6*expterm*sixdenom;
        real eedisp     = -epsilon*(delta6gam2*sixdenom + 1/delterm);
        real fgbham     = (epsilon*((-6*(delta6*(1 - expterm) + 2*gamma)*rstar5)/(gamma*gmx::square(sixterm)) +
                                    (delta6*expterm)/sixterm - (2*delta*rstardelta)/(rstar*gmx::square(delterm))))/(2.*rmin);
#ifdef DEBUG
        if (msghandler && msghandler->debug())
        {
            msghandler->writeDebug(gmx::formatString("vrep: ai %d aj %d %g vdisp: %g epsilon: %10g gamma: %10g sigma: %10g delta: %10g",
                                                     indices[0], indices[1], eerep, eedisp, epsilon, gamma, rmin, delta));
        }
#endif
        e[0]     += eerep;
        e[1]     += eedisp;
        real fbond  = fgbham*rinv;
        pairforces(fbond, dx, indices, forces);
    }
}

//! Coulomb between Gaussian (or point) charges for one pair
template <typename I, typename X, typename F>
static ACT_ALWAYS_INLINE void coulombGaussianOne(gmx_unused MsgHandler      *msghandler,
                                      const TopologyEntry        &b,
                                      const I                    &indices,
                                      const X                    &x,
                                      F                          *forces,
                                      double                     *e,
                                      const std::vector<ActAtom> &atoms,
                                      double                      epsilonr)
{
    // Global atom indices for the charges
    auto &gindices = b.atomIndices();
    auto ai       = indices[0];
    auto aj       = indices[1];
    // Get the parameters.
    auto &params= b.params();
    auto izeta  = params[coulZETA];
    auto jzeta  = params[coulZETA2];
    real qq         = ONE_4PI_EPS0*atoms[gindices[0]].charge()*atoms[gindices[1]].charge()/epsilonr;
    // Width of the combined Gaussian, as in Coulomb_GG. This depends only on
    // the parameters, which allows the compiler to compute it once for all
    // conformations.
    real zeta;
    if (izeta == 0 || jzeta == 0)
    {
        zeta = izeta + jzeta;
    }
    else
    {
        zeta = izeta*jzeta*gmx::invsqrt(izeta*izeta + jzeta*jzeta);
    }
    rvec dx;
    rvec_sub(x[ai], x[aj], dx);
    auto dr2        = iprod(dx, dx);
    if (dr2 > 0)
    {
        // As Coulomb_GG and DCoulomb_GG, but computing erf only once
        auto rinv   = gmx::invsqrt(dr2);
        auto rz     = zeta*dr2*rinv;
        auto erfrz  = std::erf(rz);
        real velec  = qq*erfrz*rinv;
        // The scalar force divided by r
        real felec  = qq*(erfrz*rinv - 2*invsqrt_pi*zeta*std::exp(-rz*rz))*rinv*rinv;
#ifdef DEBUG
        if (msghandler && msghandler->debug())
        {
            msghandler->writeDebug(gmx::formatString("vcoul ai %d aj %d %g izeta %g jzeta %g qi %g qj %g vcoul_pc %g dist %g felec %g",
                                                     gindices[0], gindices[1], velec, izeta, jzeta,
                                                     atoms[gindices[0]].charge(), atoms[gindices[1]].charge(), qq*rinv, 1/rinv, felec));
        }
#endif
        e[0] += velec;
        pairforces(felec, dx, indices, forces);
    }
    else
    {
        e[0] += qq*Nuclear_GG(0, zeta);
    }
}

/*! \brief Coulomb between point charges for one pair
 * This is the same as coulombGaussianOne with both widths zero, but the
 * branches for non-zero widths are removed by the compiler, which matters
 * for vectorization.
 */
template <typename I, typename X, typename F>
static ACT_ALWAYS_INLINE void coulombPointOne(gmx_unused MsgHandler      *msghandler,
                                              const TopologyEntry        &b,
                                              const I                    &indices,
                                              const X                    &x,
                                              F                          *forces,
                                              double                     *e,
                                              const std::vector<ActAtom> &atoms,
                                              double                      epsilonr)
{
    // Global atom indices for the charges
    auto &gindices = b.atomIndices();
    real qq        = ONE_4PI_EPS0*atoms[gindices[0]].charge()*atoms[gindices[1]].charge()/epsilonr;
    rvec dx;
    rvec_sub(x[indices[0]], x[indices[1]], dx);
    auto dr2       = iprod(dx, dx);
    auto r         = std::sqrt(dr2);
    real velec     =  qq*Nuclear_GG(r, 0.0);
    real felec     = -qq*DNuclear_GG(r, 0.0);
    e[0] += velec;
    if (dr2 > 0)
    {
        felec *= gmx::invsqrt(dr2);
        pairforces(felec, dx, indices, forces);
    }
}

//! Coulomb between Slater charges for one pair
template <typename I, typename X, typename F>
static ACT_ALWAYS_INLINE void coulombSlaterOne(gmx_unused MsgHandler      *msghandler,
                                    const TopologyEntry        &b,
                                    const I                    &indices,
                                    const X                    &x,
                                    F                          *forces,
                                    double                     *e,
                                    const std::vector<ActAtom> &atoms,
                                    double                      epsilonr)
{
    // Global atom indices for the charges and rows
    auto &gindices = b.atomIndices();
    auto ai       = indices[0];
    auto aj       = indices[1];
    // Get the parameters. We have to know their names to do this.
    auto &params= b.params();
    auto izeta  = params[coulZETA];
    auto jzeta  = params[coulZETA2];
    auto irow   = atoms[gindices[0]].row();
    auto jrow   = atoms[gindices[1]].row();
    real qq     = ONE_4PI_EPS0*atoms[gindices[0]].charge()*atoms[gindices[1]].charge()/epsilonr;
    rvec dx;
    rvec_sub(x[ai], x[aj], dx);
    auto dr2   = iprod(dx, dx);
    real r1    = std::sqrt(dr2);
    real velec =  qq*Coulomb_SS(r1, irow, jrow, izeta, jzeta);
    // The DCoulomb_SS code returns the derivative of the energy
    // so we still need a minus sign.
    // https://github.com/dspoel/ACT/issues/369
    real felec = -qq*DCoulomb_SS(r1, irow, jrow, izeta, jzeta);
#ifdef DEBUG
    if (msghandler && msghandler->debug())
    {
        auto r1 = std::sqrt(dr2);
        msghandler->writeDebug(gmx::formatString("vcoul ai %d aj %d %g fcoul %g izeta %g jzeta %g qi %g qj %g vcoul_pc %g fcoul_pc %g dist %g",
                                                 gindices[0], gindices[1], velec, felec, izeta, jzeta, atoms[gindices[0]].charge(),
                                                 atoms[gindices[1]].charge(), qq/r1, qq/dr2, r1));
    }
#endif
    e[0] += velec;
    if (dr2 > 0)
    {
        felec *= gmx::invsqrt(dr2);
        pairforces(felec, dx, indices, forces);
    }
}

/*! \brief Compute harmonic force and energy
 * \param[in]  k  The force constant
 * \param[in]  x0 The reference position/angle
 * \param[in]  x  The actual position/angle
 * \param[out] V  The energy
 * \param[out] F  The scalar force
 */
static ACT_ALWAYS_INLINE void harmonic(real k, real x0, real x, real *V, real *F)
{
    const real half = 0.5;

    real dx  = x-x0;
    real dx2 = dx*dx;

    *F  = -k*dx;
    *V  = half*k*dx2;
}

//! Harmonic bond for one bond
template <typename I, typename X, typename F>
static ACT_ALWAYS_INLINE void bondOne(gmx_unused MsgHandler *msghandler,
                           const TopologyEntry   &b,
                           const I               &indices,
                           const X               &x,
                           F                     *forces,
                           double                *e)
{
    // Get the parameters. We have to know their names to do this.
    auto &params    = b.params();
    auto bondlength = params[bondLENGTH];
    auto kb         = params[bondKB];
    auto D0         = params[bondENERGY];

    rvec dx;
    rvec_sub(x[indices[0]], x[indices[1]], dx);
    auto dr2        = iprod(dx, dx);
    auto r_1        = gmx::invsqrt(dr2);
    auto r1         = dr2*r_1;

    real vB, fbond;
    harmonic(kb, bondlength, r1, &vB, &fbond);

    fbond *= r_1;
    e[0]  += vB+D0;

    pairforces(fbond, dx, indices, forces);
}

//! Cubic bond for one bond
template <typename I, typename X, typename F>
static ACT_ALWAYS_INLINE void cubicOne(gmx_unused MsgHandler *msghandler,
                            const TopologyEntry   &b,
                            const I               &indices,
                            const X               &x,
                            F                     *forces,
                            double                *e)
{
    // Get the parameters. We have to know their names to do this.
    auto &params    = b.params();
    auto bondlength = params[cubicLENGTH];
    auto rmax       = params[cubicRMAX];
    auto kb         = params[cubicKB];
    auto De         = params[cubicDE];

    rvec dx;
    rvec_sub(x[indices[0]], x[indices[1]], dx);
    double rhi     = (2*rmax+bondlength)/3;
    double r_1, r1 = rhi;
    {
        // Do not allow the potential to go to -infinity
        auto r2 = std::min(rhi*rhi, iprod(dx, dx));
        r_1     = gmx::invsqrt(r2);
        r1      = r2*r_1;
    }
    real vB, fbond;
    vB    = kb*gmx::square(r1-bondlength)*(rmax-r1) - De;
    fbond = r_1*kb*(bondlength+2*rmax-3*r1)*(bondlength-r1);

    e[0] += vB;

    pairforces(fbond, dx, indices, forces);
}

//! Morse bond for one bond
template <typename I, typename X, typename F>
static ACT_ALWAYS_INLINE void morseOne(gmx_unused MsgHandler *msghandler,
                            const TopologyEntry   &b,
                            const I               &indices,
                            const X               &x,
                            F                     *forces,
                            double                *e)
{
    // Get the parameters. We have to know their names to do this.
    auto &params    = b.params();
    auto bondlength = params[morseLENGTH];
    auto beta       = params[morseBETA];
    auto De         = params[morseDE];
    auto D0         = params[morseD0];
    rvec dx;
    rvec_sub(x[indices[0]], x[indices[1]], dx);
    auto dr         = norm(dx);
    auto expterm    = std::exp(-beta*(dr-bondlength));
    auto expterm_1  = 1-expterm;
    //        auto vbond      = De*expterm*(expterm - 2) + D0;
    auto vbond      = De*(expterm_1*expterm_1 - 1) + D0;
    auto fbond      = 2*De*beta*expterm*(expterm-1)/dr;
    e[0]           += vbond;
#ifdef DEBUG
    if (msghandler && msghandler->debug())
    {
        msghandler->writeDebug(gmx::formatString("Morse i %d j %d dr %g vbond %g fbond %g",
                                                 indices[0], indices[1], dr, vbond, fbond));
    }
#endif
    pairforces(fbond, dx, indices, forces);
}

//! Hua bond for one bond
template <typename I, typename X, typename F>
static ACT_ALWAYS_INLINE void huaOne(gmx_unused MsgHandler *msghandler,
                          const TopologyEntry   &b,
                          const I               &indices,
                          const X               &x,
                          F                     *forces,
                          double                *e)
{
    // Get the parameters. We have to know their names to do this.
    auto &params = b.params();
    auto re      = params[huaLENGTH];
    auto hb      = params[huaB];
    auto De      = params[huaDE];
    auto c       = params[huaC];
    rvec dx;
    rvec_sub(x[indices[0]], x[indices[1]], dx);
    auto dr2        = iprod(dx, dx);
    auto dr         = std::sqrt(dr2);
    auto expterm    = std::exp(-hb*(dr-re));
    auto exp2       = (1-expterm)/(1-c*expterm);
    auto vbond      = De*(exp2*exp2-1);
    auto denom      = (c*expterm - 1);
    auto fbond      = 2*De*hb*expterm*(expterm-1)*(c-1)/(dr*denom*denom*denom);
    e[0]           += vbond;

    pairforces(fbond, dx, indices, forces);
}

//! Linear angle for one angle
template <typename I, typename X, typename F>
static ACT_ALWAYS_INLINE void linearAngleOne(gmx_unused MsgHandler *msghandler,
                                  const TopologyEntry   &aaa,
                                  const I               &indices,
                                  const X               &x,
                                  F                     *forces,
                                  double                *e)
{
    auto   &f     = *forces;
    // Get the parameters. We have to know their names to do this.
    auto &params= aaa.params();
    auto a      = params[linangA];
    auto klin   = params[linangKLIN];

    rvec r_ij, r_kj, r_ik;
    rvec_sub(x[indices[0]], x[indices[1]], r_ij);
    rvec_sub(x[indices[2]], x[indices[1]], r_kj);
    rvec_sub(r_ij, r_kj, r_ik);

    real dr2 = 0;
    real b   = 1-a;
#pragma unroll
    for (int m = 0; (m < DIM); m++)
    {
        // Minus sign due to F = -dV/dr
        real dr   = -a * r_ij[m] - b * r_kj[m];
        dr2      += dr*dr;
        real f_i  = a*klin*dr;
        real f_k  = b*klin*dr;
        real f_j  = -(f_i+f_k);
        f[indices[0]][m] += f_i;
        f[indices[1]][m] += f_j;
        f[indices[2]][m] += f_k;
    }
    e[0] += 0.5*klin*dr2;
}

/*! \brief Angle between the bonds i-j and j-k
 * This is the same as bond_angle in forcecomputerutils, but can be
 * inlined here.
 */
static ACT_ALWAYS_INLINE real bondAngle(const rvec xi,
                             const rvec xj,
                             const rvec xk,
                             rvec       r_ij,
                             rvec       r_kj,
                             real      *costh)
{
    real th;

    rvec_sub(xi, xj, r_ij);
    rvec_sub(xk, xj, r_kj);

    *costh = cos_angle(r_ij, r_kj);
    th     = std::acos(*costh);

    return th;
}

//! Harmonic angle for one angle
template <typename I, typename X, typename F>
static ACT_ALWAYS_INLINE void angleOne(gmx_unused MsgHandler *msghandler,
                            const TopologyEntry   &a,
                            const I               &indices,
                            const X               &x,
                            F                     *forces,
                            double                *e)
{
    double costh = 0;
    auto   &f    = *forces;
    // Get the parameters. We have to know their names to do this.
    auto &params    = a.params();
    auto theta0     = params[angleANGLE];
    auto ka         = params[angleKT];

    rvec r_ij, r_kj;
    auto theta = bondAngle(x[indices[0]], x[indices[1]], x[indices[2]],
                           r_ij, r_kj, &costh);

    real vA, fangle;
    harmonic(ka, theta0, theta, &vA, &fangle);
#ifdef DEBUG
    if (msghandler && msghandler->debug())
    {
        msghandler->writeDebug(gmx::formatString("Angle %g refangle %g, energy %g", theta, theta0, vA));
    }
#endif
    e[0] += vA;

    auto costh2 = gmx::square(costh);
    if (costh2 < 1)
    {
        real st, sth;
        real cik, cii, ckk;
        real nrkj2, nrij2;
        real nrkj_1, nrij_1;

        st    = fangle*gmx::invsqrt(1 - costh2);
        sth   = st*costh;

        nrij2 = iprod(r_ij, r_ij);
        nrkj2 = iprod(r_kj, r_kj);

        nrij_1 = gmx::invsqrt(nrij2);
        nrkj_1 = gmx::invsqrt(nrkj2);

        cik = st*nrij_1*nrkj_1;
        cii = sth*nrij_1*nrij_1;
        ckk = sth*nrkj_1*nrkj_1;
#pragma unroll
        for (auto m = 0; m < DIM; m++)
        {
            auto f_im    = -(cik*r_kj[m] - cii*r_ij[m]);
            auto f_km    = -(cik*r_ij[m] - ckk*r_kj[m]);
            auto f_jm    = -f_im - f_km;
            f[indices[0]][m] += f_im;
            f[indices[1]][m] += f_jm;
            f[indices[2]][m] += f_km;
        }
    }
}

//! Urey-Bradley angle for one angle
template <typename I, typename X, typename F>
static ACT_ALWAYS_INLINE void ureyBradleyOne(gmx_unused MsgHandler *msghandler,
                                  const TopologyEntry   &a,
                                  const I               &indices,
                                  const X               &x,
                                  F                     *forces,
                                  double                *e)
{
    double costh = 0;
    auto   &f    = *forces;
    const  real half = 0.5;
    // Get the parameters. We have to know their names to do this.
    auto &params    = a.params();
    auto theta0     = params[ubANGLE];
    auto ka         = params[ubKT];
    auto r13        = params[ubR13];
    auto kUB        = params[ubKUB];
    // Get the atom indices
    auto ai = indices[0];
    auto aj = indices[1];
    auto ak = indices[2];

    rvec r_ij, r_kj;
    auto theta = bondAngle(x[ai], x[aj], x[ak], r_ij, r_kj, &costh);

    // Compute deviation from the reference angle
    auto da  = theta - theta0;
    auto da2 = da*da;

    auto fangle      = -ka*da;
    e[0]            += half*ka*da2;

    // Now UB part
    rvec r_ik;
    rvec_sub(x[ai], x[ak], r_ik);
    auto rik2        = iprod(r_ik, r_ik);
    auto norm_rik    = std::sqrt(rik2);

    real vUB, fUB;
    harmonic(kUB, r13, norm_rik, &vUB, &fUB);

    e[0]   += vUB;
    fUB    *= gmx::invsqrt(rik2);

    auto costh2 = gmx::square(costh);
    if (costh2 < 1)
    {
        real st, sth;
        real cik, cii, ckk;
        real nrkj2, nrij2;
        real nrkj_1, nrij_1;

        st    = fangle*gmx::invsqrt(1 - costh2);
        sth   = st*costh;

        nrij2 = iprod(r_ij, r_ij);
        nrkj2 = iprod(r_kj, r_kj);

        nrij_1 = gmx::invsqrt(nrij2);
        nrkj_1 = gmx::invsqrt(nrkj2);

        cik = st*nrij_1*nrkj_1;
        cii = sth*nrij_1*nrij_1;
        ckk = sth*nrkj_1*nrkj_1;
#pragma unroll
        for (auto m = 0; m < DIM; m++)
        {
            real f_i = -(cik*r_kj[m] - cii*r_ij[m]);
            real f_k = -(cik*r_ij[m] - ckk*r_kj[m]);
            real f_j = -f_i - f_k;
            f[ai][m] += f_i + fUB*r_ik[m];
            f[aj][m] += f_j;
            f[ak][m] += f_k - fUB*r_ik[m];
        }
    }
}

//! Polarization for one core-shell pair
template <typename I, typename X, typename F>
static ACT_ALWAYS_INLINE void polarizationOne(gmx_unused MsgHandler      *msghandler,
                                   const TopologyEntry        &b,
                                   const I                    &indices,
                                   const X                    &x,
                                   F                          *forces,
                                   double                     *e,
                                   const std::vector<ActAtom> &atoms)
{
    const real half = 0.5;
    // Get the parameters. We have to know their names to do this.
    auto &params   = b.params();
    // Get polarizability
    auto alpha     = params[polALPHA];
    if (alpha == 0)
    {
        return;
    }
    auto rhyper    = params[polRHYPER];
    auto fchyper   = params[polFCHYPER];
    // Get shell charge
    auto q         = atoms[b.atomIndices()[1]].charge();
    double ksh     = ONE_4PI_EPS0*q*q/alpha;
    rvec dx;
    rvec_sub(x[indices[0]], x[indices[1]], dx);
    auto dr2  = iprod(dx, dx);

    if (dr2 == 0.0)
    {
        return;
    }

    auto dr         = dr2 * gmx::invsqrt(dr2);
    auto fbond      = -ksh;
    e[0]           += half*ksh*dr2;
#ifdef DEBUG
    if (msghandler && msghandler->debug())
    {
        msghandler->writeDebug(gmx::formatString("vpol %g ai %d aj %d alpha %g ksh %g dist %g",
                                                 half*ksh*dr2, indices[0],
                                                 indices[1], alpha, ksh, dr));
    }
#endif
    if (dr > rhyper && fchyper > 0)
    {
        auto ddr  = dr - rhyper;
        auto ddr2 = ddr * ddr;
        e[0]  += fchyper * ddr2 * ddr2;
        // We will multiply by the distance vector components below
        fbond -= 4 * fchyper * ddr2 * ddr/dr;
#ifdef DEBUG
        if (msghandler && msghandler->debug())
        {
            msghandler->writeDebug(gmx::formatString("Adding hyperpolarization energy correction %g kJ/mol core pos %g %g %g ddr = %g.",
                                                     fchyper * ddr2 * ddr2, x[indices[0]][XX],
                                                     x[indices[0]][YY],x[indices[0]][ZZ], ddr));
        }
#endif
    }
    pairforces(fbond, dx, indices, forces);
}

//! \brief Low-level routine to update forces from dihedrals
template <typename F>
static ACT_ALWAYS_INLINE void do_dih_fup_noshiftf(int i, int j, int k, int l, real ddphi,
                                       rvec r_ij, rvec r_kj, rvec r_kl,
                                       rvec m, rvec n,
                                       F *forces)
{
    rvec f_i, f_j, f_k, f_l;
    rvec uvec, vvec, svec;
    real iprm, iprn, nrkj, nrkj2, nrkj_1, nrkj_2;
    real a, b, p, q, toler;

    iprm  = iprod(m, m);       /*  5    */
    iprn  = iprod(n, n);       /*  5	*/
    nrkj2 = iprod(r_kj, r_kj); /*  5	*/
    toler = nrkj2*GMX_REAL_EPS;
    if ((iprm > toler) && (iprn > toler))
    {
        nrkj_1 = gmx::invsqrt(nrkj2); /* 10	*/
        nrkj_2 = nrkj_1*nrkj_1;       /*  1	*/
        nrkj   = nrkj2*nrkj_1;        /*  1	*/
        a      = -ddphi*nrkj/iprm;    /* 11	*/
        svmul(a, m, f_i);             /*  3	*/
        b     = ddphi*nrkj/iprn;      /* 11	*/
        svmul(b, n, f_l);             /*  3  */
        p     = iprod(r_ij, r_kj);    /*  5	*/
        p    *= nrkj_2;               /*  1	*/
        q     = iprod(r_kl, r_kj);    /*  5	*/
        q    *= nrkj_2;               /*  1	*/
        svmul(p, f_i, uvec);          /*  3	*/
        svmul(q, f_l, vvec);          /*  3	*/
        rvec_sub(uvec, vvec, svec);   /*  3	*/
        rvec_sub(f_i, svec, f_j);     /*  3	*/
        rvec_add(f_l, svec, f_k);     /*  3	*/
        rvec_inc((*forces)[i], f_i);          /*  3	*/
        rvec_dec((*forces)[j], f_j);          /*  3	*/
        rvec_dec((*forces)[k], f_k);          /*  3	*/
        rvec_inc((*forces)[l], f_l);          /*  3	*/
    }
}

//! \brief Low-level code to compute dihedral angles
static ACT_ALWAYS_INLINE real dih_angle(const rvec xi, const rvec xj, const rvec xk, const rvec xl,
                             rvec r_ij, rvec r_kj, rvec r_kl, rvec m, rvec n)
{
    rvec_sub(xi, xj, r_ij); /*  3        */
    rvec_sub(xk, xj, r_kj); /*  3		*/
    rvec_sub(xk, xl, r_kl); /*  3		*/

    cprod(r_ij, r_kj, m);                  /*  9        */
    cprod(r_kj, r_kl, n);                  /*  9		*/
    real phi  = gmx_angle(m, n);           /* 49 (assuming 25 for atan2) */
    real ipr  = iprod(r_ij, n);            /*  5        */
    real sign = (ipr < 0.0) ? -1.0 : 1.0;
    phi       = sign*phi;                  /*  1		*/
    /* 82 TOTAL	*/
    return phi;
}

//! Fourier dihedral for one dihedral
template <typename I, typename X, typename F>
static ACT_ALWAYS_INLINE void fourDihOne(gmx_unused MsgHandler *msghandler,
                              const TopologyEntry   &a,
                              const I               &indices,
                              const X               &x,
                              F                     *forces,
                              double                *e)
{
    // Get the parameters. We have to know their names to do this.
    auto &params  = a.params();
    // Assume we have the parameters in the correct order 0-nparam-1
    size_t nparam = params.size();
    // Get the atom indices
    auto ai       = indices[0];
    auto aj       = indices[1];
    auto ak       = indices[2];
    auto al       = indices[3];

    rvec r_ij, r_kj, r_kl, m, n;
    auto phi = M_PI + dih_angle(x[ai], x[aj], x[ak], x[al],
                                r_ij, r_kj, r_kl, m, n);

    real cos_phi = std::cos(phi);

    // Energy is
    // \sum_{j=0}^{nparam} params[j] cos^j(phi)
    // Force then becomes
    // - \sum_{j=1}^{nparam} j * params[j] cos^{j-1}(phi) sin(phi)
    real cosfac = 1;
    real v      = params[0];
    real dvdphi = 0;
    // The loop has a fixed length, such that it can be unrolled
    // and the loop over conformations can be vectorized. The potential
    // has fdihNR parameters, see potential.cpp.
#pragma unroll
    for (size_t j = 1; j < fdihNR; j++)
    {
        if (j < nparam)
        {
            dvdphi += j*params[j]*cosfac;
            cosfac *= cos_phi;
            v      += params[j]*cosfac;
        }
    }

    /* Beware of accuracy loss, cannot use 1-sqrt(cos^2) ! */
    dvdphi *= -std::sin(phi);

    do_dih_fup_noshiftf(ai, aj, ak, al, dvdphi, r_ij, r_kj, r_kl, m, n,
                        forces);
    e[0] += v;
}

//! Harmonic improper dihedral for one dihedral
template <typename I, typename X, typename F>
static ACT_ALWAYS_INLINE void improperOne(gmx_unused MsgHandler *msghandler,
                               const TopologyEntry   &a,
                               const I               &indices,
                               const X               &x,
                               F                     *forces,
                               double                *e)
{
    const  real half = 0.5;
    // Get the parameters. We have to know their names to do this.
    auto &params  = a.params();
    auto kA       = params[idihKPHI];
    // Get the atom indices
    auto ai       = indices[0];
    auto aj       = indices[1];
    auto ak       = indices[2];
    auto al       = indices[3];

    rvec r_ij, r_kj, r_kl, m, n;
    auto phi = dih_angle(x[ai], x[aj], x[ak], x[al],
                         r_ij, r_kj, r_kl, m, n);

    auto dp2 = phi*phi;

    e[0]       += half*kA*dp2;
    auto ddphi  = kA*phi;

    do_dih_fup_noshiftf(ai, aj, ak, al, ddphi, r_ij, r_kj, r_kl, m, n,
                        forces);
}

//! Proper dihedral for one dihedral
template <typename I, typename X, typename F>
static ACT_ALWAYS_INLINE void properOne(gmx_unused MsgHandler *msghandler,
                             const TopologyEntry   &a,
                             const I               &indices,
                             const X               &x,
                             F                     *forces,
                             double                *e)
{
    // Get the parameters. We have to know their names to do this.
    auto &params  = a.params();
    auto phi0     = params[pdihANGLE];
    auto kphi     = params[pdihKP];
    auto mult     = params[pdihMULT];
    // Get the atom indices
    auto ai       = indices[0];
    auto aj       = indices[1];
    auto ak       = indices[2];
    auto al       = indices[3];

    rvec r_ij, r_kj, r_kl, m, n;
    auto phi = dih_angle(x[ai], x[aj], x[ak], x[al],
                         r_ij, r_kj, r_kl, m, n);

    auto mdphi  =  mult*phi - phi0;
    auto sdphi  =  std::sin(mdphi);
    auto ddphi  = -kphi*mult*sdphi;
    auto v      = kphi*(1.0 + std::cos(mdphi));
    e[0]       += v;
#ifdef DEBUG
    if (msghandler && msghandler->debug())
    {
        msghandler->writeDebug(gmx::formatString("ai %d aj %d ak %d al %d phi %g phi0 %g kphi %g, mult %g, mdphi %g, v %g ddphi %g",
                                                 ai, aj, ak, al, phi, phi0, kphi, mult, mdphi, v, ddphi));
    }
#endif
    do_dih_fup_noshiftf(ai, aj, ak, al, ddphi, r_ij, r_kj, r_kl, m, n,
                        forces);
}

/*! \brief Two functions computing one interaction
 * The first one is used for interactions for which use returns true,
 * the second one otherwise. This allows a cheaper function to be used
 * for some interactions, deciding once per interaction rather than in
 * the loop over conformations.
 */
template <typename Use, typename One1, typename One2>
struct OneSelect
{
    //! Whether to use one1 for an interaction
    Use  use;
    //! The function used when use is true
    One1 one1;
    //! The function used otherwise
    One2 one2;
};

//! Check whether a type is a OneSelect
template <typename T>
struct IsOneSelect : std::false_type {};

//! Check whether a type is a OneSelect
template <typename Use, typename One1, typename One2>
struct IsOneSelect<OneSelect<Use, One1, One2>> : std::true_type {};

/*! \brief Call body with the function computing an interaction
 * \param[in] one   Function computing one interaction or a OneSelect
 * \param[in] entry The interaction
 * \param[in] body  Function to call with the selected function
 */
template <typename One, typename Body>
static ACT_ALWAYS_INLINE void withOne(const One           &one,
                                      const TopologyEntry &entry,
                                      Body                 body)
{
    if constexpr (IsOneSelect<One>::value)
    {
        if (one.use(entry))
        {
            body(one.one1);
        }
        else
        {
            body(one.one2);
        }
    }
    else
    {
        body(one);
    }
}

/*! \brief Compute one potential for one conformation
 * \param[in]    msghandler  For warnings and errors
 * \param[in]    entries     The atom identifiers and parameters
 * \param[in]    x           The coordinates of all particles
 * \param[inout] forces      The forces on all particles
 * \param[out]   energies    The energy per type
 * \param[in]    one         Function computing one interaction
 */
template <typename One>
static inline void computeConformation(MsgHandler                   *msghandler,
                                       const TopologyEntryVector    &entries,
                                       const std::vector<gmx::RVec> &x,
                                       std::vector<gmx::RVec>       *forces,
                                       double                        energies[MAX_ENER],
                                       One                           one)
{
    double e[MAX_ENER] = { 0, 0, 0 };
    for (const auto &b : entries)
    {
        const TopologyEntry &entry = *b.operator->();
        withOne(one, entry, [&](const auto &o)
                {
                    o(msghandler, entry, entry.atomIndices(), x, forces, e);
                });
    }
    for (int k = 0; k < MAX_ENER; k++)
    {
        energies[k] = e[k];
    }
}

/*! \brief Wrapper making a function object out of an xxxOne template
 * The name of the template is the argument of the macro.
 */
#define ONE_INTERACTION(func) \
    [](MsgHandler *mh, const TopologyEntry &b, const auto &idx, const auto &x, auto *f, double *e) \
    { func(mh, b, idx, x, f, e); }

//! As ONE_INTERACTION, for functions that also need the atoms and epsilonr
#define ONE_INTERACTION_ATOMS_EPSR(func) \
    [&atoms, epsilonr](MsgHandler *mh, const TopologyEntry &b, const auto &idx, const auto &x, auto *f, double *e) \
    { func(mh, b, idx, x, f, e, atoms, epsilonr); }

//! As ONE_INTERACTION, for functions that also need the atoms
#define ONE_INTERACTION_ATOMS(func) \
    [&atoms](MsgHandler *mh, const TopologyEntry &b, const auto &idx, const auto &x, auto *f, double *e) \
    { func(mh, b, idx, x, f, e, atoms); }


/*! \brief Return the interaction types corresponding to the energies of a potential
 * \param[in] mypot The potential
 * \return The interaction type for each energy term, at most MAX_ENER
 */
static const std::vector<InteractionType> &energyTypes(Potential mypot)
{
    // Static, since this is called for every potential in every energy calculation
    static const std::vector<InteractionType> none;
    static const std::vector<InteractionType> bonds        = { InteractionType::BONDS };
    static const std::vector<InteractionType> angles       = { InteractionType::ANGLES };
    static const std::vector<InteractionType> linearAngles = { InteractionType::LINEAR_ANGLES };
    static const std::vector<InteractionType> exchDisp     = { InteractionType::EXCHANGE, InteractionType::DISPERSION };
    static const std::vector<InteractionType> exchDispInd  = { InteractionType::EXCHANGE, InteractionType::DISPERSION, InteractionType::INDUCTION };
    static const std::vector<InteractionType> exchange     = { InteractionType::EXCHANGE };
    static const std::vector<InteractionType> elec         = { InteractionType::ELECTROSTATICS };
    static const std::vector<InteractionType> polarization = { InteractionType::POLARIZATION };
    static const std::vector<InteractionType> improper     = { InteractionType::IMPROPER_DIHEDRALS };
    static const std::vector<InteractionType> proper       = { InteractionType::PROPER_DIHEDRALS };
    static const std::vector<InteractionType> posres       = { InteractionType::POSITION_RESTRAINT };
    switch(mypot)
    {
    case Potential::HARMONIC_BONDS:
    case Potential::MORSE_BONDS:
    case Potential::HUA_BONDS:
    case Potential::CUBIC_BONDS:
        return bonds;
    case Potential::HARMONIC_ANGLES:
    case Potential::UREY_BRADLEY_ANGLES:
        return angles;
    case Potential::LINEAR_ANGLES:
        return linearAngles;
    case Potential::LJ12_6:
    case Potential::LJ8_6:
    case Potential::LJ14_7:
    case Potential::BUCKINGHAM:
    case Potential::WANG_BUCKINGHAM:
    case Potential::TANG_TOENNIES:
    case Potential::TT2b:
    case Potential::SLATER_ISA_TT:
    case Potential::GENERALIZED_BUCKINGHAM:
        return exchDisp;
    case Potential::LJ12_6_4:
        return exchDispInd;
    case Potential::BORN_MAYER:
    case Potential::SLATER_ISA:
    case Potential::MACDANIEL_SCHMIDT:
        return exchange;
    case Potential::COULOMB_POINT:
    case Potential::COULOMB_GAUSSIAN:
    case Potential::COULOMB_SLATER:
        return elec;
    case Potential::POLARIZATION:
        return polarization;
    case Potential::HARMONIC_DIHEDRALS:
        return improper;
    case Potential::PROPER_DIHEDRALS:
    case Potential::FOURIER_DIHEDRALS:
        return proper;
    case Potential::POSITION_RESTRAINT:
        return posres;
    default:
        return none;
    }
}

/*! \brief Store the energies of one potential and add them to the potential energy
 * \param[in]    mypot      The potential
 * \param[in]    coords     The conformations
 * \param[inout] cener      The energies of all conformations
 * \param[in]    energies   The energies, energy k of conformation i is stored
 *                          at energies[i*confStride + k*termStride]
 * \param[in]    confStride Stride between conformations
 * \param[in]    termStride Stride between energy terms
 */
static void storeEnergies(Potential                        mypot,
                          const std::vector<CoordsEnergy> &coords,
                          ConformationEnergies            *cener,
                          const double                    *energies,
                          size_t                           confStride,
                          size_t                           termStride)
{
    const auto &itypes = energyTypes(mypot);
    std::array<ConformationEnergies::Term, MAX_ENER> terms;
    for(size_t k = 0; k < itypes.size(); k++)
    {
        terms[k] = cener->add(itypes[k]);
    }
    auto eEpot = cener->add(InteractionType::EPOT);
    for(size_t i = 0; i < coords.size(); i++)
    {
        if (coords[i].done)
        {
            continue;
        }
        const double *myener = energies + i*confStride;
        for(size_t k = 0; k < itypes.size(); k++)
        {
            terms[k].set(i, myener[k*termStride]);
        }
        double epot = myener[0] + myener[termStride] + myener[2*termStride];
        if (!eEpot.has(i))
        {
            eEpot.set(i, 0.0);
        }
        eEpot.set(i, epot + eEpot.act[i]);
    }
}

/*! \brief Call a function for each potential
 * The function is called with the number of atoms per interaction,
 * whether the loop over conformations can be vectorized and a function
 * object computing one interaction, see ONE_INTERACTION.
 * \param[in] mypot    The potential
 * \param[in] atoms    The atoms
 * \param[in] epsilonr Relative dielectric constant
 * \param[in] func     The function
 * \return false if the potential is not implemented
 */
template <typename Func>
static bool dispatchPotential(Potential                   mypot,
                              const std::vector<ActAtom> &atoms,
                              double                      epsilonr,
                              Func                        func)
{
    switch(mypot)
    {
    case Potential::HARMONIC_BONDS:
        func.template operator()<2, true>(ONE_INTERACTION(bondOne));
        break;
    case Potential::MORSE_BONDS:
        func.template operator()<2, true>(ONE_INTERACTION(morseOne));
        break;
    case Potential::HUA_BONDS:
        func.template operator()<2, true>(ONE_INTERACTION(huaOne));
        break;
    case Potential::CUBIC_BONDS:
        func.template operator()<2, true>(ONE_INTERACTION(cubicOne));
        break;
    case Potential::HARMONIC_ANGLES:
        func.template operator()<3, true>(ONE_INTERACTION(angleOne));
        break;
    case Potential::LINEAR_ANGLES:
        func.template operator()<3, true>(ONE_INTERACTION(linearAngleOne));
        break;
    case Potential::UREY_BRADLEY_ANGLES:
        func.template operator()<3, true>(ONE_INTERACTION(ureyBradleyOne));
        break;
    case Potential::LJ12_6:
        func.template operator()<2, true>(ONE_INTERACTION(lj12_6One));
        break;
    case Potential::LJ12_6_4:
        func.template operator()<2, true>(ONE_INTERACTION(lj12_6_4One));
        break;
    case Potential::LJ8_6:
        func.template operator()<2, true>(ONE_INTERACTION(lj8_6One));
        break;
    case Potential::LJ14_7:
        func.template operator()<2, true>(ONE_INTERACTION(lj14_7One));
        break;
    case Potential::BUCKINGHAM:
        func.template operator()<2, true>(ONE_INTERACTION(buckinghamOne));
        break;
    case Potential::WANG_BUCKINGHAM:
        func.template operator()<2, true>(ONE_INTERACTION(wbhOne));
        break;
    case Potential::TANG_TOENNIES:
        func.template operator()<2, true>(ONE_INTERACTION(tangToenniesOne));
        break;
    case Potential::TT2b:
        func.template operator()<2, true>(ONE_INTERACTION(tt2bOne));
        break;
    case Potential::SLATER_ISA_TT:
        func.template operator()<2, true>(ONE_INTERACTION(slaterISATTOne));
        break;
    case Potential::GENERALIZED_BUCKINGHAM:
        func.template operator()<2, true>(ONE_INTERACTION(generalizedBhamOne));
        break;
    case Potential::BORN_MAYER:
        func.template operator()<2, true>(ONE_INTERACTION(bornMayerOne));
        break;
    case Potential::SLATER_ISA:
        func.template operator()<2, true>(ONE_INTERACTION(slaterISAOne));
        break;
    case Potential::MACDANIEL_SCHMIDT:
        func.template operator()<2, true>(ONE_INTERACTION(doubleExponentialOne));
        break;
    case Potential::COULOMB_POINT:
    case Potential::COULOMB_GAUSSIAN:
        // Point charges are cheaper, in particular when vectorized
        func.template operator()<2, true>(OneSelect { [](const TopologyEntry &b)
                                                      { return b.params()[coulZETA] == 0 && b.params()[coulZETA2] == 0; },
                                                      ONE_INTERACTION_ATOMS_EPSR(coulombPointOne),
                                                      ONE_INTERACTION_ATOMS_EPSR(coulombGaussianOne) });
        break;
    case Potential::COULOMB_SLATER:
        // The Slater integrals are not inlined, so this cannot be vectorized
        func.template operator()<2, false>(ONE_INTERACTION_ATOMS_EPSR(coulombSlaterOne));
        break;
    case Potential::POLARIZATION:
        func.template operator()<2, true>(ONE_INTERACTION_ATOMS(polarizationOne));
        break;
    case Potential::HARMONIC_DIHEDRALS:
        func.template operator()<4, true>(ONE_INTERACTION(improperOne));
        break;
    case Potential::PROPER_DIHEDRALS:
        func.template operator()<4, true>(ONE_INTERACTION(properOne));
        break;
    case Potential::FOURIER_DIHEDRALS:
        func.template operator()<4, true>(ONE_INTERACTION(fourDihOne));
        break;
    case Potential::POSITION_RESTRAINT:
        func.template operator()<1, true>(ONE_INTERACTION(fbposreOne));
        break;
    default:
        return false;
    }
    return true;
}

static inline void calcSome(MsgHandler                 *msghandler,
                            Potential                   mypot,
                            const TopologyEntryVector  &topentries,
                            const std::vector<ActAtom> &atoms,
                            std::vector<CoordsEnergy>  *coords,
                            ConformationEnergies       *cener,
                            double                      epsilonr)
{
    // Energies for each conformation, MAX_ENER per conformation
    std::vector<double> energies(coords->size()*MAX_ENER, 0.0);
    // Loop over all coordinate arrays.
    for(size_t i = 0; i < coords->size(); i++)
    {
        // This is related to shell relaxation
        if ((*coords)[i].done)
        {
            continue;
        }
        // Energies of this conformation
        double *myener    = energies.data() + i*MAX_ENER;
        // Shortcuts
        const auto &x     = (*coords)[i].x;
        // The forces are cleared in lowComputeEnergies, since all potentials add to them
        auto       &fcalc = (*coords)[i].fcalc;
        dispatchPotential(mypot, atoms, epsilonr,
                          [&]<int NATOM, bool Vectorize>(auto one)
                          {
                              computeConformation(msghandler, topentries, x, &fcalc, myener, one);
                          });
    }
    // Store the energies and sum them into the potential energy.
    storeEnergies(mypot, *coords, cener, energies.data(), MAX_ENER, 1);
}

/*! \brief Compute one potential for all conformations
 *
 * For each interaction, the coordinates and forces of its atoms in one
 * conformation are copied to small local arrays, the interaction is
 * computed by the same function as for a single conformation, and the
 * forces and energies are copied back. Since the conformations are
 * independent, the loop over them can be vectorized, provided the
 * compiler can inline everything.
 * Without forces, the forces are neither loaded nor stored, and the
 * compiler removes their computation since they are not used.
 * \tparam       NATOM     Number of atoms per interaction
 * \tparam       Vectorize Whether to ask the compiler to vectorize
 * \tparam       Forces    Whether to compute forces
 * \param[in]    entries   The atom identifiers and parameters
 * \param[inout] ca        Coordinates and forces of all conformations
 * \param[inout] energies  Energies, term k of conformation c at energies[k*stride + c]
 * \param[in]    one       Function computing one interaction
 */
template <int NATOM, bool Vectorize, bool Forces, typename One>
static void computeConformations(const TopologyEntryVector &entries,
                                 ConformationArrays        *ca,
                                 double                    *energies,
                                 One                        one)
{
    // Indices into the local arrays
    constexpr std::array<int, NATOM> lidx = []()
    {
        std::array<int, NATOM> l{};
        for (int a = 0; a < NATOM; a++)
        {
            l[a] = a;
        }
        return l;
    }();
    // The padding conformations are computed as well, which is simpler
    size_t nconf = ca->stride;
    for (const auto &b : entries)
    {
        const TopologyEntry &entry = *b.operator->();
        auto  &indices = entry.atomIndices();
        const double *xa[NATOM];
        double       *fa[NATOM];
        for (int a = 0; a < NATOM; a++)
        {
            xa[a] = ca->x.data() + DIM*indices[a]*nconf;
            fa[a] = Forces ? ca->f.data() + DIM*indices[a]*nconf : nullptr;
        }
        withOne(one, entry, [&](const auto &o)
        {
            auto lane = [&](size_t c)
            {
                std::array<gmx::RVec, NATOM> xl, fl;
                for (int a = 0; a < NATOM; a++)
                {
                    for (int m = 0; m < DIM; m++)
                    {
                        xl[a][m] = xa[a][m*nconf+c];
                        fl[a][m] = Forces ? fa[a][m*nconf+c] : 0.0;
                    }
                }
                double el[MAX_ENER];
                for (int k = 0; k < MAX_ENER; k++)
                {
                    el[k] = energies[k*nconf+c];
                }
                o(nullptr, entry, lidx, xl, &fl, el);
                if constexpr (Forces)
                {
                    for (int a = 0; a < NATOM; a++)
                    {
                        for (int m = 0; m < DIM; m++)
                        {
                            fa[a][m*nconf+c] = fl[a][m];
                        }
                    }
                }
                for (int k = 0; k < MAX_ENER; k++)
                {
                    energies[k*nconf+c] = el[k];
                }
            };
            if constexpr (Vectorize)
            {
#pragma omp simd
                for(size_t c = 0; c < nconf; c++)
                {
                    lane(c);
                }
            }
            else
            {
                for(size_t c = 0; c < nconf; c++)
                {
                    lane(c);
                }
            }
        });
    }
}

/*! \brief Compute a potential for all conformations without explicit SIMD
 * The versions with and without forces are kept in separate functions,
 * since the compiler merges their loops into loops it cannot vectorize
 * otherwise.
 * \tparam       Forces     Whether to compute forces
 * \param[in]    mypot      The potential
 * \param[in]    topentries The interactions
 * \param[in]    atoms      The atoms
 * \param[inout] ca         Coordinates and forces of all conformations
 * \param[inout] energies   Energies, term k of conformation c at energies[k*stride + c]
 * \param[in]    epsilonr   Relative dielectric constant
 * \return false if the potential is not implemented
 */
template <bool Forces>
static ACT_NOINLINE bool dispatchConformations(Potential                   mypot,
                                               const TopologyEntryVector  &topentries,
                                               const std::vector<ActAtom> &atoms,
                                               ConformationArrays         *ca,
                                               double                     *energies,
                                               double                      epsilonr)
{
    return dispatchPotential(mypot, atoms, epsilonr,
                             [&]<int NATOM, bool Vectorize>(auto one)
                             {
                                 computeConformations<NATOM, Vectorize, Forces>(topentries, ca, energies, one);
                             });
}

/*! \brief Compute a potential with conformations in the inner loop
 * \param[in]    mypot      The potential
 * \param[in]    topentries The interactions
 * \param[in]    atoms      The atoms
 * \param[in]    coords     The conformations
 * \param[inout] ca         Coordinates and forces of all conformations
 * \param[inout] cener      The energies of all conformations
 * \param[in]    epsilonr   Relative dielectric constant
 * \param[in]    forces     Whether to compute forces
 * \return false if this potential is not implemented this way
 */
static bool calcConformations(Potential                        mypot,
                              const TopologyEntryVector       &topentries,
                              const std::vector<ActAtom>      &atoms,
                              const std::vector<CoordsEnergy> &coords,
                              ConformationArrays              *ca,
                              ConformationEnergies            *cener,
                              double                           epsilonr,
                              bool                             forces)
{
    // Energy term k of conformation c is stored at k*stride + c
    auto &energies = ca->energies;
    energies.assign(MAX_ENER*ca->stride, 0.0);
    // Use explicit SIMD where available
    if (!(useExplicitSimd() &&
          computeExplicitSimd(mypot, topentries, atoms, epsilonr, ca->stride,
                              ca->x.data(), ca->f.data(), energies.data(), forces)) &&
        !(forces ?
          dispatchConformations<true>(mypot, topentries, atoms, ca, energies.data(), epsilonr) :
          dispatchConformations<false>(mypot, topentries, atoms, ca, energies.data(), epsilonr)))
    {
        return false;
    }
    storeEnergies(mypot, coords, cener, energies.data(), 1, ca->stride);
    return true;
}

/*! \brief Compute all potentials for all conformations of a topology
 * \param[in]    msghandler    For warnings and errors
 * \param[in]    pd            The force field
 * \param[in]    top           The topology
 * \param[inout] coords        The conformations, the forces are stored here
 * \param[inout] cdata         Coordinates, forces and energies of all
 *                             conformations, the energies are stored here
 * \param[in]    epsilonr      Relative dielectric constant
 * \param[in]    conformations Whether to loop over the conformations in the inner loop
 * \param[in]    forces        Whether to compute forces, if not the forces
 *                             in coords are not meaningful
 */
static inline void lowComputeEnergies(MsgHandler                *msghandler,
                                      const ForceField          *pd,
                                      const Topology            *top,
                                      std::vector<CoordsEnergy> *coords,
                                      ConformationData          *cdata,
                                      double                     epsilonr,
                                      bool                       conformations,
                                      bool                       forces)
{
    auto &cener = cdata->energies;
    cener.setNumberOfConformations(coords->size());
    // Reset the accumulated total potential energy (EPOT) for every conformation
    // before summing the contributions of the potentials. Without this,
    // repeated energy calculations (for example the many evaluations made
    // while training with the MCMC mutator) keep adding onto the value left by
    // the previous call, so the total grows without bound instead of being
    // recomputed from scratch each time.
    auto eEpot = cener.add(InteractionType::EPOT);
    for(size_t c = 0; c < coords->size(); c++)
    {
        eEpot.set(c, 0.0);
    }
    for(auto &ce : *coords)
    {
        // Also reset the done flag to make sure energies are computed once at least
        ce.done = false;
        // Clear the forces, calcSome adds the forces of each potential to them.
        // Without forces calcSome still needs the memory, but not the values.
        if (forces)
        {
            gmx::RVec fzero = { 0, 0, 0 };
            ce.fcalc.assign(ce.x.size(), fzero);
        }
        else
        {
            ce.fcalc.resize(ce.x.size());
        }
    }
    // Coordinates and forces of all conformations, for looping over
    // conformations in the inner loop. Kept between calls, since only the
    // coordinates of virtual sites and shells change.
    auto &ca = cdata->arrays;
    // With a single conformation there is nothing to vectorize.
    bool haveArrays     = (conformations && coords->size() > 1 &&
                           ca.setCoordinates(*coords, top->atoms(), forces));
    // Whether the current forces are in ca rather than in coords
    bool forcesInArrays = haveArrays;
    for(const auto &entry : top->entries())
    {
        //! \todo Can we prevent this check?
        if (entry.second.empty())
        {
            continue;
        }
        // Not all potentials are available in the force field,
        // in particular restraint potentials.
        gmx_unused Potential mypot;
        if (pd->interactionPresent(entry.first))
        {
            // Force field parameter list
            auto &ffpl = pd->findForcesConst(entry.first);
            // The function we need to do the math
            mypot = ffpl.potential();
        }
        else
        {
            mypot = defaultPotential(entry.first);
        }
        if (haveArrays)
        {
            if (forces && !forcesInArrays)
            {
                ca.loadForces(*coords);
                forcesInArrays = true;
            }
            if (calcConformations(mypot, entry.second, top->atoms(),
                                  *coords, &ca, &cener, epsilonr, forces))
            {
                continue;
            }
            if (forces)
            {
                ca.storeForces(coords);
            }
            forcesInArrays = false;
        }
        calcSome(msghandler, mypot, entry.second, top->atoms(), coords, &cener, epsilonr);
    }
    if (forces && forcesInArrays)
    {
        ca.storeForces(coords);
    }
}

void computeEnergies(MsgHandler                  *msghandler,
                     const ForceField            *pd,
                     const std::vector<ACTMol *> &allmols,
                     size_t                       what,
                     bool                         forces)
{
    double epsilonr = 1;
    auto iCoul = InteractionType::ELECTROSTATICS;
    if (pd->interactionPresent(iCoul))
    {
        auto &ffpl = pd->findForcesConst(iCoul);
        std::string epsr = "epsilonr";
        if (ffpl.optionExists(epsr))
        {
            epsilonr = std::atof(ffpl.optionValue(epsr).c_str());
        }
    }
    // Loop over molecules
    for(auto molptr : allmols)
    {
        auto &mol   = *molptr;
        auto &tops  = mol.fragmentHandler()->topologies();
        // Loop over topology entries, in principle the dimer should have everything
        // that the monomers have as well.
        //! \todo Check whether dimers and monomers have the same entries
        if (tops.size() == 1 && what > 0)
        {
            continue;
        }
        if (what >= 2)
        {
            lowComputeEnergies(msghandler, pd, mol.topology(), mol.dimerCoords(),
                               mol.conformationData(2), epsilonr, true, forces);
        }
        else
        {
            lowComputeEnergies(msghandler, pd, &tops[what], mol.monomerCoords(what),
                               mol.conformationData(what), epsilonr, true, forces);
        }
    }
}

void computePotential(MsgHandler                 *msghandler,
                      Potential                   mypot,
                      const TopologyEntryVector  &topentries,
                      const std::vector<ActAtom> &atoms,
                      std::vector<CoordsEnergy>  *coords,
                      double                      epsilonr,
                      bool                        conformations)
{
    for(auto &ce : *coords)
    {
        ce.eMap.clear();
        ce.done = false;
        gmx::RVec fzero = { 0, 0, 0 };
        ce.fcalc.assign(ce.x.size(), fzero);
    }
    ConformationData cdata;
    cdata.energies.setNumberOfConformations(coords->size());
    if (conformations && coords->size() > 1 && cdata.arrays.setCoordinates(*coords, atoms) &&
        calcConformations(mypot, topentries, atoms, *coords, &cdata.arrays, &cdata.energies, epsilonr, true))
    {
        cdata.arrays.storeForces(coords);
    }
    else
    {
        calcSome(msghandler, mypot, topentries, atoms, coords, &cdata.energies, epsilonr);
    }
    cdata.energies.store(coords);
}

} // namespace

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
#include <limits>
#include <optional>
#include <set>
#include <string>
#include <type_traits>

#include "act/basics/msg_handler.h"
#include "act/coulombintegrals/gaussian_integrals.h"
#include "act/coulombintegrals/slater_integrals.h"
#include "act/forces/conformationdata.h"
#include "act/forcefield/forcefield_parametername.h"
#include "act/forcefield/parameterrecipes.h"
#include "act/forces/forces_explicit_simd.h"
#include "act/forces/forcecomputer.h"
#include "act/forces/forces_gpu.h"
#include "act/forces/gpusystem.h"
#include "act/forces/simd_information.h"
#include "act/forces/forcecomputerutils.h"
#include "act/forces/vsitehandler.h"
#include "act/molprop/topologyentry.h"
#include "act/topology/actmol.h"
#include "act/topology/coords_energy.h"
#include "act/topology/topology.h"
#include "act/utility/timings.h"
#include "gromacs/math/vec.h"
#include "gromacs/math/units.h"
#include "gromacs/utility/basedefinitions.h"

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
 * The force is computed in the precision of the scalar force, such that
 * potentials in single precision do not convert to double precision.
 * \param[in]    fscalar The scalar force
 * \param[in]    dx      The distance vector
 * \param[in]    indices The atom indices
 * \param[inout] forces  The forces to be updated
 */
template <typename I, typename F, typename D, typename S>
static ACT_ALWAYS_INLINE void pairforces(S                fscalar,
                              const D         &dx,
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
    // The distance from the center of the restraints
    rvec dx      = { static_cast<real>(x[ai][XX] - params[fbprX0]),
                     static_cast<real>(x[ai][YY] - params[fbprY0]),
                     static_cast<real>(x[ai][ZZ] - params[fbprZ0]) };
    auto dr2     = iprod(dx, dx);
    if (dr2 > R0*R0)
    {
        auto r1   = std::sqrt(dr2);
        auto dr   = r1-R0;
        e[0]     += 0.5*K*dr*dr;
        // F = -dV/dr dx/r
        auto fpr  = -K*dr/r1;
        for (int m = 0; m < DIM; m++)
        {
            (*forces)[ai][m] += fpr*dx[m];
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

//! Whether a range of interactions holds non-bonded pairs, see PairRange
template <typename Range>
constexpr bool isPairRange = std::is_same_v<std::remove_cvref_t<Range>, PairRange>;

//! Lennard-Jones 12-6 for one pair
template <typename E, typename I, typename X, typename F>
static ACT_ALWAYS_INLINE void lj12_6One(gmx_unused MsgHandler *msghandler,
                             const E               &b,
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
template <typename E, typename I, typename X, typename F>
static ACT_ALWAYS_INLINE void lj12_6_4One(gmx_unused MsgHandler *msghandler,
                               const E               &b,
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
template <typename E, typename I, typename X, typename F>
static ACT_ALWAYS_INLINE void lj8_6One(gmx_unused MsgHandler *msghandler,
                            const E               &b,
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

//! Return x^7, in the precision of x
template <typename T>
static ACT_ALWAYS_INLINE T power7(T x)
{
    return x*gmx::power6(x);
}

//! Buffered 14-7 for one pair
template <typename E, typename I, typename X, typename F>
static ACT_ALWAYS_INLINE void lj14_7One(gmx_unused MsgHandler *msghandler,
                             const E               &b,
                             const I               &indices,
                             const X               &x,
                             F                     *forces,
                             double                *e)
{
    // Get the parameters. We have to know their names to do this. They are
    // converted to real, such that single precision computes in single
    // precision only.
    auto &params    = b.params();
    real sigma      = params[lj14_7SIGMA];
    real epsilon    = params[lj14_7EPSILON];
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
        real gamma      = params[lj14_7GAMMA];
        real delta      = params[lj14_7DELTA];
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
template <typename E, typename I, typename X, typename F>
static ACT_ALWAYS_INLINE void bornMayerOne(MsgHandler          *msghandler,
                                const E             &b,
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
template <typename E, typename I, typename X, typename F>
static ACT_ALWAYS_INLINE void slaterISAOne(MsgHandler          *msghandler,
                                const E             &b,
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
template <typename E, typename I, typename X, typename F>
static ACT_ALWAYS_INLINE void doubleExponentialOne(gmx_unused MsgHandler *msghandler,
                                        const E               &b,
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
template <typename E, typename I, typename X, typename F>
static ACT_ALWAYS_INLINE void wbhOne(gmx_unused MsgHandler *msghandler,
                          const E               &b,
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
template <typename E, typename I, typename X, typename F>
static ACT_ALWAYS_INLINE void buckinghamOne(gmx_unused MsgHandler *msghandler,
                                 const E               &b,
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
    // The sum is accumulated in fk, with terms[k] holding (br)^k/k!
    real terms[ttTerms+1];
    terms[0] = 1;
#pragma unroll
    for (int k = 0; k < ttTerms; k++)
    {
        terms[k+1] = terms[k]*bDispr/(k+1);
    }
    // For small br, 1 - exp(-br) sum_k=0^2m is the difference of nearly
    // equal numbers. There f_2m = exp(-br) sum_k>2m (br)^k/k! is used,
    // with the terms beyond ttTerms summed by Horner's rule in tail.
    real tail = ttTailCoefficients.back();
#pragma unroll
    for (int j = ttTailTerms-ttTerms-2; j >= 0; j--)
    {
        tail = tail*bDispr + ttTailCoefficients[j];
    }
    tail *= terms[ttTerms]*bDispr/(ttTerms+1);
    bool useTail = bDispr < ttTailMax;
    real fk      = 1;
    int  k       = 0;
#pragma unroll
    for (int m = 3; m <= 5; m++)
    {
        int  nn  = 2*m;
#pragma unroll
        for (; k < nn; k++)
        {
            fk += terms[k+1];
        }
        real tk = tail;
#pragma unroll
        for (int j = ttTerms; j > nn; j--)
        {
            tk += terms[j];
        }
        real damp = useTail ? ebrDisp*tk : 1-ebrDisp*fk;
        real ed   = -damp*ctt[m-3]*rinvn;
        eedisp += ed;
        fdisp  += ed*nn*rinv1;
        // The derivative of the sum is b times the sum up to 2m-1,
        // the terms up to 2m-1 cancel against the derivative of exp(-br).
        fdisp  += ctt[m-3]*rinvn*bDisp*ebrDisp*terms[nn];
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
template <typename E, typename I, typename X, typename F>
static ACT_ALWAYS_INLINE void tangToenniesOne(MsgHandler          *msghandler,
                                   const E             &b,
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
template <typename E, typename I, typename X, typename F>
static ACT_ALWAYS_INLINE void tt2bOne(MsgHandler          *msghandler,
                           const E             &b,
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
template <typename E, typename I, typename X, typename F>
static ACT_ALWAYS_INLINE void slaterISATTOne(MsgHandler          *msghandler,
                                  const E             &b,
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
template <typename E, typename I, typename X, typename F>
static ACT_ALWAYS_INLINE void generalizedBhamOne(gmx_unused MsgHandler *msghandler,
                                      const E               &b,
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
/*! \brief Width of the interaction of two Gaussian charge distributions
 * As in Coulomb_GG, a zero width means a point charge.
 * \param[in] izeta Width of the first distribution
 * \param[in] jzeta Width of the second distribution
 */
static ACT_ALWAYS_INLINE real combinedZeta(real izeta, real jzeta)
{
    if (izeta == 0 || jzeta == 0)
    {
        return izeta + jzeta;
    }
    return izeta*jzeta*gmx::invsqrt(izeta*izeta + jzeta*jzeta);
}

/*! \brief Scalar force divided by r of two Gaussian charges, without charges
 * This is (erf(x)/r - 2 zeta exp(-x^2)/sqrt(pi))/r^2 with x = zeta r, in
 * which the terms are nearly equal for small x, such that the difference
 * loses all digits in single precision. Below x^2 = smallX2 the series
 * zeta^3 (2/sqrt(pi)) (2/3 - 2x^2/5 + x^4/7 - x^6/27 + x^8/132) is used,
 * where it is more accurate than the difference.
 * \param[in] zeta  Width of the combined distribution
 * \param[in] rinv  Inverse distance
 * \param[in] x     zeta times the distance
 * \param[in] erfx  The error function of x
 * \return the scalar force divided by r
 */
static ACT_ALWAYS_INLINE real gaussianForceOverR(real zeta, real rinv, real x, real erfx)
{
    constexpr real smallX2 = GMX_DOUBLE ? 0.005 : 0.16;
    real x2 = x*x;
    if (x2 < smallX2)
    {
        return zeta*zeta*zeta*2*invsqrt_pi*(2.0/3.0 + x2*(-2.0/5.0 + x2*(1.0/7.0 + x2*(-1.0/27.0 + x2/132.0))));
    }
    return (erfx*rinv - 2*invsqrt_pi*zeta*std::exp(-x2))*rinv*rinv;
}

/*! \brief Coulomb between Gaussian charges for one pair
 * Each particle may have a second charge at the same position with its
 * own width, see ActAtom::chargeVS, such that there are up to four terms.
 * Terms without charge are skipped, except the first one.
 */
template <typename E, typename I, typename X, typename F>
static ACT_ALWAYS_INLINE void coulombGaussianOne(gmx_unused MsgHandler      *msghandler,
                                      const E                    &b,
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
    const auto &atomi = atoms[gindices[0]];
    const auto &atomj = atoms[gindices[1]];
    // The charges and widths of both particles. These depend only on the
    // atoms, which allows the compiler to compute them once for all
    // conformations.
    const real qi[2] = { static_cast<real>(atomi.charge()), static_cast<real>(atomi.chargeVS()) };
    const real qj[2] = { static_cast<real>(atomj.charge()), static_cast<real>(atomj.chargeVS()) };
    const real zi[2] = { static_cast<real>(atomi.zeta()), static_cast<real>(atomi.zetaVS()) };
    const real zj[2] = { static_cast<real>(atomj.zeta()), static_cast<real>(atomj.zetaVS()) };
    rvec dx;
    rvec_sub(x[ai], x[aj], dx);
    auto dr2        = iprod(dx, dx);
    real felec      = 0;
    for(int a = 0; a < 2; a++)
    {
        for(int c = 0; c < 2; c++)
        {
            if (a+c > 0 && qi[a]*qj[c] == 0)
            {
                continue;
            }
            // In real, such that single precision computes in single
            // precision only
            real qq   = static_cast<real>(ONE_4PI_EPS0)*qi[a]*qj[c]/static_cast<real>(epsilonr);
            // Width of the combined Gaussian, as in Coulomb_GG
            real zeta = combinedZeta(zi[a], zj[c]);
            if (zeta == 0)
            {
                // Two point charges, which only occurs for terms with a
                // second charge, since the point charge kernel is used otherwise.
                if (dr2 > 0)
                {
                    auto rinv = gmx::invsqrt(dr2);
                    e[0]     += qq*rinv;
                    felec    += qq*rinv*rinv*rinv;
                }
            }
            else if (dr2 > 0)
            {
                // As Coulomb_GG and DCoulomb_GG, but computing erf only once
                auto rinv   = gmx::invsqrt(dr2);
                auto rz     = zeta*dr2*rinv;
                auto erfrz  = std::erf(rz);
                e[0]       += qq*erfrz*rinv;
                // The scalar force divided by r
                felec      += qq*gaussianForceOverR(zeta, rinv, rz, erfrz);
            }
            else
            {
                e[0] += qq*Nuclear_GG(0, zeta);
            }
        }
    }
    if (dr2 > 0)
    {
        pairforces(felec, dx, indices, forces);
    }
}

/*! \brief Coulomb between point charges for one pair
 * This is the same as coulombGaussianOne with both widths zero, but the
 * branches for non-zero widths are removed by the compiler, which matters
 * for vectorization.
 */
template <typename E, typename I, typename X, typename F>
static ACT_ALWAYS_INLINE void coulombPointOne(gmx_unused MsgHandler      *msghandler,
                                              const E                    &b,
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

/*! \brief Coulomb between Slater charges for one pair
 * Each particle may have a second charge at the same position with its
 * own width and row, see ActAtom::chargeVS, such that there are up to four
 * terms. Terms without charge are skipped, except the first one.
 */
template <typename E, typename I, typename X, typename F>
static ACT_ALWAYS_INLINE void coulombSlaterOne(gmx_unused MsgHandler      *msghandler,
                                    const E                    &b,
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
    const auto &atomi = atoms[gindices[0]];
    const auto &atomj = atoms[gindices[1]];
    const real qi[2]   = { static_cast<real>(atomi.charge()), static_cast<real>(atomi.chargeVS()) };
    const real qj[2]   = { static_cast<real>(atomj.charge()), static_cast<real>(atomj.chargeVS()) };
    const real zi[2]   = { static_cast<real>(atomi.zeta()), static_cast<real>(atomi.zetaVS()) };
    const real zj[2]   = { static_cast<real>(atomj.zeta()), static_cast<real>(atomj.zetaVS()) };
    const int  rowi[2] = { atomi.row(), atomi.rowVS() };
    const int  rowj[2] = { atomj.row(), atomj.rowVS() };
    rvec dx;
    rvec_sub(x[ai], x[aj], dx);
    auto dr2   = iprod(dx, dx);
    real r1    = std::sqrt(dr2);
    real felec = 0;
    for(int a = 0; a < 2; a++)
    {
        for(int c = 0; c < 2; c++)
        {
            if (a+c > 0 && qi[a]*qj[c] == 0)
            {
                continue;
            }
            real qq = ONE_4PI_EPS0*qi[a]*qj[c]/epsilonr;
            e[0]   += qq*Coulomb_SS(r1, rowi[a], rowj[c], zi[a], zj[c]);
            // The DCoulomb_SS code returns the derivative of the energy
            // so we still need a minus sign.
            // https://github.com/dspoel/ACT/issues/369
            felec  -= qq*DCoulomb_SS(r1, rowi[a], rowj[c], zi[a], zj[c]);
        }
    }
    if (dr2 > 0)
    {
        felec *= gmx::invsqrt(dr2);
        pairforces(felec, dx, indices, forces);
    }
}

/* The bonded potentials are computed in double precision, also when
 * the coordinates and forces are stored in single precision. They
 * depend on differences between nearly equal numbers, such as r - r0,
 * which in single precision have an error of about GMX_FLOAT_EPS*r0,
 * giving a force error of about k*GMX_FLOAT_EPS*r0 even at the minimum.
 * The difference of two single precision coordinates is exact in double
 * precision. The explicit SIMD kernels compute in double precision too.
 */

//! Distance vector xi - xj in double precision
template <typename V>
static ACT_ALWAYS_INLINE void distanceVector(const V &xi, const V &xj, dvec dx)
{
#pragma unroll
    for (int m = 0; m < DIM; m++)
    {
        dx[m] = static_cast<double>(xi[m]) - static_cast<double>(xj[m]);
    }
}

//! Add a vector in double precision to a force
template <typename V>
static ACT_ALWAYS_INLINE void forceInc(V &&f, const dvec df)
{
#pragma unroll
    for (int m = 0; m < DIM; m++)
    {
        f[m] += df[m];
    }
}

//! Subtract a vector in double precision from a force
template <typename V>
static ACT_ALWAYS_INLINE void forceDec(V &&f, const dvec df)
{
#pragma unroll
    for (int m = 0; m < DIM; m++)
    {
        f[m] -= df[m];
    }
}

//! The cosine of the angle between two vectors, as cos_angle
static ACT_ALWAYS_INLINE double cosAngle(const dvec a, const dvec b)
{
    double ipab   = diprod(a, a)*diprod(b, b);
    double cosval = ipab > 0 ? diprod(a, b)*gmx::invsqrt(ipab) : 1;

    return std::max(-1.0, std::min(1.0, cosval));
}

/*! \brief Compute harmonic force and energy
 * \param[in]  k  The force constant
 * \param[in]  x0 The reference position/angle
 * \param[in]  x  The actual position/angle
 * \param[out] V  The energy
 * \param[out] F  The scalar force
 */
static ACT_ALWAYS_INLINE void harmonic(double k, double x0, double x, double *V, double *F)
{
    const double half = 0.5;

    double dx  = x-x0;
    double dx2 = dx*dx;

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

    dvec dx;
    distanceVector(x[indices[0]], x[indices[1]], dx);
    auto dr2        = diprod(dx, dx);
    auto r_1        = gmx::invsqrt(dr2);
    auto r1         = dr2*r_1;

    double vB, fbond;
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

    dvec dx;
    distanceVector(x[indices[0]], x[indices[1]], dx);
    double rhi     = (2*rmax+bondlength)/3;
    double r_1, r1 = rhi;
    {
        // Do not allow the potential to go to -infinity
        auto r2 = std::min<double>(rhi*rhi, diprod(dx, dx));
        r_1     = gmx::invsqrt(r2);
        r1      = r2*r_1;
    }
    double vB, fbond;
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
    dvec dx;
    distanceVector(x[indices[0]], x[indices[1]], dx);
    auto dr         = dnorm(dx);
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
    dvec dx;
    distanceVector(x[indices[0]], x[indices[1]], dx);
    auto dr2        = diprod(dx, dx);
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

    dvec r_ij, r_kj;
    distanceVector(x[indices[0]], x[indices[1]], r_ij);
    distanceVector(x[indices[2]], x[indices[1]], r_kj);

    double dr2 = 0;
    double b   = 1-a;
#pragma unroll
    for (int m = 0; (m < DIM); m++)
    {
        // Minus sign due to F = -dV/dr
        double dr   = -a * r_ij[m] - b * r_kj[m];
        dr2        += dr*dr;
        double f_i  = a*klin*dr;
        double f_k  = b*klin*dr;
        double f_j  = -(f_i+f_k);
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
template <typename V>
static ACT_ALWAYS_INLINE double bondAngle(const V &xi,
                                          const V &xj,
                                          const V &xk,
                                          dvec     r_ij,
                                          dvec     r_kj,
                                          double  *costh)
{
    double th;

    distanceVector(xi, xj, r_ij);
    distanceVector(xk, xj, r_kj);

    *costh = cosAngle(r_ij, r_kj);
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

    dvec r_ij, r_kj;
    auto theta = bondAngle(x[indices[0]], x[indices[1]], x[indices[2]],
                           r_ij, r_kj, &costh);

    double vA, fangle;
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
        double st, sth;
        double cik, cii, ckk;
        double nrkj2, nrij2;
        double nrkj_1, nrij_1;

        st    = fangle*gmx::invsqrt(1 - costh2);
        sth   = st*costh;

        nrij2 = diprod(r_ij, r_ij);
        nrkj2 = diprod(r_kj, r_kj);

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
    const  double half = 0.5;
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

    dvec r_ij, r_kj;
    auto theta = bondAngle(x[ai], x[aj], x[ak], r_ij, r_kj, &costh);

    // Compute deviation from the reference angle
    auto da  = theta - theta0;
    auto da2 = da*da;

    auto fangle      = -ka*da;
    e[0]            += half*ka*da2;

    // Now UB part
    dvec r_ik;
    distanceVector(x[ai], x[ak], r_ik);
    auto rik2        = diprod(r_ik, r_ik);
    auto norm_rik    = std::sqrt(rik2);

    double vUB, fUB;
    harmonic(kUB, r13, norm_rik, &vUB, &fUB);

    e[0]   += vUB;
    fUB    *= gmx::invsqrt(rik2);

    auto costh2 = gmx::square(costh);
    if (costh2 < 1)
    {
        double st, sth;
        double cik, cii, ckk;
        double nrkj2, nrij2;
        double nrkj_1, nrij_1;

        st    = fangle*gmx::invsqrt(1 - costh2);
        sth   = st*costh;

        nrij2 = diprod(r_ij, r_ij);
        nrkj2 = diprod(r_kj, r_kj);

        nrij_1 = gmx::invsqrt(nrij2);
        nrkj_1 = gmx::invsqrt(nrkj2);

        cik = st*nrij_1*nrkj_1;
        cii = sth*nrij_1*nrij_1;
        ckk = sth*nrkj_1*nrkj_1;
#pragma unroll
        for (auto m = 0; m < DIM; m++)
        {
            double f_i = -(cik*r_kj[m] - cii*r_ij[m]);
            double f_k = -(cik*r_ij[m] - ckk*r_kj[m]);
            double f_j = -f_i - f_k;
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
    const double half = 0.5;
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
    dvec dx;
    distanceVector(x[indices[0]], x[indices[1]], dx);
    auto dr2  = diprod(dx, dx);

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
static ACT_ALWAYS_INLINE void do_dih_fup_noshiftf(int i, int j, int k, int l, double ddphi,
                                       const dvec r_ij, const dvec r_kj, const dvec r_kl,
                                       const dvec m, const dvec n,
                                       F *forces)
{
    dvec   f_i, f_j, f_k, f_l;
    dvec   uvec, vvec, svec;
    double iprm, iprn, nrkj, nrkj2, nrkj_1, nrkj_2;
    double a, b, p, q, toler;

    iprm  = diprod(m, m);       /*  5    */
    iprn  = diprod(n, n);       /*  5	*/
    nrkj2 = diprod(r_kj, r_kj); /*  5	*/
    toler = nrkj2*GMX_REAL_EPS;
    if ((iprm > toler) && (iprn > toler))
    {
        nrkj_1 = gmx::invsqrt(nrkj2); /* 10	*/
        nrkj_2 = nrkj_1*nrkj_1;       /*  1	*/
        nrkj   = nrkj2*nrkj_1;        /*  1	*/
        a      = -ddphi*nrkj/iprm;    /* 11	*/
        dsvmul(a, m, f_i);            /*  3	*/
        b     = ddphi*nrkj/iprn;      /* 11	*/
        dsvmul(b, n, f_l);            /*  3  */
        p     = diprod(r_ij, r_kj);   /*  5	*/
        p    *= nrkj_2;               /*  1	*/
        q     = diprod(r_kl, r_kj);   /*  5	*/
        q    *= nrkj_2;               /*  1	*/
        dsvmul(p, f_i, uvec);         /*  3	*/
        dsvmul(q, f_l, vvec);         /*  3	*/
        dvec_sub(uvec, vvec, svec);   /*  3	*/
        dvec_sub(f_i, svec, f_j);     /*  3	*/
        dvec_add(f_l, svec, f_k);     /*  3	*/
        forceInc((*forces)[i], f_i);  /*  3	*/
        forceDec((*forces)[j], f_j);  /*  3	*/
        forceDec((*forces)[k], f_k);  /*  3	*/
        forceInc((*forces)[l], f_l);  /*  3	*/
    }
}

//! \brief Low-level code to compute dihedral angles
template <typename V>
static ACT_ALWAYS_INLINE double dih_angle(const V &xi, const V &xj, const V &xk, const V &xl,
                                          dvec r_ij, dvec r_kj, dvec r_kl, dvec m, dvec n)
{
    distanceVector(xi, xj, r_ij); /*  3        */
    distanceVector(xk, xj, r_kj); /*  3		*/
    distanceVector(xk, xl, r_kl); /*  3		*/

    dcprod(r_ij, r_kj, m);                       /*  9        */
    dcprod(r_kj, r_kl, n);                       /*  9		*/
    double phi  = gmx_angle_between_dvecs(m, n); /* 49 (assuming 25 for atan2) */
    double ipr  = diprod(r_ij, n);               /*  5        */
    double sign = (ipr < 0.0) ? -1.0 : 1.0;
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

    dvec r_ij, r_kj, r_kl, m, n;
    auto phi = M_PI + dih_angle(x[ai], x[aj], x[ak], x[al],
                                r_ij, r_kj, r_kl, m, n);

    double cos_phi = std::cos(phi);

    // Energy is
    // \sum_{j=0}^{nparam} params[j] cos^j(phi)
    // Force then becomes
    // - \sum_{j=1}^{nparam} j * params[j] cos^{j-1}(phi) sin(phi)
    double cosfac = 1;
    double v      = params[0];
    double dvdphi = 0;
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
    const  double half = 0.5;
    // Get the parameters. We have to know their names to do this.
    auto &params  = a.params();
    auto kA       = params[idihKPHI];
    // Get the atom indices
    auto ai       = indices[0];
    auto aj       = indices[1];
    auto ak       = indices[2];
    auto al       = indices[3];

    dvec r_ij, r_kj, r_kl, m, n;
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

    dvec r_ij, r_kj, r_kl, m, n;
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
template <typename One, typename E, typename Body>
static ACT_ALWAYS_INLINE void withOne(const One           &one,
                                      const E             &entry,
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
 * The coordinates of the particles of an interaction are relative to its
 * first particle, which is exact, see ConformationCoordinates, except for
 * interactions of one particle, such as position restraints.
 * \param[in]    msghandler  For warnings and errors
 * \param[in]    entries     The atom identifiers and parameters
 * \param[in]    x           The coordinates of all particles
 * \param[inout] forces      The forces on all particles
 * \param[out]   energies    The energy per type
 * \param[in]    one         Function computing one interaction
 */
template <typename Range, typename One, typename X, typename F>
static inline void computeConformation(MsgHandler         *msghandler,
                                       const Range        &entries,
                                       const X            &x,
                                       F                  *forces,
                                       double              energies[MAX_ENER],
                                       One                 one)
{
    double e[MAX_ENER] = { 0, 0, 0 };
    X      xe          = x;
    for (const auto &b : entries)
    {
        const auto &entry   = *b.operator->();
        const auto &indices = entry.atomIndices();
        xe.setOrigin(indices.size() > 1 ? indices[0] : -1);
        withOne(one, entry, [&](const auto &o)
                {
                    o(msghandler, entry, indices, xe, forces, e);
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
    [](MsgHandler *mh, const auto &b, const auto &idx, const auto &x, auto *f, double *e) \
    { func(mh, b, idx, x, f, e); }

//! As ONE_INTERACTION, for functions that also need the atoms and epsilonr
#define ONE_INTERACTION_ATOMS_EPSR(func) \
    [&atoms, epsilonr](MsgHandler *mh, const auto &b, const auto &idx, const auto &x, auto *f, double *e) \
    { func(mh, b, idx, x, f, e, atoms, epsilonr); }

//! As ONE_INTERACTION, for functions that also need the atoms
#define ONE_INTERACTION_ATOMS(func) \
    [&atoms](MsgHandler *mh, const auto &b, const auto &idx, const auto &x, auto *f, double *e) \
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

/*! \brief The energy terms of a potential used for an interaction type
 * A potential with one term, used for a correction, gives the energy of
 * that correction rather than its usual term, e.g. MACDANIEL_SCHMIDT for
 * INDUCTIONCORRECTION or BORN_MAYER for VDWCORRECTION.
 * \param[in] mypot The potential
 * \param[in] itype The interaction type
 * \return The energy terms
 */
static const std::vector<InteractionType> &energyTypes(Potential mypot, InteractionType itype)
{
    static const std::vector<InteractionType> inductionCorrection = { InteractionType::INDUCTIONCORRECTION };
    static const std::vector<InteractionType> vdwCorrection       = { InteractionType::VDWCORRECTION };
    const auto &types = energyTypes(mypot);
    if (types.size() == 1)
    {
        if (InteractionType::INDUCTIONCORRECTION == itype)
        {
            return inductionCorrection;
        }
        else if (InteractionType::VDWCORRECTION == itype)
        {
            return vdwCorrection;
        }
    }
    return types;
}

/*! \brief Store the energies of one term, as Term::set but without branches
 * The pointers are restricted, since otherwise each store to have, which
 * may alias anything, loads all of them again. Then the loop vectorizes if
 * the slots hold the conformations in order.
 * \tparam      nsum       Number of energies summed for the term
 * \param[in]    nslot      The number of conformations to store
 * \param[in]    confs      The conformations to store, if nullptr 0 to nslot - 1
 * \param[in]    slots      The slot of each, if nullptr conformation confs[i] is in slot i
 * \param[in]    e          The energies to sum, that of slot s at e[k][s*confStride]
 * \param[in]    confStride Stride between slots
 * \param[in]    add        Whether to add to the energies there are
 * \param[inout] act        The energies, see ConformationEnergies::Term
 * \param[inout] before     The backups
 * \param[inout] have       Whether there is an energy
 * \param[in]    backup     Whether the energy is its own backup as well, see backupEnergies
 */
template<int nsum>
static void storeTerm(size_t                     nslot,
                      const size_t *gmx_restrict confs,
                      const size_t *gmx_restrict slots,
                      const double              *e[nsum],
                      size_t                     confStride,
                      bool                       add,
                      double       *gmx_restrict act,
                      double       *gmx_restrict before,
                      char         *gmx_restrict have,
                      bool                       backup = false)
{
    const double *gmx_restrict e0 = e[0];
    const double *gmx_restrict e1 = e[nsum > 1 ? 1 : 0];
    const double *gmx_restrict e2 = e[nsum > 2 ? 2 : 0];
    for(size_t i = 0; i < nslot; i++)
    {
        size_t c    = confs ? confs[i] : i;
        size_t s    = (slots ? slots[i] : i)*confStride;
        double ei   = e0[s];
        if (nsum == 3)
        {
            ei = ei + e1[s] + e2[s];
        }
        double prev = (add && have[c]) ? act[c] : 0.0;
        double now  = prev + ei;
        before[c]   = backup ? now : (have[c] ? before[c] : 0.0);
        have[c]     = 1;
        act[c]      = now;
    }
}

/*! \brief Store the energies of one potential and add them to the potential energy
 * \param[in]    mypot      The potential
 * \param[in]    nslot      The number of conformations to store
 * \param[in]    confs      The conformations to store, if nullptr 0 to nslot - 1
 * \param[in]    slots      The slot of each, if nullptr conformation confs[i]
 *                          is in slot i
 * \param[inout] cener      The energies of all conformations
 * \param[in]    energies   The energies, energy k of slot i is stored
 *                          at energies[i*confStride + k*termStride]
 * \param[in]    confStride Stride between slots
 * \param[in]    termStride Stride between energy terms
 * \param[in]    backup     Whether the energies are their own backups as well
 */
static void storeEnergies(Potential             mypot,
                          size_t                nslot,
                          const size_t         *confs,
                          const size_t         *slots,
                          ConformationEnergies *cener,
                          const double         *energies,
                          size_t                confStride,
                          size_t                termStride,
                          bool                  backup = false)
{
    const auto &itypes = cener->termTypes() ? *cener->termTypes() : energyTypes(mypot);
    std::array<ConformationEnergies::Term, MAX_ENER> terms;
    for(size_t k = 0; k < itypes.size(); k++)
    {
        terms[k] = cener->add(itypes[k]);
    }
    auto eEpot = cener->add(InteractionType::EPOT);
    // When the entries are computed in parts, the energies are added
    bool accumulate = cener->accumulate();
    // One term at a time such that the energies are contiguous
    for(size_t k = 0; k < itypes.size(); k++)
    {
        const double *ek[1] = { energies + k*termStride };
        storeTerm<1>(nslot, confs, slots, ek, confStride, accumulate,
                     terms[k].act, terms[k].before, terms[k].have, backup);
    }
    const double *eall[3] = { energies, energies + termStride, energies + 2*termStride };
    storeTerm<3>(nslot, confs, slots, eall, confStride, true, eEpot.act, eEpot.before, eEpot.have, backup);
}

/*! \brief Whether a potential uses the charges of the atoms
 * Potentials that use ActAtom::charge must be listed here, since their
 * energies can only be used again if the charges did not change.
 * \param[in] mypot The potential
 */
static bool usesCharges(Potential mypot)
{
    switch (mypot)
    {
    case Potential::COULOMB_POINT:
    case Potential::COULOMB_GAUSSIAN:
    case Potential::COULOMB_SLATER:
    case Potential::POLARIZATION:
        return true;
    default:
        return false;
    }
}

//! Number of values per atom kept by keepCharges
static const size_t chargeValues = 5;

/*! \brief The charge and the second charge of an atom with their distributions
 * \param[in]  atom The atom
 * \param[out] q    The values
 */
static void chargeOf(const ActAtom &atom, double q[chargeValues])
{
    q[0] = atom.charge();
    q[1] = atom.zeta();
    q[2] = atom.chargeVS();
    q[3] = atom.zetaVS();
    q[4] = atom.rowVS();
}

/*! \brief Keep the charges of the atoms, if the potential uses them
 * The second charges and their distributions are kept as well.
 * \param[in]  mypot   The potential
 * \param[in]  atoms   The atoms
 * \param[out] charges The charges, empty if they are not used
 */
static void keepCharges(Potential                   mypot,
                        const std::vector<ActAtom> &atoms,
                        std::vector<double>        *charges)
{
    charges->clear();
    if (usesCharges(mypot))
    {
        for(const auto &atom : atoms)
        {
            double q[chargeValues];
            chargeOf(atom, q);
            charges->insert(charges->end(), q, q + chargeValues);
        }
    }
}

/*! \brief Whether the charges of the atoms are the ones kept
 * \param[in] mypot   The potential
 * \param[in] atoms   The atoms
 * \param[in] charges The charges from keepCharges
 * \return true if the potential does not use charges or they are the same
 */
static bool sameCharges(Potential                   mypot,
                        const std::vector<ActAtom> &atoms,
                        const std::vector<double>  &charges)
{
    if (!usesCharges(mypot))
    {
        return true;
    }
    if (charges.size() != chargeValues*atoms.size())
    {
        return false;
    }
    for(size_t i = 0; i < atoms.size(); i++)
    {
        double q[chargeValues];
        chargeOf(atoms[i], q);
        if (!std::equal(q, q + chargeValues, charges.begin() + chargeValues*i))
        {
            return false;
        }
    }
    return true;
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
//! \tparam Pairs Whether the interactions are non-bonded pairs, see PairRange, then only those potentials are available
template <bool Pairs = false, typename Func>
static bool dispatchPotential(Potential                   mypot,
                              const std::vector<ActAtom> &atoms,
                              double                      epsilonr,
                              Func                        func)
{
    switch(mypot)
    {
    case Potential::HARMONIC_BONDS:
        if constexpr (Pairs) { return false; } else { func.template operator()<2, true>(ONE_INTERACTION(bondOne)); }
        break;
    case Potential::MORSE_BONDS:
        if constexpr (Pairs) { return false; } else { func.template operator()<2, true>(ONE_INTERACTION(morseOne)); }
        break;
    case Potential::HUA_BONDS:
        if constexpr (Pairs) { return false; } else { func.template operator()<2, true>(ONE_INTERACTION(huaOne)); }
        break;
    case Potential::CUBIC_BONDS:
        if constexpr (Pairs) { return false; } else { func.template operator()<2, true>(ONE_INTERACTION(cubicOne)); }
        break;
    case Potential::HARMONIC_ANGLES:
        if constexpr (Pairs) { return false; } else { func.template operator()<3, true>(ONE_INTERACTION(angleOne)); }
        break;
    case Potential::LINEAR_ANGLES:
        if constexpr (Pairs) { return false; } else { func.template operator()<3, true>(ONE_INTERACTION(linearAngleOne)); }
        break;
    case Potential::UREY_BRADLEY_ANGLES:
        if constexpr (Pairs) { return false; } else { func.template operator()<3, true>(ONE_INTERACTION(ureyBradleyOne)); }
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
        func.template operator()<2, true>(OneSelect { [&atoms](const auto &b)
                                                      { return (atoms[b.atomIndices()[0]].zeta() == 0 &&
                                                                atoms[b.atomIndices()[1]].zeta() == 0 &&
                                                                atoms[b.atomIndices()[0]].chargeVS() == 0 &&
                                                                atoms[b.atomIndices()[1]].chargeVS() == 0); },
                                                      ONE_INTERACTION_ATOMS_EPSR(coulombPointOne),
                                                      ONE_INTERACTION_ATOMS_EPSR(coulombGaussianOne) });
        break;
    case Potential::COULOMB_SLATER:
        // The Slater integrals are not inlined, so this cannot be vectorized
        func.template operator()<2, false>(ONE_INTERACTION_ATOMS_EPSR(coulombSlaterOne));
        break;
    case Potential::POLARIZATION:
        if constexpr (Pairs) { return false; } else { func.template operator()<2, true>(ONE_INTERACTION_ATOMS(polarizationOne)); }
        break;
    case Potential::HARMONIC_DIHEDRALS:
        if constexpr (Pairs) { return false; } else { func.template operator()<4, true>(ONE_INTERACTION(improperOne)); }
        break;
    case Potential::PROPER_DIHEDRALS:
        if constexpr (Pairs) { return false; } else { func.template operator()<4, true>(ONE_INTERACTION(properOne)); }
        break;
    case Potential::FOURIER_DIHEDRALS:
        if constexpr (Pairs) { return false; } else { func.template operator()<4, true>(ONE_INTERACTION(fourDihOne)); }
        break;
    case Potential::POSITION_RESTRAINT:
        if constexpr (Pairs) { return false; } else { func.template operator()<1, true>(ONE_INTERACTION(fbposreOne)); }
        break;
    default:
        return false;
    }
    return true;
}

/*! \brief Timer for a potential
 * \param[in] mypot The potential
 * \return the index of its timer, see Timings
 */
static size_t potentialTimer(Potential mypot)
{
    // Indexed by potential, the timers are registered when first used
    static std::vector<size_t> timers;
    static const size_t        unset = std::numeric_limits<size_t>::max();
    size_t                     ip    = static_cast<size_t>(mypot);
    if (ip >= timers.size())
    {
        timers.resize(ip + 1, unset);
    }
    if (timers[ip] == unset)
    {
        timers[ip] = Timings::instance().index("Potential " + potentialToString(mypot));
    }
    return timers[ip];
}

/*! \brief Floating point operations of a potential for one conformation
 * The operations are counted in the explicit SIMD kernels in
 * forces_explicit_simd.cpp, the other potentials in the functions above.
 * Each addition, subtraction and multiplication counts as one operation,
 * a fused multiply-add as two and, as in GROMACS, a division and a square
 * root as five, since they take many cycles, such that 1/sqrt(x) counts
 * as ten. Functions count as the operations of their SIMD implementation
 * in simd_math.h: exp 30, erf with exp 90, acos 43 and atan2 40, and sin,
 * cos and log as exp. Operations on parameters only, which the compiler moves out of
 * the loop over conformations, and comparisons and selections are not
 * counted. The Slater integrals are not counted.
 * \param[in] mypot   The potential
 * \param[in] entries The interactions
 * \param[in] atoms   The atoms, for the charges
 * \param[in] forces  Whether the forces are computed
 * \return the number of operations
 */
template <typename Range>
static double interactionFlops(Potential                   mypot,
                               const Range                &entries,
                               const std::vector<ActAtom> &atoms,
                               bool                        forces)
{
    // Energy and the additional operations for the forces
    struct Cost
    {
        double energy, force;
    };
    // The distance of a pair, and spreading its force over the atoms
    const Cost pair     = { 9, 12 };
    // Some potentials leave out the interactions with parameters for which
    // they are zero, as the kernels do, these are not counted
    auto computed = [mypot](const auto &entry)
    {
        auto &params = entry.params();
        switch (mypot)
        {
        case Potential::LJ14_7:
            return params[lj14_7EPSILON] > 0 && params[lj14_7SIGMA] > 0;
        case Potential::WANG_BUCKINGHAM:
            return params[wbhEPSILON] > 0 && params[wbhGAMMA] > 0 && params[wbhSIGMA] > 0;
        case Potential::BUCKINGHAM:
            return params[bhA] > 0 && params[bhB] > 0 && params[bhC6] > 0;
        case Potential::POLARIZATION:
            return params[polALPHA] != 0;
        default:
            return true;
        }
    };
    auto constant = [&](Cost kernel, bool isPair)
    {
        double n = kernel.energy + (forces ? kernel.force : 0);
        if (isPair)
        {
            n += pair.energy + (forces ? pair.force : 0);
        }
        double ncomputed = entries.size();
        if (Potential::LJ14_7 == mypot || Potential::WANG_BUCKINGHAM == mypot ||
            Potential::BUCKINGHAM == mypot)
        {
            ncomputed = 0;
            for (const auto &b : entries)
            {
                ncomputed += computed(*b.operator->()) ? 1 : 0;
            }
        }
        return n*ncomputed;
    };
    switch (mypot)
    {
    case Potential::LJ12_6:                 return constant({ 12, 4 }, true);
    case Potential::LJ12_6_4:               return constant({ 14, 6 }, true);
    case Potential::LJ8_6:                  return constant({ 12, 4 }, true);
    case Potential::LJ14_7:                 return constant({ 39, 12 }, true);
    case Potential::BORN_MAYER:             return constant({ 44, 3 }, true);
    case Potential::SLATER_ISA:             return constant({ 49, 8 }, true);
    case Potential::MACDANIEL_SCHMIDT:      return constant({ 44, 3 }, true);
    case Potential::WANG_BUCKINGHAM:        return constant({ 58, 10 }, true);
    case Potential::BUCKINGHAM:             return constant({ 49, 5 }, true);
    case Potential::TANG_TOENNIES:          return constant({ 193, 29 }, true);
    case Potential::TT2b:                   return constant({ 193, 29 }, true);
    case Potential::SLATER_ISA_TT:          return constant({ 198, 34 }, true);
    case Potential::GENERALIZED_BUCKINGHAM: return constant({ 132, 35 }, true);
    case Potential::HARMONIC_BONDS:         return constant({ 16, 2 }, true);
    case Potential::CUBIC_BONDS:            return constant({ 18, 6 }, true);
    case Potential::MORSE_BONDS:            return constant({ 43, 8 }, true);
    case Potential::HUA_BONDS:              return constant({ 49, 13 }, true);
    case Potential::HARMONIC_ANGLES:        return constant({ 83, 71 }, false);
    case Potential::LINEAR_ANGLES:          return constant({ 23, 18 }, false);
    case Potential::UREY_BRADLEY_ANGLES:    return constant({ 97, 102 }, false);
    case Potential::HARMONIC_DIHEDRALS:     return constant({ 106, 71 }, false);
    case Potential::PROPER_DIHEDRALS:       return constant({ 132, 118 }, false);
    case Potential::POSITION_RESTRAINT:     return constant({ 17, 12 }, false);
    case Potential::COULOMB_POINT:
    case Potential::COULOMB_GAUSSIAN:
        {
            // The inverse distance, and per pair of charges, as in the
            // Coulomb kernel, see Coulomb::prepare
            double n = 0;
            for (const auto &b : entries)
            {
                const auto          &entry   = *b.operator->();
                const auto          &atomi   = atoms[entry.atomIndices()[0]];
                const auto          &atomj   = atoms[entry.atomIndices()[1]];
                const double         qi[2]   = { atomi.charge(), atomi.chargeVS() };
                const double         qj[2]   = { atomj.charge(), atomj.chargeVS() };
                const double         zi[2]   = { atomi.zeta(), atomi.zetaVS() };
                const double         zj[2]   = { atomj.zeta(), atomj.zetaVS() };
                n += pair.energy + 10 + (forces ? pair.force : 0);
                for (int a = 0; a < 2; a++)
                {
                    for (int c = 0; c < 2; c++)
                    {
                        if (a + c > 0 && qi[a]*qj[c] == 0)
                        {
                            continue;
                        }
                        bool gaussian = zi[a] != 0 || zj[c] != 0;
                        n += gaussian ? 95 : 2;
                        if (forces)
                        {
                            n += (gaussian ? 6 : 3) + (a + c > 0 ? 1 : 0);
                        }
                    }
                }
            }
            return n;
        }
    case Potential::POLARIZATION:
        {
            // More with hyperpolarization
            double n = 0;
            for (const auto &b : entries)
            {
                const auto &entry = *b.operator->();
                if (!computed(entry))
                {
                    continue;
                }
                bool hyper = entry.params()[polFCHYPER] > 0;
                n += pair.energy + 2 + (hyper ? 10 : 0);
                if (forces)
                {
                    n += pair.force + (hyper ? 8 : 0);
                }
            }
            return n;
        }
    case Potential::FOURIER_DIHEDRALS:
        {
            // Depends on the number of terms
            double n = 0;
            for (const auto &b : entries)
            {
                const auto &entry = *b.operator->();
                double nterm = std::min(entry.params().size(), static_cast<size_t>(fdihNR));
                double extra = std::max(nterm - 1, 0.0);
                n += 71 + 3*extra + (forces ? 67 + 2*extra : 0);
            }
            return n;
        }
    default:
        // Not counted, such as the Slater integrals
        return 0;
    }
}

/*! \brief Keep the energies of a potential, to use them again
 * \param[in]  energies   The energies, as passed to storeEnergies
 * \param[in]  size       Number of energies
 * \param[in]  confStride Distance between the energies of conformations
 * \param[in]  termStride Distance between the energy terms of a conformation
 * \param[out] cached     Where to keep them, may be nullptr
 */
static void keepEnergies(const double   *energies,
                         size_t          size,
                         size_t          confStride,
                         size_t          termStride,
                         CachedEnergies *cached)
{
    if (cached)
    {
        cached->energies.assign(energies, energies + size);
        cached->confStride = confStride;
        cached->termStride = termStride;
        cached->valid      = true;
    }
}

/*! \brief The coordinates of one conformation for the interaction functions
 * These use the coordinates of a particle as an rvec, as in the local arrays
 * of computeConformations. They are relative to an origin particle, from
 * the integers, such that they are exact, see coordinate.h, or absolute
 * without an origin.
 */
struct ConformationCoordinates
{
    //! The coordinates
    CoordinateVectors<const Coordinate> x;
    //! The origin particle, or -1 for absolute coordinates
    int                                 origin = -1;
    /*! \brief Set the origin
     * \param[in] o The origin particle, or -1 for absolute coordinates
     */
    void setOrigin(int o) { origin = o; }
    /*! \brief Return the coordinates of a particle
     * \param[in] a The particle
     * \return the coordinates
     */
    ACT_ALWAYS_INLINE gmx::RVec operator[](size_t a) const
    {
        return origin < 0 ? x[a].rvec() : x[a] - x[origin];
    }
};

/*! \brief Compute a potential one conformation at a time, with forces
 * The conformations in the computed slots are computed, see
 * ConformationArrays::arrange, since their forces are computed anew. The
 * energies of those for which the shells are not converged yet are stored
 * and summed into the potential energy.
 * \param[in]    msghandler The message handler, for debug output
 * \param[in]    mypot      The potential
 * \param[in]    topentries The interactions
 * \param[in]    atoms      The atoms
 * \param[inout] cdata      The conformations, the forces are added to theirs
 * \param[inout] cener      The energies of all conformations
 * \param[in]    epsilonr   Relative dielectric constant
 * \param[in]    ncompute   The number of slots to compute
 * \param[out]   cached     If not nullptr, the energies are kept here
 */
template <typename Range>
static inline void calcSome(MsgHandler                 *msghandler,
                            Potential                   mypot,
                            const Range                &topentries,
                            const std::vector<ActAtom> &atoms,
                            ConformationData           *cdata,
                            ConformationEnergies       *cener,
                            double                      epsilonr,
                            size_t                      ncompute,
                            CachedEnergies             *cached = nullptr)
{
    ScopedTimer timer(potentialTimer(mypot));
    // The conformations for which the shells are not converged, see
    // lowComputeEnergies, and their slots, the others keep their energies
    const auto         &ca = cdata->arrays;
    std::vector<size_t> confs, slots;
    for(size_t slot = 0; slot < ncompute; slot++)
    {
        size_t c = ca.conformation(slot);
        if (!cdata->confs[c].done)
        {
            confs.push_back(c);
            slots.push_back(slot);
        }
    }
    // Energies for each slot, MAX_ENER per slot
    std::vector<double> energies(ncompute*MAX_ENER, 0.0);
    for(size_t slot = 0; slot < ncompute; slot++)
    {
        size_t                  c = ca.conformation(slot);
        ConformationCoordinates x = { cdata->x(c) };
        auto                    f = cdata->f(c);
        // The forces are cleared in lowComputeEnergies, since all potentials add to them
        dispatchPotential<isPairRange<Range>>(mypot, atoms, epsilonr,
                          [&]<int NATOM, bool Vectorize>(auto one)
                          {
                              computeConformation(msghandler, topentries, x, &f,
                                                  energies.data() + slot*MAX_ENER, one);
                          });
    }
    // Store the energies and sum them into the potential energy.
    if (confs.size() == cdata->nconf() && ca.inOrder())
    {
        storeEnergies(mypot, confs.size(), nullptr, nullptr, cener, energies.data(), MAX_ENER, 1);
    }
    else
    {
        storeEnergies(mypot, confs.size(), confs.data(), slots.data(), cener, energies.data(), MAX_ENER, 1);
    }
    keepEnergies(energies.data(), energies.size(), MAX_ENER, 1, cached);
    // The forces are always computed here
    if (Timings::instance().countFlops())
    {
        Timings::instance().addFlops(ncompute*interactionFlops(mypot, topentries, atoms, true));
    }
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
template <int NATOM, bool Vectorize, bool Forces, typename Range, typename One, typename T>
static void computeConformations(const Range               &entries,
                                 ConformationArrays<T>     *ca,
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
    size_t stride = ca->stride;
    size_t nslot  = ca->nslot;
    for (const auto &b : entries)
    {
        const auto &entry = *b.operator->();
        auto  &indices = entry.atomIndices();
        const Coordinate *xa[NATOM];
        T                *fa[NATOM];
        for (int a = 0; a < NATOM; a++)
        {
            xa[a] = ca->x.data() + DIM*indices[a]*stride;
            fa[a] = Forces ? ca->f.data() + DIM*indices[a]*stride : nullptr;
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
                        // Relative to the first atom, which is exact, see
                        // coordinate.h, except for interactions of one
                        // atom, such as position restraints
                        Coordinate x0 = (NATOM > 1) ? xa[0][m*stride+c] : 0;
                        xl[a][m] = fromCoordinate(coordinateDifference(xa[a][m*stride+c], x0));
                        fl[a][m] = Forces ? fa[a][m*stride+c] : 0.0;
                    }
                }
                double el[MAX_ENER];
                for (int k = 0; k < MAX_ENER; k++)
                {
                    el[k] = energies[k*stride+c];
                }
                o(nullptr, entry, lidx, xl, &fl, el);
                if constexpr (Forces)
                {
                    for (int a = 0; a < NATOM; a++)
                    {
                        for (int m = 0; m < DIM; m++)
                        {
                            fa[a][m*stride+c] = fl[a][m];
                        }
                    }
                }
                for (int k = 0; k < MAX_ENER; k++)
                {
                    energies[k*stride+c] = el[k];
                }
            };
            if constexpr (Vectorize)
            {
#pragma omp simd
                for(size_t c = 0; c < nslot; c++)
                {
                    lane(c);
                }
            }
            else
            {
                for(size_t c = 0; c < nslot; c++)
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
 * \param[inout] ca         Coordinates and forces of all conformations, the
 *                          interactions compute them in real precision
 * \param[inout] energies   Energies, term k of conformation c at energies[k*stride + c]
 * \param[in]    epsilonr   Relative dielectric constant
 * \return false if the potential is not implemented
 */
template <bool Forces, typename Range, typename T>
static ACT_NOINLINE bool dispatchConformations(Potential                   mypot,
                                               const Range                &topentries,
                                               const std::vector<ActAtom> &atoms,
                                               ConformationArrays<T>      *ca,
                                               double                     *energies,
                                               double                      epsilonr)
{
    return dispatchPotential<isPairRange<Range>>(mypot, atoms, epsilonr,
                             [&]<int NATOM, bool Vectorize>(auto one)
                             {
                                 computeConformations<NATOM, Vectorize, Forces>(topentries, ca, energies, one);
                             });
}

/*! \brief Compute a potential with conformations in the inner loop
 * \param[in]    mypot      The potential
 * \param[in]    topentries The interactions
 * \param[in]    atoms      The atoms
 * \param[inout] ca         Coordinates and forces of the conformations
 * \param[inout] cener      The energies of all conformations
 * \param[in]    epsilonr   Relative dielectric constant
 * \param[in]    forces     Whether to compute forces
 * \param[in]    active     The conformations to store the energies of, all
 *                          if nullptr, see ConformationArrays::arrange
 * \param[in]    slots      The slot of each of active
 * \param[out]   cached     If not nullptr, the energies are kept here
 * \return false if this potential is not implemented this way
 */
template <typename Range, typename T>
static bool calcConformations(Potential                        mypot,
                              const Range                     &topentries,
                              const std::vector<ActAtom>      &atoms,
                              ConformationArrays<T>           *ca,
                              ConformationEnergies            *cener,
                              double                           epsilonr,
                              bool                             forces,
                              const std::vector<size_t>       *active,
                              const std::vector<size_t>       *slots,
                              CachedEnergies                  *cached = nullptr)
{
    ScopedTimer timer(potentialTimer(mypot));
    // Energy term k of conformation c is stored at k*stride + c
    auto &energies = ca->energies;
    energies.assign(MAX_ENER*ca->stride, 0.0);
    // Use the GPU or explicit SIMD where available
    if (!computeGpu(mypot, topentries, atoms, epsilonr, ca->stride, ca->nslot,
                    ca->x.data(), ca->f.data(), energies.data(), forces) &&
        !(useExplicitSimd() &&
          computeExplicitSimd(mypot, topentries, atoms, epsilonr, ca->stride, ca->nslot,
                              ca->x.data(), ca->f.data(), energies.data(), forces)) &&
        !(forces ?
          dispatchConformations<true>(mypot, topentries, atoms, ca, energies.data(), epsilonr) :
          dispatchConformations<false>(mypot, topentries, atoms, ca, energies.data(), epsilonr)))
    {
        return false;
    }
    if (active)
    {
        storeEnergies(mypot, active->size(), active->data(), slots->data(), cener, energies.data(), 1, ca->stride);
    }
    else
    {
        storeEnergies(mypot, ca->nconf, nullptr, nullptr, cener, energies.data(), 1, ca->stride);
    }
    keepEnergies(energies.data(), energies.size(), 1, ca->stride, cached);
    // The padding conformations are computed as well, but not counted
    if (Timings::instance().countFlops())
    {
        Timings::instance().addFlops(ca->ncomputed()*interactionFlops(mypot, topentries, atoms, forces));
    }
    return true;
}

//! Largest number of values of CachedEnergies::kernel, 32 MB
static const size_t maxCoulombKernel = 4*1024*1024;

/*! \brief Compute Coulomb energies from those of the interactions for charges of one
 * In training, often only the charges change. The energy of each
 * interaction for charges of one is computed once and kept in cached,
 * see CachedEnergies::kernel, then the energies for any charges are a
 * sum of the products of these with the charges. This needs explicit SIMD
 * and arrays in real precision, in which the kernel is kept.
 * \param[in]    mypot      The potential
 * \param[in]    topentries The interactions
 * \param[in]    atoms      The atoms
 * \param[inout] ca         Coordinates of the conformations, all in the slots
 * \param[inout] cener      The energies of all conformations
 * \param[in]    epsilonr   Relative dielectric constant
 * \param[in]    coordinateVersion ConformationData::coordinates
 * \param[inout] cached     The kernel, and the energies are kept here
 * \return false if this is not possible for the potential
 */
template <typename Range>
static bool calcCoulombKernel(Potential                   mypot,
                              const Range                &topentries,
                              const std::vector<ActAtom> &atoms,
                              ConformationArrays<real>   *ca,
                              ConformationEnergies       *cener,
                              double                      epsilonr,
                              size_t                      coordinateVersion,
                              CachedEnergies             *cached)
{
    if (!(Potential::COULOMB_POINT == mypot || Potential::COULOMB_GAUSSIAN == mypot) ||
        !cached || topentries.empty() || !useExplicitSimd())
    {
        return false;
    }
    size_t ne = topentries.size();
    size_t st = ca->stride;
    if (4*ne*st > maxCoulombKernel)
    {
        return false;
    }
    ScopedTimer timer(potentialTimer(mypot));
    auto       &kernel = cached->kernel;
    auto       &terms  = cached->kernelTerms;
    auto       &widths = cached->kernelWidths;
    const auto &x      = ca->x;
    auto       &xk     = cached->kernelX;
    // Whether kernel holds these interactions and conformations, otherwise
    // all are computed
    bool sameShape = (cached->kernelValid && cached->kernelEntries == ne &&
                      cached->kernelConformations == ca->nconf && cached->kernelStride == st &&
                      terms.size() == ne && widths.size() == 4*ne && xk.size() == x.size());
    if (!sameShape)
    {
        kernel.resize(4*ne*st);
        terms.assign(ne, 0);
        widths.assign(4*ne, 0.0);
    }
    // The particles that moved, such as virtual sites when their
    // parameters change
    bool              newCoordinates = (cached->kernelCoordinateVersion != coordinateVersion);
    std::vector<char> moved;
    if (sameShape && newCoordinates)
    {
        size_t natom = ca->natom;
        moved.resize(natom);
        for(size_t a = 0; a < natom; a++)
        {
            moved[a] = !std::equal(x.begin() + DIM*a*st, x.begin() + DIM*(a+1)*st,
                                   xk.begin() + DIM*a*st);
        }
    }
    // The interactions to compute again, with the terms of each: those
    // whose widths changed, of particles that moved, or that need terms
    // with second charges that were not computed. When the widths are
    // trained, only the interactions with the changed widths are computed.
    std::vector<unsigned char> todo(ne, 0);
    size_t                     ntodo = 0;
    for(size_t e = 0; e < ne; e++)
    {
        const auto   &ind    = topentries[e]->atomIndices();
        const auto   &ai     = atoms[ind[0]];
        const auto   &aj     = atoms[ind[1]];
        bool          vsi    = ai.chargeVS() != 0;
        bool          vsj    = aj.chargeVS() != 0;
        unsigned char needed = 1 | (vsj ? 2 : 0) | (vsi ? 4 : 0) | (vsi && vsj ? 8 : 0);
        // The widths the terms are computed from, see coulombKernelsExplicitSimd
        const double  w[4]   = { ai.zeta(), ai.zetaVS(), aj.zeta(), aj.zetaVS() };
        double       *we     = widths.data() + 4*e;
        bool          same   = (sameShape && (terms[e] & needed) == needed &&
                                std::equal(w, w + 4, we) &&
                                !(!moved.empty() && (moved[ind[0]] || moved[ind[1]])));
        if (!same)
        {
            todo[e]  = needed;
            terms[e] = needed;
            std::copy(w, w + 4, we);
            ntodo   += 1;
        }
    }
    if (ntodo > 0)
    {
        cached->kernelValid = coulombKernelsExplicitSimd(topentries, atoms, st, ca->nslot, x.data(),
                                                         todo.data(), kernel.data());
        if (!cached->kernelValid)
        {
            return false;
        }
        cached->kernelEntries       = ne;
        cached->kernelConformations = ca->nconf;
        cached->kernelStride        = st;
        // As computing the potential, where terms without charge are left
        // out, for the interactions computed
        if (Timings::instance().countFlops())
        {
            Timings::instance().addFlops(ca->nconf*interactionFlops(mypot, topentries, atoms, false)*ntodo/ne);
        }
    }
    if (newCoordinates || !sameShape)
    {
        xk.assign(x.begin(), x.end());
        cached->kernelCoordinateVersion = coordinateVersion;
    }
    // The energies for the charges of the atoms, each term as in the
    // explicit SIMD routine of the potential
    auto &energies = ca->energies;
    energies.assign(MAX_ENER*st, 0.0);
    double *en    = energies.data();
    double  flops = 0;
    for(size_t e = 0; e < ne; e++)
    {
        const auto   &ind  = topentries[e]->atomIndices();
        const auto   &ai   = atoms[ind[0]];
        const auto   &aj   = atoms[ind[1]];
        const double  qi[2] = { ai.charge(), ai.chargeVS() };
        const double  qj[2] = { aj.charge(), aj.chargeVS() };
        for(int k = 0; k < 4; k++)
        {
            double qq = ONE_4PI_EPS0*qi[k/2]*qj[k%2]/epsilonr;
            if (qq != 0 && (terms[e] & (1 << k)))
            {
                const real *kk = kernel.data() + (4*e + k)*st;
                for(size_t c = 0; c < ca->nslot; c++)
                {
                    en[c] += qq*kk[c];
                }
                flops += 2*ca->nconf;
            }
        }
    }
    storeEnergies(mypot, ca->nconf, nullptr, nullptr, cener, energies.data(), 1, st);
    keepEnergies(energies.data(), energies.size(), 1, st, cached);
    Timings::instance().addFlops(flops);
    return true;
}

/*! \brief Without arrays in real precision there is no kernel, see calcCoulombKernel
 * \return false
 */
template <typename T, typename Range>
requires (!std::is_same_v<T, real>)
static bool calcCoulombKernel(Potential, const Range &, const std::vector<ActAtom> &,
                              ConformationArrays<T> *, ConformationEnergies *, double, size_t,
                              CachedEnergies *)
{
    return false;
}

/*! \brief Use the results of the GPU for a potential, as calcConformations does
 * \param[in]    mypot    The potential
 * \param[in]    entries  The interactions
 * \param[in]    atoms    The atoms
 * \param[in]    gpuBatch The results
 * \param[in]    task     The task of gpuBatch
 * \param[inout] ca       The forces are added to
 * \param[inout] cener    The energies of all conformations
 * \param[in]    forces   Whether forces were computed
 * \param[in]    active   The conformations to store the energies of, all
 *                         if nullptr
 * \param[in]    slots    The slot of each of active
 */
template <typename Range, typename T>
static void applyGpuResults(Potential                   mypot,
                            const Range                &entries,
                            const std::vector<ActAtom> &atoms,
                            const GpuBatch             &gpuBatch,
                            size_t                      task,
                            ConformationArrays<T>      *ca,
                            ConformationEnergies       *cener,
                            bool                        forces,
                            const std::vector<size_t>  *active,
                            const std::vector<size_t>  *slots)
{
    ScopedTimer timer(potentialTimer(mypot));
    auto &energies = ca->energies;
    energies.assign(MAX_ENER*ca->stride, 0.0);
    if constexpr (std::is_same_v<T, float>)
    {
        gpuBatch.addResults(task, ca->f.data(), energies.data());
    }
    else
    {
        GMX_RELEASE_ASSERT(false, "The GPU only computes in single precision");
    }
    if (active)
    {
        storeEnergies(mypot, active->size(), active->data(), slots->data(), cener, energies.data(), 1, ca->stride);
    }
    else
    {
        storeEnergies(mypot, ca->nconf, nullptr, nullptr, cener, energies.data(), 1, ca->stride);
    }
    // The padding conformations are computed as well, but not counted
    if (Timings::instance().countFlops())
    {
        Timings::instance().addFlops(ca->ncomputed()*interactionFlops(mypot, entries, atoms, forces));
    }
}

/*! \brief Compute all potentials for all conformations of a topology
 * \param[in]    msghandler    For warnings and errors
 * \param[in]    pd            The force field
 * \param[in]    top           The topology
 * \param[inout] cdata         The conformations, with their coordinates,
 *                             forces and energies
 * \param[in]    epsilonr      Relative dielectric constant
 * \param[in]    conformations Whether to loop over the conformations in the inner loop
 * \param[in]    forces        Whether to compute forces, if not the forces
 *                             in the conformations are not meaningful
 * \param[in]    reuse         Whether the energies of an interaction type can
 *                             be taken from the previous call, if its
 *                             parameters, the coordinates and the charges
 *                             did not change, see CachedEnergies. This
 *                             requires that all changes of the coordinates
 *                             are known to cdata.
 * \param[in]    skipDone      Whether to leave out the conformations for which
 *                             the shells converged, see CoordsEnergy::done.
 *                             Their energies and forces are those of the
 *                             previous call, since their shells did not move.
 *                             The others are put in the first slots of the
 *                             arrays, see ConformationData::arrangeSlots.
 *                             If false, the done flags are reset.
 * \param[in]    intermolecular Whether to compute only the intermolecular
 *                             interactions, see Topology::separateIntermolecular.
 *                             The energies of the other interaction types are zero.
 * \param[in]    allForces     Whether the forces on all particles are needed,
 *                             otherwise only those on the shells
 * \param[inout] gpuBatch      If not nullptr, the potentials that the GPU can
 *                             compute are computed together for all molecules:
 *                             in the first call, with collect true, they are
 *                             added to gpuBatch, and nothing else is computed,
 *                             in the second call the results are used, if
 *                             gpuBatch computed them, see computeEnergies.
 *                             Both calls make the same changes to cdata
 *                             before computing, such that the arrays the
 *                             tasks point to are the same.
 * \param[in]    collect       Whether to add the tasks to gpuBatch
 * \param[inout] gpuNext       The next task of gpuBatch to use
 * \tparam       T             The floating point type of the arrays, see
 *                             ConformationData::doublePrecision
 */
template <typename T>
static void lowComputeEnergiesT(MsgHandler                *msghandler,
                                const ForceField          *pd,
                                const Topology            *top,
                                ConformationData          *cdata,
                                double                     epsilonr,
                                bool                       conformations,
                                bool                       forces,
                                bool                       reuse          = false,
                                bool                       skipDone       = false,
                                bool                       intermolecular = false,
                                bool                       allForces      = true,
                                GpuBatch                  *gpuBatch       = nullptr,
                                bool                       collect        = false,
                                size_t                    *gpuNext        = nullptr)
{
    static const size_t timer = Timings::instance().index("Energy calculation overhead");
    ScopedTimer         scoped(timer);
    auto *coords = &cdata->confs;
    auto &cener  = cdata->energies;
    cener.setNumberOfConformations(coords->size());
    // The coordinates and forces of all conformations
    auto &ra = cdata->arrays;
    GMX_RELEASE_ASSERT(ra.nconf == coords->size(), "The arrays do not hold the conformations");
    // The conformations to compute, found once such that the loops below
    // do not visit those for which the shells converged. When few shells
    // are left to relax, visiting all conformations takes longer than
    // computing the others.
    // Conformations only become done while the shells are relaxed, then
    // those of the previous call that are not done are left.
    auto &active = cdata->active;
    if (skipDone && cdata->activeKnown == coords->size())
    {
        active.erase(std::remove_if(active.begin(), active.end(),
                                    [coords](size_t c) { return (*coords)[c].done; }),
                     active.end());
    }
    else
    {
        active.clear();
        for(size_t c = 0; c < coords->size(); c++)
        {
            auto &ce = (*coords)[c];
            if (!skipDone && ce.done)
            {
                ce.done = false;
            }
            if (!ce.done)
            {
                active.push_back(c);
            }
        }
    }
    cdata->activeKnown = coords->size();
    size_t ndone = coords->size() - active.size();
    // The virtual sites requested by ConformationData::constructVsites
    VsiteHandler vsiteHandler;
    bool         constructVsites = cdata->constructVsites;
    cdata->constructVsites       = false;
    if (constructVsites)
    {
        vsiteHandler.constructPositions(top, ra.x.data(), ra.stride, ra.nconf, ra.nconf);
    }
    // Only the conformations to compute are in the slots that are computed
    cdata->arrangeSlots(active);
    if (ndone == coords->size() && ndone > 0)
    {
        return;
    }
    // The padding slots are computed as well
    ra.pad();
    // The slot of each conformation to compute, to store their energies
    std::vector<size_t> activeSlots;
    if (ndone > 0)
    {
        for(auto c : active)
        {
            activeSlots.push_back(ra.slot(c));
        }
    }
    // Energies of only some conformations cannot be used again
    reuse = reuse && ndone == 0;
    // Reset the accumulated total potential energy (EPOT) for every conformation
    // before summing the contributions of the potentials. Without this,
    // repeated energy calculations (for example the many evaluations made
    // while training with the MCMC mutator) keep adding onto the value left by
    // the previous call, so the total grows without bound instead of being
    // recomputed from scratch each time.
    // The conformations that are done keep their energy, they only get
    // one if there was none, which is zero.
    bool newEpot = !cener.find(InteractionType::EPOT).have;
    auto eEpot   = cener.add(InteractionType::EPOT);
    if (newEpot)
    {
        std::fill(eEpot.have, eEpot.have + coords->size(), 1);
    }
    // As eEpot.set(c, 0.0), without branches
    for(auto c : active)
    {
        eEpot.before[c] = eEpot.have[c] ? eEpot.before[c] : 0.0;
        eEpot.act[c]    = 0.0;
        eEpot.have[c]   = 1;
    }
    // With shells, the entries without shells do not change while the shells
    // are relaxed. Their energies and forces are computed in the first
    // calculation, in which skipDone is false, and used again in the calls
    // relaxing the shells.
    auto &fixed      = cdata->fixed;
    bool  splitShell = forces && !intermolecular && top->hasShellSplit();
    // Whether the conformations are computed in the inner loop, otherwise
    // one at a time. The arrays of another precision get copies of the
    // coordinates, see ConformationData::doublePrecision.
    bool  useArrays  = conformations && coords->size() > 1;
    auto &ca         = cdata->template arraysOf<T>();
    constexpr bool copied = !std::is_same_v<T, real>;
    if constexpr (copied)
    {
        // The computed slots are those of this precision
        ca.copyCoordinates(ra);
        ra.nslot = ca.nslot;
    }
    // The forces of the entries without shells are kept in the layout of ca
    auto &fixedForces = fixed.template forcesOfArrays<T>();
    bool  useFixed    = (splitShell && skipDone && fixed.valid &&
                         fixed.nconf == coords->size() && fixedForces.size() == ca.f.size());
    bool  keepFixed  = splitShell && !useFixed;
    if (useFixed)
    {
        // Start from the energies of the entries without shells
        for(size_t k = 0; k < fixed.types.size(); k++)
        {
            // Through pointers that the stores cannot change
            auto                       t      = cener.add(fixed.types[k]);
            const char   *gmx_restrict haveK  = fixed.have[k].data();
            const double *gmx_restrict eK     = fixed.energies[k].data();
            double       *gmx_restrict before = t.before;
            double       *gmx_restrict act    = t.act;
            char         *gmx_restrict have   = t.have;
            for(auto c : active)
            {
                before[c] = have[c] ? before[c] : 0.0;
                have[c]   = haveK[c];
                act[c]    = eK[c];
            }
        }
    }
    if (forces)
    {
        if (useFixed)
        {
            ca.loadForces(fixedForces);
        }
        else
        {
            ca.clearForces();
        }
    }
    // Whether the current forces are in ca rather than in ra, only false
    // with copied arrays, after a potential that is computed one
    // conformation at a time
    bool forcesInArrays = true;
    // Which entries to compute
    enum class Part { All, NoShells, Shells };
    // Compute the interactions of one type, given as TopologyEntryRange or PairRange
    auto computeEntries = [&](InteractionType itype, Potential mypot, const auto &entries)
    {
        // On the GPU together with the other molecules, without the
        // energies of the previous call, see gpuBatch
        bool onGpu = (std::is_same_v<T, float> && gpuBatch && useArrays && !(reuse && !forces) &&
                      GpuBatch::canCompute(mypot));
        if (collect)
        {
            if (onGpu)
            {
                gpuBatch->add(mypot, entries, top->atoms(), epsilonr, ca.stride, ca.nslot, ca.x.data());
            }
            return;
        }
        if (onGpu && gpuBatch->computed())
        {
            if (forces && !forcesInArrays)
            {
                ca.copyForces(ra);
                forcesInArrays = true;
            }
            applyGpuResults(mypot, entries, top->atoms(), *gpuBatch, (*gpuNext)++, &ca, &cener,
                            forces, ndone > 0 ? &active : nullptr, &activeSlots);
            return;
        }
        // Use the energies of the previous call if nothing changed. Forces
        // are not kept, so this is only possible if they are not needed.
        CachedEnergies *cached = nullptr;
        if (reuse && !forces)
        {
            cached = &cdata->cache[itype];
            if (cached->valid && cached->potential == mypot &&
                cached->parameterVersion == top->parameterVersion(itype) &&
                cached->coordinateVersion == cdata->coordinates.value() &&
                cached->nentries == entries.size() &&
                cached->nconf == coords->size() &&
                sameCharges(mypot, top->atoms(), cached->charges))
            {
                static const size_t timer = Timings::instance().index("Energies used again");
                ScopedTimer         scoped(timer);
                storeEnergies(mypot, coords->size(), nullptr, nullptr, &cener, cached->energies.data(),
                              cached->confStride, cached->termStride);
                return;
            }
            cached->valid             = false;
            cached->potential         = mypot;
            cached->parameterVersion  = top->parameterVersion(itype);
            cached->coordinateVersion = cdata->coordinates.value();
            cached->nentries          = entries.size();
            cached->nconf             = coords->size();
            keepCharges(mypot, top->atoms(), &cached->charges);
        }
        if (useArrays)
        {
            if (forces && !forcesInArrays)
            {
                ca.copyForces(ra);
                forcesInArrays = true;
            }
            // Without forces, only the charges may have changed
            if (cached && ca.nconf == coords->size() &&
                calcCoulombKernel(mypot, entries, top->atoms(), &ca, &cener, epsilonr,
                                  cdata->coordinates.value(), cached))
            {
                return;
            }
            if (calcConformations(mypot, entries, top->atoms(), &ca, &cener, epsilonr, forces,
                                  ndone > 0 ? &active : nullptr, &activeSlots, cached))
            {
                return;
            }
            if (forces && copied && forcesInArrays)
            {
                ca.storeForces(&ra);
            }
            forcesInArrays = !copied;
        }
        calcSome(msghandler, mypot, entries, top->atoms(), cdata, &cener, epsilonr,
                 ca.ncomputed(), cached);
    };
    auto computePart = [&](Part part)
    {
        // All types with interactions, the non-bonded pairs are stored
        // compactly, see NonBondedPairs
        for(auto itype : top->interactionTypes())
        {
            // Not all potentials are available in the force field,
            // in particular restraint potentials.
            // The function we need to do the math
            gmx_unused Potential mypot = pd->potential(itype);
            if (Potential::NONE == mypot)
            {
                mypot = defaultPotential(itype);
            }
            // Virtual sites have no energy
            if (energyTypes(mypot).empty())
            {
                continue;
            }
            // The terms depend on the interaction type for some potentials
            cener.setTermTypes(&energyTypes(mypot, itype));
            // The interactions to compute, without intramolecular ones they
            // may be none, then the energies are zero. The non-bonded pairs
            // are stored compactly, see NonBondedPairs.
            const auto &nb = top->nonBondedPairs();
            if (nb.hasType(itype))
            {
                PairRange pairs = (intermolecular ? nb.intermolecular(itype) : nb.range(itype));
                if (Part::All != part)
                {
                    pairs = nb.shellEntries(itype, Part::Shells == part);
                    if (Part::Shells == part && pairs.empty())
                    {
                        continue;
                    }
                }
                computeEntries(itype, mypot, pairs);
                continue;
            }
            const auto        &ev      = top->entry(itype);
            TopologyEntryRange entries = (intermolecular ?
                                          top->intermolecular(itype, ev) :
                                          TopologyEntryRange(ev));
            if (Part::All != part)
            {
                entries = top->shellEntries(itype, ev, Part::Shells == part);
                // Types without entries with shells keep the energies
                // of those without, the others are computed even if
                // empty, such that their energies are zero.
                if (Part::Shells == part && entries.empty())
                {
                    continue;
                }
            }
            computeEntries(itype, mypot, entries);
        }
    };
    if (collect)
    {
        // The same parts as below, only to add the tasks
        if (useFixed)
        {
            computePart(Part::Shells);
        }
        else if (keepFixed)
        {
            computePart(Part::NoShells);
            computePart(Part::Shells);
        }
        else
        {
            computePart(Part::All);
        }
        cener.setTermTypes(nullptr);
        return;
    }
    if (useFixed)
    {
        cener.setAccumulate(true);
        computePart(Part::Shells);
        cener.setAccumulate(false);
    }
    else if (keepFixed)
    {
        computePart(Part::NoShells);
        // Keep the energies and forces of the entries without shells
        if (!forcesInArrays)
        {
            ca.copyForces(ra);
            forcesInArrays = true;
        }
        // Only the forces of this precision are kept
        fixed.arrayForces.clear();
        fixed.arrayForcesDouble.clear();
        fixedForces = ca.f;
        fixed.nconf = coords->size();
        fixed.types = cener.types();
        fixed.energies.resize(fixed.types.size());
        fixed.have.resize(fixed.types.size());
        for(size_t k = 0; k < fixed.types.size(); k++)
        {
            auto t = cener.find(fixed.types[k]);
            fixed.energies[k].assign(t.act, t.act + coords->size());
            fixed.have[k].assign(t.have, t.have + coords->size());
        }
        fixed.valid = true;
        cener.setAccumulate(true);
        computePart(Part::Shells);
        cener.setAccumulate(false);
    }
    else
    {
        computePart(Part::All);
    }
    cener.setTermTypes(nullptr);
    // The van der Waals correction is part of the exchange, as in
    // ForceComputer::computeOnce
    auto vdwCorr = cener.find(InteractionType::VDWCORRECTION);
    if (vdwCorr.have)
    {
        auto exchange = cener.add(InteractionType::EXCHANGE);
        for(auto c : active)
        {
            if (vdwCorr.has(c))
            {
                exchange.set(c, (exchange.has(c) ? exchange.act[c] : 0.0) + vdwCorr.act[c]);
                vdwCorr.act[c] = 0;
            }
        }
    }
    if (forces && copied && forcesInArrays)
    {
        ca.storeForces(&ra);
    }
    // Spreading the forces on the virtual sites does not change those on
    // the shells, see ConformationData::spreadVsites
    bool spread = forces && useArrays && allForces && cdata->spreadVsites;
    if (spread)
    {
        vsiteHandler.distributeForces(top, ra.x.data(), ra.f.data(), ra.stride,
                                      ca.ncomputed(), ca.ncomputed());
    }
    // The conformations that were not computed keep the forces of the previous calculation
    cdata->vsiteForcesSpread = spread && (ndone == 0 || cdata->vsiteForcesSpread);
}

/*! \brief Compute all potentials for all conformations of a topology
 * As lowComputeEnergiesT, with the arrays of the precision that cdata asks
 * for, see ConformationData::doublePrecision.
 */
static void lowComputeEnergies(MsgHandler                *msghandler,
                               const ForceField          *pd,
                               const Topology            *top,
                               ConformationData          *cdata,
                               double                     epsilonr,
                               bool                       conformations,
                               bool                       forces,
                               bool                       reuse          = false,
                               bool                       skipDone       = false,
                               bool                       intermolecular = false,
                               bool                       allForces      = true,
                               GpuBatch                  *gpuBatch       = nullptr,
                               bool                       collect        = false,
                               size_t                    *gpuNext        = nullptr)
{
    if (cdata->useArraysDouble() && conformations && cdata->nconf() > 1)
    {
        lowComputeEnergiesT<double>(msghandler, pd, top, cdata, epsilonr, conformations, forces,
                                    reuse, skipDone, intermolecular, allForces,
                                    gpuBatch, collect, gpuNext);
    }
    else
    {
        lowComputeEnergiesT<real>(msghandler, pd, top, cdata, epsilonr, conformations, forces,
                                  reuse, skipDone, intermolecular, allForces,
                                  gpuBatch, collect, gpuNext);
    }
}

void computeConformations(MsgHandler                *msghandler,
                          const ForceField          *pd,
                          const Topology            *top,
                          ConformationData          *cdata,
                          bool                       forces,
                          bool                       skipDone)
{
    // The conformations only change if virtual sites are constructed
    // again, which ConformationData knows, or if shells move.
    bool reuse = !pd->polarizable();
    lowComputeEnergies(msghandler, pd, top, cdata, pd->epsilonr(),
                       true, forces, reuse, skipDone);
}

void computeEnergies(MsgHandler                  *msghandler,
                     const ForceField            *pd,
                     const std::vector<ACTMol *> &allmols,
                     size_t                       what,
                     bool                         forces,
                     bool                         skipDone,
                     bool                         intermolecular,
                     bool                         allForces)
{
    double epsilonr = pd->epsilonr();
    // The conformations of the molecules only change if virtual sites are
    // constructed again, which ConformationData knows, or if shells move.
    bool reuse = !pd->polarizable();
    // The molecules and their topologies to compute
    std::vector<std::pair<ACTMol *, bool>> todo;
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
        // The monomers of a dimer are not needed for intermolecular energies
        bool inter = intermolecular && tops.size() == 2 && mol.topology()->hasIntermolecular();
        if (inter && what < 2)
        {
            continue;
        }
        todo.push_back({ molptr, inter });
    }
    auto compute = [&](GpuBatch *gpuBatch, bool collect, size_t *gpuNext)
    {
        for(const auto &mi : todo)
        {
            auto &mol = *mi.first;
            if (what >= 2)
            {
                lowComputeEnergies(msghandler, pd, mol.topology(), mol.conformationData(2),
                                   epsilonr, true, forces, reuse, skipDone, mi.second, allForces,
                                   gpuBatch, collect, gpuNext);
            }
            else
            {
                lowComputeEnergies(msghandler, pd, &mol.fragmentHandler()->topologies()[what],
                                   mol.conformationData(what), epsilonr, true, forces, reuse,
                                   skipDone, false, allForces, gpuBatch, collect, gpuNext);
            }
        }
    };
    // With a GPU, the potentials it can compute are computed for all
    // molecules at once, see GpuBatch, in single precision only. Most of
    // the molecules are only collected if they can be enough work, as
    // estimated from all their interactions that the GPU could compute.
    size_t work = 0;
    if (gpu::useGpu() && !GMX_DOUBLE)
    {
        for(const auto &mi : todo)
        {
            auto       &mol   = *mi.first;
            const auto *top   = (what >= 2) ? mol.topology() : &mol.fragmentHandler()->topologies()[what];
            size_t      nconf = mol.conformationData(std::min<size_t>(what, 2))->nconf();
            for(auto itype : top->interactionTypes())
            {
                Potential mypot = pd->potential(itype);
                if (GpuBatch::canCompute(mypot))
                {
                    size_t n = top->nonBondedPairs().hasType(itype) ?
                        top->nonBondedPairs().range(itype).size() : top->entry(itype).size();
                    work += n*nconf;
                }
            }
        }
    }
    if (work >= gpuMinimumWork())
    {
        // Kept between calls, such that the memory is used again
        static thread_local GpuBatch gpuBatch;
        gpuBatch.clear();
        compute(&gpuBatch, true, nullptr);
        size_t next = 0;
        if (gpuBatch.compute(forces))
        {
            compute(&gpuBatch, false, &next);
            GMX_RELEASE_ASSERT(next == gpuBatch.size(), "Not all GPU results were used");
            return;
        }
    }
    compute(nullptr, false, nullptr);
}

#if GMX_DOUBLE

// The GPU computes in single precision only
bool relaxShellsGpu(const ForceField *, const std::vector<ACTMol *> &, size_t, const bool[2],
                    int, double, double, double, double, int, bool, size_t *nComputed)
{
    *nComputed = 0;
    return false;
}

#else

namespace
{

//! The data of a molecule for relaxing its shells on the GPU, see relaxShellsGpu
struct GpuRelaxMolecule
{
    const Topology                          *top   = nullptr;
    ConformationData                        *cdata = nullptr;
    std::vector<gpu::ForceShell>             forceShells;
    std::vector<gpu::MovingShell>            movingShells;
    //! The polarization of each moving shell, see GpuTrainingSystem::Impl::makeRowRecipes
    std::vector<const TopologyEntry *>       movingShellEntries;
    std::vector<gpu::ShellFragment>          fragments;
    std::vector<std::vector<gpu::PairEntry>> pairs;
    std::vector<Potential>                   potentials;
    std::vector<InteractionType>             itypes;
    std::vector<std::vector<double>>         energies;
    std::vector<gpu::PairBlock>              blocks;
    std::vector<char>                        done;
    std::vector<char>                        computed;
};

/*! \brief Whether a particle of the charge generation of a molecule is a virtual site
 * Only then do the coordinates of the charge generation depend on the
 * parameters of the virtual sites, since those of the particles with fixed
 * charges are not used, see GpuTrainingSystem::Impl::acmCoordinates.
 * \param[in] molptr The molecule
 * \return whether one is
 */
static bool acmHasVsites(const ACTMol *molptr)
{
    const auto *top = molptr->topology();
    if (!top->hasVsites())
    {
        return false;
    }
    std::set<int> sites;
    for (const auto &entry : top->entries())
    {
        if (isVsite(entry.first))
        {
            for (const auto &vs : entry.second)
            {
                sites.insert(vs->atomIndices().back());
            }
        }
    }
    const auto *fh = molptr->fragmentHandler();
    for (size_t f = 0; f < fh->topologies().size(); f++)
    {
        for (auto a : fh->acmParticles(f))
        {
            if (sites.count(fh->atomStart()[f] + a) > 0)
            {
                return true;
            }
        }
    }
    return false;
}

/*! \brief The shells of a monomer for the GPU, see ForceComputer::updateShells
 * \param[in]    top               The topology of the monomer
 * \param[in]    pol               Which particles are shells
 * \param[in]    offset            The first particle of the monomer
 * \param[in]    maxShellDistance2 Square of the largest distance from the parent
 * \param[inout] m                 The fragment is added
 */
void addShellFragment(const Topology   *top,
                      const PolInfo    &pol,
                      size_t            offset,
                      double            maxShellDistance2,
                      GpuRelaxMolecule *m)
{
    const auto       &atoms   = top->atoms();
    auto              isShell = pol.isShell();
    gpu::ShellFragment fr;
    fr.firstAtom        = offset;
    fr.natom            = atoms.size();
    fr.firstForceShell  = m->forceShells.size();
    fr.firstMovingShell = m->movingShells.size();
    fr.nShell           = pol.nShell();
    auto wallOf = [&](size_t i, uint32_t *parent)
    {
        bool wall = maxShellDistance2 > 0 && atoms[i].parents().size() == 1;
        *parent   = wall ? atoms[i].parents()[0] + offset : 0;
        return wall;
    };
    for (size_t i = 0; i < isShell.size(); i++)
    {
        if (isShell[i])
        {
            gpu::ForceShell fs;
            fs.shell = i + offset;
            fs.wall  = wallOf(i, &fs.parent);
            m->forceShells.push_back(fs);
        }
    }
    if (top->hasEntry(InteractionType::POLARIZATION))
    {
        for (const auto &p : top->entry(InteractionType::POLARIZATION))
        {
            size_t shellLocal = p->atomIndex(1);
            if (!isShell[shellLocal])
            {
                continue;
            }
            auto  &params = p->params();
            double alpha  = params[polALPHA];
            double qshell = atoms[shellLocal].charge();
            gpu::MovingShell ms;
            ms.shell   = shellLocal + offset;
            ms.core    = p->atomIndex(0) + offset;
            ms.wall    = wallOf(shellLocal, &ms.parent);
            ms.fcShell = (alpha > 0 && qshell != 0) ? qshell*qshell*ONE_4PI_EPS0/alpha : 0;
            ms.fcHyper = params[polFCHYPER];
            ms.rHyper  = params[polRHYPER];
            m->movingShells.push_back(ms);
            m->movingShellEntries.push_back(p.operator->());
        }
    }
    fr.nforceShell  = m->forceShells.size() - fr.firstForceShell;
    fr.nmovingShell = m->movingShells.size() - fr.firstMovingShell;
    m->fragments.push_back(fr);
}

/*! \brief The entries with shells of a topology for the GPU
 * \param[in]    pd       The force field
 * \param[in]    epsilonr Relative dielectric constant
 * \param[inout] m        The interactions are added
 * \return false if the GPU does not have a potential
 */
bool addShellEntries(const ForceField *pd, double epsilonr, GpuRelaxMolecule *m)
{
    const auto &top = *m->top;
    for (auto itype : top.interactionTypes())
    {
        Potential mypot = pd->potential(itype);
        if (Potential::NONE == mypot)
        {
            mypot = defaultPotential(itype);
        }
        if (energyTypes(mypot).empty())
        {
            continue;
        }
        // The non-bonded pairs are stored compactly, see NonBondedPairs
        const auto &nb      = top.nonBondedPairs();
        bool        isPairs = nb.hasType(itype);
        auto        entries = isPairs ? TopologyEntryRange() : top.shellEntries(itype, top.entry(itype), true);
        auto        pairs   = isPairs ? nb.shellEntries(itype, true) : PairRange();
        if (isPairs ? pairs.empty() : entries.empty())
        {
            continue;
        }
        gpu::PairKernel kernel;
        if (!gpuPairKernel(mypot, &kernel) ||
            energyTypes(mypot, itype).size() > static_cast<size_t>(gpu::pairEnergyTerms))
        {
            return false;
        }
        std::vector<gpu::PairEntry> pe;
        auto addEntries = [&](const auto &range)
        {
            pe.reserve(range.size());
            for (const auto &b : range)
            {
                gpu::PairEntry e;
                if (gpuPairEntry(kernel, *b.operator->(), top.atoms(), epsilonr, &e))
                {
                    pe.push_back(e);
                }
            }
        };
        if (isPairs)
        {
            addEntries(pairs);
        }
        else
        {
            addEntries(entries);
        }
        m->pairs.push_back(std::move(pe));
        m->potentials.push_back(mypot);
        m->itypes.push_back(itype);
        gpu::PairBlock blk;
        blk.kernel = kernel;
        m->blocks.push_back(blk);
    }
    return true;
}

} // namespace

bool relaxShellsGpu(const ForceField            *pd,
                    const std::vector<ACTMol *> &allmols,
                    size_t                       what,
                    const bool                   relaxMonomer[2],
                    int                          maxIter,
                    double                       msForceToler,
                    double                       maxShellDistance,
                    double                       minStep,
                    double                       maxStep,
                    int                          stallIterations,
                    bool                         allForces,
                    size_t                      *nComputed)
{
    *nComputed = 0;
    if (!gpu::useGpu())
    {
        return false;
    }
    double maxShellDistance2 = maxShellDistance*maxShellDistance;
    double epsilonr          = pd->epsilonr();
    // Everything is prepared before anything is changed, such that the
    // processor can relax the shells if the GPU cannot
    std::vector<GpuRelaxMolecule> mols;
    size_t                        work = 0;
    for (auto molptr : allmols)
    {
        auto  fh   = molptr->fragmentHandler();
        auto &tops = fh->topologies();
        auto &pols = fh->polInfo();
        if ((what < 2 && what >= tops.size()) || (what == 2 && tops.size() < 2))
        {
            continue;
        }
        GpuRelaxMolecule m;
        m.top   = (what < 2) ? &tops[what] : molptr->topology();
        m.cdata = molptr->conformationData(what);
        auto &ca    = m.cdata->arrays;
        auto &fixed = m.cdata->fixed;
        if (m.cdata->useArraysDouble() || !m.top->hasShellSplit() || !fixed.valid ||
            fixed.nconf != m.cdata->nconf() || fixed.arrayForces.size() != ca.f.size() ||
            ca.nconf != m.cdata->nconf() || (allForces && m.cdata->spreadVsites && m.top->hasVsites()))
        {
            return false;
        }
        if (what < 2)
        {
            addShellFragment(m.top, pols[what], 0, maxShellDistance2, &m);
        }
        else
        {
            for (size_t im = 0; im < 2; im++)
            {
                if (relaxMonomer[im])
                {
                    addShellFragment(&tops[im], pols[im], fh->atomStart()[im], maxShellDistance2, &m);
                }
            }
        }
        if (!addShellEntries(pd, epsilonr, &m))
        {
            return false;
        }
        for (const auto &pe : m.pairs)
        {
            work += pe.size()*ca.nconf;
        }
        mols.push_back(std::move(m));
    }
    // The work of all iterations, which the GPU does without the host
    if (mols.empty() || work*std::max(maxIter - 1, 1) < gpuMinimumWork())
    {
        return false;
    }
    static const size_t timer = Timings::instance().index("GPU shell relaxation");
    ScopedTimer         scoped(timer);
    std::vector<gpu::RelaxMolecule> rmols;
    for (auto &m : mols)
    {
        auto &ca = m.cdata->arrays;
        m.done.resize(ca.nconf);
        m.computed.assign(ca.nconf, 0);
        for (size_t s = 0; s < ca.nconf; s++)
        {
            m.done[s] = m.cdata->confs[ca.conformation(s)].done ? 1 : 0;
        }
        m.energies.assign(m.blocks.size(), std::vector<double>(MAX_ENER*ca.stride, 0.0));
        for (size_t b = 0; b < m.blocks.size(); b++)
        {
            m.blocks[b].pairs    = m.pairs[b].data();
            m.blocks[b].npair    = m.pairs[b].size();
            m.blocks[b].energies = m.energies[b].data();
        }
        gpu::RelaxMolecule rm;
        rm.natom        = ca.natom;
        rm.stride       = ca.stride;
        rm.nslot        = ca.nconf;
        rm.x            = ca.x.data();
        rm.f            = ca.f.data();
        rm.fixedForces  = m.cdata->fixed.arrayForces.data();
        rm.done         = m.done.data();
        rm.computed     = m.computed.data();
        rm.forceShells  = m.forceShells.data();
        rm.movingShells = m.movingShells.data();
        rm.fragments    = m.fragments.data();
        rm.nfragment    = m.fragments.size();
        rm.blocks       = m.blocks.data();
        rm.nblock       = m.blocks.size();
        rm.allForces    = allForces;
        rmols.push_back(rm);
    }
    static_assert(MAX_ENER == gpu::pairEnergyTerms, "The GPU has as many energy terms as MAX_ENER");
    gpu::RelaxSettings settings = { msForceToler, maxShellDistance2, minStep, maxStep, maxIter, stallIterations,
                                    shellWallAngleTolerance };
    bool notFinite = false;
    if (!gpu::relaxShells(rmols.data(), rmols.size(), settings, nComputed, &notFinite))
    {
        return false;
    }
    if (notFinite)
    {
        // As on the processor, see moveCoordinate
        moveCoordinate(0, std::numeric_limits<double>::quiet_NaN());
    }
    // The results as computeEnergies leaves them for the entries without
    // shells kept, see lowComputeEnergiesT
    for (auto &m : mols)
    {
        auto *cdata = m.cdata;
        auto &ca    = cdata->arrays;
        auto &cener = cdata->energies;
        auto &fixed = cdata->fixed;
        std::vector<size_t> confs, slots;
        for (size_t s = 0; s < ca.nconf; s++)
        {
            size_t c = ca.conformation(s);
            cdata->confs[c].done = m.done[s] != 0;
            if (m.computed[s])
            {
                confs.push_back(c);
                slots.push_back(s);
            }
        }
        // The active conformations are found again in the next calculation
        cdata->activeKnown = 0;
        if (confs.empty())
        {
            continue;
        }
        auto eEpot = cener.add(InteractionType::EPOT);
        for (auto c : confs)
        {
            eEpot.before[c] = eEpot.have[c] ? eEpot.before[c] : 0.0;
            eEpot.act[c]    = 0.0;
            eEpot.have[c]   = 1;
        }
        for (size_t k = 0; k < fixed.types.size(); k++)
        {
            auto t = cener.add(fixed.types[k]);
            for (auto c : confs)
            {
                t.before[c] = t.have[c] ? t.before[c] : 0.0;
                t.have[c]   = fixed.have[k][c];
                t.act[c]    = fixed.energies[k][c];
            }
        }
        cener.setAccumulate(true);
        for (size_t b = 0; b < m.blocks.size(); b++)
        {
            cener.setTermTypes(&energyTypes(m.potentials[b], m.itypes[b]));
            storeEnergies(m.potentials[b], confs.size(), confs.data(), slots.data(), &cener,
                          m.energies[b].data(), 1, ca.stride);
            if (Timings::instance().countFlops())
            {
                Timings::instance().addFlops(confs.size()*interactionFlops(m.potentials[b],
                                                                           m.top->shellEntries(m.itypes[b], m.top->entry(m.itypes[b]), true),
                                                                           m.top->atoms(), true));
            }
        }
        cener.setAccumulate(false);
        cener.setTermTypes(nullptr);
        // The van der Waals correction is part of the exchange, as in lowComputeEnergiesT
        auto vdwCorr = cener.find(InteractionType::VDWCORRECTION);
        if (vdwCorr.have)
        {
            auto exchange = cener.add(InteractionType::EXCHANGE);
            for (auto c : confs)
            {
                if (vdwCorr.has(c))
                {
                    exchange.set(c, (exchange.has(c) ? exchange.act[c] : 0.0) + vdwCorr.act[c]);
                    vdwCorr.act[c] = 0;
                }
            }
        }
    }
    return true;
}

#endif

namespace
{

//! The description of a topology and its conformations for gpu::createSystem
struct GpuSystemTopology
{
    const Topology                          *top   = nullptr;
    ConformationData                        *cdata = nullptr;
    std::vector<gpu::InteractionList>        lists;
    std::vector<std::vector<gpu::PairEntry>> pairs;
    std::vector<std::vector<gpu::BondedEntry>> bonds;
    std::vector<gpu::VsiteEntry>             vsites;
    //! For each list the row of each interaction in its parameter table, see gpu::updateParameters
    std::vector<std::vector<uint32_t>>       index;
    //! An interaction with the parameters of a row of a parameter table
    struct Row
    {
        //! A bonded interaction, or nullptr
        const TopologyEntry    *entry = nullptr;
        //! A non-bonded pair, if not entry
        std::optional<PairView> pair;
    };
    //! For each list an interaction of each row of its parameter table
    std::vector<std::vector<Row>>            rows;
    //! The kernel of each list
    std::vector<gpu::PairKernel>             pairKernels;
    std::vector<gpu::BondedKernel>           bondedKernels;
    //! The virtual sites and their types, to update their parameters
    std::vector<const TopologyEntry *>       vsiteEntries;
    std::vector<InteractionType>             vsiteTypes;
    //! The potential and interaction type of each energy array
    std::vector<Potential>                   potentials;
    std::vector<InteractionType>             itypes;

    //! Point the lists to their interactions, also after moving
    void setPointers()
    {
        for (size_t l = 0; l < lists.size(); l++)
        {
            lists[l].pairs  = pairs[l].data();
            lists[l].bonds  = bonds[l].data();
            lists[l].index  = index[l].empty() ? nullptr : index[l].data();
            lists[l].ntable = rows[l].size();
        }
    }
};

/*! \brief Describe a topology with its conformations for the GPU
 * \param[in]  pd       The force field
 * \param[in]  top      The topology
 * \param[in]  cdata    Its conformations
 * \param[out] t        The description
 * \return false if the GPU does not have a potential
 */
bool describeGpuTopology(const ForceField *pd, const Topology *top, ConformationData *cdata,
                         GpuSystemTopology *t)
{
    t->top   = top;
    t->cdata = cdata;
    double epsilonr = pd->epsilonr();
    for (auto itype : top->interactionTypes())
    {
        if (isVsite(itype))
        {
            for (const auto &b : top->entry(itype))
            {
                gpu::VsiteEntry ve;
                if (!gpuVsiteEntry(itype, *b.operator->(), &ve))
                {
                    return false;
                }
                t->vsites.push_back(ve);
                t->vsiteEntries.push_back(b.operator->());
                t->vsiteTypes.push_back(itype);
            }
            continue;
        }
        Potential mypot = pd->potential(itype);
        if (Potential::NONE == mypot)
        {
            mypot = defaultPotential(itype);
        }
        if (energyTypes(mypot).empty())
        {
            continue;
        }
        gpu::PairKernel   pk;
        gpu::BondedKernel bk;
        bool              pair   = gpuPairKernel(mypot, &pk);
        bool              bonded = !pair && gpuBondedKernel(mypot, &bk);
        if ((!pair && !bonded) ||
            energyTypes(mypot, itype).size() > static_cast<size_t>(gpu::pairEnergyTerms))
        {
            return false;
        }
        uint32_t energyIndex = t->potentials.size();
        t->potentials.push_back(mypot);
        t->itypes.push_back(itype);
        // The entries without and with shells, which are all without if
        // the topology does not split them. The non-bonded pairs are
        // stored compactly, see NonBondedPairs.
        const auto &nb = top->nonBondedPairs();
        for (int part = 0; part < 2; part++)
        {
            auto addList = [&](const auto &entries)
            {
                gpu::InteractionList il;
                il.bonded      = bonded;
                il.kernel      = pair ? static_cast<uint32_t>(pk) : static_cast<uint32_t>(bk);
                il.shells      = part == 1;
                il.energyIndex = energyIndex;
                il.nterm       = energyTypes(mypot, itype).size();
                t->pairs.emplace_back();
                t->bonds.emplace_back();
                t->index.emplace_back();
                t->rows.emplace_back();
                t->pairKernels.push_back(pk);
                t->bondedKernels.push_back(bk);
                // The interactions of the same types have the same parameters,
                // one row of the table each. Coulomb follows from the charges.
                std::map<std::string, uint32_t> rowOfKey;
                for (const auto &b : entries)
                {
                    if (!(pair && pk == gpu::PairKernel::Coulomb))
                    {
                        // A non-bonded pair has the type pair as key, which
                        // stands for the same as the identifier and the
                        // particle types of the entries
                        std::string key;
                        if constexpr (isPairRange<decltype(entries)>)
                        {
                            key = std::to_string(b.typeIndex());
                        }
                        else
                        {
                            key = b->id().id();
                            for (auto ai : b->atomIndices())
                            {
                                key += "|" + top->atoms()[ai].ffType();
                            }
                        }
                        auto rk = rowOfKey.find(key);
                        if (rk == rowOfKey.end())
                        {
                            rk = rowOfKey.emplace(key, t->rows.back().size()).first;
                            GpuSystemTopology::Row row;
                            if constexpr (isPairRange<decltype(entries)>)
                            {
                                row.pair = b;
                            }
                            else
                            {
                                row.entry = b.operator->();
                            }
                            t->rows.back().push_back(row);
                        }
                        t->index.back().push_back(rk->second);
                    }
                    if (pair)
                    {
                        // Interactions without parameters are kept, inactive,
                        // such that the lists keep their length when the
                        // parameters change, see gpu::PairEntry
                        gpu::PairEntry pe;
                        if (!gpuPairEntry(pk, *b.operator->(), top->atoms(), epsilonr, &pe))
                        {
                            pe.p[7] = 1;
                        }
                        t->pairs.back().push_back(pe);
                    }
                    else if constexpr (!isPairRange<decltype(entries)>)
                    {
                        gpu::BondedEntry be;
                        gpuBondedEntry(bk, *b.operator->(), &be);
                        t->bonds.back().push_back(be);
                    }
                }
                il.n = pair ? t->pairs.back().size() : t->bonds.back().size();
                t->lists.push_back(il);
            };
            if (nb.hasType(itype))
            {
                addList(top->hasShellSplit() ? nb.shellEntries(itype, part == 1) :
                        (part == 0 ? nb.range(itype) : PairRange()));
            }
            else
            {
                const auto &ev = top->entry(itype);
                addList(top->hasShellSplit() ?
                        top->shellEntries(itype, ev, part == 1) :
                        (part == 0 ? TopologyEntryRange(ev) : TopologyEntryRange()));
            }
        }
    }
    t->setPointers();
    return true;
}

} // namespace

bool checkGpuSystem(const ForceField *pd, const std::vector<ACTMol *> &allmols)
{
    if (!gpu::useGpu() || GMX_DOUBLE)
    {
        return false;
    }
    std::vector<GpuSystemTopology> tops;
    for (auto molptr : allmols)
    {
        auto &ftops = molptr->fragmentHandler()->topologies();
        for (size_t what = 0; what < 3; what++)
        {
            if ((what < 2 && what >= ftops.size()) || (what == 2 && ftops.size() < 2))
            {
                continue;
            }
            const Topology *top   = (what < 2) ? &ftops[what] : molptr->topology();
            auto           *cdata = molptr->conformationData(what);
            if (cdata->useArraysDouble() || cdata->arrays.nconf != cdata->nconf())
            {
                return false;
            }
            tops.emplace_back();
            if (!describeGpuTopology(pd, top, cdata, &tops.back()))
            {
                fprintf(stderr, "GPU system check: a potential is not available on the GPU\n");
                return false;
            }
        }
    }
    std::vector<gpu::SystemBlock> blocks;
    for (auto &t : tops)
    {
        auto &ca = t.cdata->arrays;
        gpu::SystemBlock sb;
        sb.natom   = ca.natom;
        sb.stride  = ca.stride;
        sb.nslot   = ca.nconf;
        sb.x       = ca.x.data();
        sb.lists   = t.lists.data();
        sb.nlist   = t.lists.size();
        sb.nenergy = t.potentials.size();
        sb.vsites  = t.vsites.data();
        sb.nvsite  = t.vsites.size();
        blocks.push_back(sb);
    }
    auto *system = gpu::createSystem(blocks.data(), blocks.size());
    if (!system)
    {
        return false;
    }
    std::vector<std::vector<double>> energies(tops.size());
    std::vector<double *>            eptr;
    for (size_t k = 0; k < tops.size(); k++)
    {
        energies[k].assign(tops[k].potentials.size()*gpu::pairEnergyTerms*blocks[k].stride, 0.0);
        eptr.push_back(energies[k].data());
    }
    gpu::computeAll(system, eptr.data());
    gpu::destroySystem(system);
    // The energy of each type of each conformation, as storeEnergies sums them
    std::map<InteractionType, std::pair<double, double>> worst;
    for (size_t k = 0; k < tops.size(); k++)
    {
        auto &t     = tops[k];
        auto &ca    = t.cdata->arrays;
        auto &cener = t.cdata->energies;
        std::map<InteractionType, std::vector<double>> gpuE;
        for (size_t e = 0; e < t.potentials.size(); e++)
        {
            const auto &types = energyTypes(t.potentials[e], t.itypes[e]);
            for (size_t term = 0; term < types.size(); term++)
            {
                auto &v = gpuE[types[term]];
                v.resize(ca.nconf, 0.0);
                for (size_t sl = 0; sl < ca.nconf; sl++)
                {
                    v[ca.conformation(sl)] += energies[k][(e*gpu::pairEnergyTerms + term)*ca.stride + sl];
                }
            }
        }
        for (const auto &g : gpuE)
        {
            auto cpu = cener.find(g.first);
            if (!cpu.have)
            {
                continue;
            }
            for (size_t c = 0; c < ca.nconf; c++)
            {
                if (!cpu.has(c))
                {
                    continue;
                }
                double diff = std::abs(g.second[c] - cpu.act[c]);
                double rel  = diff/(1 + std::abs(cpu.act[c]));
                auto  &w    = worst[g.first];
                if (rel > w.first)
                {
                    w = { rel, cpu.act[c] };
                }
            }
        }
    }
    for (const auto &w : worst)
    {
        fprintf(stderr, "GPU system check: %-22s largest relative difference %.2e (of %.6g)\n",
                interactionTypeToString(w.first).c_str(), w.second.first, w.second.second);
    }
    return true;
}

#if GMX_DOUBLE

// The GPU computes in single precision only
class GpuTrainingSystem::Impl
{
};

GpuTrainingSystem::GpuTrainingSystem() : impl_(new Impl) {}
GpuTrainingSystem::~GpuTrainingSystem() = default;
bool GpuTrainingSystem::prepare(const ForceField *, const std::vector<ACTMol *> &, double, bool, const ForceComputer *, size_t) { return false; }
size_t GpuTrainingSystem::nreplica() const { return 1; }
void GpuTrainingSystem::updateParameters(const ForceField *, const ForceComputer *) {}
void GpuTrainingSystem::stageParameters(size_t, const ForceField *, const ForceComputer *) {}
void GpuTrainingSystem::uploadParameters(const ForceField *) {}
void GpuTrainingSystem::selectReplica(size_t) {}
std::pair<double, double> GpuTrainingSystem::chargeDifference() const { return { 0, 0 }; }
void GpuTrainingSystem::copyChargesToMolecules() {}
void GpuTrainingSystem::resetShells() {}
void GpuTrainingSystem::keepShells() {}
bool GpuTrainingSystem::restoreShells(int) { return false; }
void GpuTrainingSystem::constructVsites() {}
void GpuTrainingSystem::setKeepEnergies(bool) {}
bool GpuTrainingSystem::typeEnergies(const ACTMol *, int, std::vector<GpuTypeEnergies> *) const { return false; }
bool GpuTrainingSystem::setDeviations(const GpuDeviationSetup &) { return false; }
bool GpuTrainingSystem::deviationsReady() const { return false; }
void GpuTrainingSystem::setDeviationsOnly(bool) {}
void GpuTrainingSystem::computeDeviations() {}
bool GpuTrainingSystem::deviationSums(std::vector<std::pair<double, double>> *) const { return false; }
bool GpuTrainingSystem::moleculeDeviations(size_t, GpuMoleculeDeviations *) const { return false; }
void GpuTrainingSystem::setActiveMolecules(int, const std::vector<char> &) {}
std::vector<uint64_t> GpuTrainingSystem::moleculeSignatures(const ForceField *, const ForceComputer *) { return {}; }
size_t GpuTrainingSystem::version() const { return 0; }
bool GpuTrainingSystem::recipeRows() const { return false; }
size_t GpuTrainingSystem::deviationsVersion() const { return 0; }
void GpuTrainingSystem::compute(int, bool, bool) {}
void GpuTrainingSystem::computeForces(int) {}
bool GpuTrainingSystem::haveKeptShells() const { return false; }
size_t GpuTrainingSystem::relax(int, const bool[2], const gpu::RelaxSettings &, size_t *nFailed) { *nFailed = 0; return 0; }
void GpuTrainingSystem::copyRelaxedMonomerShells() {}
void GpuTrainingSystem::saveDimers() {}
void GpuTrainingSystem::restoreDimers() {}

#else

class GpuTrainingSystem::Impl
{
public:
    //! A topology with its conformations on the GPU
    struct Block
    {
        //! The molecule and 0, 1 for its monomers or 2 for the dimer
        size_t              mol  = 0;
        int                 what = 0;
        //! The replica, see GpuTrainingSystem
        size_t              replica = 0;
        GpuSystemTopology   desc;
        //! The shells, see relaxShellsGpu
        GpuRelaxMolecule    shells;
        //! The shells with one parent and that parent
        std::vector<uint32_t> resets;
        //! The conformation of each slot and the slot of each conformation
        std::vector<size_t> confOfSlot, slotOfConf;
        //! The energies as the GPU gave them, after and before relaxing, see setKeepEnergies
        std::vector<double> rawAfter, rawBefore;
        /*! \brief The Topology::parameterVersion of the interaction type of each
         * list when its rows were staged last, 0 if they have to be made,
         * see stageRows
         */
        std::vector<size_t> rowVersion;

    };
    std::vector<Block>    blocks;
    gpu::System          *system = nullptr;
    //! Whether the shells of all blocks were kept, see keepShells
    bool                  shellsKept = false;
    //! Whether a displacement was not finite since the shells were kept
    bool                  notFinite = false;
    std::vector<ACTMol *> mols;
    std::vector<size_t>   nconf;
    double                maxShellDistance2 = 0;
    //! The molecules whose charges the GPU generates, see gpu::createAcm
    std::vector<size_t>   acmMols;
    //! The force field parameters of the charge generation, see gpu::updateAcm
    std::vector<const ForceFieldParameter *> acmValues;
    //! Whether the structures of the charge generation have virtual sites
    bool                  acmVsites = false;
    //! The molecule, replica and first particle of each molecule of the charge generation, in the order of gpu::createAcm
    struct AcmMoleculeOf
    {
        size_t mol = 0, replica = 0, firstAtom = 0;
    };
    std::vector<AcmMoleculeOf> acmMolOf;
    /*! \brief The values of the recipes of the virtual sites of the charge
     * generation, and those last sent, if the GPU constructs them, see
     * makeAcmVsites
     */
    std::vector<size_t>   acmVsiteSlots;
    std::vector<double>   acmVsiteLast;
    bool                  acmVsitesOnGpu = false;
    //! The values and coordinates last sent, and whether the charges were generated since copied
    std::vector<double>   acmLastValues, acmLastX;
    bool                  acmGenerated = false;
    //! The number of replicas and of the blocks of each, those of a replica follow each other
    size_t                nreplica = 1, nblockReplica = 0;
    //! The replica whose energies are in the molecules, see selectReplica
    size_t                resident = 0;
    //! The energies of the replicas that are not in the molecules, three per molecule
    std::vector<std::vector<ConformationEnergies>> replicaEnergies;
    //! The parameters of each replica, see stageParameters
    struct Staged
    {
        std::vector<float>            tables;
        std::vector<double>           charges;
        std::vector<gpu::MovingShell> shells;
        std::vector<gpu::VsiteEntry>  vsites;
        std::vector<double>           acmValues, acmX;
    };
    std::vector<Staged>   staged;
    //! The charges generated on the GPU, of all replicas, see copyChargesToMolecules
    std::vector<double>   acmCharges;
    //! The number of atoms of the ACM molecules of a replica
    size_t                nacmAtomReplica = 0;
    //! The replica whose generated charges are in the molecules, none if larger than nreplica
    size_t                chargesOf = std::numeric_limits<size_t>::max();

    //! The coordinates of the ACM particles, see gpu::AcmSystem::x
    std::vector<double> acmCoordinates(const ForceComputer *fc) const
    {
        std::vector<double> xa;
        for (auto k : acmMols)
        {
            auto  molptr = mols[k];
            auto  fh     = molptr->fragmentHandler();
            // The structure the processor generates the charges for, see ACTMol::GenerateAcmCharges
            std::vector<gmx::RVec> x = molptr->xOriginal();
            if (molptr->topology()->hasVsites())
            {
                Coordinates cx(x);
                fc->constructVsiteCoordinates(molptr->topology(), &cx);
                x = cx.toRVec();
            }
            for (size_t f = 0; f < fh->topologies().size(); f++)
            {
                for (auto a : fh->acmParticles(f))
                {
                    for (int d = 0; d < DIM; d++)
                    {
                        xa.push_back(x[fh->atomStart()[f] + a][d]);
                    }
                }
            }
        }
        return xa;
    }

    ~Impl()
    {
        if (system)
        {
            gpu::destroySystem(system);
        }
    }

    //! The block of each molecule, topology and replica, see find
    std::vector<int> blockIndex;
    //! Whether to keep the energies as the GPU gives them, see setKeepEnergies
    bool             keepEnergies = false;
    //! The deviations, see setDeviations: whether they are set, the setup, the number of targets and the results
    bool                devReady = false;
    GpuDeviationSetup   devSetup;
    //! The inverse temperature and the version of the conformations of each reference, which change in place
    std::vector<std::pair<double, size_t>> devRefState;
    size_t              devTargets = 0;
    bool                devComputed = false;
    //! The dimers and topologies of the deviations, and the index of each of a replica and molecule, -1 if none
    std::vector<gpu::DeviationDimer>   devDimers;
    std::vector<gpu::DeviationMonomer> devMonomers;
    std::vector<int>                   devDimerOf, devMonomerOf;
    //! The results of computeDeviations, see gpu::computeDeviations
    std::vector<double>                devDimerSums, devEpotTerms;
    //! Whether each block is computed, see setActiveMolecules
    std::vector<char>                  active;
    //! Where the parameters of each block are in those of its replica, see stageParameters
    struct Segment
    {
        size_t table = 0, charge = 0, shell = 0, vsite = 0;
        size_t ntable = 0, ncharge = 0, nshell = 0, nvsite = 0;
    };
    std::vector<Segment>               segments;
    //! Whether the parameters of all blocks of each replica were staged
    std::vector<char>                  stagedAll;
    //! The values of the charge generation of each molecule, see acmValues
    std::vector<std::vector<uint32_t>> molValues;
    //! Changes when the system is made again, and when the deviations are, see version
    size_t                             versionCount = 0, devVersionCount = 0;
    //! Whether compute keeps the energies on the GPU only, see setDeviationsOnly
    bool                devOnly = false;
    /*! \brief The recipes of the rows of the parameter tables, see makeRowRecipes
     * With them the GPU makes the rows from the values of the force field,
     * instead of stageRows. recipeRows is 0 until it is known whether they
     * can be used, then 1 if they are, -1 if not.
     */
    ParameterRecipes                   recipes;
    int                                recipeRows = 0;
    //! The values of the recipes of each replica, and whether those of each were copied, see stageParameters
    std::vector<double>                recipeValues;
    std::vector<char>                  recipeStaged;
    //! Whether stageParameters makes the rows on the processor anyway, see moleculeSignatures
    bool                               hostRows = false;

    /*! \brief Make the recipes of all rows of the parameter tables, see gpu::createRowRecipes
     * \param[in] pd The force field
     * \return whether all rows have a recipe, except those of Coulomb, which
     *         the GPU makes from the charges
     */
    bool makeRowRecipes(const ForceField *pd)
    {
        recipes = ParameterRecipes();
        std::vector<gpu::RowRecipe> rowRecipes;
        size_t tablesReplica  = staged[0].tables.size()/8;
        size_t chargesReplica = staged[0].charges.size()/4;
        for (size_t r = 1; r < nreplica; r++)
        {
            if (staged[r].tables.size()/8 != tablesReplica || staged[r].charges.size()/4 != chargesReplica)
            {
                return false;
            }
        }
        for (size_t b = 0; b < blocks.size(); b++)
        {
            const auto &bl  = blocks[b];
            const auto &d   = bl.desc;
            const auto &seg = segments[b - bl.replica*nblockReplica];
            size_t      row = bl.replica*tablesReplica + seg.table/8;
            for (size_t l = 0; l < d.lists.size(); l++)
            {
                auto itype  = d.itypes[d.lists[l].energyIndex];
                bool bonded = d.lists[l].bonded;
                for (const auto &rw : d.rows[l])
                {
                    gpu::RowRecipe rr;
                    rr.row     = row++;
                    rr.replica = static_cast<uint32_t>(bl.replica);
                    rr.bonded  = bonded ? 1 : 0;
                    if (bonded)
                    {
                        rr.kernel = static_cast<uint32_t>(d.bondedKernels[l]);
                        rr.set    = rw.entry ? recipes.directSet(pd, itype, rw.entry->id()) : -1;
                    }
                    else
                    {
                        auto kernel = d.pairKernels[l];
                        if (gpu::PairKernel::Coulomb == kernel)
                        {
                            // From the charges, see gpu::updateParameters
                            continue;
                        }
                        rr.kernel = static_cast<uint32_t>(kernel);
                        rr.set    = rw.pair ? recipes.pairSet(pd, itype, rw.pair->id()) :
                            (rw.entry ? recipes.directSet(pd, itype, rw.entry->id()) : -1);
                        if (gpu::PairKernel::Polarization == kernel)
                        {
                            // The shell is the second particle, as gpuPairEntry
                            int shell = rw.pair ? rw.pair->atomIndex(1) : rw.entry->atomIndices()[1];
                            rr.charge = static_cast<int64_t>(4*(bl.replica*chargesReplica + seg.charge/4 + shell));
                        }
                    }
                    if (rr.set < 0)
                    {
                        return false;
                    }
                    rowRecipes.push_back(rr);
                }
            }
        }
        // The moving shells and the virtual sites, in the order of stageShellsVsites
        std::vector<gpu::ShellRecipe> shellRecipes;
        std::vector<gpu::VsiteRecipe> vsiteRecipes;
        size_t shellsReplica = staged[0].shells.size();
        size_t vsitesReplica = staged[0].vsites.size();
        for (size_t b = 0; b < blocks.size(); b++)
        {
            const auto      &bl  = blocks[b];
            const auto      &d   = bl.desc;
            const auto      &seg = segments[b - bl.replica*nblockReplica];
            GpuRelaxMolecule sh;
            shellsOf(bl, &sh);
            if (sh.movingShells.size() != seg.nshell || d.vsiteEntries.size() != seg.nvsite)
            {
                return false;
            }
            for (size_t k = 0; k < sh.movingShells.size(); k++)
            {
                gpu::ShellRecipe sr;
                sr.shell   = bl.replica*shellsReplica + seg.shell + k;
                sr.replica = static_cast<uint32_t>(bl.replica);
                sr.set     = recipes.directSet(pd, InteractionType::POLARIZATION, sh.movingShellEntries[k]->id());
                sr.charge  = static_cast<int64_t>(4*(bl.replica*chargesReplica + seg.charge/4 + sh.movingShells[k].shell));
                if (sr.set < 0)
                {
                    return false;
                }
                shellRecipes.push_back(sr);
            }
            for (size_t v = 0; v < d.vsiteEntries.size(); v++)
            {
                gpu::VsiteRecipe vr;
                vr.vsite   = bl.replica*vsitesReplica + seg.vsite + v;
                vr.replica = static_cast<uint32_t>(bl.replica);
                double factor[3];
                if (!gpuVsiteParameters(d.vsiteTypes[v], *d.vsiteEntries[v], vr.index, factor))
                {
                    return false;
                }
                std::copy(factor, factor + 3, vr.factor);
                if (std::all_of(vr.index, vr.index + 3, [](int i) { return i < 0; }))
                {
                    // Without parameters
                    continue;
                }
                vr.set = recipes.directSet(pd, d.vsiteTypes[v], d.vsiteEntries[v]->id());
                if (vr.set < 0)
                {
                    return false;
                }
                vsiteRecipes.push_back(vr);
            }
        }
        gpu::createRowRecipes(system, recipes.recipes().data(), recipes.recipes().size(),
                              rowRecipes.data(), rowRecipes.size(), recipes.factors().data(),
                              recipes.values().size(), nreplica);
        gpu::createShellVsiteRecipes(system, shellRecipes.data(), shellRecipes.size(),
                                     vsiteRecipes.data(), vsiteRecipes.size());
        if (!makeAcmVsites(pd))
        {
            return false;
        }
        recipeValues.assign(nreplica*recipes.values().size(), 0.0);
        recipeStaged.assign(nreplica, 0);
        return true;
    }

    /*! \brief Let the GPU construct the virtual sites of the charge generation, see gpu::createAcmVsites
     * As acmCoordinates, from the original structures, with the parameters
     * of the recipes, in the order of the processor, see constructStrided in
     * vsitehandler.cpp.
     * \param[in] pd The force field
     * \return whether all virtual sites have a recipe
     */
    bool makeAcmVsites(const ForceField *pd)
    {
        acmVsiteSlots.clear();
        acmVsiteLast.clear();
        acmVsitesOnGpu = false;
        if (!acmVsites || acmMolOf.empty())
        {
            return true;
        }
        std::vector<int64_t>         x;
        std::vector<uint32_t>        start;
        std::vector<gpu::AcmVsite>   vsites;
        std::set<size_t>             slots;
        size_t                       nvalue = recipes.values().size();
        for (const auto &am : acmMolOf)
        {
            auto        molptr = mols[am.mol];
            const auto *top    = molptr->topology();
            Coordinates cx(molptr->xOriginal());
            const auto *c      = cx.view().data();
            GMX_RELEASE_ASSERT(3*am.firstAtom == x.size(), "The molecules of the charge generation are not in order");
            x.insert(x.end(), c, c + 3*molptr->xOriginal().size());
            start.push_back(vsites.size());
            for (const auto &entry : top->entries())
            {
                if (!isVsite(entry.first))
                {
                    continue;
                }
                for (const auto &vs : entry.second)
                {
                    gpu::AcmVsite av;
                    av.replica = static_cast<uint32_t>(am.replica);
                    double factor[3];
                    if (!gpuVsiteEntry(entry.first, *vs.operator->(), &av.entry) ||
                        !gpuVsiteParameters(entry.first, *vs.operator->(), av.index, factor))
                    {
                        return false;
                    }
                    std::copy(factor, factor + 3, av.factor);
                    if (std::any_of(av.index, av.index + 3, [](int i) { return i >= 0; }))
                    {
                        av.set = recipes.directSet(pd, entry.first, vs->id());
                        if (av.set < 0)
                        {
                            return false;
                        }
                        const auto &sr = recipes.recipes()[av.set];
                        for (int k = 0; k < sr.n; k++)
                        {
                            if (sr.slotI[k] >= 0)
                            {
                                slots.insert(am.replica*nvalue + sr.slotI[k]);
                            }
                        }
                    }
                    else
                    {
                        av.set = 0;
                    }
                    vsites.push_back(av);
                }
            }
        }
        start.push_back(vsites.size());
        gpu::createAcmVsites(system, x.data(), start.data(), vsites.data(), vsites.size());
        acmVsiteSlots.assign(slots.begin(), slots.end());
        acmVsitesOnGpu = true;
        return true;
    }

    //! Make blockIndex, after the blocks were made
    void indexBlocks()
    {
        blockIndex.assign(3*mols.size()*nreplica, -1);
        for (size_t b = 0; b < blocks.size(); b++)
        {
            blockIndex[(blocks[b].replica*mols.size() + blocks[b].mol)*3 + blocks[b].what] = b;
        }
    }

    //! The index of the block of a molecule and topology, -1 if none
    int find(size_t mol, int what, size_t replica = 0) const
    {
        return blockIndex[(replica*mols.size() + mol)*3 + what];
    }

    //! The energies of a molecule and topology of a replica, as energiesOf
    ConformationEnergies &energiesOf(size_t mol, int what, size_t replica)
    {
        if (replica == resident)
        {
            return mols[mol]->conformationData(what)->energies;
        }
        return replicaEnergies[replica][3*mol + what];
    }

    //! Whether to use each block for what, -1 for all, -2 for the monomers
    std::vector<char> use(int what) const
    {
        std::vector<char> u(blocks.size(), 0);
        for (size_t b = 0; b < blocks.size(); b++)
        {
            u[b] = ((what == -1 || blocks[b].what == what || (what == -2 && blocks[b].what < 2)) &&
                    (active.empty() || active[b])) ? 1 : 0;
        }
        return u;
    }

    //! The shells of a block, see addShellFragment
    void shellsOf(const Block &bl, GpuRelaxMolecule *shells) const
    {
        auto  molptr = mols[bl.mol];
        auto  fh     = molptr->fragmentHandler();
        auto &tops   = fh->topologies();
        auto &pols   = fh->polInfo();
        if (bl.what < 2)
        {
            addShellFragment(&tops[bl.what], pols[bl.what], 0, maxShellDistance2, shells);
        }
        else
        {
            for (size_t f = 0; f < 2; f++)
            {
                addShellFragment(&tops[f], pols[f], fh->atomStart()[f], maxShellDistance2, shells);
            }
        }
    }

    //! The energies of a block, in its molecule if its replica is resident
    ConformationEnergies &energiesOf(const Block &bl)
    {
        if (bl.replica == resident)
        {
            return mols[bl.mol]->conformationData(bl.what)->energies;
        }
        return replicaEnergies[bl.replica][3*bl.mol + bl.what];
    }

    //! The charges generated for a replica to the molecules and their fragments
    void applyCharges(size_t replica)
    {
        if (acmCharges.empty() || chargesOf == replica)
        {
            return;
        }
        chargesOf    = replica;
        size_t first = replica*nacmAtomReplica;
        // The fixed parts of the charges of the replica as it was staged,
        // the second charge and the distribution widths. The molecules have
        // those of the replica staged last, but the bounds, for instance,
        // need the total charges of this one, see ActAtom::totalCharge.
        const auto *st = replica < staged.size() && stagedAll[replica] ? &staged[replica] : nullptr;
        for (auto k : acmMols)
        {
            auto  molptr = mols[k];
            auto  fh     = molptr->fragmentHandler();
            auto &tops   = fh->topologiesPtr();
            auto  atoms  = molptr->atoms();
            for (size_t f = 0; f < tops.size(); f++)
            {
                auto   fatoms = tops[f].atomsPtr();
                size_t start  = fh->atomStart()[f];
                // Four values for each atom of the block of the fragment, see stageCharges
                const double *q4 = nullptr;
                int           b  = find(k, f);
                if (st && b >= 0 && static_cast<size_t>(b) < segments.size() &&
                    segments[b].ncharge == 4*fatoms->size())
                {
                    q4 = st->charges.data() + segments[b].charge;
                }
                for (size_t a = 0; a < fatoms->size(); a++)
                {
                    for (auto atom : { &(*fatoms)[a], &(*atoms)[start + a] })
                    {
                        atom->setCharge(acmCharges[first + start + a]);
                        if (q4)
                        {
                            atom->setChargeVS(q4[4*a + 1], q4[4*a + 3], atom->rowVS());
                            atom->setZeta(q4[4*a + 2]);
                        }
                    }
                }
            }
            first += atoms->size();
        }
    }

    /*! \brief The energies of a block to its conformation data, as lowComputeEnergiesT stores them
     * \param[in] bl       The block
     * \param[in] energies Its energies, see gpu::computeBlocksInPlace
     * \param[in] backup   Whether the energies are their backups as well, as
     *                     backupEnergies makes them of the energies before relaxing
     */
    void storeBlockEnergies(const Block &bl, const double *energies, bool backup)
    {
        auto       *cdata = mols[bl.mol]->conformationData(bl.what);
        auto       &cener = energiesOf(bl);
        size_t      nslot = bl.confOfSlot.size();
        size_t      st    = cdata->arrays.stride;
        cener.setNumberOfConformations(cdata->nconf());
        std::vector<size_t> slots(nslot);
        for (size_t s = 0; s < nslot; s++)
        {
            slots[s] = s;
        }
        auto eEpot = cener.add(InteractionType::EPOT);
        for (auto c : bl.confOfSlot)
        {
            eEpot.before[c] = eEpot.have[c] ? eEpot.before[c] : 0.0;
            eEpot.act[c]    = 0.0;
            eEpot.have[c]   = 1;
        }
        // The types stored, which are their own backups with backup
        std::vector<char> stored(static_cast<size_t>(InteractionType::EPOT) + 1, 0);
        stored[static_cast<size_t>(InteractionType::EPOT)] = 1;
        for (size_t e = 0; e < bl.desc.potentials.size(); e++)
        {
            const auto &types = energyTypes(bl.desc.potentials[e], bl.desc.itypes[e]);
            cener.setTermTypes(&types);
            storeEnergies(bl.desc.potentials[e], nslot, bl.confOfSlot.data(), slots.data(), &cener,
                          energies + e*gpu::pairEnergyTerms*st, 1, st, backup);
            for (auto t : types)
            {
                stored.resize(std::max(stored.size(), static_cast<size_t>(t) + 1), 0);
                stored[static_cast<size_t>(t)] = 1;
            }
        }
        cener.setTermTypes(nullptr);
        // The van der Waals correction is part of the exchange, as in lowComputeEnergiesT
        auto vdwCorr = cener.find(InteractionType::VDWCORRECTION);
        if (vdwCorr.have)
        {
            auto exchange = cener.add(InteractionType::EXCHANGE);
            for (auto c : bl.confOfSlot)
            {
                if (vdwCorr.has(c))
                {
                    exchange.set(c, (exchange.has(c) ? exchange.act[c] : 0.0) + vdwCorr.act[c]);
                    vdwCorr.act[c] = 0;
                }
            }
            // Backed up below, after the correction moved
            for (auto t : { InteractionType::EXCHANGE, InteractionType::VDWCORRECTION })
            {
                if (static_cast<size_t>(t) < stored.size())
                {
                    stored[static_cast<size_t>(t)] = 0;
                }
            }
        }
        if (backup)
        {
            // The other types as backupEnergies backs them up
            cener.setBackups(true);
            for (auto t : cener.types())
            {
                if (static_cast<size_t>(t) >= stored.size() || !stored[static_cast<size_t>(t)])
                {
                    auto et = cener.find(t);
                    std::copy(et.act, et.act + cener.nconf(), et.before);
                }
            }
        }
    }
};

GpuTrainingSystem::GpuTrainingSystem() : impl_(new Impl) {}

GpuTrainingSystem::~GpuTrainingSystem() = default;

bool GpuTrainingSystem::prepare(const ForceField *pd, const std::vector<ACTMol *> &allmols, double maxShellDistance,
                                bool charges, const ForceComputer *fc, size_t nreplica)
{
    auto &im = *impl_;
    nreplica  = std::max<size_t>(1, nreplica);
    bool same = im.system && im.mols == allmols && im.nreplica == nreplica;
    for (size_t k = 0; same && k < allmols.size(); k++)
    {
        same = im.nconf[k] == allmols[k]->conformationData(0)->nconf();
    }
    if (same)
    {
        return true;
    }
    if (im.system)
    {
        gpu::destroySystem(im.system);
        im.system = nullptr;
    }
    // The kept shells were those of the old system
    im.shellsKept = false;
    im.notFinite  = false;
    im.blocks.clear();
    im.acmMols.clear();
    im.acmValues.clear();
    im.acmVsites = false;
    im.acmLastValues.clear();
    im.acmLastX.clear();
    im.acmGenerated = false;
    // The energies of the other replicas are lost, those of the resident
    // one are in the molecules
    im.nreplica        = nreplica;
    im.nblockReplica   = 0;
    im.resident        = 0;
    im.replicaEnergies.clear();
    im.staged.assign(nreplica, Impl::Staged());
    im.acmCharges.clear();
    im.nacmAtomReplica = 0;
    im.chargesOf       = std::numeric_limits<size_t>::max();
    im.devReady        = false;
    im.devComputed     = false;
    im.active.clear();
    im.segments.clear();
    im.stagedAll.assign(nreplica, 0);
    // The recipes of the rows are made again for the new blocks
    im.recipeRows = 0;
    im.recipeValues.clear();
    im.recipeStaged.clear();
    im.recipes = ParameterRecipes();
    im.molValues.assign(allmols.size(), {});
    im.versionCount   += 1;
    im.mols  = allmols;
    im.nconf.clear();
    im.maxShellDistance2 = maxShellDistance*maxShellDistance;
    if (!gpu::useGpu() || !pd->polarizable())
    {
        return false;
    }
    for (size_t k = 0; k < allmols.size(); k++)
    {
        auto  molptr = allmols[k];
        auto &tops   = molptr->fragmentHandler()->topologies();
        im.nconf.push_back(molptr->conformationData(0)->nconf());
        for (int what = 0; what < 3; what++)
        {
            if ((what < 2 && what >= static_cast<int>(tops.size())) || (what == 2 && tops.size() < 2))
            {
                continue;
            }
            const Topology *top   = (what < 2) ? &tops[what] : molptr->topology();
            auto           *cdata = molptr->conformationData(what);
            auto           &ca    = cdata->arrays;
            if (cdata->useArraysDouble() || ca.nconf != cdata->nconf())
            {
                return false;
            }
            Impl::Block bl;
            bl.mol  = k;
            bl.what = what;
            if (!describeGpuTopology(pd, top, cdata, &bl.desc))
            {
                return false;
            }
            im.shellsOf(bl, &bl.shells);
            const auto &atoms = top->atoms();
            for (size_t i = 0; i < atoms.size(); i++)
            {
                if (atoms[i].pType() == ActParticle::Shell && atoms[i].parents().size() == 1)
                {
                    bl.resets.push_back(i);
                    bl.resets.push_back(atoms[i].parents()[0]);
                }
            }
            bl.confOfSlot.resize(ca.nconf);
            bl.slotOfConf.resize(ca.nconf);
            for (size_t s = 0; s < ca.nconf; s++)
            {
                bl.confOfSlot[s] = ca.conformation(s);
                bl.slotOfConf[ca.conformation(s)] = s;
            }
            im.blocks.push_back(std::move(bl));
        }
    }
    // The blocks of the other replicas, which start with the energies that
    // the molecules have now
    im.nblockReplica = im.blocks.size();
    for (size_t r = 1; r < nreplica; r++)
    {
        for (size_t b = 0; b < im.nblockReplica; b++)
        {
            Impl::Block copy = im.blocks[b];
            copy.replica     = r;
            im.blocks.push_back(std::move(copy));
        }
    }
    im.replicaEnergies.resize(nreplica);
    for (size_t r = 1; r < nreplica; r++)
    {
        for (auto molptr : allmols)
        {
            for (int what = 0; what < 3; what++)
            {
                im.replicaEnergies[r].push_back(molptr->conformationData(what)->energies);
            }
        }
    }
    im.replicaEnergies[0].resize(3*allmols.size());
    im.indexBlocks();
    std::vector<gpu::SystemBlock> sblocks;
    for (auto &bl : im.blocks)
    {
        // The pointers of the description, after the blocks were moved
        bl.desc.setPointers();
        auto &ca = im.mols[bl.mol]->conformationData(bl.what)->arrays;
        gpu::SystemBlock sb;
        sb.natom        = ca.natom;
        sb.stride       = ca.stride;
        sb.nslot        = ca.nconf;
        sb.x            = ca.x.data();
        sb.lists        = bl.desc.lists.data();
        sb.nlist        = bl.desc.lists.size();
        sb.nenergy      = bl.desc.potentials.size();
        sb.vsites       = bl.desc.vsites.data();
        sb.nvsite       = bl.desc.vsites.size();
        sb.forceShells  = bl.shells.forceShells.data();
        sb.movingShells = bl.shells.movingShells.data();
        sb.fragments    = bl.shells.fragments.data();
        sb.nfragment    = bl.shells.fragments.size();
        sb.shellParents = bl.resets.data();
        sb.nshellParent = bl.resets.size()/2;
        // The other replicas share the interactions of the first on the GPU
        if (bl.replica > 0)
        {
            sb.topologyOf = static_cast<int64_t>(sblocks.size() - bl.replica*im.nblockReplica);
        }
        sblocks.push_back(sb);
    }
    im.system = gpu::createSystem(sblocks.data(), sblocks.size());
    if (im.system && charges)
    {
        // The fragments of the molecules, the ACM particles and bonds, and
        // where their parameters are in the values
        std::vector<gpu::AcmMolecule> am;
        std::vector<gpu::AcmFragment> af;
        std::vector<uint32_t>         ap, ab, pv, bv, fv;
        std::vector<double>           bfac;
        std::vector<int32_t>          sym;
        // The values of each replica follow those of the previous one
        std::map<const ForceFieldParameter *, uint32_t> valueIndex;
        size_t replica = 0;
        size_t molecule = 0;
        auto valueOf = [&valueIndex, &im, &replica, &molecule](const ForceFieldParameter *p)
        {
            if (!p)
            {
                return gpu::noValue;
            }
            auto vi = valueIndex.find(p);
            if (vi == valueIndex.end())
            {
                vi = valueIndex.emplace(p, im.acmValues.size()).first;
                im.acmValues.push_back(p);
            }
            if (replica == 0)
            {
                im.molValues[molecule].push_back(vi->second);
            }
            return static_cast<uint32_t>(vi->second + replica*im.acmValues.size());
        };
        for (size_t kr = 0; kr < nreplica*allmols.size(); kr++)
        {
            size_t      k      = kr % allmols.size();
            replica            = kr / allmols.size();
            molecule           = k;
            if (replica > 0 && k == 0 && im.nacmAtomReplica == 0)
            {
                im.nacmAtomReplica = sym.size();
            }
            auto        molptr = allmols[k];
            const auto *fh     = molptr->fragmentHandler();
            auto       &tops   = fh->topologies();
            if (!fh->sqeOnGpu(pd) || tops.size() > 2)
            {
                continue;
            }
            gpu::AcmMolecule m;
            m.natom         = molptr->atomsConst().size();
            m.firstAtom     = sym.size();
            m.firstFragment = af.size();
            m.nfragment     = tops.size();
            int dimer       = im.find(k, 2, replica);
            m.dimerBlock    = dimer >= 0 ? static_cast<uint32_t>(dimer) : gpu::noBlock;
            for (size_t f = 0; f < tops.size(); f++)
            {
                int block = im.find(k, f, replica);
                if (block < 0)
                {
                    return false;
                }
                gpu::AcmFragment fr;
                fr.natom      = tops[f].atoms().size();
                fr.atomStart  = fh->atomStart()[f];
                fr.moleculeAtom = m.firstAtom;
                fr.block      = block;
                fr.qtotal     = fh->fragmentCharge(f);
                const auto &acm = fh->acmParticles(f);
                fr.firstAcm   = ap.size();
                fr.nacm       = acm.size();
                fr.firstBond  = ab.size()/2;
                fr.nbond      = 0;
                fr.firstFixed = fv.size();
                fr.nfixed     = 0;
                if (!acm.empty())
                {
                    ap.insert(ap.end(), acm.begin(), acm.end());
                    QgenAcm::SqeSources src;
                    fh->sqeSources(f, &src);
                    for (size_t i = 0; i < acm.size(); i++)
                    {
                        pv.push_back(valueOf(src.chi[i]));
                        pv.push_back(valueOf(src.eta[i]));
                        pv.push_back(valueOf(src.zeta[i]));
                    }
                    const auto &bonds = fh->bonds()[f];
                    for (size_t b = 0; b < bonds.size(); b++)
                    {
                        ab.push_back(bonds[b].aI());
                        ab.push_back(bonds[b].aJ());
                        bv.push_back(valueOf(src.dchi[b]));
                        bv.push_back(valueOf(src.deta[b]));
                        bfac.push_back(src.dchiFactor[b]);
                    }
                    fr.nbond = bonds.size();
                    for (auto q : src.fixedCharges)
                    {
                        fv.push_back(valueOf(q));
                    }
                    fr.nfixed = src.fixedCharges.size();
                }
                af.push_back(fr);
            }
            const auto &sc = molptr->symmetricCharges();
            for (size_t i = 0; i < m.natom; i++)
            {
                // Without symmetry information, each charge is its own
                sym.push_back(sc.size() == m.natom ? sc[i] : -1 - static_cast<int32_t>(i));
            }
            am.push_back(m);
            im.acmMolOf.push_back({ k, replica, m.firstAtom });
            if (replica == 0)
            {
                im.acmMols.push_back(k);
                im.acmVsites = im.acmVsites || acmHasVsites(molptr);
            }
        }
        if (nreplica == 1)
        {
            im.nacmAtomReplica = sym.size();
        }
        std::vector<double> xa;
        for (size_t r = 0; r < nreplica; r++)
        {
            auto xr = im.acmCoordinates(fc);
            xa.insert(xa.end(), xr.begin(), xr.end());
        }
        gpu::AcmSystem acm;
        acm.molecules      = am.data();
        acm.nmolecule      = am.size();
        acm.fragments      = af.data();
        acm.nfragment      = af.size();
        acm.acmParticles   = ap.data();
        acm.nacm           = ap.size();
        acm.x              = xa.data();
        acm.particleValues = pv.data();
        acm.bonds          = ab.data();
        acm.nbond          = ab.size()/2;
        acm.bondValues     = bv.data();
        acm.bondFactors    = bfac.data();
        acm.fixedValues    = fv.data();
        acm.nfixed         = fv.size();
        acm.nvalue         = nreplica*im.acmValues.size();
        acm.symmetry       = sym.data();
        acm.gaussian       = pd->chargeDistribution() == ChargeDistributionType::Gaussian;
        acm.coulomb        = ONE_4PI_EPS0;
        acm.electronVolt   = ELECTRONVOLT;
        gpu::createAcm(im.system, acm);
    }
    return im.system != nullptr;
}

size_t GpuTrainingSystem::nreplica() const
{
    return impl_->nreplica;
}

void GpuTrainingSystem::updateParameters(const ForceField *pd, const ForceComputer *fc)
{
    for (size_t r = 0; r < impl_->nreplica; r++)
    {
        stageParameters(r, pd, fc);
    }
    uploadParameters(pd);
}

//! The topology of a block, see GpuTrainingSystem::Impl::Block
static const Topology *blockTopology(const ACTMol *molptr, int what)
{
    auto &tops = molptr->fragmentHandler()->topologies();
    return (what < 2) ? &tops[what] : molptr->topology();
}

/*! \brief The rows of the parameter tables of one block, see GpuTrainingSystem::stageParameters
 *
 * Only the lists of which the parameters may have changed since the
 * versions in rowVersion are made, the others are left as they are in
 * tables. The parameters of the rows of a list follow from those of its
 * interaction type, which Topology::parameterVersion tells, except for the
 * polarization pairs, which depend on the charges of the shells as well.
 * \param[in]    molptr     The molecule
 * \param[in]    what       The block of the molecule, see Block
 * \param[in]    d          The description of the block
 * \param[in]    epsilonr   The relative dielectric constant
 * \param[inout] rowVersion The parameter version of each list when its rows were made
 * \param[inout] tables     The rows of the block, 8 per row
 */
static void stageRows(const ACTMol            *molptr,
                      int                      what,
                      const GpuSystemTopology &d,
                      double                   epsilonr,
                      std::vector<size_t>     *rowVersion,
                      float                   *tables)
{
    const Topology *top   = blockTopology(molptr, what);
    const auto     &atoms = top->atoms();
    rowVersion->resize(d.lists.size(), 0);
    size_t start = 0;
    for (size_t l = 0; l < d.lists.size(); l++)
    {
        const auto &rows    = d.rows[l];
        float      *p       = tables + 8*start;
        start              += rows.size();
        auto        version = top->parameterVersion(d.itypes[d.lists[l].energyIndex]);
        bool        charges = !d.lists[l].bonded && d.pairKernels[l] == gpu::PairKernel::Polarization;
        if (!charges && version == (*rowVersion)[l])
        {
            continue;
        }
        for (const auto &row : rows)
        {
            if (d.lists[l].bonded)
            {
                gpu::BondedEntry be;
                gpuBondedEntry(d.bondedKernels[l], *row.entry, &be);
                std::copy(be.p, be.p + 8, p);
            }
            else
            {
                gpu::PairEntry pe;
                bool           active = row.pair ?
                    gpuPairEntry(d.pairKernels[l], *row.pair, atoms, epsilonr, &pe) :
                    gpuPairEntry(d.pairKernels[l], *row.entry, atoms, epsilonr, &pe);
                if (!active)
                {
                    pe.p[7] = 1;
                }
                std::copy(pe.p, pe.p + 8, p);
            }
            p += 8;
        }
        (*rowVersion)[l] = version;
    }
}

//! \return the number of rows of the parameter tables of a block
static size_t tableRows(const GpuSystemTopology &d)
{
    size_t n = 0;
    for (const auto &rows : d.rows)
    {
        n += rows.size();
    }
    return n;
}

//! The charges and widths of the particles of a block, four each, see GpuTrainingSystem::stageParameters
static void stageCharges(const ACTMol *molptr, int what, double *charges)
{
    for (const auto &a : blockTopology(molptr, what)->atoms())
    {
        *charges++ = a.charge();
        *charges++ = a.chargeVS();
        *charges++ = a.zeta();
        *charges++ = a.zetaVS();
    }
}

//! The moving shells and the virtual sites of one block, see GpuTrainingSystem::stageParameters
static void stageShellsVsites(const GpuSystemTopology       &d,
                              const GpuRelaxMolecule        &sh,
                              gpu::MovingShell              *shells,
                              gpu::VsiteEntry               *vsites)
{
    // The force constants of the shells depend on their charges
    std::copy(sh.movingShells.begin(), sh.movingShells.end(), shells);
    for (size_t v = 0; v < d.vsiteEntries.size(); v++)
    {
        (void) gpuVsiteEntry(d.vsiteTypes[v], *d.vsiteEntries[v], &vsites[v]);
    }
}

void GpuTrainingSystem::stageParameters(size_t replica, const ForceField *pd, const ForceComputer *fc)
{
    // Only the parameter tables, the charges, the shells and the virtual
    // sites are sent, the parameters of the interactions follow on the GPU.
    // Each block has a segment of these, of which only those of the active
    // blocks are made again once all were made, see setActiveMolecules.
    auto  &im       = *impl_;
    double epsilonr = pd->epsilonr();
    auto  &st       = im.staged[replica];
    size_t b0       = replica*im.nblockReplica;
    // Whether the GPU makes the parameters from the recipes, see makeRowRecipes
    bool   gpuMakes = im.stagedAll[replica] && im.recipeRows == 1 && !im.hostRows;
    if (!im.stagedAll[replica])
    {
        st.tables.clear();
        st.charges.clear();
        st.shells.clear();
        st.vsites.clear();
        im.segments.resize(im.nblockReplica);
    }
    static const size_t blockT = Timings::instance().index("GPU staging of charges, shells and virtual sites");
    static const size_t acmXT  = Timings::instance().index("GPU staging of coordinates for SQE");
    {
        ScopedTimer blockTimer(blockT);
        for (size_t b = b0; b < b0 + im.nblockReplica; b++)
        {
            if (im.stagedAll[replica] && !im.active.empty() && !im.active[b])
            {
                continue;
            }
            auto &bl  = im.blocks[b];
            auto &seg = im.segments[b - b0];
            // The GPU makes the rows, the moving shells and the virtual
            // sites from the values, once it can, see makeRowRecipes
            bool onHost = !im.stagedAll[replica] || im.recipeRows != 1 || im.hostRows;
            if (!im.stagedAll[replica])
            {
                GpuRelaxMolecule sh;
                im.shellsOf(bl, &sh);
                seg.table   = st.tables.size();
                seg.charge  = st.charges.size();
                seg.shell   = st.shells.size();
                seg.vsite   = st.vsites.size();
                seg.ntable  = 8*tableRows(bl.desc);
                seg.ncharge = 4*blockTopology(im.mols[bl.mol], bl.what)->atoms().size();
                seg.nshell  = sh.movingShells.size();
                seg.nvsite  = bl.desc.vsiteEntries.size();
                st.tables.resize(seg.table + seg.ntable, 0.0f);
                st.charges.resize(seg.charge + seg.ncharge);
                st.shells.resize(seg.shell + seg.nshell);
                st.vsites.resize(seg.vsite + seg.nvsite);
                // All rows are made the first time
                bl.rowVersion.assign(bl.desc.lists.size(), 0);
                stageRows(im.mols[bl.mol], bl.what, bl.desc, epsilonr, &bl.rowVersion, st.tables.data() + seg.table);
                stageShellsVsites(bl.desc, sh, st.shells.data() + seg.shell, st.vsites.data() + seg.vsite);
            }
            else if (onHost)
            {
                // The rows in place, of the lists whose parameters changed
                stageRows(im.mols[bl.mol], bl.what, bl.desc, epsilonr, &bl.rowVersion, st.tables.data() + seg.table);
                GpuRelaxMolecule sh;
                im.shellsOf(bl, &sh);
                GMX_RELEASE_ASSERT(sh.movingShells.size() == seg.nshell && bl.desc.vsiteEntries.size() == seg.nvsite,
                                   "The parameters of a block changed in size");
                stageShellsVsites(bl.desc, sh, st.shells.data() + seg.shell, st.vsites.data() + seg.vsite);
            }
            GMX_RELEASE_ASSERT(4*blockTopology(im.mols[bl.mol], bl.what)->atoms().size() == seg.ncharge,
                               "The parameters of a block changed in size");
            stageCharges(im.mols[bl.mol], bl.what, st.charges.data() + seg.charge);
        }
    }
    im.stagedAll[replica] = 1;
    if (im.recipeRows == 1)
    {
        // The values of this replica, that the force field has now
        im.recipes.update();
        const auto &v = im.recipes.values();
        std::copy(v.begin(), v.end(), im.recipeValues.begin() + replica*v.size());
        im.recipeStaged[replica] = 1;
    }
    st.acmValues.clear();
    st.acmX.clear();
    if (!im.acmMols.empty())
    {
        // Only the values of the parameters of the charge generation,
        // and the structures if the virtual sites may have moved
        for (auto p : im.acmValues)
        {
            st.acmValues.push_back(p->value());
        }
        if (im.acmVsites && !(gpuMakes && im.acmVsitesOnGpu))
        {
            // Unless the GPU constructs them, see makeAcmVsites
            ScopedTimer acmXTimer(acmXT);
            st.acmX = im.acmCoordinates(fc);
        }
    }
}

void GpuTrainingSystem::setActiveMolecules(int replica, const std::vector<char> &active)
{
    auto &im = *impl_;
    if (im.active.size() != im.blocks.size())
    {
        im.active.assign(im.blocks.size(), 1);
    }
    for (size_t b = 0; b < im.blocks.size(); b++)
    {
        const auto &bl = im.blocks[b];
        if ((replica < 0 || bl.replica == static_cast<size_t>(replica)) && bl.mol < active.size())
        {
            im.active[b] = active[bl.mol];
        }
    }
}

std::vector<uint64_t> GpuTrainingSystem::moleculeSignatures(const ForceField *pd, const ForceComputer *fc)
{
    auto &im     = *impl_;
    auto  active = im.active;
    im.active.clear();
    // The rows on the processor, which the signatures are made of, also if
    // the GPU makes them from the recipes
    im.hostRows = true;
    stageParameters(0, pd, fc);
    im.hostRows = false;
    im.active = std::move(active);
    // FNV-1a over the bytes of the parameters of the blocks of each molecule
    auto hash = [](uint64_t h, const void *data, size_t n)
    {
        auto bytes = static_cast<const unsigned char *>(data);
        for (size_t i = 0; i < n; i++)
        {
            h = (h ^ bytes[i])*1099511628211ULL;
        }
        return h;
    };
    const auto           &st = im.staged[0];
    std::vector<uint64_t> sig(im.mols.size(), 14695981039346656037ULL);
    for (size_t b = 0; b < im.nblockReplica; b++)
    {
        const auto &seg = im.segments[b];
        auto       &h   = sig[im.blocks[b].mol];
        h = hash(h, st.tables.data() + seg.table, seg.ntable*sizeof(float));
        h = hash(h, st.charges.data() + seg.charge, seg.ncharge*sizeof(double));
        h = hash(h, st.shells.data() + seg.shell, seg.nshell*sizeof(gpu::MovingShell));
        h = hash(h, st.vsites.data() + seg.vsite, seg.nvsite*sizeof(gpu::VsiteEntry));
    }
    for (size_t k = 0; k < im.mols.size(); k++)
    {
        for (auto v : im.molValues[k])
        {
            if (v < st.acmValues.size())
            {
                sig[k] = hash(sig[k], &st.acmValues[v], sizeof(double));
            }
        }
    }
    return sig;
}

size_t GpuTrainingSystem::version() const
{
    return impl_->versionCount;
}

bool GpuTrainingSystem::recipeRows() const
{
    return 1 == impl_->recipeRows;
}

size_t GpuTrainingSystem::deviationsVersion() const
{
    return impl_->versionCount*1000003 + impl_->devVersionCount;
}

void GpuTrainingSystem::uploadParameters(const ForceField *pd)
{
    auto &im = *impl_;
    // Once all replicas were staged, whether the GPU can make the rows, the
    // moving shells and the virtual sites, see makeRowRecipes, which
    // ACT_GPU_RECIPE_ROWS=0 turns off
    if (0 == im.recipeRows &&
        std::all_of(im.stagedAll.begin(), im.stagedAll.end(), [](char c) { return c != 0; }))
    {
        const char *env = std::getenv("ACT_GPU_RECIPE_ROWS");
        bool        off = env && std::string(env) == "0";
        im.recipeRows   = (!off && im.makeRowRecipes(pd)) ? 1 : -1;
    }
    bool gpuMakes = (1 == im.recipeRows &&
                     std::all_of(im.recipeStaged.begin(), im.recipeStaged.end(), [](char c) { return c != 0; }));
    // In the order of the blocks, those of a single replica as they are
    Impl::Staged all;
    for (size_t r = 0; r < im.nreplica && im.nreplica > 1; r++)
    {
        const auto &st = im.staged[r];
        if (!gpuMakes)
        {
            all.tables.insert(all.tables.end(), st.tables.begin(), st.tables.end());
            all.shells.insert(all.shells.end(), st.shells.begin(), st.shells.end());
            all.vsites.insert(all.vsites.end(), st.vsites.begin(), st.vsites.end());
        }
        all.charges.insert(all.charges.end(), st.charges.begin(), st.charges.end());
        all.acmValues.insert(all.acmValues.end(), st.acmValues.begin(), st.acmValues.end());
        all.acmX.insert(all.acmX.end(), st.acmX.begin(), st.acmX.end());
    }
    auto &up      = im.nreplica > 1 ? all : im.staged[0];
    auto &charges = up.charges;
    auto &values  = up.acmValues;
    auto &xa      = up.acmX;
    if (!im.acmMols.empty())
    {
        // The charges of the molecules are those generated last, see
        // copyChargesToMolecules, unless the parameters changed. With
        // several replicas, the charges of the molecules are those of one
        // replica, so the charges of all are generated each time.
        bool every = im.nreplica > 1;
        if (every || values != im.acmLastValues)
        {
            gpu::updateAcm(im.system, values.data());
            im.acmLastValues = values;
            im.acmGenerated  = true;
        }
        if (gpuMakes && im.acmVsitesOnGpu)
        {
            // The GPU constructs the virtual sites again if their values changed
            std::vector<double> v;
            for (auto k : im.acmVsiteSlots)
            {
                v.push_back(im.recipeValues[k]);
            }
            if (v != im.acmVsiteLast)
            {
                gpu::acmVsitesChanged(im.system);
                im.acmVsiteLast = std::move(v);
                im.acmGenerated = true;
            }
        }
        else if (im.acmVsites && (every || xa != im.acmLastX))
        {
            gpu::updateAcmCoordinates(im.system, xa.data());
            im.acmLastX     = xa;
            im.acmGenerated = true;
        }
    }
    if (gpuMakes)
    {
        // The GPU makes the rows, the moving shells and the virtual sites from the values
        gpu::updateRecipeValues(im.system, im.recipeValues.data());
        gpu::updateParameters(im.system, nullptr, charges.data(), ONE_4PI_EPS0, pd->epsilonr());
    }
    else
    {
        gpu::updateParameters(im.system, up.tables.data(), charges.data(), ONE_4PI_EPS0, pd->epsilonr());
        gpu::updateMovingShells(im.system, up.shells.data(), up.shells.size());
        gpu::updateVsites(im.system, up.vsites.data(), up.vsites.size());
    }
}

void GpuTrainingSystem::selectReplica(size_t replica)
{
    auto &im = *impl_;
    if (replica != im.resident)
    {
        for (size_t k = 0; k < im.mols.size(); k++)
        {
            for (int what = 0; what < 3; what++)
            {
                auto &e = im.mols[k]->conformationData(what)->energies;
                std::swap(e, im.replicaEnergies[im.resident][3*k + what]);
                std::swap(e, im.replicaEnergies[replica][3*k + what]);
            }
        }
        im.resident = replica;
    }
    // Only charges that were copied after they were generated, see copyChargesToMolecules
    if (!im.acmGenerated)
    {
        im.applyCharges(replica);
    }
}

void GpuTrainingSystem::copyChargesToMolecules()
{
    auto &im = *impl_;
    if (im.acmMols.empty())
    {
        return;
    }
    if (im.acmGenerated)
    {
        im.acmGenerated = false;
        im.acmCharges.resize(im.nreplica*im.nacmAtomReplica);
        gpu::downloadAcmCharges(im.system, im.acmCharges.data());
        im.chargesOf    = std::numeric_limits<size_t>::max();
    }
    im.applyCharges(im.resident);
}

std::pair<double, double> GpuTrainingSystem::chargeDifference() const
{
    auto  &im = *impl_;
    size_t n  = 0;
    for (const auto &bl : im.blocks)
    {
        n += im.mols[bl.mol]->conformationData(bl.what)->arrays.natom;
    }
    std::vector<double> q(4*n);
    gpu::downloadCharges(im.system, q.data());
    std::pair<double, double> worst = { 0, 0 };
    size_t                    off   = 0;
    for (const auto &bl : im.blocks)
    {
        auto            molptr = im.mols[bl.mol];
        auto           &tops   = molptr->fragmentHandler()->topologies();
        const Topology *top    = (bl.what < 2) ? &tops[bl.what] : molptr->topology();
        const auto     &atoms  = top->atoms();
        // The molecules have the charges of one replica
        for (size_t i = 0; i < atoms.size() && bl.replica == im.resident; i++)
        {
            double d = std::abs(q[4*(off+i)] - atoms[i].charge());
            if (d > worst.first)
            {
                worst = { d, atoms[i].charge() };
            }
        }
        off += molptr->conformationData(bl.what)->arrays.natom;
    }
    return worst;
}

void GpuTrainingSystem::resetShells()
{
    auto u = impl_->use(-1);
    gpu::resetShells(impl_->system, u.data());
}

void GpuTrainingSystem::keepShells()
{
    auto &im = *impl_;
    im.shellsKept = !im.notFinite;
    im.notFinite  = false;
    if (im.shellsKept)
    {
        auto u = im.use(-1);
        gpu::keepShells(im.system, u.data(), false);
    }
}

bool GpuTrainingSystem::restoreShells(int what)
{
    auto &im = *impl_;
    if (!im.shellsKept)
    {
        return false;
    }
    auto u = im.use(what);
    gpu::keepShells(im.system, u.data(), true);
    return true;
}

bool GpuTrainingSystem::setDeviations(const GpuDeviationSetup &setup)
{
    auto &im = *impl_;
    // The references are recomputed in place for another temperature or
    // other conformations, see interactionReference in devcomputer.cpp
    std::vector<std::pair<double, size_t>> refState;
    for (const auto *refs : { &setup.interaction, &setup.epotReference })
    {
        for (const auto *ref : *refs)
        {
            refState.push_back(ref ? std::make_pair(ref->beta, ref->conformationVersion) : std::make_pair(0.0, size_t(0)));
        }
    }
    if (im.devReady && refState == im.devRefState &&
        setup.targets == im.devSetup.targets && setup.mols == im.devSetup.mols &&
        setup.interaction == im.devSetup.interaction && setup.interactionConfs == im.devSetup.interactionConfs &&
        setup.epot == im.devSetup.epot && setup.epotReference == im.devSetup.epotReference &&
        setup.epotConfs == im.devSetup.epotConfs && setup.lossFunction == im.devSetup.lossFunction)
    {
        return true;
    }
    im.devReady    = false;
    im.devComputed = false;
    if (!im.system || setup.mols != std::vector<const ACTMol *>(im.mols.begin(), im.mols.end()) ||
        setup.targets.size() > 8 || gpu::pairEnergyTerms != 3)
    {
        return false;
    }
    // The energy types of each block, as typeEnergies makes them
    std::vector<gpu::DeviationType> types;
    std::vector<uint32_t>           firstType;
    std::vector<uint32_t>           slotOfConf;
    std::vector<uint64_t>           firstConf;
    for (const auto &bl : im.blocks)
    {
        std::map<InteractionType, gpu::DeviationType> bt;
        auto &epot  = bt[InteractionType::EPOT];
        epot.kind   = gpu::deviationEpot;
        for (size_t e = 0; e < bl.desc.potentials.size(); e++)
        {
            const auto &itypes = energyTypes(bl.desc.potentials[e], bl.desc.itypes[e]);
            for (size_t t = 0; t < itypes.size(); t++)
            {
                auto &dt  = bt[itypes[t]];
                dt.kind   = gpu::deviationTerm;
                dt.source = e*gpu::pairEnergyTerms + t;
            }
        }
        auto vdw = bt.find(InteractionType::VDWCORRECTION);
        if (vdw != bt.end())
        {
            auto ex    = bt.find(InteractionType::EXCHANGE);
            int  exsrc = ex != bt.end() ? ex->second.source : -1;
            auto &exch   = bt[InteractionType::EXCHANGE];
            exch.kind    = gpu::deviationExchange;
            exch.source  = exsrc;
            exch.source2 = vdw->second.source;
            vdw->second.kind = gpu::deviationZero;
        }
        if (bt.size() > 32)
        {
            return false;
        }
        firstType.push_back(types.size());
        for (auto &t : bt)
        {
            t.second.type = static_cast<int32_t>(t.first);
            switch (t.first)
            {
            case InteractionType::ELECTROSTATICS:
            case InteractionType::POLARIZATION:
            case InteractionType::CHARGETRANSFER:
                t.second.role = gpu::deviationElec;
                break;
            case InteractionType::EXCHANGE:
                t.second.role = gpu::deviationExchangeRole;
                break;
            case InteractionType::INDUCTION:
                t.second.role = gpu::deviationInduction;
                break;
            case InteractionType::INDUCTIONCORRECTION:
                t.second.role = gpu::deviationInductionCorrection;
                break;
            default:
                t.second.role = 0;
            }
            types.push_back(t.second);
        }
        firstConf.push_back(slotOfConf.size());
        for (auto s : bl.slotOfConf)
        {
            slotOfConf.push_back(s);
        }
    }
    firstType.push_back(types.size());
    firstConf.push_back(slotOfConf.size());
    // The reference data of each molecule, the same for all replicas
    size_t ntarget = setup.targets.size();
    std::vector<double> weights, eqm, epotWeights, epotEqm;
    std::vector<char>   have, epotHave;
    std::vector<uint64_t> refOffset(im.mols.size(), 0), epotOffset(im.mols.size(), 0);
    for (size_t k = 0; k < im.mols.size(); k++)
    {
        refOffset[k] = weights.size();
        if (setup.interaction[k])
        {
            const auto &ref   = *setup.interaction[k];
            size_t      nconf = setup.interactionConfs[k];
            if (ref.weight.size() < nconf || ref.energy.size() != ntarget)
            {
                return false;
            }
            weights.insert(weights.end(), ref.weight.begin(), ref.weight.begin() + nconf);
            eqm.resize(weights.size()*ntarget);
            have.resize(weights.size()*ntarget);
            for (size_t j = 0; j < ntarget; j++)
            {
                for (size_t c = 0; c < nconf; c++)
                {
                    eqm[refOffset[k]*ntarget + j*nconf + c]  = ref.energy[j][c];
                    have[refOffset[k]*ntarget + j*nconf + c] = ref.have[j][c];
                }
            }
        }
        epotOffset[k] = epotWeights.size();
        if (setup.epot && setup.epotReference[k])
        {
            const auto &ref   = *setup.epotReference[k];
            size_t      nconf = setup.epotConfs[k];
            if (ref.weight.size() < nconf || ref.energy.size() != 1)
            {
                return false;
            }
            epotWeights.insert(epotWeights.end(), ref.weight.begin(), ref.weight.begin() + nconf);
            epotEqm.insert(epotEqm.end(), ref.energy[0].begin(), ref.energy[0].begin() + nconf);
            epotHave.insert(epotHave.end(), ref.have[0].begin(), ref.have[0].begin() + nconf);
        }
    }
    // The dimers and the topologies of the potential energy of each replica, in the order of the molecules
    std::vector<gpu::DeviationDimer>   dimers;
    std::vector<gpu::DeviationMonomer> monomers;
    std::vector<int32_t>               matches;
    for (size_t r = 0; r < im.nreplica; r++)
    {
        for (size_t k = 0; k < im.mols.size(); k++)
        {
            int db = im.find(k, 2, r);
            if (setup.interaction[k] && db >= 0)
            {
                gpu::DeviationDimer dm;
                dm.dimer       = db;
                dm.replica     = r;
                dm.nconf       = setup.interactionConfs[k];
                dm.refOffset   = refOffset[k];
                dm.confOffset  = dimers.empty() ? 0 : dimers.back().confOffset + dimers.back().nconf;
                dm.matchOffset = matches.size();
                for (int j = 0; j < 2; j++)
                {
                    int mb        = im.find(k, j, r);
                    dm.monomer[j] = mb >= 0 ? static_cast<uint32_t>(mb) : gpu::noBlock;
                }
                for (uint32_t i = firstType[db]; i < firstType[db + 1]; i++)
                {
                    for (int j = 0; j < 2; j++)
                    {
                        int32_t m = -1;
                        if (dm.monomer[j] != gpu::noBlock)
                        {
                            for (uint32_t i2 = firstType[dm.monomer[j]]; i2 < firstType[dm.monomer[j] + 1]; i2++)
                            {
                                if (types[i2].type == types[i].type)
                                {
                                    m = i2 - firstType[dm.monomer[j]];
                                }
                            }
                        }
                        matches.push_back(m);
                    }
                }
                dimers.push_back(dm);
            }
            int mb = im.find(k, 0, r);
            if (setup.epot && setup.epotReference[k] && mb >= 0)
            {
                gpu::DeviationMonomer mo;
                mo.block     = mb;
                mo.replica   = r;
                mo.nconf     = setup.epotConfs[k];
                mo.refOffset  = epotOffset[k];
                mo.confOffset = monomers.empty() ? 0 : monomers.back().confOffset + monomers.back().nconf;
                monomers.push_back(mo);
            }
        }
    }
    std::vector<int32_t> targets;
    for (auto t : setup.targets)
    {
        targets.push_back(static_cast<int32_t>(t));
    }
    gpu::DeviationSystem d;
    d.types          = types.data();
    d.blockFirstType = firstType.data();
    d.slotOfConf     = slotOfConf.data();
    d.blockFirstConf = firstConf.data();
    d.dimers         = dimers.data();
    d.ndimer         = dimers.size();
    d.matches        = matches.data();
    d.nmatch         = matches.size();
    d.targets        = targets.data();
    d.ntarget        = targets.size();
    d.allElecType    = static_cast<int32_t>(InteractionType::ALLELEC);
    d.exchIndType    = static_cast<int32_t>(InteractionType::EXCHIND);
    d.inductionType  = static_cast<int32_t>(InteractionType::INDUCTION);
    d.weights        = weights.data();
    d.nweight        = weights.size();
    d.eqm            = eqm.data();
    d.have           = have.data();
    d.monomers       = monomers.data();
    d.nmonomer       = monomers.size();
    d.epotWeights    = epotWeights.data();
    d.epotEqm        = epotEqm.data();
    d.epotHave       = epotHave.data();
    d.nepot          = epotWeights.size();
    d.nreplica       = im.nreplica;
    d.lossFunction   = setup.lossFunction;
    gpu::createDeviations(im.system, d);
    im.devSetup    = setup;
    im.devRefState = std::move(refState);
    // Where the results of each replica and molecule are
    im.devDimerOf.assign(im.nreplica*im.mols.size(), -1);
    im.devMonomerOf.assign(im.nreplica*im.mols.size(), -1);
    for (size_t i = 0; i < dimers.size(); i++)
    {
        im.devDimerOf[dimers[i].replica*im.mols.size() + im.blocks[dimers[i].dimer].mol] = i;
    }
    for (size_t i = 0; i < monomers.size(); i++)
    {
        im.devMonomerOf[monomers[i].replica*im.mols.size() + im.blocks[monomers[i].block].mol] = i;
    }
    im.devDimers   = std::move(dimers);
    im.devMonomers = std::move(monomers);
    im.devVersionCount += 1;
    im.devTargets = ntarget;
    im.devReady   = true;
    return true;
}

bool GpuTrainingSystem::deviationsReady() const
{
    return impl_->devReady;
}

void GpuTrainingSystem::setDeviationsOnly(bool only)
{
    impl_->devOnly = only && impl_->devReady;
}

void GpuTrainingSystem::computeDeviations()
{
    auto &im = *impl_;
    if (!im.devReady)
    {
        return;
    }
    size_t nmc = im.devMonomers.empty() ? 0 : im.devMonomers.back().confOffset + im.devMonomers.back().nconf;
    im.devDimerSums.resize(2*im.devDimers.size()*im.devTargets);
    im.devEpotTerms.resize(2*nmc);
    auto u = im.use(-1);
    gpu::computeDeviations(im.system, u.data(), im.devDimerSums.data(), im.devEpotTerms.data());
    im.devComputed = true;
}

bool GpuTrainingSystem::moleculeDeviations(size_t mol, GpuMoleculeDeviations *dev) const
{
    const auto &im = *impl_;
    if (!im.devReady || !im.devComputed || mol >= im.mols.size())
    {
        return false;
    }
    size_t rk = im.resident*im.mols.size() + mol;
    int    d  = im.devDimerOf[rk];
    dev->dimer = d >= 0;
    dev->targets.clear();
    if (d >= 0)
    {
        dev->targets.assign(im.devDimerSums.begin() + 2*d*im.devTargets,
                            im.devDimerSums.begin() + 2*(d + 1)*im.devTargets);
    }
    int m    = im.devMonomerOf[rk];
    dev->epot = m >= 0;
    dev->epotTerms.clear();
    if (m >= 0)
    {
        const auto &mo = im.devMonomers[m];
        dev->epotTerms.assign(im.devEpotTerms.begin() + 2*mo.confOffset,
                              im.devEpotTerms.begin() + 2*(mo.confOffset + mo.nconf));
    }
    return true;
}

bool GpuTrainingSystem::deviationSums(std::vector<std::pair<double, double>> *sums) const
{
    const auto &im = *impl_;
    sums->clear();
    if (!im.devReady || !im.devComputed)
    {
        return false;
    }
    // In the order of the molecules, as ForceEnergyDevComputer::calcDeviationGpu adds them
    std::vector<GpuMoleculeDeviations> devs(im.mols.size());
    for (size_t k = 0; k < im.mols.size(); k++)
    {
        (void) moleculeDeviations(k, &devs[k]);
    }
    for (size_t j = 0; j <= im.devTargets; j++)
    {
        double sw = 0, sl = 0;
        for (const auto &dev : devs)
        {
            if (j < im.devTargets && dev.dimer)
            {
                sw += dev.targets[2*j];
                sl += dev.targets[2*j + 1];
            }
            else if (j == im.devTargets && dev.epot)
            {
                for (size_t c = 0; c < dev.epotTerms.size(); c += 2)
                {
                    sw += dev.epotTerms[c];
                    sl += dev.epotTerms[c + 1];
                }
            }
        }
        sums->push_back({ sw, sl });
    }
    return true;
}

void GpuTrainingSystem::setKeepEnergies(bool keep)
{
    impl_->keepEnergies = keep;
}

bool GpuTrainingSystem::typeEnergies(const ACTMol *mol, int what, std::vector<GpuTypeEnergies> *energies) const
{
    const auto &im = *impl_;
    energies->clear();
    auto k = std::find(im.mols.begin(), im.mols.end(), mol);
    if (k == im.mols.end())
    {
        return false;
    }
    int b = im.find(k - im.mols.begin(), what, im.resident);
    if (b < 0 || im.blocks[b].rawAfter.empty() || im.blocks[b].rawBefore.empty())
    {
        return false;
    }
    const auto &bl    = im.blocks[b];
    size_t      nconf = bl.confOfSlot.size();
    size_t      st    = (*k)->conformationData(what)->arrays.stride;
    // As storeBlockEnergies: the energy of a type is that of the last
    // potential with it, the potential energy the sum of all
    std::map<InteractionType, GpuTypeEnergies> types;
    auto &epot = types[InteractionType::EPOT];
    epot.type  = InteractionType::EPOT;
    epot.act.assign(nconf, 0.0);
    epot.before.assign(nconf, 0.0);
    for (size_t e = 0; e < bl.desc.potentials.size(); e++)
    {
        const auto &itypes = energyTypes(bl.desc.potentials[e], bl.desc.itypes[e]);
        size_t      off    = e*gpu::pairEnergyTerms*st;
        for (size_t t = 0; t < itypes.size(); t++)
        {
            auto &te = types[itypes[t]];
            te.type  = itypes[t];
            te.act.resize(nconf);
            te.before.resize(nconf);
            for (size_t c = 0; c < nconf; c++)
            {
                size_t s     = off + t*st + bl.slotOfConf[c];
                te.act[c]    = 0.0 + bl.rawAfter[s];
                te.before[c] = 0.0 + bl.rawBefore[s];
            }
        }
        for (size_t c = 0; c < nconf; c++)
        {
            size_t s      = off + bl.slotOfConf[c];
            epot.act[c]    = epot.act[c] + (bl.rawAfter[s] + bl.rawAfter[s + st] + bl.rawAfter[s + 2*st]);
            epot.before[c] = epot.before[c] + (bl.rawBefore[s] + bl.rawBefore[s + st] + bl.rawBefore[s + 2*st]);
        }
    }
    // The van der Waals correction is part of the exchange
    auto vdw = types.find(InteractionType::VDWCORRECTION);
    if (vdw != types.end())
    {
        auto &exch = types[InteractionType::EXCHANGE];
        exch.type  = InteractionType::EXCHANGE;
        exch.act.resize(nconf, 0.0);
        exch.before.resize(nconf, 0.0);
        for (size_t c = 0; c < nconf; c++)
        {
            exch.act[c]          = exch.act[c] + vdw->second.act[c];
            exch.before[c]       = exch.before[c] + vdw->second.before[c];
            vdw->second.act[c]    = 0;
            vdw->second.before[c] = 0;
        }
    }
    for (auto &t : types)
    {
        energies->push_back(std::move(t.second));
    }
    return true;
}

void GpuTrainingSystem::constructVsites()
{
    gpu::constructVsites(impl_->system);
}

void GpuTrainingSystem::computeForces(int what)
{
    auto &im = *impl_;
    auto  u  = im.use(what);
    gpu::computeBlocks(im.system, u.data(), nullptr, true);
}

bool GpuTrainingSystem::haveKeptShells() const
{
    return impl_->shellsKept;
}

void GpuTrainingSystem::compute(int what, bool forces, bool backup)
{
    auto &im = *impl_;
    auto  u  = im.use(what);
    if (im.devOnly)
    {
        // The energies stay on the GPU for computeDeviations
        gpu::computeBlocksKeep(im.system, u.data(), forces, backup);
        return;
    }
    std::vector<const double *> eptr(im.blocks.size(), nullptr);
    gpu::computeBlocksInPlace(im.system, u.data(), eptr.data(), forces);
    if (im.devReady)
    {
        gpu::keepEnergies(im.system, u.data(), backup);
    }
    for (size_t b = 0; b < im.blocks.size(); b++)
    {
        if (u[b] && eptr[b])
        {
            auto &bl = im.blocks[b];
            im.storeBlockEnergies(bl, eptr[b], backup);
            if (im.keepEnergies)
            {
                size_t n   = bl.desc.potentials.size()*gpu::pairEnergyTerms*
                    im.mols[bl.mol]->conformationData(bl.what)->arrays.stride;
                auto  &raw = backup ? bl.rawBefore : bl.rawAfter;
                raw.assign(eptr[b], eptr[b] + n);
            }
        }
    }
    if (backup)
    {
        // The energies without a block are backed up as they are
        for (size_t r = 0; r < im.nreplica; r++)
        {
            for (size_t k = 0; k < im.mols.size(); k++)
            {
                for (int w = 0; w < 3; w++)
                {
                    int b = im.find(k, w, r);
                    if (b >= 0 && u[b] && eptr[b])
                    {
                        continue;
                    }
                    auto &cener = im.energiesOf(k, w, r);
                    cener.setBackups(true);
                    for (auto t : cener.types())
                    {
                        auto et = cener.find(t);
                        std::copy(et.act, et.act + cener.nconf(), et.before);
                    }
                }
            }
        }
    }
}

size_t GpuTrainingSystem::relax(int what, const bool relaxMonomer[2], const gpu::RelaxSettings &settings,
                                size_t *nFailed)
{
    auto &im = *impl_;
    auto  u  = im.use(what == 2 ? 2 : -2);
    std::vector<char> relax(2*im.blocks.size(), 0);
    for (size_t b = 0; b < im.blocks.size(); b++)
    {
        if (im.blocks[b].what < 2)
        {
            relax[2*b] = 1;
        }
        else
        {
            relax[2*b]   = relaxMonomer[0] ? 1 : 0;
            relax[2*b+1] = relaxMonomer[1] ? 1 : 0;
        }
    }
    size_t nComputed = 0;
    bool   notFinite = false;
    gpu::relaxBlocks(im.system, u.data(), relax.data(), settings, &nComputed, &notFinite, nFailed);
    if (notFinite)
    {
        im.notFinite = true;
        // As on the processor, see moveCoordinate
        moveCoordinate(0, std::numeric_limits<double>::quiet_NaN());
    }
    return nComputed;
}

void GpuTrainingSystem::copyRelaxedMonomerShells()
{
    auto &im = *impl_;
    for (size_t b = 0; b < im.blocks.size(); b++)
    {
        const auto &bl = im.blocks[b];
        if (bl.what != 2)
        {
            continue;
        }
        auto  fh     = im.mols[bl.mol]->fragmentHandler();
        auto &tops   = fh->topologies();
        auto &astart = fh->atomStart();
        for (int f = 0; f < 2; f++)
        {
            int mb = im.find(bl.mol, f, bl.replica);
            if (mb < 0)
            {
                continue;
            }
            std::vector<uint32_t> pairs;
            const auto &atoms = tops[f].atoms();
            for (size_t j = 0; j < atoms.size(); j++)
            {
                if (atoms[j].pType() == ActParticle::Shell)
                {
                    pairs.push_back(j);
                    pairs.push_back(astart[f] + j);
                }
            }
            // The same conformation in the slots of the monomer and the dimer
            std::vector<uint32_t> slots(bl.confOfSlot.size());
            for (size_t s = 0; s < slots.size(); s++)
            {
                slots[s] = im.blocks[mb].slotOfConf[bl.confOfSlot[s]];
            }
            gpu::copyCoordinates(im.system, mb, b, pairs.data(), pairs.size()/2, slots.data());
        }
    }
}

void GpuTrainingSystem::saveDimers()
{
    for (size_t b = 0; b < impl_->blocks.size(); b++)
    {
        if (impl_->blocks[b].what == 2)
        {
            gpu::saveCoordinates(impl_->system, b);
        }
    }
}

void GpuTrainingSystem::restoreDimers()
{
    for (size_t b = 0; b < impl_->blocks.size(); b++)
    {
        if (impl_->blocks[b].what == 2)
        {
            gpu::restoreCoordinates(impl_->system, b);
        }
    }
}

#endif

void computePotential(MsgHandler                 *msghandler,
                      Potential                   mypot,
                      const TopologyEntryVector  &topentries,
                      const std::vector<ActAtom> &atoms,
                      ConformationData           *cdata,
                      double                      epsilonr,
                      bool                        conformations,
                      bool                        forces)
{
    // Only the conformations are kept, the data for computing them is new
    ConformationData fresh;
    fresh.confs  = std::move(cdata->confs);
    fresh.arrays = std::move(cdata->arrays);
    *cdata       = std::move(fresh);
    auto *coords = &cdata->confs;
    for(auto &ce : *coords)
    {
        ce.done = false;
    }
    auto &ca = cdata->arrays;
    ca.clearForces();
    ca.pad();
    cdata->energies.setNumberOfConformations(coords->size());
    if (!(conformations && coords->size() > 1 &&
          calcConformations(mypot, topentries, atoms, &ca, &cdata->energies, epsilonr, forces,
                            nullptr, nullptr)))
    {
        calcSome(msghandler, mypot, topentries, atoms, cdata, &cdata->energies, epsilonr, ca.nconf);
    }
}

} // namespace

/*
 * This file is part of the GROMACS molecular simulation package.
 *
 * Copyright (c) 1991-2000, University of Groningen, The Netherlands.
 * Copyright (c) 2001-2004, The GROMACS development team.
 * Copyright (c) 2013,2014,2015,2016,2017,2018, by the GROMACS development team, led by
 * Mark Abraham, David van der Spoel, Berk Hess, and Erik Lindahl,
 * and including many others, as listed in the AUTHORS file in the
 * top-level source directory and at http://www.gromacs.org.
 * Copyright (C) 2026, the Alexandria Chemistry Toolkit developers.
 *
 * GROMACS is free software; you can redistribute it and/or
 * modify it under the terms of the GNU Lesser General Public License
 * as published by the Free Software Foundation; either version 2.1
 * of the License, or (at your option) any later version.
 *
 * GROMACS is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
 * Lesser General Public License for more details.
 *
 * You should have received a copy of the GNU Lesser General Public
 * License along with GROMACS; if not, see
 * http://www.gnu.org/licenses, or write to the Free Software Foundation,
 * Inc., 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301  USA.
 *
 * If you want to redistribute modifications to GROMACS, please
 * consider that scientific software is very special. Version
 * control is crucial - bugs must be traceable. We will be happy to
 * consider code for inclusion in the official distribution, but
 * derived work must not be called official GROMACS. Details are found
 * in the README & COPYING files - if they are missing, get the
 * official version at http://www.gromacs.org.
 *
 * To help us fund GROMACS development, we humbly ask that you cite
 * the research papers on the package. Check out http://www.gromacs.org.
 */

/*
 * Adapted for ACT by DvdS 2023-05-17
 */
#include "actpre.h"

#include "vsitehandler.h"

#include <cstdio>

#include <algorithm>
#include <array>
#include <set>
#include <vector>

#include "act/forcefield/forcefield_parametername.h"
#include "act/topology/topology.h"
#include "act/utility/timings.h"
#include "gromacs/math/functions.h"
#include "gromacs/math/vec.h"
#include "gromacs/topology/topology.h"
#include "gromacs/utility/exceptions.h"
#include "gromacs/utility/fatalerror.h"

using gmx::RVec;

/* Vsite construction routines */

static void constr_vsite1(const rvec xi, rvec x)
{
    copy_rvec(xi, x);
    /* TOTAL: 0 flops */
}

static void constr_vsite2(const rvec xi, const rvec xj, rvec x, real a)
{
    rvec dx;
    rvec_sub(xj, xi, dx);
    x[XX] = xi[XX] + a*dx[XX];
    x[YY] = xi[YY] + a*dx[YY];
    x[ZZ] = xi[ZZ] + a*dx[ZZ];
    /* TOTAL: 9 flops */
}

static void constr_vsite2fd(const rvec xi, const rvec xj, rvec x, real a)
{
    rvec dx;

    rvec_sub(xj, xi, dx);

    real factor = a*gmx::invsqrt(iprod(dx, dx));
    x[XX] = xi[XX] + factor*dx[XX];
    x[YY] = xi[YY] + factor*dx[YY];
    x[ZZ] = xi[ZZ] + factor*dx[ZZ];
    // 20 flops at least
}

static gmx_unused void constr_vsite3(const rvec xi, const rvec xj, const rvec xk, rvec x, real a, real b)
{
    real c = 1 - a - b;
    /* 2 flops */

    x[XX] = c*xj[XX] + a*xi[XX] + b*xk[XX];
    x[YY] = c*xj[YY] + a*xi[YY] + b*xk[YY];
    x[ZZ] = c*xj[ZZ] + a*xi[ZZ] + b*xk[ZZ];
    /* 15 Flops */

    /* TOTAL: 17 flops */
}

static gmx_unused void constr_vsite4(const rvec xi, const rvec xj, const rvec xk, const rvec xl,
                                     rvec x, real a, real b, real c)
{
    real d = 1 - a - b - c;
    /* 3 flops */

    for (int m = 0; m < DIM; m++)
    {
        x[m] = d*xj[m] + a*xi[m] + b*xk[m] + c*xl[m];
    }
    /* 21 Flops */

    /* TOTAL: 24 flops */
}

static void constr_vsite3FD(const rvec xi, const rvec xj, const rvec xk, rvec x, real a, real b)
{
    rvec xij, xjk, temp;
    real c;

    rvec_sub(xi, xj, xij);
    rvec_sub(xk, xj, xjk);
    /* 6 flops */

    /* temp goes from i to a point on the line jk */
    temp[XX] = xij[XX] + a*xjk[XX];
    temp[YY] = xij[YY] + a*xjk[YY];
    temp[ZZ] = xij[ZZ] + a*xjk[ZZ];
    /* 6 flops */

    c = b*gmx::invsqrt(iprod(temp, temp));
    /* 6 + 10 flops */

    x[XX] = xj[XX] + c*temp[XX];
    x[YY] = xj[YY] + c*temp[YY];
    x[ZZ] = xj[ZZ] + c*temp[ZZ];
    /* 6 Flops */

    /* TOTAL: 34 flops */
}

static void constr_vsite3OUT(const rvec xi, const rvec xj, const rvec xk, rvec x,
                             real a, real b, real c)
{
    rvec xij, xik, temp;

    /* For water-like compounds this assumes O H H as the order of atoms */
    rvec_sub(xj, xi, xij);
    rvec_sub(xk, xi, xik);
    cprod(xij, xik, temp);
    /* 15 Flops */

    x[XX] = xi[XX] + a*xij[XX] + b*xik[XX] + c*temp[XX];
    x[YY] = xi[YY] + a*xij[YY] + b*xik[YY] + c*temp[YY];
    x[ZZ] = xi[ZZ] + a*xij[ZZ] + b*xik[ZZ] + c*temp[ZZ];
    /* 18 Flops */

    /* TOTAL: 33 flops */
}

template <class Indices>
static void spread_vsite1(const Indices &ia, rvec f[])
{
    int ai = ia[0];
    int av = ia[1];
    rvec_inc(f[ai], f[av]);
    clear_rvec(f[av]);
}

template <class Indices>
static void spread_vsite2(const Indices &ia, real a,
                          rvec f[])
{
    rvec fi, fj;
    int  av, ai, aj;

    av = ia[2];
    ai = ia[0];
    aj = ia[1];

    svmul(1 - a, f[av], fi);
    svmul(    a, f[av], fj);
    /* 7 flop */

    rvec_inc(f[ai], fi);
    rvec_inc(f[aj], fj);
    /* 6 Flops */
    clear_rvec(f[av]);

    /* TOTAL: 13 flops */
}

static void spread_vsite2fd(const rvec xi,
                            const rvec xj,
                            real       a,
                            rvec       fi,
                            rvec       fj,
                            rvec       fav)
{
    rvec    dx;

    rvec_sub(xj, xi, dx);
    if (norm(dx) == 0)
    {
        GMX_THROW(gmx::InternalError(gmx::formatString("Distance between atoms is zero").c_str()));
    }
    // The virtual site is at xi + a*u, with u the unit vector from i to j,
    // see constr_vsite2fd. The derivative of u with respect to xj is
    // (1 - u u^T)/|dx|, which gives the force on j, the remainder goes to i.
    real invd  = gmx::invsqrt(iprod(dx, dx));
    real gamma = a*invd;
    rvec u;
    svmul(invd, dx, u);
    real fu = iprod(fav, u);
    rvec df;
    for(int m = 0; m < DIM; m++)
    {
        df[m] = gamma*(fav[m] - fu*u[m]);
    }
    for(int m = 0; m < DIM; m++)
    {
        fi[m] += fav[m] - df[m];
    }
    rvec_inc(fj, df);
    
    /* 6 Flops */
    clear_rvec(fav);
}

static void spread_vsite3(rvec ffi,
                          rvec ffj,
                          rvec ffk,
                          rvec fav,
                          real a, real b)
{
    rvec fi, fj, fk;
    svmul(1 - a - b, fav, fj);
    svmul(        a, fav, fi);
    svmul(        b, fav, fk);
    /* 11 flops */

    rvec_inc(ffi, fi);
    rvec_inc(ffj, fj);
    rvec_inc(ffk, fk);
    /* 9 Flops */

    clear_rvec(fav);
    /* TOTAL: 20 flops */
}

static void spread_vsite4(rvec ffi,
                          rvec ffj,
                          rvec ffk,
                          rvec ffl,
                          rvec fav,
                          real a, real b, real c)
{
    rvec fi, fj, fk, fl;
    real d = 1 - a - b - c;
    /* 3 flops */
    // Same weights as in constr_vsite4
    svmul(a, fav, fi);
    svmul(d, fav, fj);
    svmul(b, fav, fk);
    svmul(c, fav, fl);
    /* 12 flops */

    rvec_inc(ffi, fi);
    rvec_inc(ffj, fj);
    rvec_inc(ffk, fk);
    rvec_inc(ffl, fl);
    /* 12 Flops */

    clear_rvec(fav);
    /* TOTAL: 30 flops */
}

template <class Indices>
static void spread_vsite3FD(const Indices &ia, real a, real b,
                            const rvec x[], rvec f[])
{
    // Same order as in constructPositions, the virtual site is last
    int ai = ia[0];
    int aj = ia[1];
    int ak = ia[2];
    int av = ia[3];
    rvec fv;
    copy_rvec(f[av], fv);
    clear_rvec(f[av]);

    // The virtual site is at xj + b*t/|t| with t = (xi - xj) + a*(xk - xj),
    // see constr_vsite3FD. With P = 1 - t t^T/|t|^2 and c = b/|t| the
    // derivatives of the position are c P with respect to xi, a c P with
    // respect to xk, and 1 - (1 + a) c P with respect to xj.
    rvec xji, xjk, t;
    rvec_sub(x[ai], x[aj], xji);
    rvec_sub(x[ak], x[aj], xjk);
    for (int m = 0; m < DIM; m++)
    {
        t[m] = xji[m] + a*xjk[m];
    }
    real invl  = gmx::invsqrt(iprod(t, t));
    real c     = b*invl;
    real fproj = iprod(t, fv)*invl*invl;
    rvec temp;
    for (int m = 0; m < DIM; m++)
    {
        temp[m] = c*(fv[m] - fproj*t[m]);
    }
    for (int m = 0; m < DIM; m++)
    {
        f[ai][m] += temp[m];
        f[ak][m] += a*temp[m];
        f[aj][m] += fv[m] - (1 + a)*temp[m];
    }
}

template <class Indices>
static void spread_vsite3OUT(const Indices &ia, real a, real b, real c,
                             const rvec x[],
                             rvec f[])
{
    rvec    xij, xik, fv, fj, fk;
    real    cfx, cfy, cfz;
    int     av, ai, aj, ak;

    // ACT uses H O H but the code below needs
    // O H H so renumber the atoms.
    av = ia[3];
    ai = ia[1];
    aj = ia[0];
    ak = ia[2];

    rvec_sub(x[aj], x[ai], xij);
    rvec_sub(x[ak], x[ai], xik);
    /* 6 Flops */

    copy_rvec(f[av], fv);

    cfx = c*fv[XX];
    cfy = c*fv[YY];
    cfz = c*fv[ZZ];
    /* 3 Flops */

    fj[XX] = a*fv[XX]     -  xik[ZZ]*cfy +  xik[YY]*cfz;
    fj[YY] =  xik[ZZ]*cfx + a*fv[YY]     -  xik[XX]*cfz;
    fj[ZZ] = -xik[YY]*cfx +  xik[XX]*cfy + a*fv[ZZ];

    fk[XX] = b*fv[XX]     +  xij[ZZ]*cfy -  xij[YY]*cfz;
    fk[YY] = -xij[ZZ]*cfx + b*fv[YY]     +  xij[XX]*cfz;
    fk[ZZ] =  xij[YY]*cfx -  xij[XX]*cfy + b*fv[ZZ];
    /* 30 Flops */

    f[ai][XX] += fv[XX] - fj[XX] - fk[XX];
    f[ai][YY] += fv[YY] - fj[YY] - fk[YY];
    f[ai][ZZ] += fv[ZZ] - fj[ZZ] - fk[ZZ];
    rvec_inc(f[aj], fj);
    rvec_inc(f[ak], fk);
    /* 15 Flops */

    clear_rvec(f[av]);
    /* TOTAL: 54 flops */
}

namespace alexandria
{

/*! \brief Floating point operations of one virtual site, see Timings::addFlops
 * Counted from the routines above, each addition, subtraction and
 * multiplication as one operation and each division and square root as
 * five, such that 1/sqrt(x) counts as ten. The multiplication of a
 * parameter by the sign of VSITE3OUT is not counted.
 * \param[in] itype  The type of virtual site
 * \param[in] spread Whether to count distributing the force rather
 *                   than constructing the position
 * \return the number of operations
 */
static double vsiteFlops(InteractionType itype, bool spread)
{
    switch(itype)
    {
    case InteractionType::VSITE1:
        return spread ? 3 : 0;
    case InteractionType::VSITE2:
        return spread ? 13 : 9;
    case InteractionType::VSITE2FD:
        // Spreading checks the distance with a square root
        return spread ? 55 : 25;
    case InteractionType::VSITE3:
    case InteractionType::VSITE3S:
        return spread ? 20 : 17;
    case InteractionType::VSITE3FD:
        return spread ? 65 : 34;
    case InteractionType::VSITE3OUT:
    case InteractionType::VSITE3OUTS:
        return spread ? 54 : 33;
    case InteractionType::VSITE4:
    case InteractionType::VSITE4S:
    case InteractionType::VSITE4S3:
        return spread ? 27 : 24;
    default:
        return 0;
    }
}

/*! \brief Floating point operations of all virtual sites of a topology
 * \param[in] top    The topology
 * \param[in] spread Whether to count distributing the forces, see vsiteFlops
 * \return the number of operations for one conformation
 */
static double topologyVsiteFlops(const Topology *top, bool spread)
{
    double flops = 0;
    for (const auto &entry: top->entries())
    {
        if (isVsite(entry.first))
        {
            flops += vsiteFlops(entry.first, spread)*entry.second.size();
        }
    }
    return flops;
}

void VsiteHandler::constructPositions(const Topology         *top,
                                      std::vector<gmx::RVec> *coordinates) const
{
    // The operations are counted in the region that is active, since
    // timing a single conformation would take longer than the work
    Timings::instance().addFlops(topologyVsiteFlops(top, false));
    // Ugly shortcut...
    std::vector<gmx::RVec>    &x      = *coordinates;
    for (const auto &entry: top->entries())
    {
        // Only construct positions for vsites.
        if (!isVsite(entry.first))
        {
            continue;
        }
        for (const auto &vs : entry.second)
        {
            auto &atomIndices = vs->atomIndices();
            auto &params      = vs->params();
            if (params.empty())
            {
                GMX_THROW(gmx::InternalError("No parameters to generate virtual sites in topology. constructPositions should not be called before parameters are assigned."));
            }
            int ai = atomIndices[0];
            int aj = atomIndices[1];
            int ak = 0;
            int al = 0;
            switch(entry.first)
            {
            case InteractionType::VSITE1:
                constr_vsite1(x[ai], x[aj]);
                break;
            case InteractionType::VSITE2:
                ak = atomIndices[2];
                constr_vsite2(x[ai], x[aj], x[ak], params[vsite2A]);
                if (debug)
                {
                    fprintf(debug, "vsite a = %g\n", params[vsite2A]);
                }
                break;
            case InteractionType::VSITE2FD:
                ak = atomIndices[2];
                constr_vsite2fd(x[ai], x[aj], x[ak], params[vsite2fdA]);
                if (debug)
                {
                    fprintf(debug, "vsite a = %g\n", params[vsite2fdA]);
                }
                break;
            case InteractionType::VSITE3:
                ak = atomIndices[2];
                al = atomIndices[3];
                constr_vsite3(x[ai], x[aj],  x[ak], x[al],
                              params[vsite3A], params[vsite3B]);
                break;
            case InteractionType::VSITE3S:
                ak = atomIndices[2];
                al = atomIndices[3];
                constr_vsite3(x[ai], x[aj],  x[ak], x[al],
                              params[vsite3sA], params[vsite3sA]);
                break;
            case InteractionType::VSITE3FD:
                ak = atomIndices[2];
                al = atomIndices[3];
                constr_vsite3FD(x[ai], x[aj], x[ak], x[al],
                                params[vsite3fdA], params[vsite3fdB]);
                break;
            case InteractionType::VSITE3OUT:
                {
                    auto  vsite3out_vs = static_cast <const Vsite3OUT*> (vs->self());
                    ak                 = atomIndices[2];
                    al                 = atomIndices[3];
                    constr_vsite3OUT(x[aj], x[ai], x[ak], x[al],
                                     params[vsite3outA], params[vsite3outB],
                                     vsite3out_vs->sign() * params[vsite3outC]);
                    if (debug)
                    {
                        fprintf(debug, "vs3out: sign= %2d, A=%g, B=%g, C=%g\n",
                                vsite3out_vs->sign(), params[vsite3outA], params[vsite3outB],
                                vsite3out_vs->sign() * params[vsite3outC]);
                    }
                }
                break;
            case InteractionType::VSITE3OUTS:
                {
                    auto  vsite3out_vs = static_cast <const Vsite3OUT*> (vs->self());
                    ak                 = atomIndices[2];
                    al                 = atomIndices[3];
                    constr_vsite3OUT(x[aj], x[ai], x[ak], x[al],
                                     params[vsite3outsA], params[vsite3outsA],
                                     vsite3out_vs->sign() * params[vsite3outsC]);
                    if (debug)
                    {
                        fprintf(debug, "vs3out: sign= %2d, A=%g, C=%g\n",
                                vsite3out_vs->sign(), params[vsite3outsA],
                                vsite3out_vs->sign() * params[vsite3outsC]);
                    }
                }
                break;
            case InteractionType::VSITE4:
                constr_vsite4(x[atomIndices[0]], x[atomIndices[1]], x[atomIndices[2]], x[atomIndices[3]], 
                              x[atomIndices[4]], params[vsite4A], params[vsite4B], params[vsite4C]);
                break;
            case InteractionType::VSITE4S:
                constr_vsite4(x[atomIndices[0]], x[atomIndices[1]], x[atomIndices[2]], x[atomIndices[3]], 
                              x[atomIndices[4]], params[vsite4sA], params[vsite4sB], params[vsite4sB]);
                break;
            case InteractionType::VSITE4S3:
                // The three symmetric atoms have the same weight, the
                // rest is on the central atom, which comes first
                constr_vsite4(x[atomIndices[1]], x[atomIndices[0]], x[atomIndices[2]], x[atomIndices[3]], 
                              x[atomIndices[4]], params[vsite4s3A], params[vsite4s3A], params[vsite4s3A]);
                break;
            default: // throws
                GMX_THROW(gmx::InternalError(gmx::formatString("Virtual site type %s not implemented yet.", interactionTypeToString(entry.first).c_str()).c_str()));
            }
        }
    }
}

void VsiteHandler::distributeForces(const Topology               *top,
                                    const std::vector<gmx::RVec> &coords,
                                    std::vector<gmx::RVec>       *forces) const
{
    // As in constructPositions for one conformation
    Timings::instance().addFlops(topologyVsiteFlops(top, true));
    // Some shortcuts
    const auto &x    = as_rvec_array(coords.data());
    rvec       *f    = as_rvec_array((*forces).data());
    for (const auto &entry: top->entries())
    {
        if (!isVsite(entry.first))
        {
            continue;
        }
        for (const auto &vs : entry.second)
        {
            auto &atomIndices = vs->atomIndices();
            auto &params      = vs->params();

            switch(entry.first)
            {
            case InteractionType::VSITE1:
                spread_vsite1(atomIndices, f);
                break;
            case InteractionType::VSITE2:
                spread_vsite2(atomIndices, params[vsite2A], f);
                break;
            case InteractionType::VSITE2FD:
                spread_vsite2fd(x[atomIndices[0]], x[atomIndices[1]], params[vsite2fdA],
                                f[atomIndices[0]], f[atomIndices[1]], f[atomIndices[2]]);
                break;
            case InteractionType::VSITE3:
                spread_vsite3(f[atomIndices[0]], f[atomIndices[1]],
                              f[atomIndices[2]], f[atomIndices[3]],
                              params[vsite3A], params[vsite3B]);
                break;
            case InteractionType::VSITE3S:
                spread_vsite3(f[atomIndices[0]], f[atomIndices[1]],
                              f[atomIndices[2]], f[atomIndices[3]],
                              params[vsite3sA], params[vsite3sA]);
                break;
            case InteractionType::VSITE3FD:
                spread_vsite3FD(atomIndices, params[vsite3fdA], params[vsite3fdB], x, f);
                break;
            case InteractionType::VSITE3OUT:
                {
                    auto vsite3out_vs = static_cast <const Vsite3OUT*> (vs->self());
                    spread_vsite3OUT(atomIndices, params[vsite3outA], params[vsite3outB],
                                     vsite3out_vs->sign() * params[vsite3outC],
                                     x, f);
                }
                break;
            case InteractionType::VSITE3OUTS:
                {
                    auto vsite3out_vs = static_cast <const Vsite3OUT*> (vs->self());
                    spread_vsite3OUT(atomIndices, params[vsite3outsA], params[vsite3outsA], 
                                     vsite3out_vs->sign() * params[vsite3outsC],
                                     x, f);
                }
                break;
            case InteractionType::VSITE4:
                spread_vsite4(f[atomIndices[0]], f[atomIndices[1]], f[atomIndices[2]], f[atomIndices[3]], f[atomIndices[4]],
                              params[vsite4A], params[vsite4B], params[vsite4C]);
                break;
            case InteractionType::VSITE4S:
                spread_vsite4(f[atomIndices[0]], f[atomIndices[1]], f[atomIndices[2]], f[atomIndices[3]], f[atomIndices[4]],
                              params[vsite4sA], params[vsite4sB], params[vsite4sB]);
                break;
            case InteractionType::VSITE4S3:
                // The same weights as in constructing
                spread_vsite4(f[atomIndices[1]], f[atomIndices[0]], f[atomIndices[2]], f[atomIndices[3]], f[atomIndices[4]],
                              params[vsite4s3A], params[vsite4s3A], params[vsite4s3A]);
                break;
            default: // Constructing throws, so that is sufficient
                break;
            }
        }
    }
}

//! The particles of a virtual site in the local arrays of spreadSlots
static constexpr std::array<int, 5> localIndices = { 0, 1, 2, 3, 4 };

/*! \brief Construct one virtual site in many conformations
 * The coordinates of the parents in each slot, relative to the first one,
 * are copied to a local array, the site is constructed by the same routine
 * as for one conformation and its position, the last particle, is copied
 * back. The constructions do not depend on the origin, and the relative
 * coordinates are exact, see coordinate.h.
 * \tparam       n         The number of particles, known at compile time such
 *                         that the loops over them are unrolled
 * \param[in]    ai        The particles, the virtual site last
 * \param[inout] x         The coordinates, see VsiteHandler::constructPositions
 * \param[in]    stride    The stride
 * \param[in]    nslot     The number of slots
 * \param[in]    construct Constructs the site from the local coordinates
 */
template <size_t n, class Construct>
static void constructSlots(const std::vector<int> &ai,
                           Coordinate             *x,
                           size_t                  stride,
                           size_t                  nslot,
                           Construct               construct)
{
    GMX_RELEASE_ASSERT(ai.size() == n, "Unexpected number of particles of a virtual site");
    Coordinate *row[n][DIM];
    for(size_t k = 0; k < n; k++)
    {
        for(int m = 0; m < DIM; m++)
        {
            row[k][m] = x + (DIM*ai[k] + m)*stride;
        }
    }
    rvec xl[n];
    for(size_t c = 0; c < nslot; c++)
    {
        for(size_t k = 0; k+1 < n; k++)
        {
            for(int m = 0; m < DIM; m++)
            {
                xl[k][m] = fromCoordinate(coordinateDifference(row[k][m][c], row[0][m][c]));
            }
        }
        construct(xl);
        for(int m = 0; m < DIM; m++)
        {
            row[n-1][m][c] = coordinateSum(row[0][m][c], toCoordinate(xl[n-1][m]));
        }
    }
}

/*! \brief Distribute the force on one virtual site in many conformations
 * As constructSlots, the coordinates relative to the first particle and
 * the forces of the particles in each slot are copied to local arrays,
 * the force is distributed by the same
 * routine as for one conformation and the forces are copied back.
 * \tparam       n      The number of particles, as for constructSlots
 * \param[in]    ai     The particles, the virtual site last
 * \param[in]    x      The coordinates, see VsiteHandler::constructPositions
 * \param[inout] f      The forces
 * \param[in]    stride The stride
 * \param[in]    nslot  The number of slots
 * \param[in]    spread Distributes the force with the local coordinates
 *                      and forces, the particles are numbered from zero,
 *                      see localIndices
 */
template <size_t n, typename T, class Spread>
static void spreadSlots(const std::vector<int> &ai,
                        const Coordinate       *x,
                        T                      *f,
                        size_t                  stride,
                        size_t                  nslot,
                        Spread                  spread)
{
    GMX_RELEASE_ASSERT(ai.size() == n, "Unexpected number of particles of a virtual site");
    const Coordinate *xrow[n][DIM];
    T                *frow[n][DIM];
    for(size_t k = 0; k < n; k++)
    {
        for(int m = 0; m < DIM; m++)
        {
            xrow[k][m] = x + (DIM*ai[k] + m)*stride;
            frow[k][m] = f + (DIM*ai[k] + m)*stride;
        }
    }
    rvec xl[n], fl[n];
    for(size_t c = 0; c < nslot; c++)
    {
        for(size_t k = 0; k < n; k++)
        {
            for(int m = 0; m < DIM; m++)
            {
                xl[k][m] = fromCoordinate(coordinateDifference(xrow[k][m][c], xrow[0][m][c]));
                fl[k][m] = frow[k][m][c];
            }
        }
        spread(xl, fl);
        for(size_t k = 0; k < n; k++)
        {
            for(int m = 0; m < DIM; m++)
            {
                frow[k][m][c] = fl[k][m];
            }
        }
    }
}

/*! \brief Construct the virtual sites of many conformations at once
 * See VsiteHandler::constructPositions for the arguments.
 */
static void constructStrided(const Topology *top,
                             Coordinate     *x,
                             size_t          stride,
                             size_t          nslot,
                             size_t          nconf)
{
    static const size_t timer = Timings::instance().index("Virtual sites");
    const ScopedTimer   scoped(timer);
    Timings::instance().addFlops(nconf*topologyVsiteFlops(top, false));
    for (const auto &entry: top->entries())
    {
        if (!isVsite(entry.first))
        {
            continue;
        }
        for (const auto &vs : entry.second)
        {
            auto &ai     = vs->atomIndices();
            auto &params = vs->params();
            if (params.empty())
            {
                GMX_THROW(gmx::InternalError("No parameters to generate virtual sites in topology. constructPositions should not be called before parameters are assigned."));
            }
            // The same routines and arguments as for one conformation
            switch(entry.first)
            {
            case InteractionType::VSITE1:
                constructSlots<2>(ai, x, stride, nslot, [](rvec *xl)
                               { constr_vsite1(xl[0], xl[1]); });
                break;
            case InteractionType::VSITE2:
                constructSlots<3>(ai, x, stride, nslot, [&params](rvec *xl)
                               { constr_vsite2(xl[0], xl[1], xl[2], params[vsite2A]); });
                break;
            case InteractionType::VSITE2FD:
                constructSlots<3>(ai, x, stride, nslot, [&params](rvec *xl)
                               { constr_vsite2fd(xl[0], xl[1], xl[2], params[vsite2fdA]); });
                break;
            case InteractionType::VSITE3:
                constructSlots<4>(ai, x, stride, nslot, [&params](rvec *xl)
                               { constr_vsite3(xl[0], xl[1], xl[2], xl[3],
                                               params[vsite3A], params[vsite3B]); });
                break;
            case InteractionType::VSITE3S:
                constructSlots<4>(ai, x, stride, nslot, [&params](rvec *xl)
                               { constr_vsite3(xl[0], xl[1], xl[2], xl[3],
                                               params[vsite3sA], params[vsite3sA]); });
                break;
            case InteractionType::VSITE3FD:
                constructSlots<4>(ai, x, stride, nslot, [&params](rvec *xl)
                               { constr_vsite3FD(xl[0], xl[1], xl[2], xl[3],
                                                 params[vsite3fdA], params[vsite3fdB]); });
                break;
            case InteractionType::VSITE3OUT:
                {
                    double sc = static_cast<const Vsite3OUT*>(vs->self())->sign()*params[vsite3outC];
                    constructSlots<4>(ai, x, stride, nslot, [&params, sc](rvec *xl)
                                   { constr_vsite3OUT(xl[1], xl[0], xl[2], xl[3],
                                                      params[vsite3outA], params[vsite3outB], sc); });
                }
                break;
            case InteractionType::VSITE3OUTS:
                {
                    double sc = static_cast<const Vsite3OUT*>(vs->self())->sign()*params[vsite3outsC];
                    constructSlots<4>(ai, x, stride, nslot, [&params, sc](rvec *xl)
                                   { constr_vsite3OUT(xl[1], xl[0], xl[2], xl[3],
                                                      params[vsite3outsA], params[vsite3outsA], sc); });
                }
                break;
            case InteractionType::VSITE4:
                constructSlots<5>(ai, x, stride, nslot, [&params](rvec *xl)
                               { constr_vsite4(xl[0], xl[1], xl[2], xl[3], xl[4],
                                               params[vsite4A], params[vsite4B], params[vsite4C]); });
                break;
            case InteractionType::VSITE4S:
                constructSlots<5>(ai, x, stride, nslot, [&params](rvec *xl)
                               { constr_vsite4(xl[0], xl[1], xl[2], xl[3], xl[4],
                                               params[vsite4sA], params[vsite4sB], params[vsite4sB]); });
                break;
            case InteractionType::VSITE4S3:
                // The three symmetric atoms have the same weight, the
                // rest is on the central atom, which comes first
                constructSlots<5>(ai, x, stride, nslot, [&params](rvec *xl)
                               { constr_vsite4(xl[1], xl[0], xl[2], xl[3], xl[4],
                                               params[vsite4s3A], params[vsite4s3A], params[vsite4s3A]); });
                break;
            default: // throws
                GMX_THROW(gmx::InternalError(gmx::formatString("Virtual site type %s not implemented yet.", interactionTypeToString(entry.first).c_str()).c_str()));
            }
        }
    }
}

/*! \brief Distribute the forces of many conformations back to atoms
 * \tparam T The type of the coordinates and forces
 * See VsiteHandler::distributeForces for the arguments.
 */
template <typename T>
static void distributeStrided(const Topology *top,
                              const Coordinate *x,
                              T              *f,
                              size_t          stride,
                              size_t          nslot,
                              size_t          nconf)
{
    static const size_t timer = Timings::instance().index("Virtual sites");
    const ScopedTimer   scoped(timer);
    Timings::instance().addFlops(nconf*topologyVsiteFlops(top, true));
    for (const auto &entry: top->entries())
    {
        if (!isVsite(entry.first))
        {
            continue;
        }
        for (const auto &vs : entry.second)
        {
            auto &ai     = vs->atomIndices();
            auto &params = vs->params();
            const auto &ia = localIndices;
            // The same routines and arguments as for one conformation
            switch(entry.first)
            {
            case InteractionType::VSITE1:
                spreadSlots<2>(ai, x, f, stride, nslot, [&](rvec *, rvec *fl)
                            { spread_vsite1(ia, fl); });
                break;
            case InteractionType::VSITE2:
                spreadSlots<3>(ai, x, f, stride, nslot, [&](rvec *, rvec *fl)
                            { spread_vsite2(ia, params[vsite2A], fl); });
                break;
            case InteractionType::VSITE2FD:
                spreadSlots<3>(ai, x, f, stride, nslot, [&](rvec *xl, rvec *fl)
                            { spread_vsite2fd(xl[0], xl[1], params[vsite2fdA], fl[0], fl[1], fl[2]); });
                break;
            case InteractionType::VSITE3:
                spreadSlots<4>(ai, x, f, stride, nslot, [&](rvec *, rvec *fl)
                            { spread_vsite3(fl[0], fl[1], fl[2], fl[3],
                                            params[vsite3A], params[vsite3B]); });
                break;
            case InteractionType::VSITE3S:
                spreadSlots<4>(ai, x, f, stride, nslot, [&](rvec *, rvec *fl)
                            { spread_vsite3(fl[0], fl[1], fl[2], fl[3],
                                            params[vsite3sA], params[vsite3sA]); });
                break;
            case InteractionType::VSITE3FD:
                spreadSlots<4>(ai, x, f, stride, nslot, [&](rvec *xl, rvec *fl)
                            { spread_vsite3FD(ia, params[vsite3fdA], params[vsite3fdB], xl, fl); });
                break;
            case InteractionType::VSITE3OUT:
                {
                    double sc = static_cast<const Vsite3OUT*>(vs->self())->sign()*params[vsite3outC];
                    spreadSlots<4>(ai, x, f, stride, nslot, [&](rvec *xl, rvec *fl)
                                { spread_vsite3OUT(ia, params[vsite3outA], params[vsite3outB], sc,
                                                   xl, fl); });
                }
                break;
            case InteractionType::VSITE3OUTS:
                {
                    double sc = static_cast<const Vsite3OUT*>(vs->self())->sign()*params[vsite3outsC];
                    spreadSlots<4>(ai, x, f, stride, nslot, [&](rvec *xl, rvec *fl)
                                { spread_vsite3OUT(ia, params[vsite3outsA], params[vsite3outsA], sc,
                                                   xl, fl); });
                }
                break;
            case InteractionType::VSITE4:
                spreadSlots<5>(ai, x, f, stride, nslot, [&](rvec *, rvec *fl)
                            { spread_vsite4(fl[0], fl[1], fl[2], fl[3], fl[4],
                                            params[vsite4A], params[vsite4B], params[vsite4C]); });
                break;
            case InteractionType::VSITE4S:
                spreadSlots<5>(ai, x, f, stride, nslot, [&](rvec *, rvec *fl)
                            { spread_vsite4(fl[0], fl[1], fl[2], fl[3], fl[4],
                                            params[vsite4sA], params[vsite4sB], params[vsite4sB]); });
                break;
            case InteractionType::VSITE4S3:
                // The same weights as in constructing
                spreadSlots<5>(ai, x, f, stride, nslot, [&](rvec *, rvec *fl)
                            { spread_vsite4(fl[1], fl[0], fl[2], fl[3], fl[4],
                                            params[vsite4s3A], params[vsite4s3A], params[vsite4s3A]); });
                break;
            default: // Constructing throws, so that is sufficient
                break;
            }
        }
    }
}

void VsiteHandler::constructPositions(const Topology *top,
                                      Coordinate     *x,
                                      size_t          stride,
                                      size_t          nslot,
                                      size_t          nconf) const
{
    constructStrided(top, x, stride, nslot, nconf);
}

void VsiteHandler::distributeForces(const Topology   *top,
                                    const Coordinate *x,
                                    double           *f,
                                    size_t          stride,
                                    size_t          nslot,
                                    size_t          nconf) const
{
    distributeStrided(top, x, f, stride, nslot, nconf);
}

void VsiteHandler::distributeForces(const Topology   *top,
                                    const Coordinate *x,
                                    float            *f,
                                    size_t          stride,
                                    size_t          nslot,
                                    size_t          nconf) const
{
    distributeStrided(top, x, f, stride, nslot, nconf);
}

} // namespace

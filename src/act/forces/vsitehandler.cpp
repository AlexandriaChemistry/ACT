/*
 * This file is part of the GROMACS molecular simulation package.
 *
 * Copyright (c) 1991-2000, University of Groningen, The Netherlands.
 * Copyright (c) 2001-2004, The GROMACS development team.
 * Copyright (c) 2013,2014,2015,2016,2017,2018, by the GROMACS development team, led by
 * Mark Abraham, David van der Spoel, Berk Hess, and Erik Lindahl,
 * and including many others, as listed in the AUTHORS file in the
 * top-level source directory and at http://www.gromacs.org.
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
#include <set>
#include <vector>

#include "act/forcefield/forcefield_parametername.h"
#include "gromacs/math/functions.h"
#include "gromacs/math/vec.h"
#include "gromacs/pbcutil/ishift.h"
#include "gromacs/pbcutil/mshift.h"
#include "gromacs/pbcutil/pbc.h"
#include "gromacs/timing/wallcycle.h"
#include "gromacs/topology/ifunc.h"
#include "gromacs/topology/topology.h"
#include "gromacs/utility/exceptions.h"
#include "gromacs/utility/fatalerror.h"

using gmx::RVec;

static int pbc_rvec_sub(const t_pbc *pbc, const rvec xi, const rvec xj, rvec dx)
{
    if (pbc)
    {
        return pbc_dx_aiuc(pbc, xi, xj, dx);
    }
    else
    {
        rvec_sub(xi, xj, dx);
        return CENTRAL;
    }
}

/* Vsite construction routines */

static void constr_vsite1(const rvec xi, rvec x)
{
    copy_rvec(xi, x);
    /* TOTAL: 0 flops */
}

static void constr_vsite2(const rvec xi, const rvec xj, rvec x, real a, const t_pbc *pbc)
{
    real b = 1 - a;
    /* 1 flop */

    if (pbc)
    {
        rvec dx;
        pbc_dx_aiuc(pbc, xj, xi, dx);
        x[XX] = xi[XX] + a*dx[XX];
        x[YY] = xi[YY] + a*dx[YY];
        x[ZZ] = xi[ZZ] + a*dx[ZZ];
    }
    else
    {
        x[XX] = b*xi[XX] + a*xj[XX];
        x[YY] = b*xi[YY] + a*xj[YY];
        x[ZZ] = b*xi[ZZ] + a*xj[ZZ];
        /* 9 Flops */
    }

    /* TOTAL: 10 flops */
}

static void constr_vsite2fd(const rvec xi, const rvec xj, rvec x, real a, const t_pbc *pbc)
{
    rvec dx;

    if (pbc)
    {
        pbc_dx_aiuc(pbc, xj, xi, dx);
    }
    else
    {
        rvec_sub(xj, xi, dx);
    }
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

static void constr_vsite3FD(const rvec xi, const rvec xj, const rvec xk, rvec x, real a, real b,
                            const t_pbc *pbc)
{
    rvec xij, xjk, temp;
    real c;

    pbc_rvec_sub(pbc, xi, xj, xij);
    pbc_rvec_sub(pbc, xk, xj, xjk);
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

static void constr_vsite3FAD(const rvec xi, const rvec xj, const rvec xk, rvec x, real a, real b, const t_pbc *pbc)
{
    gmx_fatal(FARGS, "Fix code in constr_vsite3FAD before using");
    rvec xij, xjk, xp;
    real a1, b1, c1, invdij;

    pbc_rvec_sub(pbc, xj, xi, xij);
    pbc_rvec_sub(pbc, xk, xj, xjk);
    /* 6 flops */

    invdij = gmx::invsqrt(iprod(xij, xij));
    c1     = invdij * invdij * iprod(xij, xjk);
    xp[XX] = xjk[XX] - c1*xij[XX];
    xp[YY] = xjk[YY] - c1*xij[YY];
    xp[ZZ] = xjk[ZZ] - c1*xij[ZZ];
    a1     = a*invdij;
    b1     = b*gmx::invsqrt(iprod(xp, xp));
    /* 45 */

    x[XX] = xi[XX] + a1*xij[XX] + b1*xp[XX];
    x[YY] = xi[YY] + a1*xij[YY] + b1*xp[YY];
    x[ZZ] = xi[ZZ] + a1*xij[ZZ] + b1*xp[ZZ];
    /* 12 Flops */

    /* TOTAL: 63 flops */
}

static void constr_vsite3OUT(const rvec xi, const rvec xj, const rvec xk, rvec x,
                             real a, real b, real c, const t_pbc *pbc)
{
    rvec xij, xkj, temp;

    pbc_rvec_sub(pbc, xj, xi, xij);
    pbc_rvec_sub(pbc, xj, xk, xkj);
    cprod(xij, xkj, temp);
    /* 15 Flops */

    x[XX] = xj[XX] + a*xij[XX] + b*xkj[XX] + c*temp[XX];
    x[YY] = xj[YY] + a*xij[YY] + b*xkj[YY] + c*temp[YY];
    x[ZZ] = xj[ZZ] + a*xij[ZZ] + b*xkj[ZZ] + c*temp[ZZ];
    /* 18 Flops */

    /* TOTAL: 33 flops */
}

static gmx_unused void constr_vsite4FD(const rvec xi, const rvec xj, const rvec xk, const rvec xl, rvec x,
                                       real a, real b, real c, const t_pbc *pbc)
{
    gmx_fatal(FARGS, "Fix code in constr_vsite4FD before using");
    rvec xij, xjk, xjl, temp;
    real d;

    pbc_rvec_sub(pbc, xj, xi, xij);
    pbc_rvec_sub(pbc, xk, xj, xjk);
    pbc_rvec_sub(pbc, xl, xj, xjl);
    /* 9 flops */

    /* temp goes from i to a point on the plane jkl */
    temp[XX] = xij[XX] + a*xjk[XX] + b*xjl[XX];
    temp[YY] = xij[YY] + a*xjk[YY] + b*xjl[YY];
    temp[ZZ] = xij[ZZ] + a*xjk[ZZ] + b*xjl[ZZ];
    /* 12 flops */

    d = c*gmx::invsqrt(iprod(temp, temp));
    /* 6 + 10 flops */

    x[XX] = xi[XX] + d*temp[XX];
    x[YY] = xi[YY] + d*temp[YY];
    x[ZZ] = xi[ZZ] + d*temp[ZZ];
    /* 6 Flops */

    /* TOTAL: 43 flops */
}

static gmx_unused void constr_vsite4FDN(const rvec xi, const rvec xj, const rvec xk, const rvec xl, rvec x,
                                        real a, real b, real c, const t_pbc *pbc)
{
    gmx_fatal(FARGS, "Fix code in constr_vsite4FDN before using");
    rvec xij, xik, xil, ra, rb, rja, rjb, rm;
    real d;

    pbc_rvec_sub(pbc, xj, xi, xij);
    pbc_rvec_sub(pbc, xk, xi, xik);
    pbc_rvec_sub(pbc, xl, xi, xil);
    /* 9 flops */

    ra[XX] = a*xik[XX];
    ra[YY] = a*xik[YY];
    ra[ZZ] = a*xik[ZZ];

    rb[XX] = b*xil[XX];
    rb[YY] = b*xil[YY];
    rb[ZZ] = b*xil[ZZ];

    /* 6 flops */

    rvec_sub(ra, xij, rja);
    rvec_sub(rb, xij, rjb);
    /* 6 flops */

    cprod(rja, rjb, rm);
    /* 9 flops */

    d = c*gmx::invsqrt(norm2(rm));
    /* 5+5+1 flops */

    x[XX] = xi[XX] + d*rm[XX];
    x[YY] = xi[YY] + d*rm[YY];
    x[ZZ] = xi[ZZ] + d*rm[ZZ];
    /* 6 Flops */

    /* TOTAL: 47 flops */
}


static gmx_unused int constr_vsiten(const t_iatom *ia, const t_iparams ip[],
                                    rvec *x, const t_pbc *pbc)
{
    gmx_fatal(FARGS, "Fix code in constr_vsiten before using");
    rvec x1, dx;
    dvec dsum;
    int  n3, av, ai;
    real a;

    n3 = 3*ip[ia[0]].vsiten.n;
    av = ia[1];
    ai = ia[2];
    copy_rvec(x[ai], x1);
    clear_dvec(dsum);
    for (int i = 3; i < n3; i += 3)
    {
        ai = ia[i+2];
        a  = ip[ia[i]].vsiten.a;
        if (pbc)
        {
            pbc_dx_aiuc(pbc, x[ai], x1, dx);
        }
        else
        {
            rvec_sub(x[ai], x1, dx);
        }
        dsum[XX] += a*dx[XX];
        dsum[YY] += a*dx[YY];
        dsum[ZZ] += a*dx[ZZ];
        /* 9 Flops */
    }

    x[av][XX] = x1[XX] + dsum[XX];
    x[av][YY] = x1[YY] + dsum[YY];
    x[av][ZZ] = x1[ZZ] + dsum[ZZ];

    return n3;
}

static void spread_vsite1(const t_iatom ia[], rvec f[])
{
    t_iatom ai = ia[0];
    t_iatom av = ia[1];
    rvec_inc(f[ai], f[av]);
    clear_rvec(f[av]);
}

static void spread_vsite2(const std::vector<int> &ia, real a,
                          const rvec x[],
                          rvec f[], rvec fshift[],
                          const t_pbc *pbc, const t_graph *g)
{
    rvec    fi, fj, dx;
    t_iatom av, ai, aj;
    ivec    di;
    int     siv, sij;

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

    if (g)
    {
        ivec_sub(SHIFT_IVEC(g, ai), SHIFT_IVEC(g, av), di);
        siv = IVEC2IS(di);
        ivec_sub(SHIFT_IVEC(g, ai), SHIFT_IVEC(g, aj), di);
        sij = IVEC2IS(di);
    }
    else if (pbc)
    {
        siv = pbc_dx_aiuc(pbc, x[ai], x[av], dx);
        sij = pbc_dx_aiuc(pbc, x[ai], x[aj], dx);
    }
    else
    {
        siv = CENTRAL;
        sij = CENTRAL;
    }

    if (fshift && (siv != CENTRAL || sij != CENTRAL))
    {
        rvec_inc(fshift[siv], f[av]);
        rvec_dec(fshift[CENTRAL], fi);
        rvec_dec(fshift[sij], fj);
    }

    /* TOTAL: 13 flops */
}

static void spread_vsite2fd(const t_iatom ia[], real a,
                            const rvec x[],
                            rvec f[], const t_pbc *pbc)
{
    rvec    dx;
    t_iatom av, ai, aj;

    av = ia[3];
    ai = ia[1];
    aj = ia[2];

    if (pbc)
    {
        pbc_rvec_sub(pbc, x[aj], x[ai], dx);
    }
    else
    {
        rvec_sub(x[aj], x[ai], dx);
    }
    real gamma = a*gmx::invsqrt(iprod(dx, dx));
    rvec xis;
    for (int m = 0; m < DIM; m++)
    {
        xis[m] = x[ai][m] + gamma*dx[m];
    }
    rvec p;
    svmul(iprod(xis, f[av])/iprod(xis, xis), xis, p);
    rvec df;
    for(int m = 0; m < DIM; m++)
    {
        df[m] = gamma*(f[av][m]-p[m]);
    }
    
    rvec_inc(f[ai], df);
    rvec_dec(f[aj], df);
    
    /* 6 Flops */
    clear_rvec(f[av]);
}

static void spread_vsite3(const std::vector<int> &indices,
                          real a, real b,
                          const rvec x[], rvec f[], rvec fshift[],
                          const t_pbc *pbc, const t_graph *g)
{
    rvec    fi, fj, fk, dx;
    int     av, ai, aj, ak;
    ivec    di;
    int     siv, sij, sik;

    ai = indices[0];
    aj = indices[1];
    ak = indices[2];
    av = indices[3];

    svmul(1 - a - b, f[av], fj);
    svmul(        a, f[av], fi);
    svmul(        b, f[av], fk);
    /* 11 flops */

    rvec_inc(f[ai], fi);
    rvec_inc(f[aj], fj);
    rvec_inc(f[ak], fk);
    /* 9 Flops */

    if (g)
    {
        ivec_sub(SHIFT_IVEC(g, ai), SHIFT_IVEC(g, av), di);
        siv = IVEC2IS(di);
        ivec_sub(SHIFT_IVEC(g, ai), SHIFT_IVEC(g, aj), di);
        sij = IVEC2IS(di);
        ivec_sub(SHIFT_IVEC(g, ai), SHIFT_IVEC(g, ak), di);
        sik = IVEC2IS(di);
    }
    else if (pbc)
    {
        siv = pbc_dx_aiuc(pbc, x[ai], x[av], dx);
        sij = pbc_dx_aiuc(pbc, x[ai], x[aj], dx);
        sik = pbc_dx_aiuc(pbc, x[ai], x[ak], dx);
    }
    else
    {
        siv = CENTRAL;
        sij = CENTRAL;
        sik = CENTRAL;
    }

    if (fshift && (siv != CENTRAL || sij != CENTRAL || sik != CENTRAL))
    {
        rvec_inc(fshift[siv], f[av]);
        rvec_dec(fshift[CENTRAL], fi);
        rvec_dec(fshift[sij], fj);
        rvec_dec(fshift[sik], fk);
    }
    clear_rvec(f[av]);
    /* TOTAL: 20 flops */
}

static void spread_vsite3FD(const t_iatom ia[], real a, real b,
                            const rvec x[], rvec f[], rvec fshift[],
                            gmx_bool VirCorr, matrix dxdf,
                            const t_pbc *pbc, const t_graph *g)
{
    real    c, invl, fproj, a1;
    rvec    xvi, xij, xjk, xix, fv, temp;
    t_iatom av, ai, aj, ak;
    int     svi, sji, skj;
    ivec    di;

    av = ia[4];
    ai = ia[2];
    aj = ia[1];
    ak = ia[3];
    copy_rvec(f[av], fv);

    sji = pbc_rvec_sub(pbc, x[aj], x[ai], xij);
    skj = pbc_rvec_sub(pbc, x[ak], x[aj], xjk);
    /* 6 flops */

    /* xix goes from i to point x on the line jk */
    xix[XX] = xij[XX]+a*xjk[XX];
    xix[YY] = xij[YY]+a*xjk[YY];
    xix[ZZ] = xij[ZZ]+a*xjk[ZZ];
    /* 6 flops */

    invl = gmx::invsqrt(iprod(xix, xix));
    c    = b*invl;
    /* 4 + ?10? flops */

    fproj = iprod(xix, fv)*invl*invl; /* = (xix . f)/(xix . xix) */

    temp[XX] = c*(fv[XX]-fproj*xix[XX]);
    temp[YY] = c*(fv[YY]-fproj*xix[YY]);
    temp[ZZ] = c*(fv[ZZ]-fproj*xix[ZZ]);
    /* 16 */

    /* c is already calculated in constr_vsite3FD
       storing c somewhere will save 26 flops!     */

    a1         = 1 - a;
    f[ai][XX] += fv[XX] - temp[XX];
    f[ai][YY] += fv[YY] - temp[YY];
    f[ai][ZZ] += fv[ZZ] - temp[ZZ];
    f[aj][XX] += a1*temp[XX];
    f[aj][YY] += a1*temp[YY];
    f[aj][ZZ] += a1*temp[ZZ];
    f[ak][XX] += a*temp[XX];
    f[ak][YY] += a*temp[YY];
    f[ak][ZZ] += a*temp[ZZ];
    /* 19 Flops */

    if (g)
    {
        ivec_sub(SHIFT_IVEC(g, ia[1]), SHIFT_IVEC(g, ai), di);
        svi = IVEC2IS(di);
        ivec_sub(SHIFT_IVEC(g, aj), SHIFT_IVEC(g, ai), di);
        sji = IVEC2IS(di);
        ivec_sub(SHIFT_IVEC(g, ak), SHIFT_IVEC(g, aj), di);
        skj = IVEC2IS(di);
    }
    else if (pbc)
    {
        svi = pbc_rvec_sub(pbc, x[av], x[ai], xvi);
    }
    else
    {
        svi = CENTRAL;
    }

    if (fshift && (svi != CENTRAL || sji != CENTRAL || skj != CENTRAL))
    {
        rvec_dec(fshift[svi], fv);
        fshift[CENTRAL][XX] += fv[XX] - (1 + a)*temp[XX];
        fshift[CENTRAL][YY] += fv[YY] - (1 + a)*temp[YY];
        fshift[CENTRAL][ZZ] += fv[ZZ] - (1 + a)*temp[ZZ];
        fshift[    sji][XX] += temp[XX];
        fshift[    sji][YY] += temp[YY];
        fshift[    sji][ZZ] += temp[ZZ];
        fshift[    skj][XX] += a*temp[XX];
        fshift[    skj][YY] += a*temp[YY];
        fshift[    skj][ZZ] += a*temp[ZZ];
    }

    if (VirCorr)
    {
        /* When VirCorr=TRUE, the virial for the current forces is not
         * calculated from the redistributed forces. This means that
         * the effect of non-linear virtual site constructions on the virial
         * needs to be added separately. This contribution can be calculated
         * in many ways, but the simplest and cheapest way is to use
         * the first constructing atom ai as a reference position in space:
         * subtract (xv-xi)*fv and add (xj-xi)*fj + (xk-xi)*fk.
         */
        rvec xiv;

        pbc_rvec_sub(pbc, x[av], x[ai], xiv);

        for (int i = 0; i < DIM; i++)
        {
            for (int j = 0; j < DIM; j++)
            {
                /* As xix is a linear combination of j and k, use that here */
                dxdf[i][j] += -xiv[i]*fv[j] + xix[i]*temp[j];
            }
        }
    }

    /* TOTAL: 61 flops */
}

static void spread_vsite3FAD(const t_iatom ia[], real a, real b,
                             const rvec x[],
                             rvec f[], rvec fshift[],
                             gmx_bool VirCorr, matrix dxdf,
                             const t_pbc *pbc, const t_graph *g)
{
    rvec    xvi, xij, xjk, xperp, Fpij, Fppp, fv, f1, f2, f3;
    real    a1, b1, c1, c2, invdij, invdij2, invdp, fproj;
    t_iatom av, ai, aj, ak;
    int     svi, sji, skj, d;
    ivec    di;

    av = ia[4];
    ai = ia[1];
    aj = ia[2];
    ak = ia[3];
    copy_rvec(f[ia[1]], fv);

    sji = pbc_rvec_sub(pbc, x[aj], x[ai], xij);
    skj = pbc_rvec_sub(pbc, x[ak], x[aj], xjk);
    /* 6 flops */

    invdij    = gmx::invsqrt(iprod(xij, xij));
    invdij2   = invdij * invdij;
    c1        = iprod(xij, xjk) * invdij2;
    xperp[XX] = xjk[XX] - c1*xij[XX];
    xperp[YY] = xjk[YY] - c1*xij[YY];
    xperp[ZZ] = xjk[ZZ] - c1*xij[ZZ];
    /* xperp in plane ijk, perp. to ij */
    invdp = gmx::invsqrt(iprod(xperp, xperp));
    a1    = a*invdij;
    b1    = b*invdp;
    /* 45 flops */

    /* a1, b1 and c1 are already calculated in constr_vsite3FAD
       storing them somewhere will save 45 flops!     */

    fproj = iprod(xij, fv)*invdij2;
    svmul(fproj,                      xij,  Fpij);    /* proj. f on xij */
    svmul(iprod(xperp, fv)*invdp*invdp, xperp, Fppp); /* proj. f on xperp */
    svmul(b1*fproj,                   xperp, f3);
    /* 23 flops */

    rvec_sub(fv, Fpij, f1); /* f1 = f - Fpij */
    rvec_sub(f1, Fppp, f2); /* f2 = f - Fpij - Fppp */
    for (d = 0; (d < DIM); d++)
    {
        f1[d] *= a1;
        f2[d] *= b1;
    }
    /* 12 flops */

    c2         = 1 + c1;
    f[ai][XX] += fv[XX] - f1[XX] + c1*f2[XX] + f3[XX];
    f[ai][YY] += fv[YY] - f1[YY] + c1*f2[YY] + f3[YY];
    f[ai][ZZ] += fv[ZZ] - f1[ZZ] + c1*f2[ZZ] + f3[ZZ];
    f[aj][XX] +=          f1[XX] - c2*f2[XX] - f3[XX];
    f[aj][YY] +=          f1[YY] - c2*f2[YY] - f3[YY];
    f[aj][ZZ] +=          f1[ZZ] - c2*f2[ZZ] - f3[ZZ];
    f[ak][XX] +=                      f2[XX];
    f[ak][YY] +=                      f2[YY];
    f[ak][ZZ] +=                      f2[ZZ];
    /* 30 Flops */

    if (g)
    {
        ivec_sub(SHIFT_IVEC(g, ia[1]), SHIFT_IVEC(g, ai), di);
        svi = IVEC2IS(di);
        ivec_sub(SHIFT_IVEC(g, aj), SHIFT_IVEC(g, ai), di);
        sji = IVEC2IS(di);
        ivec_sub(SHIFT_IVEC(g, ak), SHIFT_IVEC(g, aj), di);
        skj = IVEC2IS(di);
    }
    else if (pbc)
    {
        svi = pbc_rvec_sub(pbc, x[av], x[ai], xvi);
    }
    else
    {
        svi = CENTRAL;
    }

    if (fshift && (svi != CENTRAL || sji != CENTRAL || skj != CENTRAL))
    {
        rvec_dec(fshift[svi], fv);
        fshift[CENTRAL][XX] += fv[XX] - f1[XX] - (1-c1)*f2[XX] + f3[XX];
        fshift[CENTRAL][YY] += fv[YY] - f1[YY] - (1-c1)*f2[YY] + f3[YY];
        fshift[CENTRAL][ZZ] += fv[ZZ] - f1[ZZ] - (1-c1)*f2[ZZ] + f3[ZZ];
        fshift[    sji][XX] +=          f1[XX] -    c1 *f2[XX] - f3[XX];
        fshift[    sji][YY] +=          f1[YY] -    c1 *f2[YY] - f3[YY];
        fshift[    sji][ZZ] +=          f1[ZZ] -    c1 *f2[ZZ] - f3[ZZ];
        fshift[    skj][XX] +=                          f2[XX];
        fshift[    skj][YY] +=                          f2[YY];
        fshift[    skj][ZZ] +=                          f2[ZZ];
    }

    if (VirCorr)
    {
        rvec xiv;
        int  i, j;

        pbc_rvec_sub(pbc, x[av], x[ai], xiv);

        for (i = 0; i < DIM; i++)
        {
            for (j = 0; j < DIM; j++)
            {
                /* Note that xik=xij+xjk, so we have to add xij*f2 */
                dxdf[i][j] +=
                    -xiv[i]*fv[j]
                    + xij[i]*(f1[j] + (1 - c2)*f2[j] - f3[j])
                    + xjk[i]*f2[j];
            }
        }
    }

    /* TOTAL: 113 flops */
}

static void spread_vsite3OUT(const t_iatom ia[], real a, real b, real c,
                             const rvec x[],
                             rvec f[], rvec fshift[],
                             gmx_bool VirCorr, matrix dxdf,
                             const t_pbc *pbc, const t_graph *g)
{
    rvec    xvi, xij, xik, fv, fj, fk;
    real    cfx, cfy, cfz;
    int     av, ai, aj, ak;
    ivec    di;
    int     svi, sji, ski;

    av = ia[4];
    ai = ia[1];
    aj = ia[2];
    ak = ia[3];

    sji = pbc_rvec_sub(pbc, x[aj], x[ai], xij);
    ski = pbc_rvec_sub(pbc, x[ak], x[ai], xik);
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

    if (g)
    {
        ivec_sub(SHIFT_IVEC(g, ia[1]), SHIFT_IVEC(g, ai), di);
        svi = IVEC2IS(di);
        ivec_sub(SHIFT_IVEC(g, aj), SHIFT_IVEC(g, ai), di);
        sji = IVEC2IS(di);
        ivec_sub(SHIFT_IVEC(g, ak), SHIFT_IVEC(g, ai), di);
        ski = IVEC2IS(di);
    }
    else if (pbc)
    {
        svi = pbc_rvec_sub(pbc, x[av], x[ai], xvi);
    }
    else
    {
        svi = CENTRAL;
    }

    if (fshift && (svi != CENTRAL || sji != CENTRAL || ski != CENTRAL))
    {
        rvec_dec(fshift[svi], fv);
        fshift[CENTRAL][XX] += fv[XX] - fj[XX] - fk[XX];
        fshift[CENTRAL][YY] += fv[YY] - fj[YY] - fk[YY];
        fshift[CENTRAL][ZZ] += fv[ZZ] - fj[ZZ] - fk[ZZ];
        rvec_inc(fshift[sji], fj);
        rvec_inc(fshift[ski], fk);
    }

    if (VirCorr)
    {
        rvec xiv;

        pbc_rvec_sub(pbc, x[av], x[ai], xiv);

        for (int i = 0; i < DIM; i++)
        {
            for (int j = 0; j < DIM; j++)
            {
                dxdf[i][j] += -xiv[i]*fv[j] + xij[i]*fj[j] + xik[i]*fk[j];
            }
        }
    }

    /* TOTAL: 54 flops */
}

static gmx_unused void spread_vsite4FD(const t_iatom ia[], real a, real b, real c,
                                       const rvec x[], rvec f[], rvec fshift[],
                                       gmx_bool VirCorr, matrix dxdf,
                                       const t_pbc *pbc, const t_graph *g)
{
    real    d, invl, fproj, a1;
    rvec    xvi, xij, xjk, xjl, xix, fv, temp;
    int     av, ai, aj, ak, al;
    ivec    di;
    int     svi, sji, skj, slj, m;

    av = ia[1];
    ai = ia[2];
    aj = ia[3];
    ak = ia[4];
    al = ia[5];

    sji = pbc_rvec_sub(pbc, x[aj], x[ai], xij);
    skj = pbc_rvec_sub(pbc, x[ak], x[aj], xjk);
    slj = pbc_rvec_sub(pbc, x[al], x[aj], xjl);
    /* 9 flops */

    /* xix goes from i to point x on the plane jkl */
    for (m = 0; m < DIM; m++)
    {
        xix[m] = xij[m] + a*xjk[m] + b*xjl[m];
    }
    /* 12 flops */

    invl = gmx::invsqrt(iprod(xix, xix));
    d    = c*invl;
    /* 4 + ?10? flops */

    copy_rvec(f[av], fv);

    fproj = iprod(xix, fv)*invl*invl; /* = (xix . f)/(xix . xix) */

    for (m = 0; m < DIM; m++)
    {
        temp[m] = d*(fv[m] - fproj*xix[m]);
    }
    /* 16 */

    /* c is already calculated in constr_vsite3FD
       storing c somewhere will save 35 flops!     */

    a1 = 1 - a - b;
    for (m = 0; m < DIM; m++)
    {
        f[ai][m] += fv[m] - temp[m];
        f[aj][m] += a1*temp[m];
        f[ak][m] += a*temp[m];
        f[al][m] += b*temp[m];
    }
    /* 26 Flops */

    if (g)
    {
        ivec_sub(SHIFT_IVEC(g, ia[1]), SHIFT_IVEC(g, ai), di);
        svi = IVEC2IS(di);
        ivec_sub(SHIFT_IVEC(g, aj), SHIFT_IVEC(g, ai), di);
        sji = IVEC2IS(di);
        ivec_sub(SHIFT_IVEC(g, ak), SHIFT_IVEC(g, aj), di);
        skj = IVEC2IS(di);
        ivec_sub(SHIFT_IVEC(g, al), SHIFT_IVEC(g, aj), di);
        slj = IVEC2IS(di);
    }
    else if (pbc)
    {
        svi = pbc_rvec_sub(pbc, x[av], x[ai], xvi);
    }
    else
    {
        svi = CENTRAL;
    }

    if (fshift &&
        (svi != CENTRAL || sji != CENTRAL || skj != CENTRAL || slj != CENTRAL))
    {
        rvec_dec(fshift[svi], fv);
        for (m = 0; m < DIM; m++)
        {
            fshift[CENTRAL][m] += fv[m] - (1 + a + b)*temp[m];
            fshift[    sji][m] += temp[m];
            fshift[    skj][m] += a*temp[m];
            fshift[    slj][m] += b*temp[m];
        }
    }

    if (VirCorr)
    {
        rvec xiv;
        int  i, j;

        pbc_rvec_sub(pbc, x[av], x[ai], xiv);

        for (i = 0; i < DIM; i++)
        {
            for (j = 0; j < DIM; j++)
            {
                dxdf[i][j] += -xiv[i]*fv[j] + xix[i]*temp[j];
            }
        }
    }

    /* TOTAL: 77 flops */
}


static gmx_unused void spread_vsite4FDN(const t_iatom ia[], real a, real b, real c,
                                        const rvec x[], rvec f[], rvec fshift[],
                                        gmx_bool VirCorr, matrix dxdf,
                                        const t_pbc *pbc, const t_graph *g)
{
    rvec xvi, xij, xik, xil, ra, rb, rja, rjb, rab, rm, rt;
    rvec fv, fj, fk, fl;
    real invrm, denom;
    real cfx, cfy, cfz;
    ivec di;
    int  av, ai, aj, ak, al;
    int  svi, sij, sik, sil;

    /* DEBUG: check atom indices */
    av = ia[1];
    ai = ia[2];
    aj = ia[3];
    ak = ia[4];
    al = ia[5];

    copy_rvec(f[av], fv);

    sij = pbc_rvec_sub(pbc, x[aj], x[ai], xij);
    sik = pbc_rvec_sub(pbc, x[ak], x[ai], xik);
    sil = pbc_rvec_sub(pbc, x[al], x[ai], xil);
    /* 9 flops */

    ra[XX] = a*xik[XX];
    ra[YY] = a*xik[YY];
    ra[ZZ] = a*xik[ZZ];

    rb[XX] = b*xil[XX];
    rb[YY] = b*xil[YY];
    rb[ZZ] = b*xil[ZZ];

    /* 6 flops */

    rvec_sub(ra, xij, rja);
    rvec_sub(rb, xij, rjb);
    rvec_sub(rb, ra, rab);
    /* 9 flops */

    cprod(rja, rjb, rm);
    /* 9 flops */

    invrm = gmx::invsqrt(norm2(rm));
    denom = invrm*invrm;
    /* 5+5+2 flops */

    cfx = c*invrm*fv[XX];
    cfy = c*invrm*fv[YY];
    cfz = c*invrm*fv[ZZ];
    /* 6 Flops */

    cprod(rm, rab, rt);
    /* 9 flops */

    rt[XX] *= denom;
    rt[YY] *= denom;
    rt[ZZ] *= denom;
    /* 3flops */

    fj[XX] = (        -rm[XX]*rt[XX]) * cfx + ( rab[ZZ]-rm[YY]*rt[XX]) * cfy + (-rab[YY]-rm[ZZ]*rt[XX]) * cfz;
    fj[YY] = (-rab[ZZ]-rm[XX]*rt[YY]) * cfx + (        -rm[YY]*rt[YY]) * cfy + ( rab[XX]-rm[ZZ]*rt[YY]) * cfz;
    fj[ZZ] = ( rab[YY]-rm[XX]*rt[ZZ]) * cfx + (-rab[XX]-rm[YY]*rt[ZZ]) * cfy + (        -rm[ZZ]*rt[ZZ]) * cfz;
    /* 30 flops */

    cprod(rjb, rm, rt);
    /* 9 flops */

    rt[XX] *= denom*a;
    rt[YY] *= denom*a;
    rt[ZZ] *= denom*a;
    /* 3flops */

    fk[XX] = (          -rm[XX]*rt[XX]) * cfx + (-a*rjb[ZZ]-rm[YY]*rt[XX]) * cfy + ( a*rjb[YY]-rm[ZZ]*rt[XX]) * cfz;
    fk[YY] = ( a*rjb[ZZ]-rm[XX]*rt[YY]) * cfx + (          -rm[YY]*rt[YY]) * cfy + (-a*rjb[XX]-rm[ZZ]*rt[YY]) * cfz;
    fk[ZZ] = (-a*rjb[YY]-rm[XX]*rt[ZZ]) * cfx + ( a*rjb[XX]-rm[YY]*rt[ZZ]) * cfy + (          -rm[ZZ]*rt[ZZ]) * cfz;
    /* 36 flops */

    cprod(rm, rja, rt);
    /* 9 flops */

    rt[XX] *= denom*b;
    rt[YY] *= denom*b;
    rt[ZZ] *= denom*b;
    /* 3flops */

    fl[XX] = (          -rm[XX]*rt[XX]) * cfx + ( b*rja[ZZ]-rm[YY]*rt[XX]) * cfy + (-b*rja[YY]-rm[ZZ]*rt[XX]) * cfz;
    fl[YY] = (-b*rja[ZZ]-rm[XX]*rt[YY]) * cfx + (          -rm[YY]*rt[YY]) * cfy + ( b*rja[XX]-rm[ZZ]*rt[YY]) * cfz;
    fl[ZZ] = ( b*rja[YY]-rm[XX]*rt[ZZ]) * cfx + (-b*rja[XX]-rm[YY]*rt[ZZ]) * cfy + (          -rm[ZZ]*rt[ZZ]) * cfz;
    /* 36 flops */

    f[ai][XX] += fv[XX] - fj[XX] - fk[XX] - fl[XX];
    f[ai][YY] += fv[YY] - fj[YY] - fk[YY] - fl[YY];
    f[ai][ZZ] += fv[ZZ] - fj[ZZ] - fk[ZZ] - fl[ZZ];
    rvec_inc(f[aj], fj);
    rvec_inc(f[ak], fk);
    rvec_inc(f[al], fl);
    /* 21 flops */

    if (g)
    {
        ivec_sub(SHIFT_IVEC(g, av), SHIFT_IVEC(g, ai), di);
        svi = IVEC2IS(di);
        ivec_sub(SHIFT_IVEC(g, aj), SHIFT_IVEC(g, ai), di);
        sij = IVEC2IS(di);
        ivec_sub(SHIFT_IVEC(g, ak), SHIFT_IVEC(g, ai), di);
        sik = IVEC2IS(di);
        ivec_sub(SHIFT_IVEC(g, al), SHIFT_IVEC(g, ai), di);
        sil = IVEC2IS(di);
    }
    else if (pbc)
    {
        svi = pbc_rvec_sub(pbc, x[av], x[ai], xvi);
    }
    else
    {
        svi = CENTRAL;
    }

    if (fshift && (svi != CENTRAL || sij != CENTRAL || sik != CENTRAL || sil != CENTRAL))
    {
        rvec_dec(fshift[svi], fv);
        fshift[CENTRAL][XX] += fv[XX] - fj[XX] - fk[XX] - fl[XX];
        fshift[CENTRAL][YY] += fv[YY] - fj[YY] - fk[YY] - fl[YY];
        fshift[CENTRAL][ZZ] += fv[ZZ] - fj[ZZ] - fk[ZZ] - fl[ZZ];
        rvec_inc(fshift[sij], fj);
        rvec_inc(fshift[sik], fk);
        rvec_inc(fshift[sil], fl);
    }

    if (VirCorr)
    {
        rvec xiv;
        int  i, j;

        pbc_rvec_sub(pbc, x[av], x[ai], xiv);

        for (i = 0; i < DIM; i++)
        {
            for (j = 0; j < DIM; j++)
            {
                dxdf[i][j] += -xiv[i]*fv[j] + xij[i]*fj[j] + xik[i]*fk[j] + xil[i]*fl[j];
            }
        }
    }

    /* Total: 207 flops (Yuck!) */
}


static gmx_unused int spread_vsiten(const t_iatom ia[], const t_iparams ip[],
                                    const rvec x[], rvec f[], rvec fshift[],
                                    const t_pbc *pbc, const t_graph *g)
{
    rvec xv, dx, fi;
    int  n3, av, i, ai;
    real a;
    ivec di;
    int  siv;

    n3 = 3*ip[ia[0]].vsiten.n;
    av = ia[1];
    copy_rvec(x[av], xv);

    for (i = 0; i < n3; i += 3)
    {
        ai = ia[i+2];
        if (g)
        {
            ivec_sub(SHIFT_IVEC(g, ai), SHIFT_IVEC(g, av), di);
            siv = IVEC2IS(di);
        }
        else if (pbc)
        {
            siv = pbc_dx_aiuc(pbc, x[ai], xv, dx);
        }
        else
        {
            siv = CENTRAL;
        }
        a = ip[ia[i]].vsiten.a;
        svmul(a, f[av], fi);
        rvec_inc(f[ai], fi);
        if (fshift && siv != CENTRAL)
        {
            rvec_inc(fshift[siv], fi);
            rvec_dec(fshift[CENTRAL], fi);
        }
        /* 6 Flops */
    }

    return n3;
}


namespace alexandria
{

VsiteHandler::VsiteHandler(matrix &box,
                           real    dt)
{
    set_pbc(&pbc_, -1, box);
    // TODO More checking.
    if (pbc_.ePBC != epbcNONE)
    {
        GMX_THROW(gmx::InvalidInputError(gmx::formatString("No support for periodic boundara conditions").c_str()));
    }
    dt_ = dt;
}

void VsiteHandler::constructPositions(const Topology          *top,
                                      std::vector<gmx::RVec>  *coordinates,
                                      const gmx_unused matrix &box)
{
    // Ugly shortcut...
    std::vector<gmx::RVec>    &x      = *coordinates;
    std::set<InteractionType>  vsites = {
        InteractionType::VSITE1,
        InteractionType::VSITE2,
        InteractionType::VSITE2FD,
        InteractionType::VSITE3,
        InteractionType::VSITE3S,
        InteractionType::VSITE3FD,
        InteractionType::VSITE3FAD,
        InteractionType::VSITE3OUT,
        InteractionType::VSITE3OUTS
    };
    for (const auto &entry: top->entries())
    {
        // Only construct positions for vsites.
        if (vsites.find(entry.first) == vsites.end())
        {
            continue;
        }
        for (const auto &vs : entry.second)
        {
            auto &atomIndices = vs->atomIndices();
            auto &params      = vs->params();
            if (params.empty())
            {
                if (debug)
                {
                    fprintf(debug, "WARNING: No parameters to generate virtual sites in topology (yet).\n");
                }
                continue;
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
                constr_vsite2(x[ai], x[aj], x[ak], params[vsite2A], &pbc_);
                if (debug)
                {
                    fprintf(debug, "vsite a = %g\n", params[vsite2A]);
                }
                break;
            case InteractionType::VSITE2FD:
                ak = atomIndices[2];
                constr_vsite2fd(x[ai], x[aj], x[ak], params[vsite2A], &pbc_);
                if (debug)
                {
                    fprintf(debug, "vsite a = %g\n", params[vsite2A]);
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
                                params[vsite3fdA], params[vsite3fdB], &pbc_);
                break;
            case InteractionType::VSITE3FAD:
                ak = atomIndices[2];
                al = atomIndices[3];
                constr_vsite3FAD(x[ai], x[aj], x[ak], x[al],
                                 params[vsite3fadA], params[vsite3fadB], &pbc_);
                break;

            case InteractionType::VSITE3OUT:
                {
                    auto  vsite3out_vs = static_cast <const Vsite3OUT*> (vs->self());
                    ak                 = atomIndices[2];
                    al                 = atomIndices[3];
                    constr_vsite3OUT(x[ai], x[aj], x[ak], x[al],
                                     params[vsite3outA], params[vsite3outB],
                                     vsite3out_vs->sign() * params[vsite3outC], &pbc_);
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
                    constr_vsite3OUT(x[ai], x[aj], x[ak], x[al],
                                     params[vsite3outsA], params[vsite3outsA],
                                     vsite3out_vs->sign() * params[vsite3outsC], &pbc_);
                    if (debug)
                    {
                        fprintf(debug, "vs3out: sign= %2d, A=%g, C=%g\n",
                                vsite3out_vs->sign(), params[vsite3outsA],
                                vsite3out_vs->sign() * params[vsite3outsC]);
                    }
                }
                break;

#ifdef LATER
            case F_VSITE4FD:
                aj = ia[3];
                ak = ia[4];
                al = ia[5];
                b1 = ip[tp].vsite.b;
                c1 = ip[tp].vsite.c;
                constr_vsite4FD(x[ai], x[aj], x[ak], x[al], x[avsite], a1, b1, c1,
                                pbc_null2);
                break;
            case F_VSITE4FDN:
                aj = ia[3];
                ak = ia[4];
                al = ia[5];
                b1 = ip[tp].vsite.b;
                c1 = ip[tp].vsite.c;
                constr_vsite4FDN(x[ai], x[aj], x[ak], x[al], x[avsite], a1, b1, c1,
                                 pbc_null2);
                break;
            case F_VSITEN:
                inc = constr_vsiten(ia, ip, x, pbc_null2);
                break;
#endif
            default: // throws
                GMX_THROW(gmx::InternalError(gmx::formatString("Virtual site type %s not implemented yet.", interactionTypeToString(entry.first).c_str()).c_str()));
            }
        }
    }
}

void VsiteHandler::distributeForces(const Topology               *top,
                                    const std::vector<gmx::RVec> &coords,
                                    std::vector<gmx::RVec>       *forces,
                                    const gmx_unused matrix      &box)
{
    bool     VirCorr = false;
    matrix   dxdf    = { { 0 } };
    rvec    *fshift  = nullptr;
    t_graph *g       = nullptr;
    // Some shortcuts
    const auto &x    = as_rvec_array(coords.data());
    rvec       *f    = as_rvec_array((*forces).data());
    for (const auto &entry: top->entries())
    {
        for (const auto &vs : entry.second)
        {
            auto &atomIndices = vs->atomIndices();
            auto &params      = vs->params();
            // Ugly hack to minimize change in underlying gromacs code
            std::vector<t_iatom> ia = { -1 };
            for(auto &ai : atomIndices)
            {
                ia.push_back(ai);
            }
            switch(entry.first)
            {
            case InteractionType::VSITE1:
                spread_vsite1(ia.data(), f);
                break;
            case InteractionType::VSITE2:
                spread_vsite2(atomIndices, params[vsite2A], x, f, fshift, &pbc_, g);
                break;
            case InteractionType::VSITE2FD:
                spread_vsite2fd(ia.data(), params[vsite2A], x, f, &pbc_);
                break;
            case InteractionType::VSITE3:
                spread_vsite3(atomIndices, params[vsite3A], params[vsite3B], x, f, fshift, &pbc_, g);
                break;
            case InteractionType::VSITE3S:
                spread_vsite3(atomIndices, params[vsite3sA], params[vsite3sA], x, f, fshift, &pbc_, g);
                //spread_vsite3s(atomIndices, params[vsite3sA], f);
                break;
            case InteractionType::VSITE3FD:
                spread_vsite3FD(ia.data(), params[vsite3A], params[vsite3B], x, f, fshift,
                                VirCorr, dxdf, &pbc_, g);
                break;
            case InteractionType::VSITE3FAD:
                spread_vsite3FAD(ia.data(), params[vsite3fadA], params[vsite3fadB], x, f, fshift, VirCorr, dxdf, &pbc_, g);
                break;
            case InteractionType::VSITE3OUT:
                {
                    auto vsite3out_vs = static_cast <const Vsite3OUT*> (vs->self());
                    spread_vsite3OUT(ia.data(), params[vsite3outA], params[vsite3outB],
                                     vsite3out_vs->sign() * params[vsite3outC],
                                     x, f, fshift, VirCorr, dxdf, &pbc_, g);
                }
                break;
            case InteractionType::VSITE3OUTS:
                {
                    auto vsite3out_vs = static_cast <const Vsite3OUT*> (vs->self());
                    spread_vsite3OUT(ia.data(), params[vsite3outsA], params[vsite3outsA], 
                                     vsite3out_vs->sign() * params[vsite3outsC],
                                     x, f, fshift, VirCorr, dxdf, &pbc_, g);
                }
                break;
                //case F_VSITE4FD:
                //b1 = ip[tp].vsite.b;
                //c1 = ip[tp].vsite.c;
                //spread_vsite4FD(ia, a1, b1, c1, x, f, fshift, VirCorr, dxdf, &pbc_2, g);
                //break;
                //case F_VSITE4FDN:
                //b1 = ip[tp].vsite.b;
                //c1 = ip[tp].vsite.c;
                //spread_vsite4FDN(ia, a1, b1, c1, x, f, fshift, VirCorr, dxdf, &pbc_2, g);
                //break;
                //case F_VSITEN:
                //      inc = spread_vsiten(ia, ip, x, f, fshift, &pbc_2, g);
                //break;
            default: // check
                break;
            }
        }
    }
}

} // namespace

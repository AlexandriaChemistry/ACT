/*
 * This file is part of the GROMACS molecular simulation package.
 *
 * Copyright (c) 1991-2000, University of Groningen, The Netherlands.
 * Copyright (c) 2001-2004, The GROMACS development team.
 * Copyright (c) 2012,2014,2015,2016,2017 by the GROMACS development team.
 * Copyright (c) 2018,2019,2020, by the GROMACS development team, led by
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
/* Adopted for ACT structure by DvdS 2023-05 */

#include "actpre.h"

#include "princ.h"

#include <cmath>

#include "gromacs/linearalgebra/nrjac.h"
#include "gromacs/math/functions.h"
#include "gromacs/math/vec.h"

#define NDIM 4

namespace alexandria
{

void principal_comp(const std::vector<int>       &index,
                    const std::vector<real>      &mass,
                    const std::vector<gmx::RVec> &x, 
                    matrix                       *trans,
                    gmx::RVec                    *inertia)
{
    int      i, j, ai, m, nrot;
    real     mm, rx, ry, rz;
    double   dd[NDIM], tvec[NDIM];
    double **inten, **ev;
#ifdef DEBUG
    real e[NDIM];
#endif
    real temp;

    inten = new double *[NDIM];
    ev    = new double *[NDIM];
    for (int i = 0; (i < NDIM); i++)
    {
        inten[i] = new double[NDIM];
        ev[i]    = new double[NDIM];
        dd[i] = 0.0;
#ifdef DEBUG
        e[i] = 0.0;
#endif
    }

    for (int i = 0; (i < NDIM); i++)
    {
        for (int m = 0; (m < NDIM); m++)
        {
            inten[i][m] = 0;
        }
    }
    for (size_t i = 0; (i < index.size()); i++)
    {
        ai = index[i];
        mm = mass[ai];
        rx = x[ai][XX];
        ry = x[ai][YY];
        rz = x[ai][ZZ];
        inten[0][0] += mm * (gmx::square(ry) + gmx::square(rz));
        inten[1][1] += mm * (gmx::square(rx) + gmx::square(rz));
        inten[2][2] += mm * (gmx::square(rx) + gmx::square(ry));
        inten[1][0] -= mm * (ry * rx);
        inten[2][0] -= mm * (rx * rz);
        inten[2][1] -= mm * (rz * ry);
    }
    inten[0][1] = inten[1][0];
    inten[0][2] = inten[2][0];
    inten[1][2] = inten[2][1];
#ifdef DEBUG
    ptrans("initial", inten, dd, e);
#endif

    for (i = 0; (i < DIM); i++)
    {
        for (m = 0; (m < DIM); m++)
        {
            (*trans)[i][m] = inten[i][m];
        }
    }

    /* Call numerical recipe routines */
    jacobi(inten, 3, dd, ev, &nrot);
#ifdef DEBUG
    ptrans("jacobi", ev, dd, e);
#endif

    /* Sort eigenvalues in ascending order */
#define SWAPPER(i)                               \
    if (std::abs(dd[(i) + 1]) < std::abs(dd[i])) \
    {                                            \
        temp = dd[i];                            \
        for (j = 0; (j < NDIM); j++)             \
        {                                        \
            tvec[j] = ev[j][i];                  \
        }                                        \
        dd[i] = dd[(i) + 1];                     \
        for (j = 0; (j < NDIM); j++)             \
        {                                        \
            ev[j][i] = ev[j][(i) + 1];           \
        }                                        \
        dd[(i) + 1] = temp;                      \
        for (j = 0; (j < NDIM); j++)             \
        {                                        \
            ev[j][(i) + 1] = tvec[j];            \
        }                                        \
    }
    SWAPPER(0)
    SWAPPER(1)
    SWAPPER(0)
#ifdef DEBUG
    ptrans("swap", ev, dd, e);
    t_trans(*trans, dd, ev);
#endif

    for (i = 0; (i < DIM); i++)
    {
        (*inertia)[i] = dd[i];
        for (m = 0; (m < DIM); m++)
        {
            (*trans)[i][m] = ev[m][i];
        }
    }

    for (i = 0; (i < NDIM); i++)
    {
      delete[] inten[i];
      delete[] ev[i];
    }
    delete[] inten;
    delete[] ev;
}

void rotate_atoms(const std::vector<int> &index,
                  std::vector<gmx::RVec> *x,
                  const matrix            trans)
{
    size_t n = x->size();
    if (!index.empty())
    {
        n = index.size();
    }
    n = index.size();
    for (size_t i = 0; i < n; i++)
    {
        size_t ii    = i;
        if (!index.empty())
        {
            ii = index[i];
        }
        real xt   = (*x)[ii][XX];
        real yt   = (*x)[ii][YY];
        real zt   = (*x)[ii][ZZ];
        (*x)[ii][XX] = trans[XX][XX] * xt + trans[XX][YY] * yt + trans[XX][ZZ] * zt;
        (*x)[ii][YY] = trans[YY][XX] * xt + trans[YY][YY] * yt + trans[YY][ZZ] * zt;
        (*x)[ii][ZZ] = trans[ZZ][XX] * xt + trans[ZZ][YY] * yt + trans[ZZ][ZZ] * zt;
    }
}

real calc_xcm(const std::vector<gmx::RVec> &x,
              const std::vector<int>       &index,
              const std::vector<ActAtom>   &atoms,
              gmx::RVec                    *xcm,
              bool                          bQ)
{
    clear_rvec(*xcm);
    real   tm = 0;
    size_t n  = x.size();
    if (!index.empty())
    {
        n = index.size();
    }
    for (size_t i = 0; i < n; i++)
    {
        size_t ii    = i;
        if (!index.empty())
        {
            ii = index[i];
        }
        real m0 = 1;
        if (!atoms.empty())
        {
            if (bQ)
            {
                m0 = std::abs(atoms[ii].charge());
            }
            else
            {
                m0 = atoms[ii].mass();
            }
        }
        tm += m0;
        for (int m = 0; (m < DIM); m++)
        {
            (*xcm)[m] += m0 * x[ii][m];
        }
    }
    for (int m = 0; (m < DIM); m++)
    {
        (*xcm)[m] /= tm;
    }

    return tm;
}

real sub_xcm(std::vector<gmx::RVec>       *x,
             const std::vector<int>       &index,
             const std::vector<ActAtom>   &atoms,
             gmx::RVec                    *xcm,
             bool                          bQ)
{
    real tm  = calc_xcm(*x, index, atoms, xcm, bQ);
    size_t n = x->size();
    if (!index.empty())
    {
        n = index.size();
    }
    n = index.size();
    for (size_t i = 0; i < n; i++)
    {
        size_t ii    = i;
        if (!index.empty())
        {
            ii = index[i];
        }
        rvec_dec((*x)[ii], *xcm);
    }
    return tm;
}

void add_xcm(std::vector<gmx::RVec>       *x,
             const std::vector<int>       &index,
             gmx::RVec                    &xcm)
{
    size_t n  = x->size();
    if (!index.empty())
    {
        n = index.size();
    }
    n = index.size();
    for (size_t i = 0; i < n; i++)
    {
        size_t ii    = i;
        if (!index.empty())
        {
            ii = index[i];
        }
        rvec_inc((*x)[ii], xcm);
    }
}

} // namespace

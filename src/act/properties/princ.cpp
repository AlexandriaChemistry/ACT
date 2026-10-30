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

#include "act/topology/actatom.h"
#include "act/utility/rotation.h"
#include "gromacs/math/functions.h"
#include "gromacs/math/vec.h"

namespace alexandria
{

void principal_comp(const std::vector<int>       &index,
                    const std::vector<real>      &mass,
                    const std::vector<gmx::RVec> &x, 
                    matrix                       *trans,
                    gmx::RVec                    *inertia)
{
    // The inertia tensor
    double inten[DIM][DIM] = { { 0 } };
    for (size_t i = 0; (i < index.size()); i++)
    {
        int  ai = index[i];
        real mm = mass[ai];
        real rx = x[ai][XX];
        real ry = x[ai][YY];
        real rz = x[ai][ZZ];
        inten[0][0] += mm * (gmx::square(ry) + gmx::square(rz));
        inten[1][1] += mm * (gmx::square(rx) + gmx::square(rz));
        inten[2][2] += mm * (gmx::square(rx) + gmx::square(ry));
        inten[1][0] -= mm * (ry * rx);
        inten[2][0] -= mm * (rx * rz);
        inten[2][1] -= mm * (rz * ry);
    }
    // The moments of inertia in increasing order and the principal axes,
    // which form a rotation.
    double dd[DIM], ev[DIM][DIM];
    diagonalizeSymmetric(inten, dd, ev);
    for (int i = 0; (i < DIM); i++)
    {
        (*inertia)[i] = dd[i];
        for (int m = 0; (m < DIM); m++)
        {
            (*trans)[i][m] = ev[i][m];
        }
    }
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
                m0 = std::abs(atoms[ii].totalCharge());
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
             const gmx::RVec              &xcm)
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

/*
 * This source file is part of the Alexandria Chemistry Toolkit.
 *
 * Copyright (C) 2014-2022
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
#include "forcecomputerutils.h"

#include <cmath>

#include "gromacs/math/vec.h"
    
real bond_angle(const rvec xi,
                const rvec xj,
                const rvec xk,
                rvec r_ij,
                rvec r_kj,
                real *costh)
/* Return value is the angle between the bonds i-j and j-k */
{
    real th;

    rvec_sub(xi, xj, r_ij);
    rvec_sub(xk, xj, r_kj);

    *costh = cos_angle(r_ij, r_kj);
    th     = std::acos(*costh);

    return th;
}

real dih_angle(const rvec xi,
               const rvec xj,
               const rvec xk,
               const rvec xl,
               rvec r_ij,
               rvec r_kj,
               rvec r_kl,
               rvec m,
               rvec n)
{
    rvec_sub(xi, xj, r_ij);
    rvec_sub(xk, xj, r_kj);
    rvec_sub(xk, xl, r_kl);

    cprod(r_ij, r_kj, m);
    cprod(r_kj, r_kl, n);
    real phi  = gmx_angle(m, n);
    real ipr  = iprod(r_ij, n);
    real sign = (ipr < 0.0) ? -1.0 : 1.0;

    return sign*phi;
}




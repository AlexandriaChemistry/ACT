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
#include "gromacs/math/vectypes.h"

/*! \brief Calculate bond-angle. 
 * \param[in] xi The first coordinate i
 * \param[in] xj The second coordinate j
 * \param[in] xk The third coordinate k
 * \param[out] r_ij Distance vector r_j - r_i
 * \param[out] r_kj Distance vector r_j - r_k
 * \param[out] costh The cosine between the bond vectors
 * \return the bond angle
 */
real bond_angle(const rvec  xi, 
                const rvec  xj, 
                const rvec  xk,
                rvec        r_ij, 
                rvec        r_kj, 
                real       *costh);

/*! \brief Calculate dihedral-angle.
 * \param[in] xi The first coordinate i
 * \param[in] xj The second coordinate j
 * \param[in] xk The third coordinate k
 * \param[in] xl The third coordinate l
 * \param[out] r_ij Distance vector r_j - r_i
 * \param[out] r_kj Distance vector r_j - r_k
 * \param[out] r_kl Distance vector r_l - r_k
 * \param[out] m    Normal to plane i-j-k
 * \param[out] n    Normal to plane j-k-l
 * \return the dihedral angle
 */
real dih_angle(const rvec xi,
               const rvec xj,
               const rvec xk,
               const rvec xl,
               rvec       r_ij,
               rvec       r_kj,
               rvec       r_kl,
               rvec       m,
               rvec       n);


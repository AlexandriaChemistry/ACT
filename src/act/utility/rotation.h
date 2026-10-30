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
 * Diagonalization of symmetric 3x3 matrices and fitting of structures,
 * using the Eigen library.
 * \author David van der Spoel <david.vanderspoel@icm.uu.se>
 */
#ifndef ACT_UTILITY_ROTATION_H
#define ACT_UTILITY_ROTATION_H

#include <vector>

#include "gromacs/math/vectypes.h"
#include "gromacs/utility/real.h"

namespace alexandria
{

/*! \brief Diagonalize a symmetric 3x3 matrix
 * The eigenvalues are sorted in increasing order. The eigenvectors are
 * normalized, the largest component of the first two is positive, and
 * the third is their cross product, such that they form a rotation
 * rather than a reflection.
 * \param[in]  A    The matrix, only the lower triangle is used
 * \param[out] eval The eigenvalues
 * \param[out] evec The eigenvectors, evec[k] belongs to eval[k]
 */
void diagonalizeSymmetric(const double A[DIM][DIM],
                          double       eval[DIM],
                          double       evec[DIM][DIM]);

/*! \brief Rotate coordinates to fit a reference structure
 * The rotation minimizes the weighted root mean square deviation, see
 * W. Kabsch, Acta Cryst. A 32 (1976) 922. Both structures should have
 * their (weighted) center in the origin. A reflection is never used,
 * also not for planar structures.
 * \param[in]    weight The weight of each particle, zero for particles
 *                      that should not be used in the fit
 * \param[in]    xref   The reference coordinates
 * \param[inout] x      The coordinates to rotate
 */
void do_fit(const std::vector<real>      &weight,
            const std::vector<gmx::RVec> &xref,
            std::vector<gmx::RVec>       *x);

/*! \brief Compute the rotation matrix of the fit, see do_fit
 * \param[in]  weight The weight of each particle
 * \param[in]  xref   The reference coordinates
 * \param[in]  x      The coordinates to rotate
 * \param[out] R      The rotation, x is fitted by R x
 */
void calcFitRotation(const std::vector<real>      &weight,
                     const std::vector<gmx::RVec> &xref,
                     const std::vector<gmx::RVec> &x,
                     double                        R[DIM][DIM]);

} // namespace alexandria

#endif

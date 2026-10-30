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
 * Implements the diagonalization of symmetric 3x3 matrices and the
 * fitting of structures.
 * \author David van der Spoel <david.vanderspoel@icm.uu.se>
 */
#include "actpre.h"

#include "rotation.h"

#include <cmath>

#include <Eigen/Dense>

#include "gromacs/utility/gmxassert.h"

namespace alexandria
{

void diagonalizeSymmetric(const double A[DIM][DIM],
                          double       eval[DIM],
                          double       evec[DIM][DIM])
{
    Eigen::Matrix3d M;
    for (int m = 0; m < DIM; m++)
    {
        for (int n = 0; n <= m; n++)
        {
            M(m, n) = A[m][n];
            M(n, m) = A[m][n];
        }
    }
    // The eigenvalues are in increasing order, the eigenvectors in the columns
    Eigen::SelfAdjointEigenSolver<Eigen::Matrix3d> solver(M);
    GMX_RELEASE_ASSERT(solver.info() == Eigen::Success, "Diagonalization failed");
    Eigen::Matrix3d v = solver.eigenvectors();
    // A unique choice of the signs, the largest component positive
    for (int k = 0; k < DIM-1; k++)
    {
        Eigen::Index imax;
        v.col(k).cwiseAbs().maxCoeff(&imax);
        if (v(imax, k) < 0)
        {
            v.col(k) *= -1;
        }
    }
    // The third axis completes a right-handed system
    v.col(DIM-1) = v.col(0).cross(v.col(1));
    for (int k = 0; k < DIM; k++)
    {
        eval[k] = solver.eigenvalues()(k);
        for (int m = 0; m < DIM; m++)
        {
            evec[k][m] = v(m, k);
        }
    }
}

void calcFitRotation(const std::vector<real>      &weight,
                     const std::vector<gmx::RVec> &xref,
                     const std::vector<gmx::RVec> &x,
                     double                        R[DIM][DIM])
{
    GMX_RELEASE_ASSERT(weight.size() == xref.size() && x.size() == xref.size(),
                       "Number of particles does not match");
    // The weighted correlation matrix of the coordinates
    Eigen::Matrix3d H = Eigen::Matrix3d::Zero();
    for (size_t i = 0; i < x.size(); i++)
    {
        if (weight[i] != 0)
        {
            Eigen::Vector3d a(x[i][XX], x[i][YY], x[i][ZZ]);
            Eigen::Vector3d b(xref[i][XX], xref[i][YY], xref[i][ZZ]);
            H += weight[i] * a * b.transpose();
        }
    }
    Eigen::JacobiSVD<Eigen::Matrix3d, Eigen::ComputeFullU | Eigen::ComputeFullV> svd(H);
    Eigen::Matrix3d U = svd.matrixU();
    Eigen::Matrix3d V = svd.matrixV();
    // Change the sign of the axis with the smallest singular value if the
    // rotation would be a reflection
    Eigen::Matrix3d D = Eigen::Matrix3d::Identity();
    if ((V * U.transpose()).determinant() < 0)
    {
        D(DIM-1, DIM-1) = -1;
    }
    Eigen::Matrix3d rot = V * D * U.transpose();
    for (int m = 0; m < DIM; m++)
    {
        for (int n = 0; n < DIM; n++)
        {
            R[m][n] = rot(m, n);
        }
    }
}

void do_fit(const std::vector<real>      &weight,
            const std::vector<gmx::RVec> &xref,
            std::vector<gmx::RVec>       *x)
{
    double R[DIM][DIM];
    calcFitRotation(weight, xref, *x, R);
    for (auto &xi : *x)
    {
        gmx::RVec xold = xi;
        for (int m = 0; m < DIM; m++)
        {
            xi[m] = R[m][XX]*xold[XX] + R[m][YY]*xold[YY] + R[m][ZZ]*xold[ZZ];
        }
    }
}

} // namespace alexandria

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
 * Tests of the diagonalization of symmetric matrices and the fitting of
 * structures.
 * \author David van der Spoel <david.vanderspoel@icm.uu.se>
 */

#include "actpre.h"

#include <cmath>
#include <random>
#include <vector>

#include <gtest/gtest.h>

#include "act/utility/rotation.h"
#include "gromacs/math/vec.h"

namespace alexandria
{

namespace
{

//! \return the determinant of a 3x3 matrix
double determinant(const double R[DIM][DIM])
{
    return (R[0][0]*(R[1][1]*R[2][2] - R[1][2]*R[2][1]) -
            R[0][1]*(R[1][0]*R[2][2] - R[1][2]*R[2][0]) +
            R[0][2]*(R[1][0]*R[2][1] - R[1][1]*R[2][0]));
}

//! A symmetric matrix with entries from a seeded generator
void randomSymmetric(int seed, double A[DIM][DIM])
{
    std::mt19937                           gen(seed);
    std::uniform_real_distribution<double> dis(-2.0, 2.0);
    for (int m = 0; m < DIM; m++)
    {
        for (int n = 0; n <= m; n++)
        {
            A[m][n] = A[n][m] = dis(gen);
        }
    }
}

// Test: The eigenvalues and eigenvectors are correct, sorted, orthonormal
// and form a rotation, and the eigenvalues give the invariants of the matrix
TEST(DiagonalizeSymmetricTest, EigenvaluesAndEigenvectors)
{
    for (int seed : { 1, 2, 3, 1993 })
    {
        SCOPED_TRACE(seed);
        double A[DIM][DIM], eval[DIM], evec[DIM][DIM];
        randomSymmetric(seed, A);
        diagonalizeSymmetric(A, eval, evec);
        for (int k = 0; k < DIM; k++)
        {
            if (k > 0)
            {
                EXPECT_LE(eval[k-1], eval[k]);
            }
            for (int m = 0; m < DIM; m++)
            {
                double Av = 0;
                for (int n = 0; n < DIM; n++)
                {
                    Av += A[m][n]*evec[k][n];
                }
                EXPECT_NEAR(eval[k]*evec[k][m], Av, 1e-12) << k << " " << m;
            }
            for (int l = 0; l < DIM; l++)
            {
                double dot = 0;
                for (int m = 0; m < DIM; m++)
                {
                    dot += evec[k][m]*evec[l][m];
                }
                EXPECT_NEAR(k == l ? 1 : 0, dot, 1e-12) << k << " " << l;
            }
        }
        EXPECT_NEAR(1, determinant(evec), 1e-12);
        // The largest component of the first two axes is positive
        for (int k = 0; k < DIM-1; k++)
        {
            int imax = 0;
            for (int m = 1; m < DIM; m++)
            {
                if (std::abs(evec[k][m]) > std::abs(evec[k][imax]))
                {
                    imax = m;
                }
            }
            EXPECT_GT(evec[k][imax], 0) << k;
        }
        // The trace, the sum of the principal minors and the determinant
        double minors = (A[0][0]*A[1][1] - A[0][1]*A[1][0] +
                         A[0][0]*A[2][2] - A[0][2]*A[2][0] +
                         A[1][1]*A[2][2] - A[1][2]*A[2][1]);
        EXPECT_NEAR(A[0][0] + A[1][1] + A[2][2], eval[0] + eval[1] + eval[2], 1e-12);
        EXPECT_NEAR(minors, eval[0]*eval[1] + eval[0]*eval[2] + eval[1]*eval[2], 1e-12);
        EXPECT_NEAR(determinant(A), eval[0]*eval[1]*eval[2], 1e-12);
    }
}

// Test: Degenerate eigenvalues, for a sphere-like inertia tensor
TEST(DiagonalizeSymmetricTest, Degenerate)
{
    double A[DIM][DIM] = { { 2, 0, 0 }, { 0, 2, 0 }, { 0, 0, 1 } };
    double eval[DIM], evec[DIM][DIM];
    diagonalizeSymmetric(A, eval, evec);
    EXPECT_NEAR(1, eval[0], 1e-12);
    EXPECT_NEAR(2, eval[1], 1e-12);
    EXPECT_NEAR(2, eval[2], 1e-12);
    EXPECT_NEAR(1, std::abs(evec[0][ZZ]), 1e-12);
    EXPECT_NEAR(1, determinant(evec), 1e-12);
}

//! Coordinates centered at the weighted origin
std::vector<gmx::RVec> makeStructure(const std::vector<real> &w, int seed, bool planar)
{
    std::mt19937                         gen(seed);
    std::uniform_real_distribution<real> dis(-0.5, 0.5);
    std::vector<gmx::RVec>               x;
    gmx::RVec                            xc = { 0, 0, 0 };
    real                                 wt = 0;
    for (size_t i = 0; i < w.size(); i++)
    {
        gmx::RVec xi = { dis(gen), dis(gen), planar ? 0 : dis(gen) };
        x.push_back(xi);
        for (int m = 0; m < DIM; m++)
        {
            xc[m] += w[i]*xi[m];
        }
        wt += w[i];
    }
    for (auto &xi : x)
    {
        for (int m = 0; m < DIM; m++)
        {
            xi[m] -= xc[m]/wt;
        }
    }
    return x;
}

//! Rotate coordinates by a rotation around a, then b, then c
std::vector<gmx::RVec> rotate(const std::vector<gmx::RVec> &x, double a, double b, double c)
{
    double Rx[DIM][DIM] = { { 1, 0, 0 }, { 0, std::cos(a), -std::sin(a) }, { 0, std::sin(a), std::cos(a) } };
    double Ry[DIM][DIM] = { { std::cos(b), 0, std::sin(b) }, { 0, 1, 0 }, { -std::sin(b), 0, std::cos(b) } };
    double Rz[DIM][DIM] = { { std::cos(c), -std::sin(c), 0 }, { std::sin(c), std::cos(c), 0 }, { 0, 0, 1 } };
    std::vector<gmx::RVec> y = x;
    for (auto R : { Rx, Ry, Rz })
    {
        for (auto &yi : y)
        {
            gmx::RVec t = yi;
            for (int m = 0; m < DIM; m++)
            {
                yi[m] = R[m][XX]*t[XX] + R[m][YY]*t[YY] + R[m][ZZ]*t[ZZ];
            }
        }
    }
    return y;
}

//! \return the weighted root mean square deviation
double rmsd(const std::vector<real> &w, const std::vector<gmx::RVec> &x, const std::vector<gmx::RVec> &y)
{
    double msd = 0, wt = 0;
    for (size_t i = 0; i < x.size(); i++)
    {
        gmx::RVec dx;
        rvec_sub(x[i], y[i], dx);
        msd += w[i]*iprod(dx, dx);
        wt  += w[i];
    }
    return std::sqrt(msd/wt);
}

// Test: A rotated structure is fitted back onto the original
TEST(DoFitTest, RecoversRotation)
{
    std::vector<real> w = { 12, 1, 1, 16, 14, 1, 32 };
    for (bool planar : { false, true })
    {
        SCOPED_TRACE(planar ? "planar" : "3D");
        auto xref = makeStructure(w, 1332, planar);
        auto x    = rotate(xref, 0.3, -1.2, 2.5);
        EXPECT_GT(rmsd(w, xref, x), 0.1);
        double R[DIM][DIM];
        calcFitRotation(w, xref, x, R);
        // A rotation, also for planar structures
        EXPECT_NEAR(1, determinant(R), 1e-10);
        do_fit(w, xref, &x);
        EXPECT_NEAR(0, rmsd(w, xref, x), 1e-6);
    }
}

// Test: The fit of a noisy structure is optimal, and does not depend on
// particles that have no weight
TEST(DoFitTest, OptimalWithNoiseAndZeroWeights)
{
    std::vector<real> w = { 12, 1, 0, 16, 14, 1, 0, 32, 1 };
    auto xref = makeStructure(w, 1477, false);
    auto x    = rotate(xref, -0.7, 0.4, 1.1);
    std::mt19937                         gen(5);
    std::uniform_real_distribution<real> noise(-0.02, 0.02);
    for (auto &xi : x)
    {
        for (int m = 0; m < DIM; m++)
        {
            xi[m] += noise(gen);
        }
    }
    // Particles without weight do not change the fit
    double R1[DIM][DIM], R2[DIM][DIM];
    calcFitRotation(w, xref, x, R1);
    x[2] = { 5, 5, 5 };
    x[6] = { -3, 2, 7 };
    calcFitRotation(w, xref, x, R2);
    for (int m = 0; m < DIM; m++)
    {
        for (int n = 0; n < DIM; n++)
        {
            EXPECT_NEAR(R1[m][n], R2[m][n], 1e-10) << m << " " << n;
        }
    }
    do_fit(w, xref, &x);
    double best = rmsd(w, xref, x);
    // Within the noise
    EXPECT_LT(best, 0.03);
    // Any further rotation makes the fit worse
    double delta = 1e-3;
    for (int axis = 0; axis < DIM; axis++)
    {
        for (double sign : { -1.0, 1.0 })
        {
            double a = axis == 0 ? sign*delta : 0;
            double b = axis == 1 ? sign*delta : 0;
            double c = axis == 2 ? sign*delta : 0;
            EXPECT_GT(rmsd(w, xref, rotate(x, a, b, c)), best) << axis << " " << sign;
        }
    }
}

} // namespace

} // namespace alexandria

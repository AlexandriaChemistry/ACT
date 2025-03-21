/*
 * This source file is part of the Alexandria program.
 *
 * Copyright (C) 2014-2018,2021
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
 * Implements part of the alexandria program.
 * \author Mohammad Mehdi Ghahremanpour <mohammad.ghahremanpour@icm.uu.se>
 * \author David van der Spoel <david.vanderspoel@icm.uu.se>
 */
 
 
#include "../regression.h"

#include <cmath>
#include <cstdlib>

#include <gtest/gtest.h>

#include "testutils/cmdlinetest.h"
#include "testutils/refdata.h"
#include "testutils/testasserts.h"
#include "testutils/testfilemanager.h"

//! Class to test the regression analysis code
class RegressionTest : public gmx::test::CommandLineTestBase
{
    protected:
        //! Checking data structure
        gmx::test::TestReferenceChecker checker_;

        //! Init set tolecrance
        RegressionTest () : checker_(this->rootChecker())
        {
            auto tolerance = gmx::test::relativeToleranceAsFloatingPoint(1.0, 5e-2);
            checker_.setDefaultTolerance(tolerance);
        }

        /*! \brief Do the actual test for solving A x = b
         * \param[in] a A matrix
         * \param[in] b Right-hand side
         * \return error code from the solver.
         */
        int testRegression(MatrixWrapper &a, std::vector<double> b)
        {
            std::vector<double> x(a.nColumn());
            int info = a.solve(b, &x);
            if (info == 0)
            {
                checker_.checkSequence(x.begin(), x.end(), "solution");
            }
            return info;
        }
};

TEST_F (RegressionTest, Solve_A_x_is_B_2)
{
#define NCOL 3
#define NROW 3
    double        aa[NROW][NCOL] =
    {
        { 3, 0, 0 },
        { 0, 2, 0 },
        { 0, 0, 1 }
    };
    MatrixWrapper a(NCOL, NROW);
    for (int i = 0; i < NROW; i++)
    {
        for (int j = 0; j < NCOL; j++)
        {
            a.set(j, i, aa[i][j]);
        }
    }
    std::vector<double> b({ 3, 4, 5 });

    // Diagonal matrix, should give as answer ( 1, 2, 5 )
    EXPECT_FALSE(testRegression(a, b));
#undef NCOL
#undef NROW
}

TEST_F (RegressionTest, Solve_A_x_is_B_3)
{
#define NCOL 2
#define NROW 4
    double        aa[NROW][NCOL] =
    {
        { 3, 0 },
        { 0, 2 },
        { 1, 0 },
        { 1, 2 }
    };
    MatrixWrapper a(NCOL, NROW);
    for (int i = 0; i < NROW; i++)
    {
        for (int j = 0; j <NCOL; j++)
        {
            a.set(j, i, aa[i][j]);
        }
    }
    std::vector<double> b({ 6, 2, 2, 4 });
    // Answer should be ( 2, 1 )
    EXPECT_FALSE(testRegression(a, b));
#undef NCOL
#undef NROW
}

TEST_F (RegressionTest, Solve_A_x_is_B_4)
{
#define NCOL 1
#define NROW 2
    double        aa[NROW][NCOL] =
    {
        { 3 },
        { 2 }
    };
    MatrixWrapper a(NCOL, NROW);
    for (int i = 0; i < NROW; i++)
    {
        for (int j = 0; j <NCOL; j++)
        {
            a.set(j, i, aa[i][j]);
        }
    }
    std::vector<double> b({ 6, 4 });
    // Answer should be ( 2 )
    EXPECT_FALSE(testRegression(a, b));
#undef NCOL
#undef NROW
}

TEST_F (RegressionTest, Solve_A_x_is_B_5)
{
#define NCOL 2
#define NROW 2
    double        aa[NROW][NCOL] =
    {
        { 3, 0 },
        { 0, 2 }
    };
    MatrixWrapper a(NCOL, NROW);
    for (int i = 0; i < NROW; i++)
    {
        for (int j = 0; j <NCOL; j++)
        {
            a.set(j, i, aa[i][j]);
        }
    }
    std::vector<double> b({ 5, -5 });
    // Answer should be ( 1.6666, -2.5)
    EXPECT_FALSE(testRegression(a, b));
#undef NCOL
#undef NROW
}

TEST_F (RegressionTest, Solve_A_x_is_B_6)
{
#define NCOL 2
#define NROW 3
    double        aa[NROW][NCOL] =
    {
        { 1, 0 },
        { 0, 1 },
        { 1, 2 }
    };
    MatrixWrapper a(NCOL, NROW);
    for (int i = 0; i < NROW; i++)
    {
        for (int j = 0; j <NCOL; j++)
        {
            a.set(j, i, aa[i][j]);
        }
    }
    std::vector<double> b({ 1, 2, 5 });
    // Answer should be ( 1, 2 )
    EXPECT_FALSE(testRegression(a, b));
#undef NCOL
#undef NROW
}

TEST_F (RegressionTest, Fail_A_x_is_B)
{
#define NCOL 3
#define NROW 3
    double        aa[NROW][NCOL] =
    {
        { 1, 0, 0 },
        { 0, 0, 0 },
        { 0, 0, 0 }
    };
    MatrixWrapper a(NCOL, NROW);
    for (int i = 0; i < NROW; i++)
    {
        for (int j = 0; j <NCOL; j++)
        {
            a.set(j, i, aa[i][j]);
        }
    }
    std::vector<double> b({ 1, -1, 1 });
    // Answer should be ( 1, 2 )
    EXPECT_FALSE(testRegression(a, b));
#undef NCOL
#undef NROW
}


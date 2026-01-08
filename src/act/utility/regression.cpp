/*
 * This source file is part of the Alexandria Chemistry Toolkit.
 *
 * Copyright (C) 2014-2026
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

#include "regression.h"

#include <cmath>
#include <cstdio>
#include <cstring>
#include <iostream>

#include <vector>

#pragma GCC diagnostic push
#pragma GCC diagnostic ignored "-Wpragmas"
#pragma GCC diagnostic ignored "-Wdeprecated-copy-dtor"
#pragma GCC diagnostic ignored "-Wdeprecated-redundant-constexpr-static-def"
#include <Eigen/Dense>
#include <Eigen/SVD> 
#pragma GCC diagnostic pop

#include "gromacs/math/matrix.h"
#include "gromacs/math/vec.h"
#include "gromacs/utility/exceptions.h"
#include "gromacs/utility/stringutil.h"

/*! \brief Linear regression with many columns
 *
 * \param[in]  rhs  The right hand side of the equation
 * \param[in]  a    The matrix
 * \param[out] x    The solution of the problem
 * \return 0 if all is fine, the problematic row number otherwise
 */
static int multi_regression1(const std::vector<double>  *rhs,
                             double                    **a, 
                             std::vector<double>        *x)
{
    int nrow  = rhs->size();
    int ncol  = x->size();
   
    Eigen::MatrixXd Aeigen(nrow, ncol);
    Eigen::VectorXd Beigen(nrow);
    for(int i = 0; i < nrow; i++)
    {
        Beigen(i) = (*rhs)[i];
        for(int j = 0; j < ncol; j++)
        {
            Aeigen(i, j) = a[j][i];
        }
    }
    
    Eigen::ColPivHouseholderQR<Eigen::MatrixXd> dec(Aeigen);
    Eigen::VectorXd Xeigen = dec.solve(Beigen);
    for(int i = 0; i < ncol; i++)
    {
        (*x)[i] = Xeigen(i);
    }
    return 0;
}

MatrixWrapper::MatrixWrapper(int ncolumn, int nrow)
: ncol_(ncolumn), nrow_(nrow)
{
    a_ = alloc_matrix(ncolumn, nrow);
}

MatrixWrapper::MatrixWrapper(const std::vector<double> &flat, const int number, const char order)
: MatrixWrapper(order == 'C' ? flat.size()/number : number, order == 'C' ? number : flat.size()/number)
{
    GMX_RELEASE_ASSERT(order=='C' || order=='R', "Invalid flattening order. Please use 'C' or 'R'");
    if (order == 'C')
    {
        for (int j = 0; j < ncol_; j++)
        {
            for (int i = 0; i < nrow_; i++)
            {
                a_[j][i] = flat[j*number + i];
            }
        }
    }
    else  // order == 'R'
    {
        for (int j = 0; j < ncol_; j++)
        {
            for (int i = 0; i < nrow_; i++)
            {
                a_[j][i] = flat[i*number + j];
            }
        }
    }
}

MatrixWrapper::MatrixWrapper(const MatrixWrapper &source)
: MatrixWrapper(source.ncol_, source.nrow_)
{
    std::copy(&source.a_[0][0], &source.a_[0][0] + source.ncol_ * source.nrow_, &a_[0][0]);
}

MatrixWrapper::MatrixWrapper(MatrixWrapper &&source)
: a_(source.a_), ncol_(source.ncol_), nrow_(source.nrow_)
{
    source.a_ = nullptr;
}

MatrixWrapper &MatrixWrapper::operator=(const MatrixWrapper &rhs)
{
    if (&rhs != this)
    {
        MatrixWrapper tmp(rhs);
        double **d = a_;
        a_ = tmp.a_;
        tmp.a_ = d;  // When tmp goes out of scope it will free the memory
    }
    return *this;
}

MatrixWrapper &MatrixWrapper::operator=(MatrixWrapper &&rhs)
{
    MatrixWrapper tmp(std::move(rhs));
    double **d = a_;
    a_ = tmp.a_;
    tmp.a_ = d;
    return *this;
}

MatrixWrapper::~MatrixWrapper()
{
    if (a_ != nullptr)
    {
        free_matrix(a_);
    }
}

std::vector<double> MatrixWrapper::flatten(const char order) const
{
    GMX_RELEASE_ASSERT(order=='C' || order=='R', "Invalid flattening order. Please use 'C' or 'R'");
    std::vector<double> vec(ncol_*nrow_);
    if (order == 'C')
    {
        for (int j = 0; j < ncol_; j++)
        {
            for (int i = 0; i < nrow_; i++)
            {
                vec[j*nrow_+i] = a_[j][i];
            }
        }
    }
    else {  // order == 'R'
        for (int i = 0; i < nrow_; i++)
        {
            for (int j = 0; j < ncol_; j++)
            {
                vec[i*ncol_+j] = a_[j][i];
            }
        }
    }
    return vec;
}

bool MatrixWrapper::isSymmetric(double tolerance)
{
    GMX_RELEASE_ASSERT(
        nrow_ == ncol_,
        gmx::formatString("Matrix with dimensions (%d,%d) is not square!", nrow_, ncol_).c_str()
    );
    for (int j = 0; j < ncol_; j++)
    {
        for (int i = j+1; i < nrow_; i++)
        {
            double sum = a_[j][i] + a_[i][j];
            if (std::abs(sum) > tolerance &&
                std::abs(a_[j][i] - a_[i][j])/sum > tolerance)
            {
                return false;
            }
        }
    }
    return true;
}

void MatrixWrapper::averageTriangle()
{
    GMX_RELEASE_ASSERT(
        nrow_ == ncol_,
        gmx::formatString("Matrix with dimensions (%d,%d) is not square!", nrow_, ncol_).c_str()
    );
    for (int j = 0; j < ncol_; j++)
    {
        for (int i = j+1; i < nrow_; i++)
        {
            a_[j][i] = a_[i][j] = (a_[j][i]+a_[i][j])/2;
        }
    }
}

int MatrixWrapper::solve(std::vector<double> rhs, std::vector<double> *solution)
{
    return multi_regression1(&rhs, a_, solution);
}

std::string MatrixWrapper::toString() const
{
    const int FLOAT_SIZE = 13;
    std::string str;
    for (int i = 0; i < nrow_; i++)
    {
        str.append("[ ");
        for (int j = 0; j < ncol_; j++)
        {
            str.append(gmx::formatString("%-*g ", FLOAT_SIZE, a_[j][i]));
        }
        str.append("]\n");
    }
    return str;
}

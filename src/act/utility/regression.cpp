/*
 * This source file is part of the Alexandria Chemistry Toolkit.
 *
 * Copyright (C) 2014-2020
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

#include <Eigen/Dense>
#include <Eigen/SVD> 

#include "gromacs/linearalgebra/matrix.h"
#include "gromacs/math/vec.h"
#include "gromacs/utility/exceptions.h"
#include "gromacs/utility/stringutil.h"

extern "C"
{
#ifdef OLD
void dgelsd_(int* m, int* n, int* nrhs, double* a, int* lda,
             const double* b, int* ldb, double* s, double* rcond, int* rank,
             double* work, int* lwork, int* iwork, int* info );

void dgels_(const char* trans, int* m, int* n, int* nrhs, double* a, int* lda,
            const double* b, int* ldb, double* work, int* lwork, int* info );

void dgesvd_(const char* jobu, const char* jobvt, int* m, int* n, double* a,
             int* lda, double* s, double* u, int* ldu, double* vt, int* ldvt,
             double* work, int* lwork, int* info );
#endif
}

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

/*! \brief Linear regression with many columns
 *
 * \param[in]  rhs  The right hand side of the equation
 * \param[in]  a    The matrix
 * \param[out] x    The solution of the problem
 * \return 0 if all is fine, the problematic row number otherwise
 */
#ifdef OLD
static int multi_regression2(const std::vector<double>  *rhs,
                             double                    **a, 
                             std::vector<double>        *x)
{
    /* Query and allocate the optimal workspace */
    int                 nrow  = rhs->size();
    int                 ncol  = x->size();
    int                 lwork = -1;
    int                 lda   = std::max(1, nrow);
    int                 ldb   = std::max(1, std::max(nrow, ncol));
    int                 nrhs  = 1;
    int                 rank;
    double              rcond = -1.0;
    double              wkopt;
    std::vector<double> s;
    s.resize(nrow);
    // Compute length of integer array iwork according to
    // https://software.intel.com/sites/products/documentation/doclib/mkl_sa/11/mkl_lapack_examples/dgelsd_ex.c.htm
    int              smlsiz = 25;
    int              nlvl   = std::max(0L, std::lround(std::log2(std::min(nrow, ncol)/smlsiz + 1) ) + 1);
    int              liwork = 3*std::min(ncol, nrow)*nlvl + 11*std::min(nrow, ncol);
    std::vector<int> iwork;
    iwork.resize(liwork);
    int              info;
    bool             bDgelsd = true;
    if (bDgelsd)
    {
        dgelsd_ (&nrow, &ncol, &nrhs, a[0], &lda, rhs->data(), &ldb, s.data(),
                 &rcond, &rank, &wkopt, &lwork,
                 iwork.data(), &info );
    }
    else
    {
        dgels_ ("No transpose", &nrow, &ncol, &nrhs, a[0], &lda, rhs->data(), &ldb,
                &wkopt, &lwork, &info );
    }
    lwork = (int)wkopt;
    std::vector<double> work;
    work.resize(lwork);
    GMX_RELEASE_ASSERT(rhs->size() - nrow == 0, gmx::formatString("rhs.size = %ld nrow = %d", rhs->size(), nrow).c_str());
    /* Solve the equations A*X = B */
    if (bDgelsd)
    {
        dgelsd_ (&nrow, &ncol, &nrhs, a[0], &lda, rhs->data(), &ldb, s.data(),
                 &rcond, &rank, work.data(), &lwork,
                 iwork.data(), &info );
    }
    else
    {
        dgels_ ("No transpose", &nrow, &ncol, &nrhs, a[0], &lda, rhs->data(), &ldb,
                work.data(), &lwork, &info );
    }
    /* Check for convergence */
    if (info > 0)
    {
        fprintf(stderr, "The algorithm computing SVD failed to converge and the least squares solution could not be computed. Info = %d, nrow = %d, ncol = %d", info, nrow, ncol);
    }
    else
    {
        GMX_RELEASE_ASSERT(x->size() - ncol == 0, 
                           gmx::formatString("x has the wrong size (%zu) instead of %d", x->size(), ncol).c_str());
        for (int i = 0; i < ncol; i++)
        {
            (*x)[i] = (*rhs)[i];
        }
    }
    return info;
}
#endif

static void tensor2matrix(tensor c, double **a)
{
    for (int i = 0; i < DIM; i++)
    {
        for (int j = 0; j < DIM; j++)
        {
            a[i][j] = c[i][j];
        }
    }

}


static void rowwise2tensor(double **a, tensor c)
{
    c[XX][XX] = a[XX][XX];
    c[XX][YY] = a[XX][YY];
    c[XX][ZZ] = a[XX][ZZ];
    c[YY][XX] = a[YY][XX];
    c[YY][YY] = a[YY][YY];
    c[YY][ZZ] = a[YY][ZZ];
    c[ZZ][XX] = a[ZZ][XX];
    c[ZZ][YY] = a[ZZ][YY];
    c[ZZ][ZZ] = a[ZZ][ZZ];
}

static void columnwise2tensor(double **a, tensor c)
{
    c[XX][XX] = a[XX][ZZ];
    c[XX][YY] = a[XX][YY];
    c[XX][ZZ] = a[XX][XX];
    c[YY][XX] = a[ZZ][ZZ];
    c[YY][YY] = a[ZZ][YY];
    c[YY][ZZ] = a[ZZ][XX];
    c[ZZ][XX] = a[YY][ZZ];
    c[ZZ][YY] = a[YY][YY];
    c[ZZ][ZZ] = a[YY][XX];
}

#ifdef OLD
static void SVD(int*    m,    int*     n,     double** a,
                int*    lda,  double** s,     double** u,
                int*    ldu,  double** vt,    int*     ldvt)
{
    int                 lwork;
    int                 info;
    double              wkopt;
    std::vector<double> work;

    lwork = -1;
    dgesvd_("All", "All", m, n, a[0], lda, s[0], u[0], ldu, vt[0], ldvt, &wkopt, &lwork, &info);
    lwork = (int)wkopt;
    work.resize(lwork);
    dgesvd_("All", "All", m, n, a[0], lda, s[0], u[0], ldu, vt[0], ldvt, work.data(), &lwork, &info);
    if (info > 0)
    {
        char buf[256];
        snprintf(buf, sizeof(buf), "SVD algorithm failed to converge. Info = %d.", info);
        GMX_THROW(gmx::InvalidInputError(buf));
    }
}

void kabsch_rotation(tensor P, tensor Q, tensor rotated_P)
{
    int       m      = DIM;
    int       n      = DIM;
    int       lda    = m;
    int       ldu    = m;
    int       ldvt   = n;

    double  **A      = alloc_matrix(m, n);
    double  **s      = alloc_matrix(m, n);
    double  **u      = alloc_matrix(m, n);
    double  **vt     = alloc_matrix(m, n);

    tensor    C, U, UT, V, VT;
    tensor    Pt, R;


    /*transpose of p (Pt)*/
    transpose(P, Pt);

    /*Covariance matrix (A)*/
    mmul(Pt, Q, C);
    tensor2matrix(C, A);

    /* SVD returns:
       u  -> Left singular vectors (stored columnwise)
       s  -> Singular values
       vt -> Transpose of the Right singular vectors (stored rowwise)
     */
    SVD(&m, &n, A, &lda, s, u, &ldu, vt, &ldvt);

    columnwise2tensor(u, U);
    rowwise2tensor(vt, VT);

    transpose(U, UT);
    transpose(VT, V);

    /*Check for correction*/
    if ((det(UT) * det(V)) < 0)
    {
        UT[ZZ][ZZ] = -UT[ZZ][ZZ];
    }

    /*Rotation matrix (R)*/
    mmul(UT, V, R);

    /*Rotate*/
    mmul(P, R, rotated_P);

    free_matrix(A);
    free_matrix(s);
    free_matrix(u);
    free_matrix(vt);
}
#endif

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

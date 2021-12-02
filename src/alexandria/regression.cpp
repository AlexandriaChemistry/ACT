/*
 * This source file is part of the Alexandria Chemistry Toolkit.
 *
 * Copyright (C) 2014-2020
 *
 * Developers:
 *             Mohammad Mehdi Ghahremanpour,
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

#include <vector>

#include "gromacs/linearalgebra/matrix.h"
#include "gromacs/math/vec.h"
#include "gromacs/utility/exceptions.h"
#include "gromacs/utility/stringutil.h"

extern "C"
{
void dgelsd_(int* m, int* n, int* nrhs, double* a, int* lda,
             double* b, int* ldb, double* s, double* rcond, int* rank,
             double* work, int* lwork, int* iwork, int* info );

void dgels_(const char* trans, int* m, int* n, int* nrhs, double* a, int* lda,
            double* b, int* ldb, double* work, int* lwork, int* info );

void dgesvd_(const char* jobu, const char* jobvt, int* m, int* n, double* a,
             int* lda, double* s, double* u, int* ldu, double* vt, int* ldvt,
             double* work, int* lwork, int* info );

}

/*! \brief Linear regression with many columns
 *
 * \param[in]  rhs  The right hand side of the equation
 * \param[in]  a    The matrix
 * \param[out] x    The solution of the problem
 * \return 0 if all is fine, the problematic row number otherwise
 */
static int multi_regression2(std::vector<double> *rhs,
                             double              **a, 
                             std::vector<double> *x)
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
        fprintf(stderr, "The algorithm computing SVD failed to converge and the least squares solution could not be computed. Info = %d.", info);
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

MatrixWrapper::MatrixWrapper(int ncolumn, int nrow)
{
    a_       = alloc_matrix(ncolumn, nrow);
}

MatrixWrapper::~MatrixWrapper()
{
    free_matrix(a_);
}

int MatrixWrapper::solve(std::vector<double> rhs, std::vector<double> *solution)
{
    return multi_regression2(&rhs, a_, solution);
}

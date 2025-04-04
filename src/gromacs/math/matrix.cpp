/*
 * This file is part of the GROMACS molecular simulation package.
 *
 * Copyright (c) 1991-2000, University of Groningen, The Netherlands.
 * Copyright (c) 2001-2008, The GROMACS development team.
 * Copyright (c) 2012,2013,2014,2015,2017, by the GROMACS development team, led by
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
#include "actpre.h"

#include "matrix.h"

#include "config.h"

#include <stdio.h>

#include "gromacs/utility/smalloc.h"

double **alloc_matrix(int n, int m)
{
    double **ptr;
    int      i;

    if (n <= 0 || m <= 0)
    {
        return nullptr;
    }
    /* There's always time for more pointer arithmetic! */
    /* This is necessary in order to be able to work with LAPACK */
    snew(ptr, n);
    snew(ptr[0], n*m);
    for (i = 1; (i < n); i++)
    {
        ptr[i] = ptr[i-1]+m;
    }
    return ptr;
}

void free_matrix(double **a)
{
    sfree(a[0]);
    sfree(a);
}

#define DEBUG_MATRIX 1
void matrix_multiply(FILE *fp, int n, int m, double **x, double **y, double **z)
{
    int i, j, k;

#ifdef DEBUG_MATRIX
    if (fp)
    {
        fprintf(fp, "Multiplying %d x %d matrix with a %d x %d matrix\n",
                n, m, m, n);
    }
    if (fp)
    {
        for (i = 0; (i < n); i++)
        {
            for (j = 0; (j < m); j++)
            {
                fprintf(fp, " %7g", x[i][j]);
            }
            fprintf(fp, "\n");
        }
    }
#endif
    for (i = 0; (i < m); i++)
    {
        for (j = 0; (j < m); j++)
        {
            z[i][j] = 0;
            for (k = 0; (k < n); k++)
            {
                z[i][j] += x[k][i]*y[j][k];
            }
        }
    }
}


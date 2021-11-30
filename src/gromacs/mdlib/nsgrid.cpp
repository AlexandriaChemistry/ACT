/*
 * This file is part of the GROMACS molecular simulation package.
 *
 * Copyright (c) 1991-2000, University of Groningen, The Netherlands.
 * Copyright (c) 2001-2004, The GROMACS development team.
 * Copyright (c) 2013,2014,2015,2017,2018, by the GROMACS development team, led by
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
/* This file is completely threadsafe - keep it that way! */
#include "actpre.h"

#include "nsgrid.h"

#include <cmath>
#include <cstdio>
#include <cstdlib>

#include <algorithm>

#include "gromacs/fileio/pdbio.h"
#include "gromacs/gmxlib/network.h"
#include "gromacs/math/vec.h"
#include "gromacs/mdtypes/forcerec.h"
#include "gromacs/pbcutil/pbc.h"
#include "gromacs/utility/fatalerror.h"
#include "gromacs/utility/futil.h"
#include "gromacs/utility/smalloc.h"

/***********************************
 *         Grid Routines
 ***********************************/

static const char *range_warn =
    "Explanation: During neighborsearching, we assign each particle to a grid\n"
    "based on its coordinates. If your system contains collisions or parameter\n"
    "errors that give particles very high velocities you might end up with some\n"
    "coordinates being +-Infinity or NaN (not-a-number). Obviously, we cannot\n"
    "put these on a grid, so this is usually where we detect those errors.\n"
    "Make sure your system is properly energy-minimized and that the potential\n"
    "energy seems reasonable before trying again.";

static void calc_x_av_stddev(int n, rvec *x, rvec av, rvec stddev)
{
    dvec s1, s2;
    int  i, d;

    clear_dvec(s1);
    clear_dvec(s2);

    for (i = 0; i < n; i++)
    {
        for (d = 0; d < DIM; d++)
        {
            s1[d] += x[i][d];
            s2[d] += x[i][d]*x[i][d];
        }
    }

    dsvmul(1.0/n, s1, s1);
    dsvmul(1.0/n, s2, s2);

    for (d = 0; d < DIM; d++)
    {
        av[d]     = s1[d];
        stddev[d] = std::sqrt(s2[d] - s1[d]*s1[d]);
    }
}

static void get_nsgrid_boundaries_vac(real av, real stddev,
                                      real *bound0, real *bound1,
                                      real *bdens0, real *bdens1)
{
    /* Set the grid to 2 times the standard deviation of
     * the charge group centers in both directions.
     * For a uniform bounded distribution the width is sqrt(3)*stddev,
     * so all charge groups fall within the width.
     * For a sphere stddev is r/sqrt(5): 99.2% falls within the width.
     * For a Gaussian distribution 98% fall within the width.
     */
    *bound0 = av - NSGRID_STDDEV_FAC*stddev;
    *bound1 = av + NSGRID_STDDEV_FAC*stddev;

    *bdens0 = av - GRID_STDDEV_FAC*stddev;
    *bdens1 = av + GRID_STDDEV_FAC*stddev;
}

void get_nsgrid_boundaries(int nboundeddim, matrix box,
                           rvec *gr0, rvec *gr1,
                           int ncg, rvec *cgcm,
                           rvec grid_x0, rvec grid_x1,
                           real *grid_density)
{
    rvec av, stddev;
    real vol, bdens0, bdens1;
    int  d;

    if (nboundeddim < DIM)
    {
        calc_x_av_stddev(ncg, cgcm, av, stddev);
    }

    vol = 1;
    for (d = 0; d < DIM; d++)
    {
        if (d < nboundeddim)
        {
            grid_x0[d] = (gr0 != nullptr ? (*gr0)[d] : 0);
            grid_x1[d] = (gr1 != nullptr ? (*gr1)[d] : box[d][d]);
            vol       *= (grid_x1[d] - grid_x0[d]);
        }
        else
        {
            get_nsgrid_boundaries_vac(av[d], stddev[d],
                                      &grid_x0[d], &grid_x1[d],
                                      &bdens0, &bdens1);
            vol *= (bdens1 - bdens0);
        }

        if (debug)
        {
            fprintf(debug, "Set grid boundaries dim %d: %f %f\n",
                    d, grid_x0[d], grid_x1[d]);
        }
    }

    *grid_density = ncg/vol;
}

static void set_grid_sizes(rvec izones_x0, rvec izones_x1, real rlist,
                           t_grid *grid,
                           real grid_density)
{
    int      i;
    gmx_bool bDD, bDDRect;
    rvec     izones_size;
    real     inv_r_ideal, size;

    for (i = 0; (i < DIM); i++)
    {
        if (debug)
        {
            fprintf(debug,
                    "set_grid_sizes, i-zone bounds for dim %d: %6.3f %6.3f\n",
                    i, izones_x0[i], izones_x1[i]);
        }
        izones_size[i] = izones_x1[i] - izones_x0[i];
    }

    /* Use the ideal number of cg's per cell to set the ideal cell size */
    inv_r_ideal = std::cbrt(grid_density/grid->ncg_ideal);
    if (rlist > 0 && inv_r_ideal*rlist < 1)
    {
        inv_r_ideal = 1/rlist;
    }
    if (debug)
    {
        fprintf(debug, "CG density %f ideal ns cell size %f\n",
                grid_density, 1/inv_r_ideal);
    }

    clear_rvec(grid->cell_offset);
    for (i = 0; (i < DIM); i++)
    {
        /* Initial settings, for DD might change below */
        grid->cell_offset[i] = izones_x0[i];
        size                 = izones_size[i];

        bDD = false;
        if (!bDD)
        {
            bDDRect = FALSE;
        }
        if (!bDDRect)
        {
            /* No DD or the box is triclinic is this direction:
             * we will use the normal grid ns that checks all cells
             * that are within cut-off distance of the i-particle.
             */
            grid->n[i] = gmx::roundToInt(size*inv_r_ideal);
            if (grid->n[i] < 2)
            {
                grid->n[i] = 2;
            }
            grid->cell_size[i] = size/grid->n[i];
            grid->ncpddc[i]    = 0;
        }
        if (debug)
        {
            fprintf(debug, "grid dim %d size %d x %f: %f - %f\n",
                    i, grid->n[i], grid->cell_size[i],
                    grid->cell_offset[i],
                    grid->cell_offset[i]+grid->n[i]*grid->cell_size[i]);
        }
    }

    if (debug)
    {
        fprintf(debug, "CG ncg ideal %d, actual density %.1f\n",
                grid->ncg_ideal, grid_density*grid->cell_size[XX]*grid->cell_size[YY]*grid->cell_size[ZZ]);
    }
}

t_grid *init_grid(FILE *fplog, t_forcerec *fr)
{
    char   *ptr;
    t_grid *grid;

    snew(grid, 1);

    grid->npbcdim = ePBC2npbcdim(fr->ePBC);

    if (fr->ePBC == epbcXY && fr->nwall == 2)
    {
        grid->nboundeddim = 3;
    }
    else
    {
        grid->nboundeddim = grid->npbcdim;
    }

    if (debug)
    {
        fprintf(debug, "The coordinates are bounded in %d dimensions\n",
                grid->nboundeddim);
    }

    /* The ideal number of cg's per ns grid cell seems to be 10 */
    grid->ncg_ideal = 10;
    ptr             = getenv("GMX_NSCELL_NCG");
    if (ptr)
    {
        sscanf(ptr, "%d", &grid->ncg_ideal);
        if (fplog)
        {
            fprintf(fplog, "Set ncg_ideal to %d\n", grid->ncg_ideal);
        }
        if (grid->ncg_ideal <= 0)
        {
            gmx_fatal(FARGS, "The number of cg's per cell should be > 0");
        }
    }
    if (debug)
    {
        fprintf(debug, "Set ncg_ideal to %d\n", grid->ncg_ideal);
    }

    return grid;
}

void done_grid(t_grid *grid)
{
    if (grid == nullptr)
    {
        return;
    }
    grid->nr      = 0;
    clear_ivec(grid->n);
    grid->ncells  = 0;
    sfree(grid->cell_index);
    sfree(grid->a);
    sfree(grid->index);
    sfree(grid->nra);
    grid->cells_nalloc = 0;
    sfree(grid->dcx2);
    sfree(grid->dcy2);
    sfree(grid->dcz2);
    grid->dc_nalloc = 0;

    if (debug)
    {
        fprintf(debug, "Successfully freed memory for grid pointers.");
    }
    sfree(grid);
}

int xyz2ci_(int nry, int nrz, int x, int y, int z)
/* Return the cell index */
{
    return (nry*nrz*x+nrz*y+z);
}

void ci2xyz(t_grid *grid, int i, int *x, int *y, int *z)
/* Return x,y and z from the cell index */
{
    int ci;

    range_check_mesg(i, 0, grid->nr, range_warn);

    ci  = grid->cell_index[i];
    *x  = ci / (grid->n[YY]*grid->n[ZZ]);
    ci -= (*x)*grid->n[YY]*grid->n[ZZ];
    *y  = ci / grid->n[ZZ];
    ci -= (*y)*grid->n[ZZ];
    *z  = ci;
}

static int ci_not_used(const ivec n)
{
    /* Return an improbable value */
    return xyz2ci(n[YY], n[ZZ], 3*n[XX], 3*n[YY], 3*n[ZZ]);
}

static void set_grid_ncg(t_grid *grid, int ncg)
{
    int nr_old, i;

    grid->nr = ncg;
    if (grid->nr+1 > grid->nr_alloc)
    {
        nr_old         = grid->nr_alloc;
        grid->nr_alloc = over_alloc_dd(grid->nr) + 1;
        srenew(grid->cell_index, grid->nr_alloc);
        for (i = nr_old; i < grid->nr_alloc; i++)
        {
            grid->cell_index[i] = 0;
        }
        srenew(grid->a, grid->nr_alloc);
    }
}

void grid_first(FILE *fplog, t_grid *grid,
                rvec izones_x0, rvec izones_x1,
                real rlist, real grid_density)
{
    int    i, m;

    set_grid_sizes(izones_x0, izones_x1, rlist, grid, grid_density);

    grid->ncells = grid->n[XX]*grid->n[YY]*grid->n[ZZ];

    if (grid->ncells+1 > grid->cells_nalloc)
    {
        /* Allocate double the size so we have some headroom */
        grid->cells_nalloc = 2*grid->ncells;
        srenew(grid->nra,  grid->cells_nalloc+1);
        srenew(grid->index, grid->cells_nalloc+1);

        if (fplog)
        {
            fprintf(fplog, "Grid: %d x %d x %d cells\n",
                    grid->n[XX], grid->n[YY], grid->n[ZZ]);
        }
    }

    m = std::max(grid->n[XX], std::max(grid->n[YY], grid->n[ZZ]));
    if (m > grid->dc_nalloc)
    {
        /* Allocate with double the initial size for box scaling */
        grid->dc_nalloc = 2*m;
        srenew(grid->dcx2, grid->dc_nalloc);
        srenew(grid->dcy2, grid->dc_nalloc);
        srenew(grid->dcz2, grid->dc_nalloc);
    }

    grid->nr = 0;
    for (i = 0; (i < grid->ncells); i++)
    {
        grid->nra[i] = 0;
    }
}

static void calc_bor(int cg0, int cg1, int ncg, int CG0[2], int CG1[2])
{
    if (cg1 > ncg)
    {
        CG0[0] = cg0;
        CG1[0] = ncg;
        CG0[1] = 0;
        CG1[1] = cg1-ncg;
    }
    else
    {
        CG0[0] = cg0;
        CG1[0] = cg1;
        CG0[1] = 0;
        CG1[1] = 0;
    }
    if (debug)
    {
        int m;

        fprintf(debug, "calc_bor: cg0=%d, cg1=%d, ncg=%d\n", cg0, cg1, ncg);
        for (m = 0; (m < 2); m++)
        {
            fprintf(debug, "CG0[%d]=%d, CG1[%d]=%d\n", m, CG0[m], m, CG1[m]);
        }
    }

}

void calc_elemnr(t_grid *grid, int cg0, int cg1, int ncg)
{
    int     CG0[2], CG1[2];
    int    *cell_index = grid->cell_index;
    int    *nra        = grid->nra;
    int     i, m, ncells;
    int     ci, not_used;

    ncells = grid->ncells;
    if (ncells <= 0)
    {
        gmx_fatal(FARGS, "Number of grid cells is zero. Probably the system and box collapsed.\n");
    }

    not_used = ci_not_used(grid->n);

    calc_bor(cg0, cg1, ncg, CG0, CG1);
    for (m = 0; (m < 2); m++)
    {
        for (i = CG0[m]; (i < CG1[m]); i++)
        {
            ci = cell_index[i];
            if (ci != not_used)
            {
                range_check_mesg(ci, 0, ncells, range_warn);
                nra[ci]++;
            }
        }
    }
}

void calc_ptrs(t_grid *grid)
{
    int *index = grid->index;
    int *nra   = grid->nra;
    int  ix, iy, iz, ci, nr;
    int  nnra, ncells;

    ncells = grid->ncells;
    if (ncells <= 0)
    {
        gmx_fatal(FARGS, "Number of grid cells is zero. Probably the system and box collapsed.\n");
    }

    ci = nr = 0;
    for (ix = 0; (ix < grid->n[XX]); ix++)
    {
        for (iy = 0; (iy < grid->n[YY]); iy++)
        {
            for (iz = 0; (iz < grid->n[ZZ]); iz++, ci++)
            {
                range_check_mesg(ci, 0, ncells, range_warn);
                index[ci] = nr;
                nnra      = nra[ci];
                nr       += nnra;
                nra[ci]   = 0;
            }
        }
    }
}

void grid_last(t_grid *grid, int cg0, int cg1, int ncg)
{
    int     CG0[2], CG1[2];
    int     i, m;
    int     ci, not_used, ind, ncells;
    int    *cell_index = grid->cell_index;
    int    *nra        = grid->nra;
    int    *index      = grid->index;
    int    *a          = grid->a;

    ncells = grid->ncells;
    if (ncells <= 0)
    {
        gmx_fatal(FARGS, "Number of grid cells is zero. Probably the system and box collapsed.\n");
    }

    not_used = ci_not_used(grid->n);

    calc_bor(cg0, cg1, ncg, CG0, CG1);
    for (m = 0; (m < 2); m++)
    {
        for (i = CG0[m]; (i < CG1[m]); i++)
        {
            ci     = cell_index[i];
            if (ci != not_used)
            {
                range_check_mesg(ci, 0, ncells, range_warn);
                ind    = index[ci]+nra[ci]++;
                range_check_mesg(ind, 0, grid->nr, range_warn);
                a[ind] = i;
            }
        }
    }
}

void fill_grid(gmx_domdec_zones_t *dd_zones,
               t_grid *grid, int ncg_tot,
               int cg0, int cg1, rvec cg_cm[])
{
    int       *cell_index;
    int        nry, nrz;
    rvec       n_box, offset;
    int        cg, d;
    ivec       ind;

    if (cg0 == -1)
    {
        /* We have already filled the grid up to grid->ncg,
         * continue from there.
         */
        cg0 = grid->nr;
    }

    set_grid_ncg(grid, ncg_tot);

    cell_index = grid->cell_index;

    /* Initiate cell borders */
    nry = grid->n[YY];
    nrz = grid->n[ZZ];
    for (d = 0; d < DIM; d++)
    {
        if (grid->cell_size[d] > 0)
        {
            n_box[d] = 1/grid->cell_size[d];
        }
        else
        {
            n_box[d] = 0;
        }
    }
    copy_rvec(grid->cell_offset, offset);

    if (debug)
    {
        fprintf(debug, "Filling grid from %d to %d\n", cg0, cg1);
    }

    if (dd_zones == nullptr)
    {
        for (cg = cg0; cg < cg1; cg++)
        {
            for (d = 0; d < DIM; d++)
            {
                ind[d] = static_cast<int>((cg_cm[cg][d] - offset[d])*n_box[d]);
                /* With pbc we should be done here.
                 * Without pbc cg's outside the grid
                 * should be assigned to the closest grid cell.
                 */
                if (ind[d] < 0)
                {
                    ind[d] = 0;
                }
                else if (ind[d] >= grid->n[d])
                {
                    ind[d] = grid->n[d] - 1;
                }
            }
            cell_index[cg] = xyz2ci(nry, nrz, ind[XX], ind[YY], ind[ZZ]);
        }
    }

}

void check_grid(t_grid *grid)
{
    int ix, iy, iz, ci, cci, nra;

    if (grid->ncells <= 0)
    {
        gmx_fatal(FARGS, "Number of grid cells is zero. Probably the system and box collapsed.\n");
    }

    ci  = 0;
    cci = 0;
    for (ix = 0; (ix < grid->n[XX]); ix++)
    {
        for (iy = 0; (iy < grid->n[YY]); iy++)
        {
            for (iz = 0; (iz < grid->n[ZZ]); iz++, ci++)
            {
                if (ci > 0)
                {
                    nra = grid->index[ci]-grid->index[cci];
                    if (nra != grid->nra[cci])
                    {
                        gmx_fatal(FARGS, "nra=%d, grid->nra=%d, cci=%d",
                                  nra, grid->nra[cci], cci);
                    }
                }
                cci = xyz2ci(grid->n[YY], grid->n[ZZ], ix, iy, iz);
                range_check_mesg(cci, 0, grid->ncells, range_warn);

                if (cci != ci)
                {
                    gmx_fatal(FARGS, "ci = %d, cci = %d", ci, cci);
                }
            }
        }
    }
}

void print_grid(FILE *log, t_grid *grid)
{
    int i, nra, index;
    int ix, iy, iz, ci;

    fprintf(log, "nr:        %d\n", grid->nr);
    fprintf(log, "nrx:       %d\n", grid->n[XX]);
    fprintf(log, "nry:       %d\n", grid->n[YY]);
    fprintf(log, "nrz:       %d\n", grid->n[ZZ]);
    fprintf(log, "ncg_ideal: %d\n", grid->ncg_ideal);
    fprintf(log, "    i  cell_index\n");
    for (i = 0; (i < grid->nr); i++)
    {
        fprintf(log, "%5d  %5d\n", i, grid->cell_index[i]);
    }
    fprintf(log, "cells\n");
    fprintf(log, " ix iy iz   nr  index  cgs...\n");
    ci = 0;
    for (ix = 0; (ix < grid->n[XX]); ix++)
    {
        for (iy = 0; (iy < grid->n[YY]); iy++)
        {
            for (iz = 0; (iz < grid->n[ZZ]); iz++, ci++)
            {
                index = grid->index[ci];
                nra   = grid->nra[ci];
                fprintf(log, "%3d%3d%3d%5d%5d", ix, iy, iz, nra, index);
                for (i = 0; (i < nra); i++)
                {
                    fprintf(log, "%5d", grid->a[index+i]);
                }
                fprintf(log, "\n");
            }
        }
    }
    fflush(log);
}

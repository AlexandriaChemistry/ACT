/*
 * This file is part of the GROMACS molecular simulation package.
 *
 * Copyright (c) 1991-2000, University of Groningen, The Netherlands.
 * Copyright (c) 2001-2004, The GROMACS development team.
 * Copyright (c) 2012,2014,2015,2017,2018, by the GROMACS development team, led by
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

#include "statistics.h"

#include <cmath>

#include <map>

#include "gromacs/math/functions.h"
#include "gromacs/math/vec.h"
#include "gromacs/utility/fatalerror.h"
#include "gromacs/utility/real.h"

static int gmx_dnint(double x)
{
    return gmx::roundToInt(x);
}

eStats gmx_stats::add_point_ydy(double y, double dy)
{
    return add_point(x.size(), y, 0, dy);
}

eStats gmx_stats::add_point(double xx, double yy,
                            double dxx, double dyy)
{
    x.push_back(xx);
    y.push_back(yy);
    dx.push_back(dxx);
    dy.push_back(dyy);
    computed = false;

    return eStats::OK;
}

eStats gmx_stats::get_point(real *xx, real *yy,
                            real *dxx, real *dyy, real level)
{
    eStats     ok;
    int        outlier;
    real       rmsd, r;

    if ((ok = gmx_stats::get_rmsd(&rmsd)) != eStats::OK)
    {
        return ok;
    }
    outlier = 0;
    while ((outlier == 0) && (np_c < x.size()))
    {
        r       = std::abs(x[np_c] - y[np_c]);
        outlier = static_cast<int>(r > rmsd*level);
        if (outlier)
        {
            if (nullptr != xx)
            {
                *xx  = x[np_c];
            }
            if (nullptr != yy)
            {
                *yy  = y[np_c];
            }
            if (nullptr != dxx)
            {
                *dxx = dx[np_c];
            }
            if (nullptr != dyy)
            {
                *dyy = dy[np_c];
            }
        }
        np_c++;

        if (outlier)
        {
            return eStats::OK;
        }
    }

    np_c = 0;

    return eStats::NO_POINTS;
}

eStats gmx_stats::add_points(int n, real *xx, real *yy,
                             real *dxx, real *dyy)
{
    for (int i = 0; (i < n); i++)
    {
        eStats ok = add_point(xx[i], yy[i],
                              (nullptr != dxx) ? dxx[i] : 0,
                              (nullptr != dyy) ? dyy[i] : 0);
        if (ok != eStats::OK)
        {
            return ok;
        }
    }
    return eStats::OK;
}

eStats gmx_stats::compute(int weight)
{
    double yy, yx, xx, sx, sy, chi2, chi2aa, d2;
    double ssxx, ssyy, ssxy;
    double w, wtot, yx_nw, sy_nw, sx_nw, yy_nw, xx_nw, dx2, dy2;

    int    N = x.size();

    if (!computed)
    {
        if (N < 1)
        {
            return eStats::NO_POINTS;
        }

        xx   = xx_nw = 0;
        yy   = yy_nw = 0;
        yx   = yx_nw = 0;
        sx   = sx_nw = 0;
        sy   = sy_nw = 0;
        wtot = 0;
        d2   = 0;
        mae  = 0, mse = 0;
        for (int i = 0; (i < N); i++)
        {
            double dd = y[i]-x[i];
            d2 += gmx::square(dd);
            
            mae += fabs(dd);
            mse += dd;
            if ((dy[i] != 0.0) && (weight == elsqWEIGHT_Y))
            {
                w = 1/gmx::square(dy[i]);
            }
            else
            {
                w = 1;
            }

            wtot  += w;

            xx    += w*gmx::square(x[i]);
            xx_nw += gmx::square(x[i]);

            yy    += w*gmx::square(y[i]);
            yy_nw += gmx::square(y[i]);

            yx    += w*y[i]*x[i];
            yx_nw += y[i]*x[i];

            sx    += w*x[i];
            sx_nw += x[i];

            sy    += w*y[i];
            sy_nw += y[i];
        }

        /* Compute average, sigma and error */
        mae        = mae/N;
        mse        = mse/N; 
        aver       = sy_nw/N;

        double dd = yy_nw/N - gmx::square(sy_nw/N);
        if (dd > 0)
        {
            sigma_aver = std::sqrt(dd);
        }
        else if (N == 1)
        {
            sigma_aver = dy[0];
        }
        else
        {
            sigma_aver = 0;
        }
        error      = sigma_aver/std::sqrt(static_cast<double>(N));

        /* Compute RMSD between x and y */
        rmsd = std::sqrt(d2/N);

        /* Correlation coefficient for data */
        yx_nw       /= N;
        xx_nw       /= N;
        yy_nw       /= N;
        sx_nw       /= N;
        sy_nw       /= N;
        ssxx         = N*(xx_nw - gmx::square(sx_nw));
        ssyy         = N*(yy_nw - gmx::square(sy_nw));
        ssxy         = N*(yx_nw - (sx_nw*sy_nw));
        Rdata = std::sqrt(gmx::square(ssxy)/(ssxx*ssyy));

        /* Compute straight line through datapoints, either with intercept
           zero (result in aa) or with intercept variable (results in a
           and b) */
        yx = yx/wtot;
        xx = xx/wtot;
        sx = sx/wtot;
        sy = sy/wtot;

        aa = (yx/xx);
        a  = (yx-sx*sy)/(xx-sx*sx);
        b  = (sy)-(a)*(sx);

        /* Compute chi2, deviation from a line y = ax+b. Also compute
           chi2aa which returns the deviation from a line y = ax. */
        chi2   = 0;
        chi2aa = 0;
        for (int i = 0; (i < N); i++)
        {
            real ddy = 1;
            if (dy[i] > 0)
            {
                ddy = dy[i];
            }
            chi2aa += gmx::square((y[i]-(aa*x[i]))/ddy);
            chi2   += gmx::square((y[i]-(a*x[i]+b))/ddy);
        }
        if (N > 2)
        {
            chi2   = std::sqrt(chi2/(N-2));
            chi2aa = std::sqrt(chi2aa/(N-2));

            /* Look up equations! */
            dx2            = (xx-sx*sx);
            dy2            = (yy-sy*sy);
            sigma_a = std::sqrt(chi2/((N-2)*dx2));
            sigma_b = sigma_a*std::sqrt(xx);
            Rfit    = std::abs(ssxy)/std::sqrt(ssxx*ssyy);
            Rfitaa  = Rfit; //aa*std::sqrt(dx2/dy2);
        }
        else
        {
            chi2    = 0;
            chi2aa  = 0;
            sigma_a = 0;
            sigma_b = 0;
            Rfit    = 0;
            Rfitaa  = 0;
        }

        computed = true;
    }

    return eStats::OK;
}

eStats gmx_stats::get_ab(int weight,
                         real *aa, real *bb, real *daa, real *dbb,
                         real *chi2a, real *Rfita)
{
    eStats     ok;

    if ((ok = compute(weight)) != eStats::OK)
    {
        return ok;
    }
    if (nullptr != aa)
    {
        *aa    = a;
    }
    if (nullptr != bb)
    {
        *bb    = b;
    }
    if (nullptr != daa)
    {
        *daa   = sigma_a;
    }
    if (nullptr != dbb)
    {
        *dbb   = sigma_b;
    }
    if (nullptr != chi2a)
    {
        *chi2a = chi2;
    }
    if (nullptr != Rfita)
    {
        *Rfita = Rfit;
    }

    return eStats::OK;
}

eStats gmx_stats::get_a(int weight, real *aa, real *daa,
                        real *chi2a, real *Rfita)
{
    eStats     ok;

    if ((ok = compute(weight)) != eStats::OK)
    {
        return ok;
    }
    if (nullptr != aa)
    {
        *aa    = a;
    }
    if (nullptr != daa)
    {
        *daa   = sigma_a;
    }
    if (nullptr != chi2a)
    {
        *chi2a = chi2;
    }
    if (nullptr != Rfita)
    {
        *Rfita = Rfit;
    }

    return eStats::OK;
}

eStats gmx_stats::get_average(real *avera)
{
    eStats     ok;

    if ((ok = compute(elsqWEIGHT_NONE)) != eStats::OK)
    {
        return ok;
    }

    *avera = aver;

    return eStats::OK;
}

eStats gmx_stats::get_mse_mae(real *msea, real *maea)
{
    eStats     ok;

    if ((ok = compute(elsqWEIGHT_NONE)) != eStats::OK)
    {
        return ok;
    }

    if (NULL != msea)
    {
        *msea = mse;
    }
    if (NULL != maea)
    {
        *maea = mae;
    }

    return eStats::OK;
}

eStats gmx_stats::get_ase(real *avera, real *sigmaa, real *errora)
{
    eStats     ok;

    if ((ok = compute(elsqWEIGHT_NONE)) != eStats::OK)
    {
        return ok;
    }

    if (nullptr != avera)
    {
        *avera  = aver;
    }
    if (nullptr != sigmaa)
    {
        *sigmaa = sigma_aver;
    }
    if (nullptr != errora)
    {
        *errora = error;
    }

    return eStats::OK;
}

eStats gmx_stats::get_sigma(real *sigmaa)
{
    eStats     ok;

    if ((ok = compute(elsqWEIGHT_NONE)) != eStats::OK)
    {
        return ok;
    }

    *sigmaa = sigma_aver;

    return eStats::OK;
}

eStats gmx_stats::get_error(real *errors)
{
    eStats     ok;

    if ((ok = compute(elsqWEIGHT_NONE)) != eStats::OK)
    {
        return ok;
    }

    *errors = error;

    return eStats::OK;
}

eStats gmx_stats::get_corr_coeff(real *Ra)
{
    eStats     ok;

    if ((ok = gmx_stats::compute(elsqWEIGHT_NONE)) != eStats::OK)
    {
        return ok;
    }

    *Ra = Rdata;

    return eStats::OK;
}

eStats gmx_stats::get_rmsd(real *rmsda)
{
    eStats        ok;

    if ((ok = compute(elsqWEIGHT_NONE)) != eStats::OK)
    {
        return ok;
    }

    *rmsda = rmsd;

    return eStats::OK;
}

eStats gmx_stats::dump_xy(FILE *fp)
{
    for (size_t i = 0; (i < x.size()); i++)
    {
        fprintf(fp, "%12g  %12g  %12g  %12g\n", x[i], y[i],
                dx[i], dy[i]);
    }

    return eStats::OK;
}

eStats gmx_stats::remove_outliers(double level)
{
    int        iter  = 1, done = 0;
    eStats     ok;
    real       rmsd, r;

    while ((x.size() >= 10) && !done)
    {
        if ((ok = get_rmsd(&rmsd)) != eStats::OK)
        {
            return ok;
        }
        done = 1;
        size_t i;
        for (i = 0; (i < x.size()); )
        {
            r = std::abs(x[i]-y[i]);
            if (r > level*rmsd)
            {
                fprintf(stderr, "Removing outlier, iter = %d, rmsd = %g, x = %g, y = %g\n",
                        iter, rmsd, x[i], y[i]);
                if (i < x.size()-1)
                {
                    x[i]  = x[x.size()-1];
                    y[i]  = y[x.size()-1];
                    dx[i] = dx[x.size()-1];
                    dy[i] = dy[x.size()-1];
                }
                done = 0;
            }
            else
            {
                i++;
            }
        }
        x.resize(i);
        y.resize(i);
        dx.resize(i);
        dy.resize(i);
        iter++;
    }

    return eStats::OK;
}

eStats gmx_stats::make_histogram(real binwidth, int *nb,
                                eHisto ehisto, int normalized,
                                std::vector<double> *xx,
                                std::vector<double> *yy)
{
    int        index = 0, nbins = *nb;
    std::vector<int> nindex;
    double     minx, maxx, maxy, miny, delta, dd, minh;

    if (((binwidth <= 0) && (nbins <= 0)) ||
        ((binwidth > 0) && (nbins > 0)))
    {
        return eStats::INVALID_INPUT;
    }
    if (x.size() <= 1)
    {
        return eStats::NO_POINTS;
    }
    minx = maxx = x[0];
    miny = maxy = y[0];
    for (size_t i = 1; (i < x.size()); i++)
    {
        miny = (y[i] < miny) ? y[i] : miny;
        maxy = (y[i] > maxy) ? y[i] : maxy;
        minx = (x[i] < minx) ? x[i] : minx;
        maxx = (x[i] > maxx) ? x[i] : maxx;
    }
    if (ehisto == eHisto::X)
    {
        delta = maxx-minx;
        minh  = minx;
    }
    else if (ehisto == eHisto::Y)
    {
        delta = maxy-miny;
        minh  = miny;
    }
    else
    {
        return eStats::INVALID_INPUT;
    }

    if (binwidth == 0)
    {
        binwidth = (delta)/nbins;
    }
    else
    {
        nbins = std::min(1, gmx_dnint((delta)/binwidth + 0.5));
    }
    if (nbins == 0)
    {
        return eStats::NO_POINTS;
    }
    xx->resize(nbins, 0.0);
    nindex.resize(nbins, -1);
    for (int i = 0; (i < nbins); i++)
    {
        (*xx)[i] = minh + binwidth*(i+0.5);
    }
    if (normalized == 0)
    {
        dd = 1;
    }
    else
    {
        dd = 1.0/(binwidth*x.size());
    }

    yy->resize(nbins, 0.0);
    for (size_t i = 0; (i < x.size()); i++)
    {
        if (ehisto == eHisto::Y)
        {
            index = static_cast<int>((y[i]-miny)/binwidth);
        }
        else if (ehisto == eHisto::X)
        {
            index = static_cast<int>((x[i]-minx)/binwidth);
        }
        if (index < 0)
        {
            index = 0;
        }
        if (index > nbins-1)
        {
            index = nbins-1;
        }
        (*yy)[index] += dd;
        nindex[index]++;
    }
    if (*nb == 0)
    {
        *nb = nbins;
    }

    return eStats::OK;
}

std::map<eStats,  const char *> stats_error =
{
    { eStats::OK, "All well in STATS land" },
    { eStats::NO_POINTS, "Not enough points" },
    { eStats::NO_MEMORY, "Not enough memory" },
    { eStats::ERROR,     "Unknown error" },
    { eStats::INVALID_INPUT, "Invalid histogram input" },
    { eStats::NOT_IMPLEMENTED, "Not implemented yet" }
};

const char *gmx_stats_message(eStats estats)
{
    return stats_error[estats];
}

/* Old convenience functions, should be merged with the core
   statistics above. */
eStats lsq_y_ax(int n, real x[], real y[], real *a)
{
    gmx_stats lsq;
    eStats    ok;
    real      da, chi2, Rfit;

    lsq.add_points(n, x, y, nullptr, nullptr);
    ok = lsq.get_a(elsqWEIGHT_NONE, a, &da, &chi2, &Rfit);

    return ok;
}

static eStats low_lsq_y_ax_b(int n, const real *xr, const double *xd, real yr[],
                          real *a, real *b, real *r, real *chi2)
{
    gmx_stats lsq;
    eStats    ok;

    for (int i = 0; (i < n); i++)
    {
        double pt;

        if (xd != nullptr)
        {
            pt = xd[i];
        }
        else if (xr != nullptr)
        {
            pt = xr[i];
        }
        else
        {
            gmx_incons("Either xd or xr has to be non-NULL in low_lsq_y_ax_b()");
        }

        if ((ok = lsq.add_point(pt, yr[i], 0, 0)) != eStats::OK)
        {
            return ok;
        }
    }
    ok = lsq.get_ab(elsqWEIGHT_NONE, a, b, nullptr, nullptr, chi2, r);

    return ok;
}

eStats lsq_y_ax_b(int n, real x[], real y[], real *a, real *b, real *r, real *chi2)
{
    return low_lsq_y_ax_b(n, x, nullptr, y, a, b, r, chi2);
}

eStats lsq_y_ax_b_xdouble(int n, double x[], real y[], real *a, real *b,
                       real *r, real *chi2)
{
    return low_lsq_y_ax_b(n, nullptr, x, y, a, b, r, chi2);
}

eStats lsq_y_ax_b_error(int n, real x[], real y[], real dy[],
                     real *a, real *b, real *da, real *db,
                     real *r, real *chi2)
{
    gmx_stats lsq;
    eStats    ok;

    for (int i = 0; (i < n); i++)
    {
        ok = lsq.add_point(x[i], y[i], 0, dy[i]);
        if (ok != eStats::OK)
        {
            return ok;
        }
    }
    ok = lsq.get_ab(elsqWEIGHT_Y, a, b, da, db, chi2, r);

    return ok;
}

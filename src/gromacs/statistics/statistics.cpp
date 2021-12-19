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
    return add_point(x_.size(), y, 0, dy);
}

eStats gmx_stats::add_point(double xx, double yy,
                            double dxx, double dyy)
{
    x_.push_back(xx);
    y_.push_back(yy);
    dx_.push_back(dxx);
    dy_.push_back(dyy);
    computed_ = false;

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
    while ((outlier == 0) && (np_c_ < x_.size()))
    {
        r       = std::abs(x_[np_c_] - y_[np_c_]);
        outlier = static_cast<int>(r > rmsd*level);
        if (outlier)
        {
            if (nullptr != xx)
            {
                *xx  = x_[np_c_];
            }
            if (nullptr != yy)
            {
                *yy  = y_[np_c_];
            }
            if (nullptr != dxx)
            {
                *dxx = dx_[np_c_];
            }
            if (nullptr != dyy)
            {
                *dyy = dy_[np_c_];
            }
        }
        np_c_++;

        if (outlier)
        {
            return eStats::OK;
        }
    }

    np_c_ = 0;

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
    double yy, yx, xx, sx, sy, d2;
    double ssxx, ssyy, ssxy;
    double w, wtot, yx_nw, sy_nw, sx_nw, yy_nw, xx_nw, dx2, dy2;

    int    N = x_.size();

    if (!computed_)
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
        mae_ = 0, mse_ = 0;
        for (int i = 0; (i < N); i++)
        {
            double dd = y_[i]-x_[i];
            d2 += gmx::square(dd);
            
            mae_ += fabs(dd);
            mse_ += dd;
            if ((dy_[i] != 0.0) && (weight == elsqWEIGHT_Y))
            {
                w = 1/gmx::square(dy_[i]);
            }
            else
            {
                w = 1;
            }

            wtot  += w;

            xx    += w*gmx::square(x_[i]);
            xx_nw += gmx::square(x_[i]);

            yy    += w*gmx::square(y_[i]);
            yy_nw += gmx::square(y_[i]);

            yx    += w*y_[i]*x_[i];
            yx_nw += y_[i]*x_[i];

            sx    += w*x_[i];
            sx_nw += x_[i];

            sy    += w*y_[i];
            sy_nw += y_[i];
        }

        /* Compute average, sigma and error */
        mae_        = mae_/N;
        mse_        = mse_/N; 
        aver_       = sy_nw/N;

        double dd = yy_nw/N - gmx::square(sy_nw/N);
        if (dd > 0)
        {
            sigma_aver_ = std::sqrt(dd);
        }
        else if (N == 1)
        {
            sigma_aver_ = dy_[0];
        }
        else
        {
            sigma_aver_ = 0;
        }
        error_      = sigma_aver_/std::sqrt(static_cast<double>(N));

        /* Compute RMSD between x and y */
        rmsd_ = std::sqrt(d2/N);

        /* Correlation coefficient for data */
        yx_nw       /= N;
        xx_nw       /= N;
        yy_nw       /= N;
        sx_nw       /= N;
        sy_nw       /= N;
        ssxx         = N*(xx_nw - gmx::square(sx_nw));
        ssyy         = N*(yy_nw - gmx::square(sy_nw));
        ssxy         = N*(yx_nw - (sx_nw*sy_nw));
        Rdata_       = std::sqrt(gmx::square(ssxy)/(ssxx*ssyy));

        /* Compute straight line through datapoints, either with intercept
           zero (result in aa) or with intercept variable (results in a
           and b) */
        yx = yx/wtot;
        xx = xx/wtot;
        sx = sx/wtot;
        sy = sy/wtot;

        aa_ = (yx/xx);
        a_  = (yx-sx*sy)/(xx-sx*sx);
        b_  = (sy)-(a_)*(sx);

        /* Compute chi2, deviation from a line y = ax+b. Also compute
           chi2aa which returns the deviation from a line y = ax. */
        chi2_   = 0;
        chi2aa_ = 0;
        for (int i = 0; (i < N); i++)
        {
            real ddy = 1;
            if (dy_[i] > 0)
            {
                ddy = dy_[i];
            }
            chi2aa_ += gmx::square((y_[i]-(aa_*x_[i]))/ddy);
            chi2_   += gmx::square((y_[i]-(a_*x_[i]+b_))/ddy);
        }
        if (N > 2)
        {
            chi2_   = std::sqrt(chi2_/(N-2));
            chi2aa_ = std::sqrt(chi2aa_/(N-2));

            /* Look up equations! */
            dx2            = (xx-sx*sx);
            dy2            = (yy-sy*sy);
            sigma_a_ = std::sqrt(chi2_/((N-2)*dx2));
            sigma_b_ = sigma_a_*std::sqrt(xx);
            Rfit_    = std::abs(ssxy)/std::sqrt(ssxx*ssyy);
            Rfitaa_  = Rfit_; //aa*std::sqrt(dx2/dy2);
        }
        else
        {
            chi2_    = 0;
            chi2aa_  = 0;
            sigma_a_ = 0;
            sigma_b_ = 0;
            Rfit_    = 0;
            Rfitaa_  = 0;
        }

        computed_ = true;
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
        *aa    = a_;
    }
    if (nullptr != bb)
    {
        *bb    = b_;
    }
    if (nullptr != daa)
    {
        *daa   = sigma_a_;
    }
    if (nullptr != dbb)
    {
        *dbb   = sigma_b_;
    }
    if (nullptr != chi2a)
    {
        *chi2a = chi2_;
    }
    if (nullptr != Rfita)
    {
        *Rfita = Rfit_;
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
        *aa    = a_;
    }
    if (nullptr != daa)
    {
        *daa   = sigma_a_;
    }
    if (nullptr != chi2a)
    {
        *chi2a = chi2_;
    }
    if (nullptr != Rfita)
    {
        *Rfita = Rfit_;
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

    *avera = aver_;

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
        *msea = mse_;
    }
    if (NULL != maea)
    {
        *maea = mae_;
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
        *avera  = aver_;
    }
    if (nullptr != sigmaa)
    {
        *sigmaa = sigma_aver_;
    }
    if (nullptr != errora)
    {
        *errora = error_;
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

    *sigmaa = sigma_aver_;

    return eStats::OK;
}

eStats gmx_stats::get_error(real *errors)
{
    eStats     ok;

    if ((ok = compute(elsqWEIGHT_NONE)) != eStats::OK)
    {
        return ok;
    }

    *errors = error_;

    return eStats::OK;
}

eStats gmx_stats::get_corr_coeff(real *Ra)
{
    eStats     ok;

    if ((ok = gmx_stats::compute(elsqWEIGHT_NONE)) != eStats::OK)
    {
        return ok;
    }

    *Ra = Rdata_;

    return eStats::OK;
}

eStats gmx_stats::get_rmsd(real *rmsda)
{
    eStats        ok;

    if ((ok = compute(elsqWEIGHT_NONE)) != eStats::OK)
    {
        return ok;
    }

    *rmsda = rmsd_;

    return eStats::OK;
}

eStats gmx_stats::dump_xy(FILE *fp)
{
    for (size_t i = 0; (i < x_.size()); i++)
    {
        fprintf(fp, "%12g  %12g  %12g  %12g\n", x_[i], y_[i],
                dx_[i], dy_[i]);
    }

    return eStats::OK;
}

eStats gmx_stats::remove_outliers(double level)
{
    int        iter  = 1, done = 0;
    eStats     ok;
    real       rmsd, r;

    while ((x_.size() >= 10) && !done)
    {
        if ((ok = get_rmsd(&rmsd)) != eStats::OK)
        {
            return ok;
        }
        done = 1;
        size_t i;
        for (i = 0; (i < x_.size()); )
        {
            r = std::abs(x_[i]-y_[i]);
            if (r > level*rmsd)
            {
                fprintf(stderr, "Removing outlier, iter = %d, rmsd = %g, x = %g, y = %g\n",
                        iter, rmsd, x_[i], y_[i]);
                if (i < x_.size()-1)
                {
                    x_[i]  = x_[x_.size()-1];
                    y_[i]  = y_[x_.size()-1];
                    dx_[i] = dx_[x_.size()-1];
                    dy_[i] = dy_[x_.size()-1];
                }
                done = 0;
            }
            else
            {
                i++;
            }
        }
        x_.resize(i);
        y_.resize(i);
        dx_.resize(i);
        dy_.resize(i);
        iter++;
    }

    return eStats::OK;
}

eStats gmx_stats::make_histogram(real binwidth, int *nb,
                                 eHisto ehisto, bool normalized,
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
    if (x_.size() <= 1)
    {
        return eStats::NO_POINTS;
    }
    minx = maxx = x_[0];
    miny = maxy = y_[0];
    for (size_t i = 1; (i < x_.size()); i++)
    {
        miny = std::min(y_[i], miny);
        maxy = std::max(y_[i], maxy);
        minx = std::min(x_[i], minx);
        maxx = std::max(x_[i], maxx);
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
        nbins = std::max(1, gmx_dnint((delta)/binwidth + 0.5));
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
    if (!normalized)
    {
        dd = 1;
    }
    else
    {
        dd = 1.0/(binwidth*x_.size());
    }

    yy->resize(nbins, 0.0);
    for (size_t i = 0; (i < x_.size()); i++)
    {
        if (ehisto == eHisto::Y)
        {
            index = static_cast<int>((y_[i]-miny)/binwidth);
        }
        else if (ehisto == eHisto::X)
        {
            index = static_cast<int>((x_[i]-minx)/binwidth);
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

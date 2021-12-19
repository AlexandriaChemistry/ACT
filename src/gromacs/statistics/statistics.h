/*
 * This file is part of the GROMACS molecular simulation package.
 *
 * Copyright (c) 1991-2000, University of Groningen, The Netherlands.
 * Copyright (c) 2001-2008, The GROMACS development team.
 * Copyright (c) 2010,2014,2015, by the GROMACS development team, led by
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
/*! \libinternal \file
 * \brief
 * Declares simple statistics toolbox
 *
 * \authors David van der Spoel <david.vanderspoel@icm.uu.se>
 * \inlibraryapi
 */
#ifndef GMX_STATISTICS_H
#define GMX_STATISTICS_H

#include <cstdio>

#include <vector>

#include "gromacs/utility/real.h"

//! Error codes returned by the routines
enum class eStats 
{
    //! All well
    OK,
    //! Not enough points
    NO_POINTS,
    //! Not enough memory
    NO_MEMORY,
    //! Unknown error
    ERROR,
    //! Invalid user input
    INVALID_INPUT,
    //! Something not yet implementd
    NOT_IMPLEMENTED
};

//! Enum for statistical weights
enum {
    elsqWEIGHT_NONE, elsqWEIGHT_X, elsqWEIGHT_Y,
    elsqWEIGHT_XY, elsqWEIGHT_NR
};

//! Enum determining which coordinate to histogram
enum class eHisto {
    //! Histogram of X values
    X, 
    //! Histogram of Y values
    Y
};

class gmx_stats
{
private:
    std::vector<double> x_;
    std::vector<double> y_;
    std::vector<double> dx_;
    std::vector<double> dy_;
    double              aa_         = 0;
    double              a_          = 0;
    double              b_          = 0;
    double              sigma_a_    = 0;
    double              sigma_b_    = 0;
    double              aver_       = 0;
    double              sigma_aver_ = 0;
    double              error_      = 0;
    double              rmsd_       = 0;
    double              Rdata_      = 0;
    double              Rfit_       = 0;
    double              Rfitaa_     = 0;
    double              chi2_       = 0;
    double              chi2aa_     = 0;
    double              mse_        = 0;
    double              mae_        = 0;
    bool                computed_   = false;
    size_t              np_c_       = 0;

    eStats compute(int weight);
    
public:    
    //! Constructor
    gmx_stats() {}
    
    /*! \brief
     * Remove outliers from a straight line, where level in units of
     * sigma. Level needs to be larger than one obviously.
     * \param[in] stats The data structure
     * \param[in] level The sigma level
     * \return error code
     */
    eStats remove_outliers(double level);

    /*! \brief
     * Add a point to the data set with just y values.
     * \param[in] stats The data structure
     * \param[in] y   The y value
     * \param[in] dy  The error in the y value
     * \return error code
     */
    eStats add_point_ydy(double y, double dy);
                        
    /*! \brief
     * Add a point to the data set
     * \param[in] stats The data structure
     * \param[in] x   The x value
     * \param[in] y   The y value
     * \param[in] dx  The error in the x value
     * \param[in] dy  The error in the y value
     * \return error code
     */
    eStats add_point(double x, double y,
                     double dx, double dy);

    /*! \brief
     * Add a series of datapoints at once. The arrays dx and dy may
     * be NULL in that case zero uncertainties will be assumed.
     *
     * \param[in] stats The data structure
     * \param[in] n   Number of points
     * \param[in] x   The array of x values
     * \param[in] y   The array of y values
     * \param[in] dx  The error in the x value
     * \param[in] dy  The error in the y value
     * \return error code
     */
    eStats add_points(int n, real *x, real *y,
                      real *dx, real *dy);

    const std::vector<double> getX() const { return x_; }
    const std::vector<double> getY() const { return y_; }
    /*! \brief
     * Delivers data points from the statistics.
     *
     * Should be used in a while loop. Variables for either
     * pointer may be NULL, in which case the routine can be used as an
     * expensive point counter.
     * Return the data points one by one. Return eStats::OK while there are
     * more points, and returns eStats::NOPOINTS when the last point has
     * been returned.
     * If level > 0 then the outliers outside level*sigma are reported
     * only.
     * \param[in] stats The data structure
     * \param[out] x    An x value
     * \param[out] y    An y value
     * \param[out] dx   The error in the x value
     * \param[out] dy   The error in the y value
     * \param[in]  level sigma level (see above)
     * \return error code
     */
    eStats get_point(real *x, real *y,
                     real *dx, real *dy, real level);

    /*! \brief
     * Fit the data to y = ax + b, possibly weighted, if uncertainties
     * have been input. da and db may be NULL.
     * \param[in] stats The data structure
     * \param[in] weight type of weighting
     * \param[out] a slope
     * \param[out] b intercept
     * \param[out] da sigma in a
     * \param[out] db sigma in b
     * \param[out] chi2 normalized quality of fit
     * \param[out] Rfit correlation coefficient
     * \return error code
     */
    eStats get_ab(int weight,
                  real *a, real *b,
                  real *da, real *db, real *chi2, real *Rfit);
    
    /*! \brief
     * Fit the data to y = ax, possibly weighted, if uncertainties have
     * have been input. da and db may be NULL.
     * \param[in] stats The data structure
     * \param[in] weight type of weighting
     * \param[out] a slope
     * \param[out] da sigma in a
     * \param[out] chi2 normalized quality of fit
     * \param[out] Rfit correlation coefficient
     * \return error code
     */
    eStats get_a(int weight,
                 real *a, real *da, real *chi2, real *Rfit);

    /*! \brief
     * Get the correlation coefficient.
     * \param[in]  stats The data structure
     * \param[out] R the correlation coefficient between the data (x and y) as input to the structure.
     * \return error code
     */
    eStats get_corr_coeff(real *R);

    /*! \brief
     * Get the root mean square deviation.
     * \param[in]  stats The data structure
     * \param[out] rmsd  the root mean square deviation between x and y values.
     * \return error code
     */
    eStats get_rmsd(real *rmsd);

    /*! \brief
     * Get the number of points.
     * \return the number
     */
    size_t get_npoints() const { return x_.size(); };
    
    /*! \brief
     * Computes and returns the average value.
     * \param[in]  stats The data structure
     * \param[out] aver  Average value
     * \return error code
     */
    eStats get_average(real *aver);
    
    /*! \brief
     * Computes and returns the standard deviation.
     * \param[in]  stats The data structure
     * \param[out] sigma  Standard deviation
     * \return error code
     */
    eStats get_sigma(real *sigma);
    
    /*! \brief
     * Computes and returns the standard error.
     * \param[in]  stats The data structure
     * \param[out] error Standard error
     * \return error code
     */
    eStats get_error(real *error);

    /*! \brief
     * Pointers may be null, in which case no assignment will be done.
     * \param[in]  stats The data structure
     * \param[out] aver  Average value
     * \param[out] sigma  Standard deviation
     * \param[out] error Standard error
     * \return error code
     */
    eStats get_ase(real *aver, real *sigma, real *error);

    /*! \brief
     * Return mean signed error and mean absolute error
     * \param[in]  stats The data structure
     * \param[out] mse  mean signed error
     * \param[out] mae  mean absolute
     * \return error code
     */
    eStats get_mse_mae(real *mse, real *mae);

    /*! \brief
     * Dump the x, y, dx, dy data to a text file
     * \param[in]  stats The data structure
     * \param[in] fp  File pointer
     * \return error code
     */
    eStats dump_xy(FILE *fp);
    
    /*! \brief
     * Make a histogram of the data present.
     *
     * Uses either binwidth to
     * determine the number of bins, or nbins to determine the binwidth,
     * therefore one of these should be zero, but not the other. If *nbins = 0
     * the number of bins will be returned in this variable. ehisto should be one of
     * ehistoX or ehistoY. If
     * normalized, the integral of the histogram will be
     * normalized to one. The output is in two arrays, *x and *y, to which
     * you should pass a pointer. Memory for the arrays will be allocated
     * as needed. Function returns one of the eStats codes.
     * \param[in] binwidth For the histogram
     * \param[in] nbins    Number of bins
     * \param[in] ehisto   Type (see enum above)
     * \param[in] normalized see above
     * \param[out] x see above
     * \param[out] y see above
     * \return error code
     */
    eStats make_histogram(real binwidth, int *nbins,
                          eHisto ehisto, bool normalized,
                          std::vector<double> *x,
                          std::vector<double> *y);
};

/*! \brief
 * Return message belonging to error code
 * \param[in] eStats error code
 */
const char *gmx_stats_message(eStats estats);

/****************************************************
 * Some statistics utilities for convenience: useful when a complete data
 * set is available already from another source, e.g. an xvg file.
 ****************************************************/
/*! \brief
 * Fit a straight line y=ax thru the n data points x, y, return the
 * slope in *a.
 * \param[in] n number of points
 * \param[in] x data points x
 * \param[in] y data point y
 * \param[out] a slope
 * \return error code
 */
eStats lsq_y_ax(int n, real x[], real y[], real *a);

/*! \brief
 * Fit a straight line y=ax+b thru the n data points x, y.
 * \param[in] n number of points
 * \param[in] x data points x
 * \param[in] y data point y
 * \param[out] a slope
 * \param[out] b intercept
 * \param[out] r correlation coefficient
 * \param[out] chi2 quality of fit
 * \return error code
 */
eStats lsq_y_ax_b(int n, real x[], real y[], real *a, real *b, real *r,
                  real *chi2);

/*! \copydoc lsq_y_ax_b
 */
eStats lsq_y_ax_b_xdouble(int n, double x[], real y[],
                          real *a, real *b, real *r, real *chi2);

/*! \brief
 * Fit a straight line y=ax+b thru the n data points x, y.
 * \param[in] n number of points
 * \param[in] x data points x
 * \param[in] y data point y
 * \param[in] dy uncertainty in data point y
 * \param[out] a slope
 * \param[out] b intercept
 * \param[out] da error in slope
 * \param[out] db error in intercept
 * \param[out] r correlation coefficient
 * \param[out] chi2 quality of fit
 * \return error code
 */
eStats lsq_y_ax_b_error(int n, real x[], real y[], real dy[],
                        real *a, real *b, real *da, real *db,
                        real *r, real *chi2);

#endif

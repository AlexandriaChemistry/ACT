/*
 * This source file is part of the Alexandria Chemistry Toolkit.
 *
 * Copyright (C) 2014-2021
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

#ifndef ALEXANDRIA_BAYES_H
#define ALEXANDRIA_BAYES_H

#include <functional>
#include <random>
#include <vector>

#include "gromacs/commandline/pargs.h"
#include "gromacs/fileio/oenv.h"
#include "gromacs/fileio/xvgr.h"
#include "gromacs/math/functions.h"
#include "gromacs/math/units.h"
#include "gromacs/utility/fatalerror.h"
#include "gromacs/utility/real.h"

#include "molselect.h"
#include "confighandler.h"

namespace alexandria
{

//! How to perform the calculation of deviations (chi-squared)
enum class CalcDev {
    //! Do it in parallel
    Parallel = 1,
    //! Do it on the master only
    Master = 2,
    //! Do the final one only (typically on the master)
    Final = 3
};

class Sensitivity
{
private:
    //! \brief parameter values used
    std::vector<double> p_;
    //! \brief chi2 values obtains
    std::vector<double> chi2_;
    //! \brief Constants for the parabola fitting
    double a_ = 0, b_ = 0, c_ = 0;
public:
    //! \brief Constructor
    Sensitivity() {}

    /*! \brief
     * Add a point
     * \param[in] p    The parameter value
     * \param[in] chi2 The chi-squared value
     */
    void add(double p, double chi2)
    {
        p_.push_back(p);
        chi2_.push_back(chi2);
    }
    /*! \brief
     * Compute the fit to the curve
     * \param[in] fp File pointer for debugging output
     */
    void computeForceConstants(FILE *fp);

    //! Return the constants after computation
    double a() const { return a_; }
    double b() const { return b_; }
    double c() const { return c_; }

    /*! \brief Print output
     * \param[in] fp    File pointer for output
     * \param[in] label Label for identifying the parameter
     */
    void print(FILE *fp, const std::string &label);
};

}  //namespace alexandria

#endif //ALEXANDRIA_BAYES_H

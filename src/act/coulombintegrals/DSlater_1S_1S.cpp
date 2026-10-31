/*
 * This file is part of the GROMACS molecular simulation package.
 *
 * Copyright (c) 2014,2015, by the GROMACS development team, led by
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

#include "slater_low.h"

#if HAVE_LIBCLN
cl_R DSlater_1S_1S(cl_R r, cl_R xi, cl_R xj)
{
    cl_R S, rxi, rxj;

    rxi = rxj = S = ZERO;
    rxi = r*xi;
    rxj = r*xj;
    if (xi == xj)
    {
        if (r == 0LL)
        {
            S = 0LL

            ;
        }
        else
        {
            S = -(-33LL*xi + 48LL*exp(2LL*r*xi)*xi - 36LL*r*Power(xi, 2LL) -

                  12LL*Power(r, 2LL)*Power(xi, 3LL))/(24LL*exp(2LL*r*xi)*r) +

                (-24LL + 24LL*exp(2LL*r*xi) - 33LL*r*xi - 18LL*Power(r, 2LL)*Power(xi, 2LL) -

                 4LL*Power(r, 3LL)*Power(xi, 3LL))/(24LL*exp(2LL*r*xi)*Power(r, 2LL)) +

                (xi*(-24LL + 24LL*exp(2LL*r*xi) - 33LL*r*xi - 18LL*Power(r, 2LL)*Power(xi, 2LL) -

                     4LL*Power(r, 3LL)*Power(xi, 3LL)))/(12LL*exp(2LL*r*xi)*r)

            ;
        }

    }
    else
    {
        if (r == 0LL)
        {
            S = 0LL

            ;
        }
        else
        {
            S = (exp(2LL*r*(xi + xj))*Power(Power(xi, 2LL) - Power(xj, 2LL), 3LL) +

                 exp(2LL*r*xj)*Power(xj, 4LL)*

                 (-3LL*Power(xi, 2LL) - r*Power(xi, 3LL) + Power(xj, 2LL) + r*xi*Power(xj, 2LL)) -

                 exp(2LL*r*xi)*Power(xi, 4LL)*

                 (Power(xi, 2LL)*(1LL + r*xj) - Power(xj, 2LL)*(3LL + r*xj)))/

                (exp(2LL*r*(xi + xj))*Power(r, 2LL)*Power(xi - xj, 3LL)*Power(xi + xj, 3LL)) +

                (2LL*(exp(2LL*r*(xi + xj))*Power(Power(xi, 2LL) - Power(xj, 2LL), 3LL) +

                      exp(2LL*r*xj)*Power(xj, 4LL)*

                      (-3LL*Power(xi, 2LL) - r*Power(xi, 3LL) + Power(xj, 2LL) + r*xi*Power(xj, 2LL)) -

                      exp(2LL*r*xi)*Power(xi, 4LL)*

                      (Power(xi, 2LL)*(1LL + r*xj) - Power(xj, 2LL)*(3LL + r*xj))))/

                (exp(2LL*r*(xi + xj))*r*Power(xi - xj, 3LL)*Power(xi + xj, 2LL)) -

                (2LL*exp(2LL*r*(xi + xj))*(xi + xj)*Power(Power(xi, 2LL) - Power(xj, 2LL), 3LL) +

                 exp(2LL*r*xj)*Power(xj, 4LL)*(-Power(xi, 3LL) + xi*Power(xj, 2LL)) +

                 2LL*exp(2LL*r*xj)*Power(xj, 5LL)*

                 (-3LL*Power(xi, 2LL) - r*Power(xi, 3LL) + Power(xj, 2LL) + r*xi*Power(xj, 2LL)) -

                 exp(2LL*r*xi)*Power(xi, 4LL)*(Power(xi, 2LL)*xj - Power(xj, 3LL)) -

                 2LL*exp(2LL*r*xi)*Power(xi, 5LL)*

                 (Power(xi, 2LL)*(1LL + r*xj) - Power(xj, 2LL)*(3LL + r*xj)))/

                (exp(2LL*r*(xi + xj))*r*Power(xi - xj, 3LL)*Power(xi + xj, 3LL))

            ;
        }

    }
    return S;
}

#else
#include <cmath>

static double pow3(double x)
{
    return x*x*x;
}
static double pow4(double x)
{
    return x*x*x*x;
}

/*! \brief The derivative of the integral with the formulas for equal and different widths
 * \param[in] r  The distance
 * \param[in] xi The width of the first charge
 * \param[in] xj The width of the second charge
 * \return the derivative of the integral
 */
static double DSlater_1S_1SFormula(double r, double xi, double xj)
{
    double S, rxi, rxj;

    rxi = rxj = S = 0;
    rxi = r*xi;
    rxj = r*xj;
    if (xi == xj)
    {
        if (r == 0)
        {
            S = 0

            ;
        }
        else
        {
            double exp2rxi = exp(2*rxi);
            S = -(-33*xi + 48*exp2rxi*xi - 36*rxi*xi -

                  12*r*r*pow3(xi))/(24*exp2rxi*r) +

                (-24 + 24*exp2rxi - 33*rxi - 18*rxi*rxi -

                 4*pow3(rxi))/(24*exp2rxi*r*r) +

                (xi*(-24 + 24*exp2rxi - 33*rxi - 18*rxi*rxi -

                     4*pow3(rxi)))/(12*exp2rxi*r)

            ;
        }

    }
    else
    {
        if (r == 0)
        {
            S = 0

            ;
        }
        else
        {
            double exp2rxixj = exp(2*r*(xi + xj));
            double exp2rxj   = exp(2*rxj);
            double exp2rxi   = exp(2*rxi);
            S = (exp2rxixj*pow3(xi*xi - xj*xj) +

                 exp2rxj*pow4(xj)*

                 (-3*xi*xi - r*pow3(xi) + xj*xj + rxi*xj*xj) -

                 exp2rxi*pow4(xi)*

                 (xi*xi*(1 + rxj) - xj*xj*(3 + rxj)))/

                (exp2rxixj*r*r*pow3(xi - xj)*pow3(xi + xj)) +

                (2*(exp2rxixj*pow3(xi*xi - xj*xj) +

                      exp2rxj*pow4(xj)*

                      (-3*xi*xi - r*pow3(xi) + xj*xj + rxi*xj*xj) -

                      exp2rxi*pow4(xi)*

                      (xi*xi*(1 + rxj) - pow(xj, 2)*(3 + rxj))))/

                (exp2rxixj*r*pow3(xi - xj)*pow(xi + xj, 2)) -

                (2*exp2rxixj*(xi + xj)*pow3(xi*xi - xj*xj) +

                 exp2rxj*pow4(xj)*(-pow3(xi) + xi*xj*xj) +

                 2*exp2rxj*pow4(xj)*xj*

                 (-3*xi*xi - r*pow3(xi) + xj*xj + rxi*xj*xj) -

                 exp2rxi*pow4(xi)*(xi*xi*xj - pow3(xj)) -

                 2*exp2rxi*pow4(xi)*xi*

                 (xi*xi*(1 + rxj) - xj*xj*(3 + rxj)))/

                (exp2rxixj*r*pow3(xi - xj)*pow3(xi + xj))

            ;
        }

    }
    return S;
}

//! Relative difference of the widths below which the formula for different widths is not used
static const double closeWidths = 3e-2;

double DSlater_1S_1S(double r, double xi, double xj)
{
    // For close widths the formula for different widths loses its precision,
    // since it divides by the cube of the difference of the widths. As the
    // derivative of the integral is symmetric in the widths, it differs from that of the mean width
    // by a term of second order in their difference, which is scaled from
    // the formula for different widths at a relative difference closeWidths.
    // The formula for different widths is less precise for the derivative
    // than for the integral, therefore closeWidths is larger here. This is
    // accurate to about 1e-7 and continuous at closeWidths.
    double xm = 0.5*(xi + xj);
    double d  = (xm > 0) ? std::abs(xi - xj)/xm : 0;
    if (r == 0 || xi == xj || d >= closeWidths)
    {
        return DSlater_1S_1SFormula(r, xi, xj);
    }
    double s0 = DSlater_1S_1SFormula(r, xm, xm);
    double s1 = DSlater_1S_1SFormula(r, xm*(1 - closeWidths/2), xm*(1 + closeWidths/2));
    double f  = d/closeWidths;
    return s0 + (s1 - s0)*f*f;
}

#endif

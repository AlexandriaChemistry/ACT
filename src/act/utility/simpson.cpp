/*
 * This source file is part of the Alexandria Chemistry Toolkit.
 *
 * Copyright (C) 2026
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
 * \author David van der Spoel <david.vanderspoel@icm.uu.se>
 */

#include "simpson.h"

#include <cmath>

#include "act/basics/msg_handler.h"
#include "gromacs/utility/stringutil.h"

double simpsonIntegrate(alexandria::MsgHandler    *msghandler,
                        bool                       spherical,
                        const std::vector<double> &x,
                        const std::vector<double> &y,
                        bool                       zeroPadding)
{
    if (x.size() != y.size())
    {
        msghandler->msg(alexandria::ACTStatus::Error,
                        gmx::formatString("x and y vectors of different size (%zu vs. %zu)", x.size(), y.size()));
    }
    if (x.size() % 2 == 0 && !zeroPadding)
    {
        msghandler->msg(alexandria::ACTStatus::Error,
                        gmx::formatString("vectors have even number of data points (%zu) and it is not allowed to add a point", x.size()));
    }
    if (x.size() < 2)
    {
        msghandler->msg(alexandria::ACTStatus::Warning,
                        gmx::formatString("vectors have just %zu data points, return zero", x.size()));
        return 0.0;
    }
    double toler    = 1e-4;
    double dx       = x[1] - x[0];
    double dx_3     = dx/3.0;
    double integral = y[0]*dx_3;
    if (spherical)
    {
        integral *= 4*M_PI*x[0]*x[0];
    }
    for(size_t i = 1; i < x.size()-1; i++)
    {
        double dx2 = x[i+1]-x[i];
        if ((dx-dx2)/(dx+dx2) > toler)
        {
            msghandler->msg(alexandria::ACTStatus::Error,
                            gmx::formatString("x vectors is not uniformly space. Got %g first, and now %g", dx, dx2));
        }
        double yy = y[i];
        if (spherical)
        {
            yy *= 4*M_PI*x[i]*x[i];
        }
        if (i % 2 == 1)
        {
            integral += 4*yy*dx_3;
        }
        else
        {
            integral += 2*yy*dx_3;
        }
    }
    if (x.size() % 2 == 1)
    {
        // Last point!
        double yy = y[x.size()-1];
        double xx = x[x.size()-1]+dx;
        if (spherical)
        {
            yy *= 4*M_PI*xx*xx;
        }
        integral += yy*dx_3;
    }
    // If we use zeroPadding we assume that the last entry in y is zero,
    // so we do not have to add anything.
    return integral;
}


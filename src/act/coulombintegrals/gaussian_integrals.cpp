/*
 * This source file is part of the Alexandria Chemistry Toolkit.
 *
 * Copyright (C) 2014-2022
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

#include "actpre.h"

#include <math.h>
#include <stdio.h>

#include "coulombintegrals.h"

static double sqr(double x)
{
    return x*x;
}

double Nuclear_GG(double r, double zeta)
{
    /* This routine may be called with zeta 0.
     * In that case it is a simple 1/r interaction.
     */
    if (zeta == 0)
    {
        if (r == 0)
        {
            return 0;
        }
        else
        {
            return 1.0/r;
        }
    }
    else if (r == 0)
    {
        return 2.0*zeta/sqrt(M_PI);
    }
    else
    {
        return erf(zeta*r)/r;
    }
}

double Coulomb_GG(double r, double zi, double zj)
{
    double zeff = 0;

    /* This routine may be called with one or both zeta 0.
     * In that case it is either a Nuclear_GG or simple 1/r interaction.
     */
    if ((zi != 0) || (zj != 0))
    {
        if (zi == 0)
        {
            zeff = zj;
        }
        else if (zj == 0)
        {
            zeff = zi;
        }
        else
        {
            zeff = zi*zj/sqrt(sqr(zi)+sqr(zj));
        }
    }
    
    return Nuclear_GG(r, zeff);
}

double DNuclear_GG(double r, double z)
{
    if (r == 0)
    {
        return 0;
    }
    if (z == 0)
    {
        return 1.0/sqr(r);
    }
    else
    {
        double rz = r * z;
        
        return  (-1)*((2.0/sqrt(M_PI))*exp(-sqr(rz))*(z/r) - erf(rz)/sqr(r));
    }
}

double DCoulomb_GG(double r, double zi, double zj)
{
    double zeff = 0;

    if ((zi != 0) || (zj != 0))
    {
        if (zi == 0)
        {
            zeff = zj;
        }
        else if (zj == 0)
        {
            zeff = zi;
        }
        else
        {
            zeff = zi*zj/sqrt(sqr(zi)+sqr(zj));
        }
    }
    
    return DNuclear_GG(r, zeff);
}

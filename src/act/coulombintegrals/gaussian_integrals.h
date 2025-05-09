/*
 * This source file is part of the Alexandria Chemistry Toolkit.
 *
 * Copyright (C) 2014-2024
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

/*! \file
 * \brief
 * Provides routines for computing Coulomb integrals analytically.
 * The Slater based code depends on the optional CLN library for arbitrary
 * precision arithmetic.
 *
 * \author David van der Spoel <david.vanderspoel@icm.uu.se>
 * \inpublicapi
 * \ingroup module_alexandria
 */
#ifndef _COULOMBINTEGRALS_H
#define _COULOMBINTEGRALS_H


#include "gromacs/math/functions.h"

const double invsqrt_pi = 1.0/std::sqrt(M_PI);

static inline gmx_unused double sqr(double x)
{
    return x*x;
}

/*! \brief 
 * Compute the Coulomb overlap integral between a Gaussian distributed charge
 * and a point charge. The width xi may be zero.
 * No units are introduced or converted but nm are assumed.
 * No factor 1/4 pi epsilon0 is introduced.
 *
 * \param[in] r    distance in nm
 * \param[in] zeta Gaussian width in 1/nm 
 * \return    Integral value
 */    
static inline gmx_unused double Nuclear_GG(double r, double zeta)
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
        return 2.0*zeta*invsqrt_pi;
    }
    else
    {
        return erf(zeta*r)/r;
    }
}

/*! \brief 
 * Compute the Coulomb overlap integral between two Gaussian distributed charges.
 * The widths zi and zj may both be zero.
 * No units are introduced or converted but nm are assumed.
 * No factor 1/4 pi epsilon0 is introduced.
 *
 * \param[in] r  distance in nm
 * \param[in] zi Gaussian width in 1/nm 
 * \param[in] zj Gaussian width in 1/nm 
 * \return    Integral value
 */    
static inline gmx_unused double Coulomb_GG(double r, double zi, double zj)
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
            zeff = zi*zj*gmx::invsqrt(sqr(zi)+sqr(zj));
        }
    }
    
    return Nuclear_GG(r, zeff);
}

/*! \brief 
 * Compute the derivative of the Coulomb overlap integral between a Gaussian 
 * distributed charge and a point charge. The width zeta may be zero.
 * No units are introduced or converted but nm are assumed.
 * No factor 1/4 pi epsilon0 is introduced.
 *
 * \param[in] r    distance in nm
 * \param[in] zeta Gaussian width in 1/nm 
 * \return    Integral value
 */    
static inline gmx_unused double DNuclear_GG(double r, double zeta)
{
    if (r == 0)
    {
        return 0;
    }
    if (zeta == 0)
    {
        return -1.0/sqr(r);
    }
    else
    {
        double rz = r * zeta;
        
        return (2.0*invsqrt_pi)*exp(-sqr(rz))*(zeta/r) - erf(rz)/sqr(r);
    }
}

/*! \brief 
 * Compute the derivative of the Coulomb overlap integral between two Gaussian 
 * distributed charges with respect to the distance.
 * The widths zi and zj may both be zero.
 * No units are introduced or converted but nm are assumed.
 * No factor 1/4 pi epsilon0 is introduced.
 *
 * \param[in] r  distance in nm
 * \param[in] zi Gaussian width in 1/nm 
 * \param[in] zj Gaussian width in 1/nm 
 * \return    Integral value
 */    
static inline gmx_unused double DCoulomb_GG(double r, double zi, double zj)
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
            zeff = zi*zj*gmx::invsqrt(sqr(zi)+sqr(zj));
        }
    }
    
    return DNuclear_GG(r, zeff);
}

#endif

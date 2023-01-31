/*
 * This source file is part of the Alexandria Chemistry Toolkit.
 *
 * Copyright (C) 2014-2023
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
#ifndef SLATERINTEGRALS_H
#define slaterINTEGRALS_H

/*! \brief 
 * Compute the Slater overlap integral between two Slater distributed charges.
 * The widths xi and xj may both be zero. The Slater function types i and j
 * should be in the range 1 .. SLATER_MAX
 * No units are introduced or converted but nm are assumed.
 * No factor 1/4 pi epsilon0 is introduced.
 *
 * \param[in] r  distance in nm
 * \param[in] i  Slater function type
 * \param[in] j  Slater function type
 * \param[in] xi gaussian width in 1/nm 
 * \param[in] xj gaussian width in 1/nm 
 * \return    Integral value
 */    
double Coulomb_SS(double r,int i,int j,double xi,double xj);

/*! \brief 
 * Compute the Slater overlap integral between a Slater distributed charge and a
 * point charge. The width xi may be zero. The Slater function type i should be in the
 * range 1 .. SLATER_MAX
 * No units are introduced or converted but nm are assumed.
 * No factor 1/4 pi epsilon0 is introduced.
 *
 * \param[in] r  distance in nm
 * \param[in] i  Slater function type
 * \param[in] xi gaussian width in 1/nm 
 * \return    Integral value
 */    
double Nuclear_SS(double r,int i,double xi);

/*! \brief 
 * Compute the derivative of the Slater overlap integral between two Slater 
 * distributed charges with respect to r.
 * The widths xi and xj may both be zero. The Slater function types i and j
 * should be in the range 1 .. SLATER_MAX
 *
 * \param[in] r  distance in nm
 * \param[in] i  Slater function type
 * \param[in] j  Slater function type
 * \param[in] xi gaussian width in 1/nm 
 * \param[in] xj gaussian width in 1/nm 
 * \return    Integral value
 */    
double DCoulomb_SS(double r,int i, int j, double xi, double xj);

/*! \brief 
 * Compute the derivative of the Slater overlap integral between a Slater 
 * distributed charge and a point charge with respect to r.
 * The width xi may be zero. The Slater function type i should be in the
 * range 1 .. SLATER_MAX
 * No units are introduced or converted but nm are assumed.
 * No factor 1/4 pi epsilon0 is introduced.
 *
 * \param[in] r  distance in nm
 * \param[in] i  Slater function type
 * \param[in] xi gaussian width in 1/nm 
 * \return    Integral value
 */    
double DNuclear_SS(double r,int i, double xi);

#endif

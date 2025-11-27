/*
 * This source file is part of the Alexandria Chemistry Toolkit.
 *
 * Copyright (C) 2014-2025
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
 * \ingroup module_alexandria
 */
#ifndef COULOMB_GAUSSIAN_H
#define COULOMB_GAUSSIAN_H

#include "act/coulombintegrals/gaussian_integrals.h"

/*! \brief Compute Coulomb interaction energy and force
 * \param[in] qq     Product of charges
 * \param[in] izeta  Distribution width for first particle
 * \param[in] jzeta  Distribution width for first particle
 * \param[in] r      Distance
 * \param[out] velec Energy
 * \param[out] felec Force
 */
static void coulomb_gaussian(real qq, real izeta, real jzeta,
                             real r, real *velec, real *felec)
{
    *velec =  qq*Coulomb_GG(r, izeta, jzeta);
    *felec = -qq*DCoulomb_GG(r, izeta, jzeta);
}

#endif

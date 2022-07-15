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

#ifndef GAUSSIAN_INTEGRALS_H
#define GAUSSIAN_INTEGRALS_H


/*! \brief Compute Coulomb interaction between Gaussian distributed charges
 * \param[in] r  Distance
 * \param[in] xi Distribution width for atom i (1/distance unit)
 * \param[in] xj Distribution width for atom j (1/distance unit)
 * \return The interaction energy (1/4 pi epsilon_0 is not included).
 */
double Coulomb_GG(double r, double xi, double xj);

/*! \brief Compute interaction between Point charge and Gaussian distributed charge
 * \param[in] r  Distance
 * \param[in] xi Distribution width for atom i (1/distance unit)
 * \return The interaction energy (1/4 pi epsilon_0 is not included).
 */
double Nuclear_GG(double r, double xi);

/*! \brief Compute Coulomb force between Gaussian distributed charges
 * \param[in] r  Distance
 * \param[in] xi Distribution width for atom i (1/distance unit)
 * \param[in] xj Distribution width for atom j (1/distance unit)
 * \return The force (1/4 pi epsilon_0 is not included).
 */
double DCoulomb_GG(double r, double xi, double xj);

/*! \brief Compute Coulomb force between Point charge and Gaussian distributed charge
 * \param[in] r  Distance
 * \param[in] xi Distribution width for atom i (1/distance unit)
 * \return The force (1/4 pi epsilon_0 is not included).
 */
double DNuclear_GG(double r, double xi);

#endif

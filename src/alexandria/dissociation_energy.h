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

#ifndef DISSOCIATION_ENERGY_H
#define DISSOCIATION_ENERGY_H

/*! \brief Compute the dissociation energies for all the bonds.
 *
 * Given all the bonds and the enthalpies of formation of all
 * molecules, we can approximate the dissociation enthalpy (D0 in the
 * Morse potential by least squares fitting the D0 to reproduce the
 * molecular energy (Delta H formation of molecule - Delta H formation of
 * the atoms). This is a crude approximation since all other energy
 * terms in the force field are ignored, however the dissociation
 * energy is the largest contribution to the molecular energy.
 * \param[in]    fplog         File pointer to write to
 * \param[inout] pd            The force field to update
 * \param[in]    mymol         The molecules
 * \param[in]    nDissociation The number of dissociation energies to compute
 */
void getDissociationEnergy(FILE                     *fplog,
                           Poldata                  *pd,
                           const std::vector<MyMol> &molset,
                           int                       nDissociation);

#endif

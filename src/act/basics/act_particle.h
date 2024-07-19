/*
 * This source file is part of the Alexandria Chemistry Toolkit.
 *
 * Copyright (C) 2024
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
#ifndef ACT_PARTICLE_H
#define ACT_PARTICLE_H

#include <string>

namespace alexandria
{
//! ActParticle types
enum class ActParticle
{
    //! Atom
    Atom,
    //! Virtual Site
    Vsite,
    //! Shell
    Shell,
    //! Sigma hole
    SigmaHole
};

/*! \brief
 * Convert ActParticle to string.
 * \param[in] apType The ActParticle type
 * \return The corresponding string
 */
const std::string &actParticleToString(ActParticle apType);

/*! \brief
 * Convert string to ActParticle
 * \param[in]  name   Name of the ActParticle
 * \param[out] apType The corresponding ActParticle type
 * \return whether or not a corresponding type was found
 */
bool stringToActParticle(const std::string &name, ActParticle *apType);

} // namespace

#endif

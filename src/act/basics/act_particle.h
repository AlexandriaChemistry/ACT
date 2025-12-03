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
//! ActParticle types used throughout the toolkit.
/*! 
 * ActParticle is an enumeration of particle "types" recognized
 * by various parts of the Alexandria Chemistry Toolkit (ACT).
 *
 * The enumeration values and their typical meanings:
 * - Atom : A physical atom (holds element, formal charge, etc.).
 * - Vsite : A virtual site (off-atom interaction site). Note the
 *           canonical string representation is "VSite" (capital S).
 * - Shell : A shell particle used in shell-model representations.
 * - SigmaHole : A sigma hole site (commonly used to model anisotropic
 *               electron density on certain atoms, e.g. halogens).
 *
 * These enum values are intended to be compact identifiers that
 * can be converted to/from their canonical string representations
 * via the helper functions declared below.
 */
enum class ActParticle
{
    //! Atom
    Atom,
    //! Virtual Site (canonical string: "VSite")
    Vsite,
    //! Shell particle
    Shell,
    //! Sigma hole particle (canonical string: "SigmaHole")
    SigmaHole
};

/*! \brief Convert ActParticle to its canonical string representation.
 *
 * The canonical strings currently used by ACT are:
 *  - ActParticle::Atom -> "Atom"
 *  - ActParticle::Vsite -> "VSite"
 *  - ActParticle::Shell -> "Shell"
 *  - ActParticle::SigmaHole -> "SigmaHole"
 *
 * \param[in] apType The ActParticle enum value to convert.
 * \return A const reference to the canonical string for the enum value.
 *
 * \note The returned reference refers to internal storage managed by
 *       the implementation (a static map). The reference remains valid
 *       for the lifetime of the program. Do not attempt to modify the
 *       returned string.
 *
 * \warning This function assumes apType is one of the defined enum
 *          values. Behavior for out-of-range/invalid enum values is
 *          unspecified.
 *
 * \example
 * \code
 * ActParticle p = ActParticle::Vsite;
 * std::string s = actParticleToString(p); // s == "VSite"
 * \endcode
 */
const std::string &actParticleToString(ActParticle apType);

/*! \brief Parse a canonical ActParticle string into the enum value.
 *
 * The function performs an exact (case-sensitive) comparison against
 * the canonical strings documented above. If a matching canonical name
 * is found, the corresponding enum value is written to \p apType and
 * the function returns true.
 *
 * \param[in]  name   Input string containing the canonical particle name.
 *                    Comparisons are case-sensitive; for example
 *                    "VSite" matches, but "vsite" does not.
 * \param[out] apType Pointer to an ActParticle variable that will be set
 *                    to the corresponding enum value when a match is found.
 * \return true if a corresponding enum value was found and written to \p apType,
 *         false otherwise.
 *
 * \note To obtain the canonical strings, use actParticleToString().
 *
 * \example
 * \code
 * ActParticle p;
 * bool ok = stringToActParticle("SigmaHole", &p); // ok == true, p == ActParticle::SigmaHole
 * ok = stringToActParticle("sigmahole", &p);     // ok == false (case mismatch)
 * \endcode
 */
bool stringToActParticle(const std::string &name, ActParticle *apType);

} // namespace alexandria

#endif
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
#ifndef POTENTIAL_H
#define POTENTIAL_H

#include <string>

#include "act/basics/chargemodel.h"

namespace alexandria
{
    enum class Potential
    { 
        NONE,
            LJ8_6, LJ12_6, LJ14_7,
            GENERALIZED_BUCKINGHAM, WANG_BUCKINGHAM,
            EXPONENTIAL,
            COULOMB_POINT, COULOMB_GAUSSIAN, COULOMB_SLATER,
            HARMONIC_BONDS, MORSE_BONDS, CUBIC_BONDS,
            HARMONIC_ANGLES, LINEAR_ANGLES, UREY_BRADLEY_ANGLES,
            HARMONIC_DIHEDRALS, FOURIER_DIHEDRALS, PROPER_DIHEDRALS,
            POLARIZATION,
            VSITE2, VSITE2FD, VSITE3, VSITE3FD, VSITE3FAD, VSITE3OUT, VSITE3OUTS
            };

    /*! \brief Convert potential into string
     * \param[in] p The potential type
     * \return A string
     */
    const std::string &potentialToString(Potential p);

    /*!\brief Convert a string to a Potential type
     * \param[in]  pname The string
     * \param[out] p     The potential
     * \return true if successfull, false if the string did not match a known potential
     */
    bool stringToPotential(const std::string &pname, Potential *p);

    /*! \brief Backward compatability routine
     * \param[in] p  The potential
     * \return The gromacs type or -1 if not found
     */

    int potentialToGromacsType(Potential p);
    /*! \brief Convenience function
     * \param[in] p  The potential
     * \return The gromacs name or nulltpr if not found
     */
    const char *potentialToGromacsString(Potential p);

    /*! \brief Convenience function
     * \param[in] c The charge type
     * \return Corresponding Coulomb potential
     */
    Potential chargeTypeToPotential(ChargeType c);

    /*! \brief Convenience function
     * \param[in] p The potential
     * \return Corresponding charge type
     */
    ChargeType potentialToChargeType(Potential p);
} // namespace alexandria

#endif

/*
 * This source file is part of the Alexandria Chemistry Toolkit.
 *
 * Copyright (C) 2024,2025
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

//! \brief Interaction potentials implemented in the ACT
enum class Potential
    { 
        NONE,
        LJ8_6, LJ12_6, LJ14_7, LJ12_6_4,
        GENERALIZED_BUCKINGHAM, WANG_BUCKINGHAM, SLATER_ISA, SLATER_ISA_TT,
        BORN_MAYER, MACDANIEL_SCHMIDT, BUCKINGHAM, TANG_TOENNIES, TT2b,
        COULOMB_POINT, COULOMB_GAUSSIAN, COULOMB_SLATER,
        HARMONIC_BONDS, MORSE_BONDS, CUBIC_BONDS, HUA_BONDS,
        HARMONIC_ANGLES, LINEAR_ANGLES, UREY_BRADLEY_ANGLES,
        HARMONIC_DIHEDRALS, FOURIER_DIHEDRALS, PROPER_DIHEDRALS,
        POLARIZATION,
        VSITE1, VSITE2, VSITE2FD, VSITE3, VSITE3S, VSITE3FD, VSITE3FAD, VSITE3OUT, VSITE3OUTS, VSITE4, VSITE4S, VSITE4S3
    };

//! \brief Structure containing information about Potential functions
typedef struct {
    //! \brief String corresponding to potential function
    const std::string         name;
    //! \brief Corresponding type in GROMACS (if present)
    int                       ftype;
    //! \brief List of parameter names
    std::vector<const char *> param;
    //! \brief Energy expression used to populate OpenMM force field files
    const std::string         energy;
    //! \brief A prefactor for implementing combination rules in OpenMM force field files
    const std::string         prefactor;
} PotentialProperties;

/*! \brief Map from Potential to it's PotentialProperties
 */
extern std::map<Potential, PotentialProperties> potprops;

/*! \brief Convert potential into string
 * \param[in] p The potential type
 * \return A string
 */
const std::string &potentialToString(Potential p);

/*! \brief Convert potential into energy expression
 * \param[in] p The potential type
 * \return A string
 */
const std::string &potentialToEnergy(Potential p);

/*! \brief Convert potential into prefactors for the energy expression
 * \param[in] p The potential type
 * \return A string
 */
const std::string &potentialToPreFactor(Potential p);

/*! \brief Get parameter names for potential
 * \param[in] p The potential type
 * \return A vector of strings
 */
const std::vector<const char *> potentialToParameterName(Potential p);

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
Potential chargeDistributionTypeToPotential(ChargeDistributionType c);

/*! \brief Convenience function
 * \param[in] p The potential
 * \return Corresponding charge type
 */
ChargeDistributionType potentialToChargeDistributionType(Potential p);

} // namespace alexandria

#endif

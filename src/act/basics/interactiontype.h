/*
 * This source file is part of the Alexandria Chemistry Toolkit.
 *
 * Copyright (C) 2014-2020 
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
#ifndef INTERACTIONTYPE_H
#define INTERACTIONTYPE_H

#include <string>

namespace alexandria
{
//! Interaction type
enum class InteractionType
{
    //! Chemical bonds
    BONDS,
    //! Angles between two bonds
    ANGLES,
    //! Atom triplets i,j,k that are in a linear conformation
    LINEAR_ANGLES,
    //! Atom quadruplets i,j,k,l. The dihedral angle is the angle between planes made up by atoms i,j,k and j,k,l.
    PROPER_DIHEDRALS,
    //! Out-of-plane dihedral angles.
    IMPROPER_DIHEDRALS,
    //! Van der Waals interaction between atoms
    VDW,
    //! Charge interactions
    COULOMB,
    //! Polarization interaction between an atom (core) or a virtual site and its directly connected shell.
    POLARIZATION,
    //! Potential energy, sum over the above
    EPOT,
    //! Constrained bonds, that are not allowed to fluctuate.
    CONSTR,
    //! Virtual interaction sites along the bond between two atoms.
    VSITE2,
    //! Virtual interaction sites determined by three atoms in a plane
    VSITE3FAD,
    //! Virtual interaction sites determined by three atoms out of the plane
    VSITE3OUT,
    //! Bond hardness and electronegativity needed to use the split charge equilibration algorithm
    BONDCORRECTIONS,
    //! Correction to the electronegativity difference between two atoms connected by a bond
    ELECTRONEGATIVITYEQUALIZATION,
    //! Abusing this structure for something that is a charge rathers than an interaction type.
    CHARGE
};

/*! \brief
 * Convert interaction type to string.
 * \param[in] iType The interaction type
 * \return The corresponding string
 */
const std::string &interactionTypeToString(InteractionType iType);

/*! \brief
 * Convert interaction type to descriptive string rather than 
 * what is force field files.
 * \param[in] iType The interaction type
 * \return The corresponding string
 */
const std::string &interactionTypeToDescription(InteractionType iType);

/*! \brief
 * Convert string to interaction type.
 * \param[in] name Name of the interaction
 * \return The corresponding interaction type
 * \throws if there is no corresponding interaction type
 */
InteractionType stringToInteractionType(const std::string &name);

/*! \brief
 * Return number of atoms involved with this interaction type
 * \param[in] iType The interactionType
 * \return number of atoms typically 1-4.
 */
int interactionTypeToNatoms(InteractionType iType);

} // namespace

#endif

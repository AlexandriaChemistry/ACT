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
#ifndef FUNCTIONTYPE_H
#define FUNCTIONTYPE_H

#include <string>

namespace alexandria
{

enum class FunctionType
{
 // Potential functions for describing FunctionType::BOND
 
 //! Harmonic function
 HARMONIC,
    
 //! Morse function
 MORSE,

 // Potential functions for describing FunctionType::ANGLE
 
 //!
 UREY_BRADLEY,

 //!
 COSINE,


 // Potential functions for describing FunctionType::LINEAR_ANGLES
 LINEAR_ANGLE,


 // Potential functions for describing FunctionType::IMPROPER_DIHEDRAL
 IDIHS,

 // Potential functions for describing FunctionType::PROPER_DIHEDRAL

 //!
 PDIHS,

 //!
 FOURDIHS,

 // Potential functions for describing FunctionType::COULOMB
 COUL_SR,
 
 //! Dirac delta function
 DIRAC_DELTA,

 //! Gaussain distribution function
 GAUSSIAN,

 //! Slater distribution function
 SLATER,

 // Potential functions for describing FunctionType::VDW

 //! Lennard-Jones function
 LENNARD_JONES,

 //! Buckingham function
 BUCKINGHAM,

 //! Wang-Buckingham function
 WANG_BUCKINGHAM
};


/*! \brief
 * Convert function type to string.
 * \param[in] fType The function type
 * \return The corresponding string
 */
const std::string &functionTypeToString(FunctionType fType);

/*! \brief
 * Convert function type to descriptive string rather than 
 * what is force field files.
 * \param[in] fType The function type
 * \return The corresponding string
 */
const std::string &functionTypeToDescription(FunctionType fType);

/*! \brief
 * Convert string to function type.
 * \param[in] name Name of the function
 * \return The corresponding function type
 * \throws if there is no corresponding function type
 */
FunctionType stringToFunctionType(const std::string &name);

/*! \brief
 * Return number of parameters s involved with this function type
 * \param[in] fType The functionType
 * \return number of parameters.
 */
int functionTypeToNparams(FunctionType fType);


} // namespace alexandria

#endif


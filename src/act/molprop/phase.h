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
 
#ifndef PHASE_H
#define PHASE_H
 
#include <string>
    
namespace alexandria
{

//! Enum to describe the phase corresponding to a property
enum class ePhase {
    //! The gas phase
    GAS,
    //! The liquid pahse
    LIQUID,
    //! The solid phase
    SOLID,
    //! Plasma!
    PLASMA
};

/*! \brief
 * Yield string corresponding to phase
 * 
 * \param[in] ep The phase enum
 * \return string corresponding to phase
 */
const std::string &phase2string(ePhase ep);

/*! \brief
 * Yield phase corresponding to string
 * 
 * \param[in] phase String corresponding to phase
 * \return The phase enum
 */
ePhase string2phase(const std::string &phase);

} // namespace alexandria

#endif

/*
 * This source file is part of the Alexandria Chemistry Toolkit.
 *
 * Copyright (C) 2020 
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
 * \author David van der Spoel <david.vanderspoel@icm.uu.se>
 */

#ifndef MUTABILITY_H
#define MUTABILITY_H

#include <string>

namespace alexandria
{

//! \brief Enum determining whether a parameter can be changed
enum class Mutability
{
    //! Parameter that cannot be changed
    Fixed,
    //! Parameter that is dependent on another parameter and should not be changed independently
    Dependent,
    //! Parameter (charge) that should be changed by the Alexandria Charge Method algorithms only
    ACM,
    //! Parameter that can be modified within bounds
    Bounded,
    //! Parameter that can be modified without restrictions
    Free
};

/*! \brief Return a string corresponding to the mutability
 */
const std::string &mutabilityName(Mutability mutability);

/*! \brief Lookup a string and return mutability value
 *
 * \param[in]  name       String
 * \param[out] mutability Point to Mutability encoded in the string
 * \return true if the name could be converted to a Mutability, false otherwise
 */
bool nameToMutability(const std::string &name, Mutability *mutability);

} // namespace alexandria

#endif

/*
 * This source file is part of the Alexandria Chemistry Toolkit.
 *
 * Copyright (C) 2023
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
#ifndef ACT_BASIC_ATOMPROPS_H
#define ACT_BASIC_ATOMPROPS_H

#include <map>
#include <string>

namespace alexandria
{

class AtomProp
{
private:
    //! Element name in full
    std::string name_;
    //! Atomic number
    int         atomnumber_;
    //! Atomic mass
    double      mass_;
    //! Formal charge
    int        charge_;
    //! Multiplicity at lowest energy
    int        mult_;
public:
    //! Constructor
    AtomProp(std::string name, int atomnumber, double mass, int charge, int mult) :
        name_(name), atomnumber_(atomnumber), mass_(mass), charge_(charge), mult_(mult) {};
    
    //! \return full name
    const std::string &name() const { return name_; }
    
    //! \return the atomic number
    int atomnumber() const { return atomnumber_; }
    
    //! \return the atomic mass
    double mass() const { return mass_; }

    //! \return the formal charge
    int charge() const { return charge_; }

    //! \return the multiplicity at lowest energy
    int mult() const { return mult_; }
};

/*! \brief Return a table of atom properties with the element as the key.
 */ 
std::map<std::string, AtomProp> readAtomProps();

} // namespace alexandria

#endif

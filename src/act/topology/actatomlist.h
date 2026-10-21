/*
 * This source file is part of the Alexandria Chemistry Toolkit.
 *
 * Copyright (C) 2021-2026
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
 * \author Mohammad Mehdi Ghahremanpour <mohammad.ghahremanpour@icm.uu.se>
 * \author David van der Spoel <david.vanderspoel@icm.uu.se>
 */
#ifndef ACT_ACTATOMLIST_H
#define ACT_ACTATOMLIST_H

#include <list>

#include "act/topology/actatom.h"
#include "gromacs/math/vectypes.h"

namespace alexandria
{

class ActAtomListItem
{
private:
    //! The atom information
    ActAtom   atom_;
    //! The original index
    size_t    index_ = 0;
    //! Atomic coordinates
    gmx::RVec x_     = { 0, 0, 0 };
public:
    //! constructor
    ActAtomListItem(const ActAtom   atom,
                    size_t          index,
                    const gmx::RVec x) : atom_(atom), index_(index), x_(x) {}

    //! \return the atom
    const ActAtom atom() const { return atom_; }

    //! \return pointer to the atom for editing
    ActAtom *atomPtr() { return &atom_; }

    //! \return the index
    size_t index() const { return index_; }

    //! \return the coordinates
    const gmx::RVec x() const { return x_; }

    //! Equal index
    bool operator==(size_t index) const { return index == index_; }
};

//! Shortcut for passing this around.
typedef std::list<ActAtomListItem> AtomList;

} // namespace alexandria

#endif

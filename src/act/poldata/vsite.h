﻿/*
 * This source file is part of the Alexandria Chemistry Toolkit.
 *
 * Copyright (C) 2014-2020
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
 
#ifndef ACT_VSITE_H
#define ACT_VSITE_H

#include "gromacs/mdtypes/commrec.h"

#include <string>
#include <vector>
 
#include "act/utility/communicationrecord.h"

namespace alexandria
{

//! Class to distguish virtual site types
enum class VsiteType
{
    //! Linear virtual site (gv_linear) dependent on two atoms
    LINEAR,
    //! Planar virtual site (gv_planar) dependent on three atoms
    PLANAR,
    //! Planar virtual site (gv_ringplanar) dependent on three atoms in a ring
    RING_PLANAR,
    //! Planar virtual site (gv_inplane) dependent on three atoms
    IN_PLANE,
    //! Out-of-plane virtual site (gv_outplane) dependent on three atoms
    OUT_OF_PLANE,
    //! All possible virtual sites
    ALL         
};

/*! \brief
 * Convert virtual site type to string
 * \param[in] vType The virtual site type
 * \return VsiteType name
 */
const char *vsiteType2string(VsiteType vType);

/*! \brief
 * Convert string to virtual site type
 * \param[in] string The string to convert
 * \returns The vsite type
 * \throws if the string is not recognized
 */
VsiteType string2vsiteType(const char *string);

class Vsite
{
    public:

        Vsite () {};

        Vsite(const std::string &atype,
              const std::string &type,
              int                number,
              double             distance,
              double             angle,
              int                ncontrolatoms);

        /*! \brief
         * Return the atom type to which the vsites are connected.
         */
        const std::string &atype() const { return atype_; }

        /*! \brief
         * Return the type pf virtual site.
         */
        const VsiteType &type() const { return type_; }

        /*! \brief
         * Return the distance between the atom and the virtual site.
         */
        double distance() const {return distance_; }

        /*! \brief
         * Return the angle between the atom and the virtual site.
         */
        double angle() const {return angle_; }

        /*! \brief
         * Return the number of virtual sites connected to the atom.
         */
        int nvsite() const {return number_; }

        /*! \brief
         * Return the number of atoms needed to locate the vsite.
         */
        int ncontrolatoms() const {return ncontrolatoms_; }

        CommunicationStatus Send(const CommunicationRecord *cr, int dest);

        CommunicationStatus Bcast(const CommunicationRecord *cr);

        CommunicationStatus Receive(const CommunicationRecord *cr, int src);


    private:
        std::string atype_;
        VsiteType   type_;
        int         number_;
        double      distance_;
        double      angle_;
        int         ncontrolatoms_;
};

using VsiteIterator      = typename std::vector<Vsite>::iterator;
using VsiteConstIterator = typename std::vector<Vsite>::const_iterator;

} // namespace alexandria

#endif

/*
 * This source file is part of the Alexandria Chemistry Toolkit.
 *
 * Copyright (C) 2014-2025
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
#ifndef SYMCHARGES_H
#define SYMCHARGES_H

#include <string>

#include "act/utility/communicationrecord.h"

namespace alexandria
{

/*! \brief Rule for symmetrizing charges
 * Atoms of element attached_ bonded to an atom of element central_
 * get equal charges if there are exactly numattach_ of them,
 * e.g. the hydrogens in a methyl group.
 */
class Symcharges
{
    public:

        //! Default constructor
        Symcharges () {}

        /*! \brief Constructor
         * \param[in] central   Element of the central atom
         * \param[in] attached  Element of the attached atoms
         * \param[in] numattach Number of attached atoms
         */
        Symcharges(const std::string &central,
                   const std::string &attached,
                   int                numattach);

        //! \return the element of the central atom
        const std::string &getCentral() const { return central_; }

        //! \return the element of the attached atoms
        const std::string &getAttached() const { return attached_; }

        //! \return the number of attached atoms
        int getNumattach() const { return numattach_; }

        /*! \brief Send the contents to another processor
         * \param[in] cr   Communication record
         * \param[in] dest Destination processor
         * \return the communication status
         */
        CommunicationStatus Send(const CommunicationRecord *cr, int dest);

        /*! \brief Broadcast the contents to other processors
         * \param[in] cr   Communication record
         * \param[in] root Root processor of the broadcast (unused)
         * \param[in] comm MPI communicator
         * \return the communication status
         */
        CommunicationStatus BroadCast(const CommunicationRecord *cr, int root, MPI_Comm comm);

        /*! \brief Receive the contents from another processor
         * \param[in] cr  Communication record
         * \param[in] src Source processor
         * \return the communication status
         */
        CommunicationStatus Receive(const CommunicationRecord *cr, int src);

    private:
        //! Element of the central atom
        std::string central_;
        //! Element of the attached atoms
        std::string attached_;
        //! Number of attached atoms
        int         numattach_;
};

//!\brief  Shortcut for iterator over Symcharges
using SymchargesIterator      = typename std::vector<Symcharges>::iterator;
//!\brief  Shortcut for const_iterator over Symcharges
using SymchargesConstIterator = typename std::vector<Symcharges>::const_iterator;

/*! \brief Utility to make command line information about combrules
 * \param[inout] crinfo Array of strings to be edited
 */
void add_comb_rule_info(std::vector<const char *> *crinfo);

} // namespace aleaxndria

#endif

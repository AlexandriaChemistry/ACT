/*
 * This source file is part of the Alexandria Chemistry Toolkit.
 *
 * Copyright (C) 2020-2024
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

#ifndef IDENTIFIER_H
#define IDENTIFIER_H

#include <string>
#include <vector>

#include "act/utility/communicationrecord.h"
#include "interactiontype.h"

namespace alexandria
{
    //! \brief enum describing whether and how atom order swapping is allowed in interactions
    enum class CanSwap {
        //! The order of atoms in an interaction cannot be swapped
        No,
        //! The order of atoms in an interaction can be swapped by reversing the order
        Yes,
        //! The order of atoms can be swapped specifically for improper dihedrals
        Idih,
        //! Special case for linear angles
        Linear,
        //! Special case for vsite2
        Vsite2,
        //! Special case for vsite3
        Vsite3,
    };

/*! \brief Convert string to CanSwap
 * \param[in] str The string
 * \return  the enum value
 * \throws a gmx::InvalidInputError in case of an unknown string
 * TODO: return a boolean instead of crashing
 */
CanSwap stringToCanSwap(const std::string &str);

/*! \brief Convert CanSwap to string
 * \param[in]  canSwap The enum value
 * \return The string
 */
const std::string &canSwapToString(CanSwap canSwap);

/*! \brief Class to hold an identifier
 */
class Identifier
{
 public:
    /*! \brief Empty constructor
     */
    Identifier() {}

    /*! \brief Constructor
     * Will extract atoms and bondorders from the id string. If the string
     * is incorrect the program will be terminated with a message.
     * \param[in] iType   The interaction type for this identifier
     * \param[in] id      String containing atom/bond names and bond orders
     * \param[in] canSwap Can the order of the atoms be swapped
     * \throws with a gmx::InvalidInputError if the input is not consistent
     */
    Identifier(InteractionType    iType,
               const std::string &id,
               CanSwap            canSwap);

    /*! \brief Constructor
     * \param[in] atoms      Vector containing atom/bond names
     * \param[in] bondorders The bond order in case of bonds
     * \param[in] canSwap    Can the order of the atoms be swapped
     */
    Identifier(const std::vector<std::string> &atoms,
               const std::vector<double>      &bondorders,
               CanSwap                         canSwap);

    /*! \brief Constructor for a single atom
     * \param[in] atom  Atom name
     */
    Identifier(const std::string &atom);

    /*! \brief Return standard identifier string
     * \throws if there are no ids
     */
    const std::string &id() const;

    //! \return whether a swapped id exists
    bool haveSwapped() const { return ids_.size() > 1; }

    /*! \brief Return swapped identifier string
     * \throws if there are no ids
     */
    const std::string &swapped() const;

    //! \return swappability of this identifier
    CanSwap canSwap() const { return canSwap_; }

    //! \brief Return the bond orders
    const std::vector<double> &bondOrders() const { return bondOrders_; }

    /*! \brief Comparison operator
     *
     * Either both identifiers are swappable or neither.
     * \param a An identifier
     * \param b Another identifier
     * \return Whether a < b
     */
    friend bool operator<(const Identifier &a, const Identifier &b)
    {
        // TODO check implementation
        if (a.canSwap() < b.canSwap())
        {
            return true;
        }
        else if (a.canSwap() > b.canSwap())
        {
            return false;
        }
        else if (a.atoms_.size() < b.atoms_.size())
        {
            return true;
        }
        else if (a.atoms_.size() > b.atoms_.size())
        {
            return false;
        }
        else if (a.ids_[0] == b.ids_[0])
        {
            // Identical!
            return false;
        }
        else if (a.canSwap() != CanSwap::No)
        {
            for(auto &id : a.ids_)
            {
                if (id == b.id())
                {
                    // Identical as well
                    return false;
                }
            }
        }
        return a.ids_[0] < b.ids_[0];
    }

    /*! \brief Comparison operator
     *
     * Either both identifiers are swappable or neither.
     * \param a An identifier
     * \param b Another identifier
     * \return Whether a == b
     */
    friend bool operator==(const Identifier &a, const Identifier &b)
    {
        return !(a < b) && !(b < a);
    }

    //! \brief Return the atoms, or rather the components of the id
    const std::vector<std::string> &atoms() const { return atoms_; }

    /*! \brief Send the contents to another processor
     * \param[in] cr   Communication data structure
     * \param[in] dest Processor id to send the data to
     */
    CommunicationStatus Send(const CommunicationRecord *cr, int dest)  const;

    /*! \brief Broadcast content to and from other processors
     * \param[in] cr   Communication data structure
     * \param[in] root The root of this communication
     * \param[in] comm MPI communicator
     */
    CommunicationStatus BroadCast(const CommunicationRecord *cr,
                                  int                        root,
                                  MPI_Comm                   comm);

    /*! \brief Receive contents from another processor
     * \param[in] cr  Communication data structure
     * \param[in] src Processor id to receive the data from
     */
    CommunicationStatus Receive(const CommunicationRecord *cr, int src);

 private:
    //! Whether we can swap and if so, how
    CanSwap                  canSwap_ = CanSwap::No;
    //! The id in short of these parameter, with alternate ids
    std::vector<std::string> ids_;
    //! The bond orders
    std::vector<double>      bondOrders_;
    //! The atoms
    std::vector<std::string> atoms_;

    //! Correct the order of atoms and created permuted ids
    void orderAtoms();

    void update();
};

} // namespace alexandria

#endif

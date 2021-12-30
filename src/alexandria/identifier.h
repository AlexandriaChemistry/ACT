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

#ifndef IDENTIFIER_H
#define IDENTIFIER_H

#include <string>
#include <vector>

#include "gromacs/mdtypes/commrec.h"

#include "communication.h"
#include "interactiontype.h"

namespace alexandria
{

enum class CanSwap {
    //! The order of atoms in an interaction cannot be swapped
    No,
    //! The order of atoms in an interaction can be swapped
    Yes
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

    /*!
     * Copy constructor
     * \param[in] other   reference Identifier object
     */
    Identifier(const Identifier &other)
    : id_(other.id()), swappedId_(other.swappedId()),
      bondOrders_(other.bondOrders()), atoms_(other.atoms()) {}

    /*! \brief Simple constructor
     * \param[in] atoms   Vector containing atom/bond names
     * \param[in] canSwap Can the order of the atoms be swapped
     */
    Identifier(const std::vector<std::string> &atoms,
               CanSwap                         canSwap);

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
     * \param[in] atoms     Vector containing atom/bond names
     * \param[in] bondorder The bond order in case of bonds
     * \param[in] canSwap   Can the order of the atoms be swapped
     */
    Identifier(const std::vector<std::string> &atoms,
               const std::vector<double>      &bondorders,
               CanSwap                         canSwap);

    /*! \brief Constructor for a single atom
     * \param[in] atom  Atom name
     */
    Identifier(const std::string &atom);

    //! \brief Return standard identifier string
    const std::string &id() const { return id_; }

    //! \brief Return swapped identifier string
    const std::string &swappedId() const { return swappedId_; }

    //! \brief Return the bond orders
    const std::vector<double> &bondOrders() const { return bondOrders_; }

    /*! \brief Comparison operator
     *
     * Either both identifiers are swappable or neither.
     * \param a An identifier
     * \param b Another identifier
     * \return Whether a < b
     */
    friend bool operator<(const Identifier &a, const Identifier &b);

    /*! \brief Comparison operator
     *
     * Either both identifiers are swappable or neither.
     * \param a An identifier
     * \param b Another identifier
     * \return Whether a == b
     */
    friend bool operator==(const Identifier &a, const Identifier &b);

    //! \brief Return the atoms
    const std::vector<std::string> &atoms() const { return atoms_; }

    /*! \brief Send the contents to another processor
     * \param[in] cr   Communication data structure
     * \param[in] dest Processor id to send the data to
     */
    CommunicationStatus Send(const t_commrec *cr, int dest)  const;

    /*! \brief Receive contents from another processor
     * \param[in] cr  Communication data structure
     * \param[in] src Processor id to receive the data from
     */
    CommunicationStatus Receive(const t_commrec *cr, int src);

 private:

    //! The id in short of these parameter
    std::string              id_;
    //! The swapped id
    std::string              swappedId_;
    //! The bond orders
    std::vector<double>      bondOrders_;
    //! The atoms
    std::vector<std::string> atoms_;

    //! Correct the order of atoms
    void orderAtoms();

    /*! Create the swapped version of the id if needed
     * \param[in] canSwap What we are allowed to do in this case
     */
    void createSwapped(CanSwap canSwap);
};

} // namespace alexandria

#endif

/*
 * This source file is part of the Alexandria Chemistry Toolkit.
 *
 * Copyright (C) 2022
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

#ifndef FRAGMENT_H
#define FRAGMENT_H

#include <vector>

#include "act/utility/communicationrecord.h"

namespace alexandria
{

/*! \brief
 * Contains fragment properties
 *
 * \inpublicapi
 * \ingroup module_alexandria
 */
class Fragment
{
 private:
    //! Mass of this fragment
    double           mass_         = 0.0;
    //! Charge of this fragment
    int              charge_       = 0;
    //! Multiplicity of this fragment
    int              multiplicity_ = 1;
    //! The atom indices starting at zero in this fragment
    std::vector<int> atoms_;
    //! String containing the atom indices
    std::string      atomString_;
    
    void makeAtomString();
 public:
    //! Default constructor
    Fragment() {}
    
    /*! \brief Constructor with initiation
     */
    Fragment(double                  mass,
             int                     charge,
             int                     multiplicity,
             const std::vector<int> &atoms) : mass_(mass), charge_(charge),
        multiplicity_(multiplicity), atoms_(atoms)
    {
        makeAtomString();
    }
       
    //! Return the mass
    double mass() const { return mass_; }
    
    /*! Set the mass
     * \param[in] m The new mass
     */
    void setMass(double m) { mass_ = m; }
    
    //! Return the charge
    int charge() const { return charge_; }
    
    //! Return the multiplicity
    int multiplicity() const { return multiplicity_; }
    
    //! Set the atoms
    void setAtoms(const std::vector<int> &atoms)
    {
        atoms_ = atoms;
        makeAtomString();
    }
    
    //! Return the list of atoms (numbers starting from zero)
    const std::vector<int> &atoms() const { return atoms_; }
    
    //! Return a string containing the atom numbers
    const std::string &atomString() const { return atomString_; }
    
    /*! \brief
     * Sends this object over an MPI connection
     *
     * \param[in] cr   GROMACS data structure for MPI communication
     * \param[in] dest Destination processor
     * \return the CommunicationStatus of the operation
     */
    CommunicationStatus Send(const CommunicationRecord *cr,
                             int                        dest) const;
    
    /*! \brief
     * Receives this object over an MPI connection
     *
     * \param[in] cr  GROMACS data structure for MPI communication
     * \param[in] src Source processor
     * \return the CommunicationStatus of the operation
     */
    CommunicationStatus Receive(const CommunicationRecord *cr,
                                int                        src);

};

} // namespace alexandria

#endif

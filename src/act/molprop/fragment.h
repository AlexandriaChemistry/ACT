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

#include <cstdio>
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
    //! Identifier for this fragment
    std::string      id_;
    //! Mass of this fragment
    double           mass_           = 0.0;
    //! Charge of this fragment
    int              charge_         = 0;
    //! Multiplicity of this fragment
    int              multiplicity_   = 1;
    //! Symmetry number of this fragment
    int              symmetryNumber_ = 1;
    //! The formula of this fragment
    std::string      formula_;
    //! The formula of this fragment formatted for LaTeX
    std::string      texform_;
    //! The atom indices starting at 0 in this fragment
    std::vector<int> atoms_;
    //! String containing the atom indices, starting from 1 for storing in XML files.
    std::string      atomString_;
    
    //! Convert atom numbers to a string
    void makeAtomString();
    
    //! Convert formula to a LaTeX compatible one.
    void makeTexFormula();

 public:
    //! Default constructor
    Fragment() {}
    
    /*! \brief Constructor with initiation
     */
    Fragment(const std::string      &id,
             double                  mass,
             int                     charge,
             int                     multiplicity,
             int                     symmetryNumber,
             const std::string      &formula,
             const std::vector<int> &atoms) : 
        id_(id), mass_(mass), charge_(charge),
        multiplicity_(multiplicity),
        formula_(formula), atoms_(atoms)
    {
        makeAtomString();
        makeTexFormula();
        symmetryNumber_ = std::max(1, symmetryNumber);
    }
       
    //! Return the id
    const std::string &id() const { return id_; }
    
    /*! Set the identifier
     * \param[in] id The new identifier
     */
    void setId(const std::string &id) { id_ = id; }
    
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
    
    //! Return the symmetry number
    int symmetryNumber() const { return symmetryNumber_; } 
    
    //! Return the formula
    const std::string &formula() const { return formula_; }

    //! Return the formula
    void setFormula(const std::string &formula)
    { 
        formula_ = formula;
        makeTexFormula();
    }
    
    //! Return the formula for LaTeX
    const std::string &texFormula() const { return texform_; }

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
    
    /*! Write the content of this fragment to a file
     * \param[in] fp The file pointer, if nullptr function will do nothing
     */
    void dump(FILE *fp) const;
    
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

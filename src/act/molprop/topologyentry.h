/*
 * This source file is part of the Alexandria Chemistry Toolkit.
 *
 * Copyright (C) 2021-2023
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

#ifndef ACT_MOLPROP_TOPOLOGYENTRY_H
#define ACT_MOLPROP_TOPOLOGYENTRY_H

#include <cstdio>
#include <any>
#include <map>
#include <vector>

#include "act/basics/basecontainer.h"
#include "act/basics/identifier.h"
#include "act/forcefield/particletype.h"
#include "act/forcefield/forcefield.h"

namespace alexandria
{

/*! \brief
 * Base class for topology components
 */
class TopologyEntry
{
private:
    //! Atom indices
    std::vector<int>                            indices_;
    //! Parameters belonging to this interaction
    std::vector<double>                         params_;
    //! The bond orders
    std::vector<double>                         bondOrder_;
    //! The force field identifier belonging to this entry
    Identifier                                  id_;
    //! The gromacs topology index (in idef.ffparams)
    int                                         gromacsType_ = -1;
public:
    //! Default empty constructor
    TopologyEntry() {}

    //! Return myself or an inherited class in const form
    TopologyEntry *self() { return this; }
    
    const TopologyEntry *self() const { return this; }
        
    //! Return myself or an inherited class in non-const form
    TopologyEntry *selfPtr() { return this; }
    
    /*! \brief Add one atom index
     * \param[ai] The atom index
     */
    void addAtom(int ai) { indices_.push_back(ai); }

    /*! Renumber atoms if needed
     * If shells are added to a molecule, additional particles
     * are inserted and atoms making up bonds have to be
     * renumbered.
     * \param[in] renumber The mapping array
     */
    void renumberAtoms(const std::vector<int> &renumber);

    /*! Returns an atom id
     * \param[in] ai The index number
     */
    int atomIndex(size_t ai) const
    {
        if (ai >= indices_.size())
        {
            GMX_THROW(gmx::InternalError(gmx::formatString("Atom index %ld out of range, should be <= %ld", 
                                                           ai, indices_.size()).c_str()));
        }
        return indices_[ai];
    }

    //! \return atom index array
    const std::vector<int> &atomIndices() const { return indices_; }

    //! \return the bondorders
    const std::vector<double> &bondOrders() const { return bondOrder_; }

    //! Add a bond order
    void addBondOrder(double bo) { bondOrder_.push_back(bo); }

    /*! Set a specific bond order
     * \param[in] ai The index in the array
     * \param[in] bo The new bond order
     * \throws if ai is outside the range of the vector length
     */
    void setBondOrder(size_t ai, double bo);

    double bondOrder(size_t ai) const { return bondOrder_[ai]; }

    double bondOrder() const { return bondOrder_[0]; }

    //! Return the gromacs type
    int gromacsType() const { return gromacsType_; }

    //! Set the gromacs type
    void setGromacsType(int gromacsType) { gromacsType_ = gromacsType; }

    /*! \brief Set the identifier(s)
     * \param[in] id    The identifier(s)
     */
    void setId(const Identifier &id) { id_ = id; }
    
    /*! \brief Set the parameters
     * \param[in] param The parameter(s)
     */
    void setParams(const std::vector<double> &param) { params_ = param; }

    const std::vector<double> &params() const { return params_; }
    
    //! Return my identifier(s)
    const Identifier &id() const { return id_; }

    /*! \brief Check whether this entry has nAtom atoms
     * \param[in] nAtom the expected number of atoms
     * \throws if the number of atoms does not match the expectation
     */
    void check(size_t nAtom) const;
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
     * \param[in] cr   Data structure for MPI communication
     * \param[in] root The MPI root
     * \param[in] comm MPI communicator
     * \return the CommunicationStatus of the operation
     */
    CommunicationStatus BroadCast(const CommunicationRecord *cr,
                                  int                        root,
                                  MPI_Comm                   comm);

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

//! Vector of TopologyEntry items or inherited variants
typedef std::vector<BaseContainer<TopologyEntry>> TopologyEntryVector;

/*! \brief
 * Atom pair in a molecule
 *
 * \inpublicapi
 * \ingroup module_alexandria
 */
class AtomPair : public TopologyEntry
{
 public:
    //! Default constructor
    AtomPair() {}
    
    //! Constructor setting the ids of the atoms
    AtomPair(int ai, int aj) { Set(ai, aj); }

    //! Sets the ids of the atoms
    void Set(int ai, int aj)
    {
        addAtom(ai);
        addAtom(aj);
        // TODO: This is a fake bond order.
        addBondOrder(1);
    }
    
    //! Returns the ids of the atoms and the bondorder
    void get(int *ai, int *aj) const;
    
    //! Returns the first atom id
    int aI() const
    {
        return atomIndex(0);
    }
      
     //! Returns the second atom id
    int aJ() const
    {
        return atomIndex(1);
    }
    
    //! Return an AtomPair with the order of atoms swapped
    AtomPair swap() const;
    
    /*! \brief Return whether two AtomPairs are the same
     * \param[in] other The other AtomPair
     * \return true if they are the same
     */
    bool operator==(const AtomPair &other) const;
    
    /*! \brief Return whether one AtomPair is smaller than the other
     * \param[in] other The other AtomPair
     * \return true if this pair is smaller
     */
    bool operator<(const AtomPair &other) const;
};

/*! \brief
 * Chemical bond in a molecule with associated bond order.
 *
 * The chemical bonds in a molecule are stored here along with the bond order
 * which together can be used for determining the atom types in a force field
 * calculation. In the present implementation this data is generated by OpenBabel
 * and stored in molprop files. 
 *
 * \inpublicapi
 * \ingroup module_alexandria
 */
class Bond : public TopologyEntry
{
 public:
    //! Default constructor
    Bond() {}
    
    //! Constructor setting the ids of the atoms and the bondorder
    Bond(int ai, int aj, double bondorder) { Set(ai, aj, bondorder); }

     //! Sets the ids of the atoms and the bondorder
    void Set(int ai, int aj, double bondorder)
    {
        addAtom(ai);
        addAtom(aj);
        addBondOrder(bondorder);
    }
    
    //! Returns the ids of the atoms and the bondorder
    void get(int *ai, int *aj, double *bondorder) const;
    
    //! Returns the first atom id
    int aI() const
    {
        return atomIndex(0);
    }
      
     //! Returns the second atom id
    int aJ() const
    {
        return atomIndex(1);
    }
    
    //! Return a Bond with the order of atoms swapped
    Bond swap() const;
 
    /*! \brief Return whether two Bonds are the same
     * \param[in] other The other bond
     * \return true if they are the same
     */
    bool operator==(const Bond &other) const;
};
//! Iterator over Bond items
using BondIterator      = typename std::vector<Bond>::iterator;
//! Const iterator over Bond items
using BondConstIterator = typename std::vector<Bond>::const_iterator;

/*! \brief
 * Angle between two bonds i-j and j-k in a molecule.
 * Atom j of the first bond should be identical to atom i
 * of the second bond.
 *
 * \inpublicapi
 * \ingroup module_alexandria
 */
class Angle : public TopologyEntry
{
 private:
    //! The bonds
    Bond             b_[2];
    //! Whether this is a linear angle
    bool             isLinear_ = false;
 public:
    //! Default constructor
    Angle() {}
    
    /*! Angle constructor with initial data
     * \param[in] bij  First bond
     * \param[in] bjk  Second bond
     */
    Angle(const Bond bij, const Bond bjk);

    /*! \brief Set the two bonds
     * \param[in] bij First bond
     * \param[in] bjk Second bond
     */
    void setBonds(const Bond &bij, const Bond &bjk);

    //! Return the first bond
    const Bond &bij() const { return b_[0]; }
    
    //! Return the second bond
    const Bond &bjk() const { return b_[1]; }
    
    //! Return whether this is a linear angle
    bool isLinear() const { return isLinear_; }

    //! Set whether this is a linear angle
    void setLinear(bool linear) { isLinear_ = linear; }

    //! @copydoc Bond::renumberAtoms
    void renumberAtoms(const std::vector<int> &renumber);
};

/*! \brief
 * Improper dihedral defined by one central atom i and three bonds
 * bij, bik, bil in a molecule in a planar arrangement.
 * Atom i should be the first atom in all bonds.
 *
 * \inpublicapi
 * \ingroup module_alexandria
 */
class Improper : public TopologyEntry
{
 private:
    //! The bonds
    Bond             b_[3];
 public:
    //! Default constructor
    Improper() {}
    
    /*! Improper dihedral constructor with initial data
     * \param[in] bij  First bond
     * \param[in] bik  Second bond
     * \param[in] bil  Third bond
     */
    Improper(Bond bij, Bond bik, Bond bil);

    //! Return the first bond
    const Bond &bij() const { return b_[0]; }
    
    //! Return the second bond
    const Bond &bik() const { return b_[1]; }
    
    //! Return the third bond
    const Bond &bil() const { return b_[2]; }
    
    //! @copydoc Bond::renumberAtoms
    void renumberAtoms(const std::vector<int> &renumber);
};

/*! \brief
 * Proper dihedral defined by three bonds in a row 
 * bij, bjk, bkl in a molecule.
 *
 * \inpublicapi
 * \ingroup module_alexandria
 */
class Proper : public TopologyEntry
{
 private:
    //! The bonds
    Bond             b_[3];
 public:
    //! Default constructor
    Proper() {}
    
    /*! Improper dihedral constructor with initial data
     * \param[in] bij  First bond
     * \param[in] bjk  Second bond
     * \param[in] bkl  Third bond
     */
    Proper(Bond bij, Bond bjk, Bond bkl);

    //! Return the first bond
    const Bond &bij() const { return b_[0]; }
    
    //! Return the second bond
    const Bond &bjk() const { return b_[1]; }
    
    //! Return the third bond
    const Bond &bkl() const { return b_[2]; }
    
    //! @copydoc Bond::renumberAtoms
    void renumberAtoms(const std::vector<int> &renumber);
};

/*! \brief
 * Virtual site on a linear bond.
 *
 * \inpublicapi
 * \ingroup module_alexandria
 */
class Vsite2 : public TopologyEntry
{
 public:
    //! Default constructor
    Vsite2() {}
    
    //! Constructor setting the ids of the atoms and the bondorder
    Vsite2(int ai, int aj, int vs)
    {
        addAtom(ai);
        addAtom(aj);
        addAtom(vs);
    }
    
    //! Returns the ids of the atoms and the bondorder
    void get(int *ai, int *aj, int *vs) const;
    
    //! Returns the first atom id
    int aI() const
    {
        return atomIndex(0);
    }
      
     //! Returns the second atom id
    int aJ() const
    {
        return atomIndex(1);
    }
    
    //! \return virtual site ID
    int vs() const
    {
        return atomIndex(2);
    }
    //! Return a Vsite2 with the order of atoms swapped
    Vsite2 swap() const;
 
    /*! \brief Return whether two Vsite2 instances are the same
     * \param[in] other The other Vsite2
     * \return true if they are the same
     */
    bool operator==(const Vsite2 &other) const;
};

/*! \brief
 * Virtual site based on 3 particles
 *
 * \inpublicapi
 * \ingroup module_alexandria
 */
class Vsite3 : public TopologyEntry
{
 public:
    //! Default constructor
    Vsite3() {}
    
    //! Constructor setting the ids of the atoms and the bondorder
    Vsite3(int ai, int aj,  int ak, int vs)
    {
        addAtom(ai);
        addAtom(aj);
        addAtom(ak);
        addAtom(vs);
    }
    
    //! Returns the ids of the atoms and the bondorder
    void get(int *ai, int *aj, int *ak, int *vs) const;
    
    //! Returns the first atom id
    int aI() const
    {
        return atomIndex(0);
    }
      
     //! Returns the second atom id
    int aJ() const
    {
        return atomIndex(1);
    }
    
     //! Returns the third atom id
    int aK() const
    {
        return atomIndex(2);
    }
    
    //! \return virtual site ID
    int vs() const
    {
        return atomIndex(3);
    }
    //! Return a Vsite3 with the order of atoms swapped
    Vsite3 swap() const;
 
    /*! \brief Return whether two Vsite3 instances are the same
     * \param[in] other The other Vsite3
     * \return true if they are the same
     */
    bool operator==(const Vsite3 &other) const;
};

} // namespace alexandria

#endif

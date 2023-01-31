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

#ifndef ACT_TOPOLOGY_H
#define ACT_TOPOLOGY_H

#include <cstdio>
#include <map>
#include <vector>

#include "act/alexandria/actmol_low.h"
#include "act/basics/identifier.h"
#include "act/forcefield/particletype.h"
#include "act/forcefield/forcefield.h"
#include "act/forcefield/forcefield_parameterlist.h"
#include "act/utility/communicationrecord.h"
#include "gromacs/gmxpreprocess/grompp-impl.h"
#include "gromacs/gpu_utils/hostallocator.h"

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
            GMX_THROW(gmx::InternalError("Atom index out of range"));
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
    Angle(Bond bij, Bond bjk);

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

class ActAtom
{
private:
    //! My identifier
    Identifier       id_;
    //! The atom name
    std::string      name_;
    //! The chemical element
    std::string      elem_;
    //! The atom type in the force field
    std::string      ffType_;
    //! The particle type
    int              pType_;
    //! My shell particles, if any
    std::vector<int> shells_;
    //! My atom particle, if any,or -1
    int              core_ = -1;
    //! The atomic number
    int              atomicNumber_;
    //! The mass
    double           mass_;
    //! The charge
    double           charge_;
    //! Residue number
    int              residueNumber_ = 1;
public:
    ActAtom(const std::string &name,
            const std::string &elem,
            const std::string &ffType,
            int                pType,
            int                atomicNumber,
            double             newmass,
            double             newcharge) :
        id_({ name }), name_(name), elem_(elem), ffType_(ffType), pType_(pType), atomicNumber_(atomicNumber), mass_(newmass), charge_(newcharge)
    {}
    
    ActAtom(const ParticleType &pt) :
        id_({ pt.id() }), name_(pt.id().id() ), elem_(pt.element()), ffType_( pt.id().id() ),
        pType_( pt.gmxParticleType()), atomicNumber_(pt.atomnumber()), mass_(pt.mass()), charge_(pt.charge())
    {}
    
    //! \return Identifier
    const Identifier &id() const { return id_; }

    //! \return the name
    const std::string &name() const { return name_; }
    
    //! \return the element
    const std::string &element() const { return elem_; }
    
    //! \return the ffType
    const std::string &ffType() const { return ffType_; }
    
    //! \return the particle type
    int pType() const { return pType_; }

    //! \return the list of shells for this particle
    const std::vector<int> &shells() const { return shells_; }

    /*! \brief Add a shell particle index
     * \param[in] index The index
     */
    void addShell(int index) { shells_.push_back(index); }
    
    /*! \brief Set the core particle index
     * \param[in] index The index
     */
    void setCore(int index) { core_ = index; }

    /*! \brief Return the core particle connected to this shell
     * \return the core particle or -1 if there is none
     */
    int core() const { return core_; }

    //! \return the atomic number
    int atomicNumber() const { return atomicNumber_; }

    //! \return the mass
    double mass() const { return mass_; }
    
    //! \return the charge
    double charge() const { return charge_; }
    
    /*! set the charge
     * \param[in] charge The new value
     */
    void setCharge(double newcharge) { charge_ = newcharge; }

    //! return the residue number
    int residueNumber() const { return residueNumber_; }
};

class Topology
{
private:
    //! Map from InteractionType to topology element pointers
    std::map<InteractionType, std::vector<TopologyEntry *> >  entries_;
    //! Non bonded exclusions, array is length of number of atoms
    std::vector<std::vector<int> >                            exclusions_;
    //! List of atoms
    std::vector<ActAtom>                                      atoms_;
    //! The (residue) name
    std::string                                               name_;
 public:
    Topology()
    {
        name_.assign("MOL");
    }

    /*! Constructor
     * This code copies relevant structures from the outside world
     * \param[in] bonds The bonds connecting this molecule.
     */
    Topology(const std::vector<Bond> &bonds);

    //! Return the name
    const std::string &name() const { return name_; }

    /*! Find a bond between two atoms
     * \param[in] ai The first atom
     * \param[in] aj The second atom
     * \return pointer to the Bond
     * \throws if not found
     */
    const Bond *findBond(int ai, int aj) const;

    /*! \brief Add one bond to the list
     * \param[in] bond The bond to add
     */
    void addBond(const Bond &bond);

    /*! \brief Initiate or update the atoms data.
     * Must be called every time the data changes (e.g. charges).
     * \param[in] atoms Gromacs atoms structure
     */
    void setAtoms(const t_atoms *atoms);

    /*! \brief Copy shell info to atoms
     * Must be called after setAtoms
     */
    void shellsToAtoms();

    //! Add an atom to the topology
    void addAtom(const ActAtom &atom) { atoms_.push_back(atom); }
    
    /*! \brief Build the topology internals
     * Calls the functions to generate angles, impropers and dihedrals (if flag set).
     * Will also generate non-bonded atom pairs and exclusions.
     * \param[in] pd             The force field structure
     * \param[in] x              The atomic coordinates
     * \param[in] LinearAngleMin Minimum angle to be considered linear (degrees)
     * \param[in] PlanarAngleMax Maximum angle to be considered planar (degrees)
     * \param[in] missing        How to treat missing parameters
     */
    void build(const ForceField                *pd,
               const std::vector<gmx::RVec> &x,
               double                        LinearAngleMin,
               double                        PlanarAngleMax,
               missingParameters             missing);

    //! \return the vector of atoms
    const std::vector<ActAtom> &atoms() const { return atoms_; }
    
    //! \return the vector of atoms for editing
    std::vector<ActAtom> *atomsPtr() { return &atoms_; }
    
    //! \return the mass of the compound
    double mass() const;
    
    /*! \brief Find a topology entry matching the inputs/
     * \param[in] itype     The InteractionType
     * \param[in] aindex    The atom indices
     * \param[in] bondOrder The array of bond orders
     * \param[in] cs        Whether or not the order of the atoms can be swapped
     * \return the entry
     */
    const TopologyEntry *findTopologyEntry(InteractionType            itype,
                                           const std::vector<int>    &aindex,
                                           const std::vector<double> &bondOrder,
                                           CanSwap                    cs) const;

    /*! Generate the angles
     * To generate angles we need the coordinates to check whether
     * there is a linear geometry.
     * \param[in] x              The atomic coordinates
     * \param[in] LinearAngleMin Minimum angle to be considered linear (degrees)
     */
    void makeAngles(const std::vector<gmx::RVec> &x,
                    double                        LinearAngleMin);

    /*! Generate the impropers
     * To generate impropers we need the coordinates to check whether
     * there is a planar geometry.
     * \param[in] x              The atomic coordinates
     * \param[in] PlanarAngleMax Maximum angle to be considered planar (degrees)
      */
    void makeImpropers(const std::vector<gmx::RVec> &x,
                       double                        PlanarAngleMax);

    /*! \brief Add a custom list of interactions
     * \param[in] itype The interaction type (should not yet exist)
     * \param[in] vec   The new interactions
     */
    void addEntry(InteractionType                     itype,
                  const std::vector<TopologyEntry *> &entry);

    /*! \brief Generate exclusiones
     * \param[in]  nrexcl    The number of exclusions to generate (max 2)
     * \param[in]  nratom    The number of atoms in the system
     */
    void generateExclusions(int nrexcl,
                            int nratoms);
    
    //! \return the number of atoms                    
    size_t nAtoms() const { return atoms_.size(); }
        
    /*! \brief Convert to gromacs exclusiones
     */
    t_excls *gromacsExclusions();
    /*! Generate the non-bonded pair list based on just the atoms
     * \param[in] natoms The number of atoms
     */
    void makePairs(int natoms);

    /*! Add shell pairs
     * Based on the atom pair list, add interaction between atoms and shell
     * and between pairs of shells. This should be called after makePairs and
     * generateExclusions have been called (and in that order) and after shells
     * have been added.
     */
    void addShellPairs();
    
    /*! Generate the proper dihedrals
     */
    void makePropers();

    /*! Generate virtual sites for bonds
     */
    void makeVsite2s(const ForceFieldParameterList &vsite2);

    //! @copydoc Bond::renumberAtoms
    void renumberAtoms(const std::vector<int> &renumber);

    /*! \brief Check whether an interaction type is present
     * \param[in] itype the looked after type
     * \return true if found
     */
    bool hasEntry(InteractionType itype) const { return entries_.find(itype) != entries_.end(); }

    /*! \brief Return the desired entry or throw
     * \param[in] itype The desired interaction type
     * \return a vector of pointers
     * \throws if not found. If you do not want that to happen, check first!
     */
    const std::vector<TopologyEntry *> &entry(InteractionType itype) const;

    //! \return the whole map of parameters
    std::map<InteractionType, std::vector<TopologyEntry *> > *entries() { return &entries_; }

    //! \return the whole map of parameters, const style
    const std::map<InteractionType, std::vector<TopologyEntry *> > &entries() const { return entries_; }
    
    /*! \brief Print structure to a file
     * \param[in] fp The file pointer
     */
    void dump(FILE *fp) const;

    /*! \brief Fill in the parameters in the topology entries.
     * Must be called repeatedly during optimizations of energy.
     * \param[in] pd The force field structure
     */
    void fillParameters(const ForceField *pd);
         
    /*! \brief Add identifiers to interactions
     * \param[in] pd The force field structure
     */                              
    void setIdentifiers(const ForceField *pd);
};

} // namespace alexandria

#endif

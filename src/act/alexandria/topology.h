/*
 * This source file is part of the Alexandria Chemistry Toolkit.
 *
 * Copyright (C) 2021-2024
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
#include <list>
#include <map>
#include <vector>

#include "act/alexandria/actmol_low.h"
#include "act/basics/identifier.h"
#include "act/forcefield/particletype.h"
#include "act/forcefield/forcefield.h"
#include "act/forcefield/forcefield_parameterlist.h"
#include "act/molprop/fragment.h"
#include "act/molprop/molprop.h"
#include "act/molprop/topologyentry.h"
#include "act/utility/communicationrecord.h"

namespace alexandria
{

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
    std::vector<int> cores_;
    //! The atomic number
    int              atomicNumber_;
    //! The mass
    double           mass_;
    //! The charge
    double           charge_;
    //! Residue number. If at -1 it means it has not been set.
    int              residueNumber_ = -1;
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

    //! \return the list of shells for this particle for editing
    std::vector<int> *shellsPtr() { return &shells_; }

    /*! \brief Add a shell particle index
     * \param[in] index The index
     */
    void addShell(int index) { shells_.push_back(index); }

    /*! \brief Set the core particle index
     * \param[in] index The index
     */
    void addCore(int index) { cores_.push_back(index); }

    /*! \brief Return the core particle connected to this shell
     * \return the core particle or -1 if there is none
     */
    const std::vector<int> &cores() const { return cores_; }

    /*! \brief Return the core particle connected to this shell
     * \return the core particle or -1 if there is none
     */
    std::vector<int> *coresPtr() { return &cores_; }

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

    /*! set the residue number
     * \param[in] resnr the new number
     */
    void setResidueNumber(int resnr) { residueNumber_ = resnr; }
};

class ActAtomListItem;

//! Shortcut for passing this around.
typedef std::list<ActAtomListItem> AtomList;

class Topology
{
private:
    //! Map from InteractionType to topology element pointers
    std::map<InteractionType, TopologyEntryVector>            entries_;
    //! Non bonded exclusions, array is length of number of atoms
    std::map<InteractionType, std::vector<std::vector<int>>>  exclusions_;
    //! List of atoms
    std::vector<ActAtom>                                      atoms_;
    //! The molecule (compound) name
    std::string                                               moleculeName_;
    //! The residue names. If just one, this will likely be the same as the moleculeName_
    std::vector<std::string>                                  residueNames_;
    //! The real atoms, that is, not the vsites or shells
    std::vector<int>                                          realAtoms_;

    /*! Generate virtual sites for bonds.
     * \param[in]    pd       The force field
     * \param[inout] atomList The atom and coordinate list that may be extended
     * \return map containing number of entries added for each interaction type
     */
    std::map<InteractionType, size_t> makeVsite2s(const ForceField *pd,
                                                  AtomList         *atomList);

    /*! \brief Generate three body vsites
     * \param[in]    pd       The force field
     * \param[inout] atomList The linked list of particles
     * \return map containing number of entries added for each interaction type
     */
    std::map<InteractionType, size_t> makeVsite3s(const ForceField *pd,
                                                  AtomList         *atomList);

    /*! \brief Add identifiers to interactions
     * \param[in] pd The force field structure
     * \param[in] itype The interaction type for which to do this
     */
    void setEntryIdentifiers(const ForceField *pd,
                             InteractionType   itype);

    /*! Add polarizabilities to the topology if needed
     * \param[in]    pd       Force field
     * \param[inout] atomList The atoms and coordinates will be updated as well
     */
    void addShells(const ForceField *pd,
                   AtomList         *atomList);

    /*! \brief Find a topology entry matching the inputs if it exists
     * \param[in] itype     The InteractionType
     * \param[in] aindex    The atom indices
     * \param[in] bondOrder The array of bond orders
     * \param[in] cs        Whether or not the order of the atoms can be swapped
     * \return The entry you were looking for or nullptr
     */
    const TopologyEntry *findTopologyEntry(const TopologyEntryVector &entries,
                                           const std::vector<int>    &aindex,
                                           const std::vector<double> &bondOrder,
                                           CanSwap                    cs) const;

    /*! \brief Remove interactions that should not be there in the first place.
     * If two atoms are in the exclusion list, their shells should be as well.
     * This routine checks for such interactions (Coulomb or VDW) and removes them.
     * For vsites, they inherit the exclusions from their constructing atoms.
     * \param[inout] pairs      The pair list
     * \param[in]    exclusions The exclusions
     */
    void fixExclusions(TopologyEntryVector                 *pairs,
                       const std::vector<std::vector<int>> &exclusions);

    /*! \brief Generate exclusions and remove pairs
     * \param[inout] pairs The list of all atom pairs
     * \param[in] nrexcl   The number of exclusions to generate for Coulomb interactions (max 2)
     * \return list of exclusions
     */
    std::vector<std::vector<int>> generateExclusions(TopologyEntryVector *pairs,
                                                     int                  nrexcl);

 public:
    Topology()
    {
    }

    /*! Constructor
     * This code copies relevant structures from the outside world
     * \param[in] bonds The bonds connecting this molecule.
     */
    Topology(const std::vector<Bond> &bonds);

    //! Return the name
    const std::string &moleculeName() const { return moleculeName_; }

    /*! Find a bond between two atoms
     * \param[in] ai The first atom
     * \param[in] aj The second atom
     * \return Copy of the Bond
     * \throws if not found
     */
    Bond findBond(int ai, int aj) const;

    /*! \brief Add one bond to the list
     * \param[in] bond The bond to add
     */
    void addBond(const Bond &bond);

    //! Add an atom to the topology
    void addAtom(const ActAtom &atom) { atoms_.push_back(atom); }

    //! Debugging stuff
    void dumpPairlist(FILE *fp, InteractionType itype) const;

    //! \return the array of real atoms
    const std::vector<int> &realAtoms() const { return realAtoms_; }

    /*! Generate atoms in the topology from an experiment (QM calc)
     * \param[in] pd  The force field
     * \param[in] mol The molecule
     * \param[out] x The coordinate
     * \returns status variable
     */
    immStatus GenerateAtoms(const ForceField       *pd,
                            const MolProp          *mol,
                            std::vector<gmx::RVec> *x);
    /*! \brief Build the topology internals
     * Calls the functions to generate angles, impropers and dihedrals (if flag set).
     * Will also generate non-bonded atom pairs and exclusions. Shells will be generated
     * if they are present in the force field.
     * \param[in] pd             The force field structure
     * \param[inout] x           The atomic coordinates
     * \param[in] LinearAngleMin Minimum angle to be considered linear (degrees)
     * \param[in] PlanarAngleMax Maximum angle to be considered planar (degrees)
     * \param[in] missing        How to treat missing parameters
     */
    void build(const ForceField       *pd,
               std::vector<gmx::RVec> *x,
               double                  LinearAngleMin,
               double                  PlanarAngleMax,
               missingParameters       missing);

    /*! \brief Fill in the parameters in the topology entries.
     * Must be called repeatedly during optimizations of energy.
     * \param[in] pd The force field structure
     */
    void fillParameters(const ForceField *pd);
    //! \return the vector of atoms
    const std::vector<ActAtom> &atoms() const { return atoms_; }

    //! \return the residue names
    const std::vector<std::string> residueNames() const { return residueNames_; }

    /*! \brief Add residue information
     * \param[in] residueNumber The original number
     * \param[in] residueName   The new number
     */
    void addResidue(int                residueNumber,
                    const std::string &residueName);

    //! \return the vector of atoms for editing
    std::vector<ActAtom> *atomsPtr() { return &atoms_; }

    //! \return the mass of the compound
    double mass() const;

    /*! Generate the angles
     * To generate angles we need the coordinates to check whether
     * there is a linear geometry.
     * \param[in] pd             The force field
     * \param[in] x              The atomic coordinates
     * \param[in] LinearAngleMin Minimum angle to be considered linear (degrees)
     */
    void makeAngles(const ForceField             *pd,
                    const std::vector<gmx::RVec> &x,
                    double                        LinearAngleMin);

    /*! Generate the impropers
     * To generate impropers we need the coordinates to check whether
     * there is a planar geometry.
     * \param[in] pd             The force field
     * \param[in] x              The atomic coordinates
     * \param[in] PlanarAngleMax Maximum angle to be considered planar (degrees)
      */
    void makeImpropers(const ForceField             *pd,
                       const std::vector<gmx::RVec> &x,
                       double                        PlanarAngleMax);

    /*! Generate the proper dihedrals
     * \param[in] pd The force field
     */
    void makePropers(const ForceField *pd);

    /*! \brief Add a custom list of interactions
     * \param[in] itype The interaction type (should not yet exist)
     * \param[in] vec   The new interactions
     */
    void addEntry(InteractionType            itype,
                  const TopologyEntryVector &entry);

    //! \return the number of atoms
    size_t nAtoms() const { return atoms_.size(); }

    /*! Generate the non-bonded pair list based on just the atoms
     * \param[in] pd     Force field needed to set identifiers.
     * \param[in] natoms The number of atoms
     */
    void makePairs(const ForceField *pd,
                   InteractionType   itype);

    /*! Add shell pairs
     * Based on the atom pair list, add interaction between atoms and shell
     * and between pairs of shells. This should be called after makePairs and
     * generateExclusions have been called (and in that order) and after shells
     * have been added.
     */
    void addShellPairs();

    //! @copydoc Bond::renumberAtoms
    void renumberAtoms(const std::vector<int> &renumber);

    /*! \brief Check whether an interaction type is present
     * \param[in] itype the looked after type
     * \return true if found
     */
    bool hasEntry(InteractionType itype) const { return entries_.find(itype) != entries_.end(); }

    //! \return whether there are any virtual sites
    bool hasVsites() const;

    /*! \brief Return the desired entry or throw
     * \param[in] itype The desired interaction type
     * \return a vector of pointers
     * \throws if not found. If you do not want that to happen, check first!
     */
    const TopologyEntryVector &entry(InteractionType itype) const;

    //! \return the whole map of parameters
    std::map<InteractionType, TopologyEntryVector> *entries() { return &entries_; }

    //! \return the whole map of parameters, const style
    const std::map<InteractionType, TopologyEntryVector> &entries() const { return entries_; }

    /*! \brief Determine whether an exclusion type is present
     * \param[in] itype the InteractionType (COULOMB or VANDERWAALS)
     * \return the result
     */
    bool hasExclusions(InteractionType itype) const
    {
        return exclusions_.find(itype) != exclusions_.end();
    }

    /*! \brief Return exclusions for interaction type
     * \param[in] itype the InteractionType (COULOMB or VANDERWAALS)
     * \return the exclusions.
     * \throws if itype is not found 
     */
    const std::vector<std::vector<int>> &exclusions(InteractionType itype) const;

    /*! \brief Print structure to a file
     * \param[in] fp The file pointer
     */
    void dump(FILE *fp) const;

};

} // namespace alexandria

#endif

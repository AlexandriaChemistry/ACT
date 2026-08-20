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

#ifndef ACT_ACTATOM_H
#define ACT_ACTATOM_H

#include <list>
#include <vector>

#include "act/basics/act_particle.h"
#include "act/basics/identifier.h"
#include "act/forcefield/particletype.h"

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
    //! Row in the periodic table (for Slater electrostatics)
    int              row_ = 0;
    //! The atom type in the force field
    std::string      ffType_;
    //! The particle type
    ActParticle      apType_;
    //! My shell particles, if any
    std::vector<int> shells_;
    //! My vsite particles, if any
    std::vector<int> vsites_;
    /*! My "parent" particles, if any. 
     * Both VSite and Shell particles have "parents"
     * If shells are connected to vsites, the parent maybe a vsite,
     * but typically these will be atoms.
     */
    std::vector<int> parents_;
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
            ActParticle        apType,
            int                atomicNumber,
            double             newmass,
            double             newcharge,
            int                row) :
        id_({ name }), name_(name), elem_(elem), row_(row), ffType_(ffType), apType_(apType), atomicNumber_(atomicNumber), mass_(newmass), charge_(newcharge)
    {}

    ActAtom(const ParticleType &pt) :
        id_({ pt.id() }), name_(pt.id().id() ), elem_(pt.element()), ffType_( pt.id().id() ),
        apType_( pt.apType()), atomicNumber_(pt.atomnumber()), mass_(pt.mass()), charge_(pt.charge())
    {}

    //! \return Identifier
    const Identifier &id() const { return id_; }

    //! \return the name
    const std::string &name() const { return name_; }

    //! \return the element
    const std::string &element() const { return elem_; }

    //! \return the row
    int row() const { return row_; }

    //! \return the ffType
    const std::string &ffType() const { return ffType_; }

    //! \return the particle type
    ActParticle pType() const { return apType_; }

    //! \return the list of shells for this particle
    const std::vector<int> &shells() const { return shells_; }

    //! \return the list of vsites for this particle
    const std::vector<int> &vsites() const { return vsites_; }

    //! \return the list of shells for this particle for editing
    std::vector<int> *shellsPtr() { return &shells_; }

    /*! \brief Add a shell particle index
     * \param[in] index The index
     */
    void addShell(int index) { shells_.push_back(index); }

    /*! \brief Add a vsite particle index
     * \param[in] index The index
     */
    void addVsite(int index) { vsites_.push_back(index); }

    /*! \brief Set the core particle index
     * \param[in] index The index
     */
    void addCore(int index) { parents_.push_back(index); }

    /*! \brief Return the core particle connected to this shell
     * \return the core particle or -1 if there is none
     */
    const std::vector<int> &parents() const { return parents_; }

    /*! \brief Return the core particle connected to this shell
     * \return the core particle or -1 if there is none
     */
    std::vector<int> *parentsPtr() { return &parents_; }

    //! \return the atomic number
    int atomicNumber() const { return atomicNumber_; }

    //! \return the mass
    double mass() const { return mass_; }

    //! \return the charge
    double charge() const { return charge_; }

    /*! set the charge
     * \param[in] newcharge The new value
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


} // namespace alexandria

#endif

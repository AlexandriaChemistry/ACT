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

#include <vector>

#include "act/basics/act_particle.h"
#include "act/basics/identifier.h"
#include "act/forcefield/particletype.h"

namespace alexandria
{

/*! \brief A particle in a molecule: an atom, a shell or a virtual site,
 * with its properties and its relations to other particles.
 */
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
    //! The inverse mass
    double           invmass_ = 0;
    //! The charge
    double           charge_;
    //! Second charge at the same position, see ParticleType::chargeVS
    double           chargeVS_ = 0;
    //! Width of the distribution of the second charge
    double           zetaVS_   = 0;
    //! Row of the second charge, for Slater distributions
    int              rowVS_    = 0;
    //! Residue number. If at -1 it means it has not been set.
    int              residueNumber_ = -1;
public:
    /*! \brief Constructor
     * \param[in] name         The atom name, also used as identifier
     * \param[in] elem         The chemical element
     * \param[in] ffType       The atom type in the force field
     * \param[in] apType       The particle type
     * \param[in] atomicNumber The atomic number
     * \param[in] newmass      The mass
     * \param[in] newcharge    The charge
     * \param[in] row          Row in the periodic table
     */
    ActAtom(const std::string &name,
            const std::string &elem,
            const std::string &ffType,
            ActParticle        apType,
            int                atomicNumber,
            double             newmass,
            double             newcharge,
            int                row) :
        id_({ name }), name_(name), elem_(elem), row_(row), ffType_(ffType), apType_(apType), atomicNumber_(atomicNumber), mass_(newmass), charge_(newcharge)
    {
        if (mass_ > 0)
        {
            invmass_ = 1.0/mass_;
        }
    }

    /*! \brief Constructor copying the properties of a force field particle type
     * \param[in] pt The particle type
     */
    ActAtom(const ParticleType &pt) :
        id_({ pt.id() }), name_(pt.id().id() ), elem_(pt.element()), ffType_( pt.id().id() ),
        apType_( pt.apType()), atomicNumber_(pt.atomnumber()), mass_(pt.mass()), charge_(pt.charge()),
        chargeVS_(pt.chargeVS()), zetaVS_(pt.zetaVS()), rowVS_(pt.rowVS())
    {
        if (pt.mass() > 0)
        {
            invmass_ = 1.0/pt.mass();
        }
    }

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

    //! \return the inverse mass, 0 if the mass is not positive
    double invmass() const { return invmass_; }

    //! \return the charge
    double charge() const { return charge_; }

    /*! set the charge
     * \param[in] newcharge The new value
     */
    void setCharge(double newcharge) { charge_ = newcharge; }

    /*! \return the second charge, at the same position but with its own
     * distribution, see ParticleType::chargeVS
     */
    double chargeVS() const { return chargeVS_; }

    //! \return the width of the distribution of the second charge
    double zetaVS() const { return zetaVS_; }

    //! \return the row of the second charge
    int rowVS() const { return rowVS_; }

    //! \return the sum of the charge and the second charge
    double totalCharge() const { return charge_ + chargeVS_; }

    /*! \brief Set the second charge
     * \param[in] chargeVS The charge
     * \param[in] zetaVS   The width of its distribution
     * \param[in] rowVS    Its row
     */
    void setChargeVS(double chargeVS, double zetaVS, int rowVS)
    {
        chargeVS_ = chargeVS;
        zetaVS_   = zetaVS;
        rowVS_    = rowVS;
    }

    //! \return the residue number, -1 if not set
    int residueNumber() const { return residueNumber_; }

    /*! set the residue number
     * \param[in] resnr the new number
     */
    void setResidueNumber(int resnr) { residueNumber_ = resnr; }
};

} // namespace alexandria

#endif

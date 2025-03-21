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
#ifndef QTYPE_H
#define QTYPE_H

#include <map>
#include <string>

#include "act/molprop/molpropobservable.h"
#include "gromacs/math/vec.h"
#include "gromacs/math/vectypes.h"

namespace alexandria
{

class ActAtom;
class ForceField;
class ForceComputer;
class QgenResp;
class Topology;
 
/*! \brief Enumerated type to differentiate the charge types 
 * and properties derived from the charges.
 */
enum class qType { 
    //! Alexandria Charge Model derived property
    Calc,
    //! Electrostatic potential derived point charges
    ESP,
    //! Restrained electrostatic potential method
    RESP,
    //! AM1 Bond-charge correction method
    BCC,
    //! Mulliken charges
    Mulliken,
    //! Hirshfeld charges
    Hirshfeld,
    //! CM5 charges
    CM5,
    //! Gasteiger charges
    Gasteiger,
    //! Electronic properties straight from DFT or QC calcs
    Elec,
    //! Alexandria charge model charges
    ACM
};

/*! \brief return string corresponding to charge type
 */
const std::string &qTypeName(qType qt);

/*! \brief convert string to qtype
 * \param[in] type The string
 * \return a qType.
 * \throws if not found
 */
qType stringToQtype(const std::string &type);

/*! \brief Return a complete map of qTypes and their names
 */
const std::map<qType, std::string> &qTypes();

/*! Class to hold electrostatic properties.
 * To compare the properties of different models we have this class
 * to hold electrostatic moments and charges, and indeed a grid 
 * structure to hold the electrostatic potential.
 */
class QtypeProps
{
 private:
    //! Identity
    qType                  qtype_      = qType::Calc;
    //! Electrostatic moments
    std::map<MolPropObservable, std::vector<double> > multipoles_;
    //! Polarizability tensor
    bool                   hasAlpha_   = false;
    tensor                 alpha_      = { { 0 } };
    //! Polarizability anisotropy
    double                 anisotropy_ = 0;
    //! Polarizability isotropic value
    double                 isotropy_   = 0;
    //! The coordinates
    std::vector<gmx::RVec> x_;
    //! Center of charge
    rvec                   coc_        = { 0 };
    //! Atomic charges
    std::vector<double>    q_;
    //! Atomic numbers
    std::vector<int>       atomNumber_;
    //! Reset all the calc moments to zero
    void resetMoments();
    //! Recompute the center of charge
    void computeCoC();
    
    
 public:
    //! Empty constructor
    QtypeProps() {}
    /*! \brief Constructor
     * \param[in] qtype  My own identity
     * \param[in] atoms  The atoms in this compound
     * \param[in] coords and their coordinates
     */
    QtypeProps(qType                         qtype,
               const std::vector<ActAtom>   &atoms,
               const std::vector<gmx::RVec> &coords);

    //! Return my qType
    qType qtype() const { return qtype_; }

    /*! Set my qType
     * \param[in] qtype The type
     */
    void setQtype(qType qtype) { qtype_ = qtype; }

    //! Initialize variables for multipoles
    void initializeMoments();

    /*! \brief Set charges and coordinates.
     *
     * \param[in] q The charges
     * \param[in] x The coordinates
     */
    void setQandX(const std::vector<double>    &q,
                  const std::vector<gmx::RVec> &x);
    
    /*! \brief Set charges and coordinates.
     *
     * \param[in] atoms The atoms
     * \param[in] x     The coordinates
     */
    void setQandX(const std::vector<ActAtom>   &atoms,
                  const std::vector<gmx::RVec> &x);
    
    /*! \brief Set charges.
     *
     * \param[in] q The charges
     */
    void setQ(const std::vector<double> &q);
    /*! \brief Set charges.
     *
     * \param[in] atoms An ActAtoms vector
     */
    void setQ(const std::vector<ActAtom> &atoms);
    /*! \brief Set coordinates.
     *
     * \param[in] x The coordinates
     */
    void setX(const std::vector<gmx::RVec> &x);
    
    //! Return the coordinates
    const std::vector<gmx::RVec> &x() const { return x_; }
    
    /*! \brief Store center of charge
     * \param[in] coc The center of charge
     */
    void setCenterOfCharge(const rvec &coc) { copy_rvec(coc, coc_); }
    
    /*! \brief Compute electric moments and store them.
     */
    void calcMoments();
    
    /*! \brief Return dipole for charge type qt.
     */
    double dipole() const;
    
    //! Check whether a polarizability is present
    bool hasPolarizability() const { return hasAlpha_; }
    
    /*! \brief Return polarizability tensor
     */
    const tensor &polarizabilityTensor() const { return alpha_; }

    /*! \brief Set the polarizability tensor
     * \param[in] alpha The tensor
     */
    void setPolarizabilityTensor(const tensor &alpha);

    /*! \brief Compute the polarizability tensor for the system.
     * The result is stored internally
     * \param[in] pd        The force field
     * \param[in] top       The molecular topology
     * \param[in] forceComp The force computer    
     */
    void calcPolarizability(const ForceField    *pd,
                            const Topology      *top,
                            const ForceComputer *forceComp);
    
    /*! \brief Return isotropic polarizability
     */
    double isotropicPolarizability() const { return isotropy_; }
    
    /*! \brief Return anisotropic polarizability
     */
    double anisotropicPolarizability() const { return anisotropy_; }
   
    /*! \brief Check whether a multipole observable is present
     * \param[in] mpo The type of multipole
     * \return true if found
     */
    bool hasMultipole(MolPropObservable mpo) const;
    
    /*! \brief Return multipole
     * \param[in] mpo The type of multipole
     * \return The multipole values
     * \throws with invalid input
     */
    const std::vector<double> getMultipole(MolPropObservable mpo) const;

    /*! \brief Set quadrupole tensor
     * \param[in] quad The quadrupole tensor
     */
    void setMultipole(MolPropObservable mpo, const std::vector<double> &mult);

    /*! \brief Return charges
     * \return The atomic charges
     */
    const std::vector<double> &charge() const { return q_; };

};
 
} // namespace alexandria

#endif

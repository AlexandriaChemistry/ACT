/*
 * This source file is part of the Alexandria Chemistry Toolkit.
 *
 * Copyright (C) 2021
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
 * \author Mohammad Mehdi Ghahremanpour <mohammad.ghahremanpour@icm.uu.se>
 * \author David van der Spoel <david.vanderspoel@icm.uu.se>
 */
#ifndef QTYPE_H
#define QTYPE_H

#include <map>
#include <string>

#include "gromacs/gpu_utils/hostallocator.h"
#include "gromacs/math/vec.h"
#include "gromacs/math/vectypes.h"

struct t_atoms;

namespace alexandria
{

class QgenResp;

/*! \brief Enumerated type to differentiate the charge types 
 * and properties derived from the charges.
 */
enum class qType { 
    //! Alexandria Charge Model derived property
    Calc,
    //! Electrostatic potential derived point charges
    ESP,
    //! Mulliken charges
    Mulliken,
    //! Hirshfeld charges
    Hirshfeld,
    //! CM5 charges
    CM5,
    //! Gasteiger charges
    Gasteiger,
    //! Electronic properties straight from DFT or QC calcs
    Elec
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

/*! Class to hold electrostatic properties
 * To compare the properties of different models we have this class
 * to hold electrostatic moments and charges, and indeed a grid 
 * structure to hold the electrostatic potential.
 */
class QtypeProps
{
 private:
    //! Identity
    qType                  qtype_;
    //! Electrostatic moments
    rvec                   mu_         = { 0 };
    tensor                 quadrupole_ = {{ 0 }};
    //! Norm of the dipole
    double                 dipole_     = 0;
    //! The coordinates
    gmx::HostVector<gmx::RVec> x_;
    //! Center of charge
    rvec                   coc_        = {0, 0, 0};
    //! Atomic charges
    std::vector<double>    q_;
    //! Resp calculation structure
    QgenResp               *QgenResp_   = nullptr;
 public:
    /*! \brief Constructor
     * \param[in] qtype  My own identity
     */
    QtypeProps(qType qtype) : qtype_(qtype)  {}
    
    /*! \brief Set charges and coordinates.
     *
     * \param[in] q The charges
     * \param[in] x The coordinates
     */
    void setQandX(const std::vector<double>        &q,
                  const gmx::HostVector<gmx::RVec> &x);
    
    /*! \brief Set charges.
     *
     * \param[in] atoms An t_atoms structure
     */
    void setQ(const t_atoms *atoms);
    
    /*! \brief Set charges.
     *
     * \param[in] q The charges
     */
    void setQ(const std::vector<double> &q);
    
    /*! \brief Set coordinates.
     *
     * \param[in] x The coordinates
     */
    void setX(const gmx::HostVector<gmx::RVec> &x);
    
    //! Return the coordinates
    const gmx::HostVector<gmx::RVec> &x() const { return x_; }
    
    /*! \brief Store center of charge
     * \param[in] coc The center of charge
     */
    void setCenterOfCharge(const rvec &coc) { copy_rvec(coc, coc_); }
    
    /*! \brief Compute electric moments and store them.
     */
    void calcMoments();
    
    /*! \brief Return dipole vector
     * \return The dipole vector
     */
    const rvec &mu() const { return mu_; };
    
    /*! \brief Set dipole vector
     * \param[in] mu The dipole vector
     */
    void setMu(const rvec mu);
    
    /*! \brief
     * Return computed dipole for charge type qt.
     */
    double dipole() const { return dipole_; }

    /*! \brief Return quadrupole
     * \return The quadrupole
     */
    const tensor &quad() const { return quadrupole_; };

    /*! \brief Set quadrupole tensor
     * \param[in] quad The quadrupole tensor
     */
    void setQuadrupole(const tensor quad);
    
    /*! \brief Return charges
     * \return The atomic charges
     */
    const std::vector<double> &charge() const { return q_; };

    /*! \brief Return internal structure
     * \return the QgenResp_ data structure
     */
    QgenResp *qgenResp();
    
    /*! \brief Return internal structure
     * \return the QgenResp_ data structure
     */
    const QgenResp *qgenRespConst();

    //! Copy back the charges from Resp
    void copyRespQ();
};
 
} // namespace alexandria

#endif

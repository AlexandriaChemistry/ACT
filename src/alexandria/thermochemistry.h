/*
 * This file is part of the GROMACS molecular simulation package.
 *
 * Copyright (c) 2017,2018,2019,2020, by the GROMACS development team, led by
 * Mark Abraham, David van der Spoel, Berk Hess, and Erik Lindahl,
 * and including many others, as listed in the AUTHORS file in the
 * top-level source directory and at http://www.gromacs.org.
 *
 * GROMACS is free software; you can redistribute it and/or
 * modify it under the terms of the GNU Lesser General Public License
 * as published by the Free Software Foundation; either version 2.1
 * of the License, or (at your option) any later version.
 *
 * GROMACS is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
 * Lesser General Public License for more details.
 *
 * You should have received a copy of the GNU Lesser General Public
 * License along with GROMACS; if not, see
 * http://www.gnu.org/licenses, or write to the Free Software Foundation,
 * Inc., 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301  USA.
 *
 * If you want to redistribute modifications to GROMACS, please
 * consider that scientific software is very special. Version
 * control is crucial - bugs must be traceable. We will be happy to
 * consider code for inclusion in the official distribution, but
 * derived work must not be called official GROMACS. Details are found
 * in the README & COPYING files - if they are missing, get the
 * official version at http://www.gromacs.org.
 *
 * To help us fund GROMACS development, we humbly ask that you cite
 * the research papers on the package. Check out http://www.gromacs.org.
 */
/*! \internal \file
 * \brief
 * Code for computing entropy and heat capacity from eigenvalues
 *
 * \author David van der Spoel <david.vanderspoel@icm.uu.se>
 */
#ifndef ALEXANDRIA_THERMOCHEMISTRY_H
#define ALEXANDRIA_THERMOCHEMISTRY_H

#include <map>
#include <vector>

#include "alexandria/mymol.h"
#include "gromacs/math/units.h"
#include "gromacs/math/vec.h"
#include "gromacs/utility/basedefinitions.h"

namespace alexandria
{

enum class TCComponent {
    Translation, Rotation, Vibration, Total
};

/*! \brief Get name corresponding to TCComponent
 * \param[in] tcc The ThermoChemistry component
 * \return string corresponding to tcc
 */
const std::string &tcComponentName(TCComponent tc);

/*! \return map of TCCs and strings
 */
const std::map<TCComponent, std::string> &tccmap();

class ThermoChemistry
{
private:
    double zpe_     = 0;
    double dhForm_  = 0;
    std::map<TCComponent, double> cv_;
    std::map<TCComponent, double> S0_;
    std::map<TCComponent, double> E_;

    /*! \brief Compute zero point energy from an array of eigenvalues.
     *
     * This routine computes the zero point energy.
     *
     * \param[in] frequencies  The frequencies
     * \param[in] scale_factor Factor to scale frequencies by before computing cv
     * \return The zero point energy (kJ/mol)
     */
    double zeroPointEnergy(const std::vector<double> &frequencies,
                           double                     scale_factor);
    
    
    /*! \brief Compute properties due to vibrational motion
     *
     * Results for cv, S0 and E are stored in internal variables.
     * \param[in] frequencies  The vibrational frequencies
     * \param[in] temperature  Temperature (K)
     * \param[in] scale_factor Factor to scale frequencies by before computing properties
     */
    void calcVibrationalProperties(const std::vector<double> &frequencies,
                                   double                     temperature,
                                   double                     scale_factor);
                                       
    /*! \brief Compute entropy due to translational motion
     *
     * Following the equations in J. W. Ochterski,
     * Thermochemistry in Gaussian, Gaussian, Inc., 2000
     * Pitssburg PA
     *
     * \param[in] mass         Molecular mass (Dalton)
     * \param[in] temperature  Temperature (K)
     * \param[in] pressure     Pressure (bar) at which to compute
     * \returns The translational entropy (J/mol K)
     */
    double translationalEntropy(double mass,
                                double temperature,
                                double pressure);

    /*! \brief Compute entropy due to rotational motion
     *
     * Following the equations in J. W. Ochterski,
     * Thermochemistry in Gaussian, Gaussian, Inc., 2000
     * Pitssburg PA
     *
     * \param[in] temperature  Temperature (K)
     * \param[in] natom        Number of atoms
     * \param[in] linear       TRUE if this is a linear molecule
     * \param[in] theta        The principal moments of inertia (unit of Energy)
     * \param[in] sigma_r      Symmetry factor, should be >= 1
     * \returns The rotational entropy (J/mol K)
     */
    double rotationalEntropy(double      temperature,
                             int         natom,
                             bool        linear,
                             const rvec &theta,
                             double      sigma_r);
    
public:
    /*! \brief Constructor that computes everything at once
     * TODO add comments
     */
    ThermoChemistry(const MyMol               *mymol,
                    const std::vector<double> &frequencies,
                    double                     temperature,
                    double                     pressure,
                    double                     scale_factor);

    //! return the zero point energy
    double ZPE() const { return zpe_; }

    //! return the heat capacity
    double cv(TCComponent tcc) const { return cv_.find(tcc)->second; }
    
    //! return the internal energy
    double Einternal(TCComponent tcc) const { return E_.find(tcc)->second; }
    
    //! return the standard entropy
    double S0(TCComponent tcc) const { return S0_.find(tcc)->second; }
    
    //! return the enthalpy of formation
    double DHform() const { return dhForm_; }
};

}

#endif

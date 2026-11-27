/*
 * This file is part of the GROMACS molecular simulation package.
 *
 * Copyright (c) 2017,2018,2019,2020,2021,2023 by the GROMACS development team, led by
 * Mark Abraham, David van der Spoel, Berk Hess, and Erik Lindahl,
 * and including many others, as listed in the AUTHORS file in the
 * top-level source directory and at http://www.gromacs.org.
 * Copyright (C) 2026, the Alexandria Chemistry Toolkit developers.
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
#include "actpre.h"

#include <algorithm>
#include <cmath>
#include <cstdio>
#include <map>
#include <string>
#include <vector>

#include "act/basics/act_particle.h"
#include "act/basics/atomization_energy.h"
#include "act/basics/msg_handler.h"
#include "act/molprop/molpropobservable.h"
#include "act/properties/princ.h"
#include "act/properties/symmetrynumber.h"
#include "act/topology/actmol.h"
#include "act/topology/actmol_low.h"
#include "act/topology/topology.h"
#include "act/utility/units.h"
#include "gromacs/math/functions.h"
#include "gromacs/math/units.h"
#include "gromacs/math/vec.h"
#include "gromacs/math/vectypes.h"
#include "gromacs/utility/gmxassert.h"
#include "gromacs/utility/real.h"
#include "gromacs/utility/stringutil.h"
#include "thermochemistry.h"

#define universalGasConstant RGAS

namespace alexandria
{

//! \brief Map from enum to string
std::map<TCComponent, std::string> tccMap = {
    { TCComponent::Translation, "Translation" },
    { TCComponent::Rotation,    "Rotation"    },
    { TCComponent::Vibration,   "Vibration"   },
    { TCComponent::Total,       "Total"       }
};

const std::string &tcComponentName(TCComponent tc)
{
    return tccMap[tc];
}

const std::map<TCComponent, std::string> &tccmap()
{
    return tccMap;
}

double ThermoChemistry::zeroPointEnergy(const std::vector<double> &frequencies,
                                        double                     scale_factor)
{
    // Convert frequency (ps^-1) to energy (kJ/mol)
    double factor = PLANCK;
    double zpe    = 0;
    for (const auto &omega : frequencies)
    {
        // Imaginary modes are stored as negative frequencies and,
        // like in calcVibrationalProperties, do not contribute.
        if (omega > 0)
        {
            zpe += 0.5 * factor * scale_factor * omega;
        }
    }
    return zpe;
}

void ThermoChemistry::calcVibrationalProperties(MsgHandler                *msghandler,
                                                const std::vector<double> &frequencies,
                                                double                     temperature,
                                                double                     scale_factor)
{
    for (size_t i = 0; i < frequencies.size(); i++)
    {
        if (frequencies[i] > 0)
        {
            // Scaled frequency
            double sf    = scale_factor*frequencies[i];
            // Vibrational temperature
            // https://en.wikipedia.org/wiki/Vibrational_temperature
            double theta = PLANCK*sf/BOLTZ;
            if (msghandler && msghandler->debug())
            {
                msghandler->writeDebug(gmx::formatString("scaled freq %10g (1/ps) %10g (1/cm) theta %10g\n", sf, 
                                                         convertFromGromacs(sf, mpo_unit2(MolPropObservable::FREQUENCY)),
                                                         theta));
            }
            // The zero point energy is part of the internal energy at
            // all temperatures, including zero.
            E_[TCComponent::Vibration] += 0.5 * BOLTZ * theta;
            if (temperature > 0)
            {
                double tt = theta / temperature;
                // For larger numbers the thermal contributions are nil.
                if (tt < 100)
                {
                    E_[TCComponent::Vibration]  += BOLTZ * theta / std::expm1(tt);
                    cv_[TCComponent::Vibration] += KILO * BOLTZ * std::exp(tt) * gmx::square(tt / std::expm1(tt));
                    S0_[TCComponent::Vibration] += KILO * BOLTZ *  (tt/ std::expm1(tt) - std::log1p(-std::exp(-tt)));
                }
            }
        }
    }
}

double ThermoChemistry::translationalEntropy(double mass,
                                             double temperature,
                                             double pressure)
{
    if (temperature <= 0)
    {
        // The equation for the ideal gas does not hold here
        return 0;
    }
    double kT = BOLTZ * temperature;

    GMX_RELEASE_ASSERT(mass > 0, "Molecular mass should be larger than zero");
    GMX_RELEASE_ASSERT(pressure > 0, "Pressure should be larger than zero");
    // Convert bar to Pascal
    double P   = pressure * 1e5;
    double qT  = (std::pow(2 * M_PI * mass * kT / gmx::square(PLANCK), 1.5) * (kT / P)
                  * (1e30 / AVOGADRO));
    return universalGasConstant * (2.5 + std::log(qT));
}

double ThermoChemistry::rotationalEntropy(double      temperature,
                                          int         natom,
                                          bool        linear,
                                          const rvec &theta,
                                          double      sigma_r)
{
    if (temperature <= 0 || natom <= 1)
    {
        // The equation for the rigid rotor does not hold at zero
        // temperature, and an atom does not rotate.
        return 0;
    }
    GMX_RELEASE_ASSERT(sigma_r > 0, "Symmetry factor should be larger than zero");
    if (linear)
    {
        GMX_RELEASE_ASSERT(theta[0] > 0, "Theta should be larger than zero");
        return universalGasConstant * (1 + std::log(temperature / (sigma_r * theta[0])));
    }
    double Q = theta[XX] * theta[YY] * theta[ZZ];
    GMX_RELEASE_ASSERT(Q > 0, "Q should be larger than zero");
    double qR = std::sqrt(M_PI * std::pow(temperature, 3) / Q) / sigma_r;
    return universalGasConstant * (1.5 + std::log(qR));
}

//! \brief Compute rotational constants
static void calcTheta(MsgHandler                   *msghandler,
                      const ACTMol                 *actmol,
                      const std::vector<gmx::RVec> &coords,
                      rvec                          theta)
{
    gmx::RVec              xcm;
    const auto            &atoms = actmol->topology()->atoms();
    std::vector<int>       index;
    std::vector<gmx::RVec> xx;
    for(size_t i = 0; i < atoms.size(); i++)
    {
        if (atoms[i].pType() == ActParticle::Atom)
        {
            index.push_back(i);
        }
        rvec zzz;
        for(int m = 0; m < DIM; m++)
        {
            zzz[m] = coords[i][m];
        }
        xx.emplace_back(zzz);
    }
    
    (void) calc_xcm(xx, index, atoms, &xcm, false);
    std::vector<gmx::RVec> x_com;
    x_com.resize(atoms.size());
    clear_rvec(theta);
    std::vector<real> mass;
    for (size_t i = 0; i < atoms.size(); i++)
    {
        copy_rvec(coords[i], x_com[i]);
        mass.push_back(atoms[i].mass());
    }
    (void)sub_xcm(&x_com, index, atoms, &xcm, false);

    gmx::RVec inertia;
    matrix    trans;
    principal_comp(index, mass, x_com, &trans, &inertia);
    // (kJ/mol ps)^2/(Dalton nm^2 kJ/mol K) =
    // c_kilo kg m^2 ps^2/(s^2 mol g/mol nm^2 K) =
    // c_kilo^2 10^18 / 10^24 K = 1/K
    double rot_const = gmx::square(PLANCK) / (8 * gmx::square(M_PI) * BOLTZ);
    // Rotational temperature (1/K)
    if (actmol->linearMolecule())
    {
        // For linear molecules the first element of the inertia
        // vector is zero.
        theta[0] = rot_const / inertia[1];
    }
    else
    {
        for (int m = 0; m < DIM; m++)
        {
            theta[m] = rot_const / inertia[m];
        }
    }
    if (msghandler && msghandler->debug())
    {
        msghandler->writeDebug(gmx::formatString("Rotational temperatures (Kelvin)%12.5f%12.5f%12.5f\n", theta[XX], theta[YY], theta[ZZ]));
    }
}

/*! \brief Determine the rotational symmetry number
 * The symmetry number is computed from the structure. If that fails, the
 * value stored for the compound (e.g. from the database) is used.
 * \param[in] msghandler For warnings, may be nullptr
 * \param[in] actmol     The molecule
 * \param[in] coords     The coordinates
 * \param[in] warn       Whether to warn about differences
 * \return The symmetry number
 */
static int symmetryNumberFromStructure(MsgHandler                   *msghandler,
                                       const ACTMol                 *actmol,
                                       const std::vector<gmx::RVec> &coords,
                                       bool                          warn)
{
    int stored = 1;
    if (actmol->fragments().size() == 1)
    {
        stored = std::max(1, actmol->fragments()[0].symmetryNumber());
    }
    // Only real atoms, not shells or virtual sites
    std::vector<int>       atomicNumbers;
    std::vector<double>    masses;
    std::vector<gmx::RVec> x;
    const auto            &atoms = actmol->topology()->atoms();
    for(size_t i = 0; i < atoms.size(); i++)
    {
        if (atoms[i].pType() == ActParticle::Atom)
        {
            atomicNumbers.push_back(atoms[i].atomicNumber());
            masses.push_back(atoms[i].mass());
            x.push_back(coords[i]);
        }
    }
    int computed = rotationalSymmetryNumber(atomicNumbers, masses, x);
    if (computed == 0)
    {
        if (warn && msghandler)
        {
            msghandler->msg(ACTStatus::Warning,
                            gmx::formatString("Could not compute the rotational symmetry number of %s, using %d",
                                              actmol->getMolname().c_str(), stored));
        }
        return stored;
    }
    if (computed != stored && warn && msghandler)
    {
        msghandler->msg(ACTStatus::Warning,
                        gmx::formatString("Rotational symmetry number of %s computed from the structure is %d, the stored value is %d. Using %d.",
                                          actmol->getMolname().c_str(), computed, stored, computed));
    }
    return computed;
}

ThermoChemistry::ThermoChemistry(MsgHandler                   *msghandler,
                                 const ACTMol                 *actmol,
                                 const std::vector<gmx::RVec> &coords,
                                 const AtomizationEnergy      &atomenergy,
                                 const std::vector<double>    &frequencies,
                                 double                        epot,
                                 double                        temperature,
                                 double                        pressure,
                                 double                        scale_factor)
{
    // Shells and virtual sites are not atoms
    int natom = 0;
    for(const auto &a : actmol->atomsConst())
    {
        if (a.pType() == ActParticle::Atom)
        {
            natom += 1;
        }
    }
    // An atom does not rotate
    bool rotates = natom > 1;
    rvec theta   = { 0, 0, 0 };
    if (rotates)
    {
        calcTheta(msghandler, actmol, coords, theta);
    }
    symmetryNumber_ = symmetryNumberFromStructure(msghandler, actmol, coords, temperature > 0);
    double sigma_r  = symmetryNumber_;
    zpe_ = zeroPointEnergy(frequencies, scale_factor);
    for(const auto &tc: tccmap())
    {
        cv_[tc.first] = 0;
        S0_[tc.first] = 0;
        E_[tc.first]  = 0;
    }
    // The classical heat capacities of translation and rotation do
    // not hold at zero temperature, where there is no entropy either.
    if (temperature > 0)
    {
        double cvRot = actmol->linearMolecule() ? RGAS : 1.5*RGAS;
        cv_[TCComponent::Translation] = 1.5*RGAS;
        E_[TCComponent::Translation]  = 1.5*RGAS*temperature/KILO;
        if (rotates)
        {
            cv_[TCComponent::Rotation] = cvRot;
            E_[TCComponent::Rotation]  = cvRot*temperature/KILO;
        }
    }
    S0_[TCComponent::Translation]  = translationalEntropy(actmol->totalMass(), temperature, pressure);
    S0_[TCComponent::Rotation]     = rotationalEntropy(temperature, natom, actmol->linearMolecule(),
                                                       theta, sigma_r);
    // Compute componenents of thermochemistry
    calcVibrationalProperties(msghandler, frequencies, temperature, scale_factor);
    // and sum them
    for(const auto &tc: tccmap())
    {
        if (tc.first != TCComponent::Total)
        {
            cv_[TCComponent::Total] += cv_[tc.first];
            S0_[TCComponent::Total] += S0_[tc.first];
            E_[TCComponent::Total]  += E_[tc.first];
        }
    }
    // Using the equations from Ochterski2000a. The energy of the atoms
    // is zero, that is the potential energy has to be relative to the
    // separated atoms.
    double sumAtomicEps0 = 0;
    double D0M           = sumAtomicEps0 - epot - zpe_;
    // Enthalpy of formation of the atoms at zero K
    double sumAtomicH0   = computeAtomizationEnergy(msghandler, actmol->atomsConst(), atomenergy, 0);
    double dhF0          = sumAtomicH0 - D0M;
    // What the enthalpy of the compound gains from zero K to this
    // temperature. The internal energy includes the zero point energy.
    double dHcompound    = E_[TCComponent::Total] - zpe_ + BOLTZ*temperature;
    // And the same for the elements in their standard states
    double dHelements    = computeElementEnthalpy(msghandler, actmol->atomsConst(), atomenergy, temperature);
    dhForm_ = dhF0 + dHcompound - dHelements;
}

}

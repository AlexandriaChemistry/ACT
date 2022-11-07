/*
 * This file is part of the GROMACS molecular simulation package.
 *
 * Copyright (c) 2017,2018,2019,2020,2021, by the GROMACS development team, led by
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
#include "actpre.h"

#include "thermochemistry.h"

#include <map>

#include <cmath>
#include <cstdio>

#include "act/utility/units.h"
#include "alexandria/mymol_low.h"
#include "alexandria/princ.h"
#include "gromacs/math/units.h"
#include "gromacs/utility/arrayref.h"
#include "gromacs/utility/fatalerror.h"
#include "gromacs/utility/gmxassert.h"

#define universalGasConstant RGAS

namespace alexandria
{

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
        zpe += 0.5 * factor * scale_factor * omega;
    }
    return zpe;
}

void ThermoChemistry::calcVibrationalProperties(const std::vector<double> &frequencies,
                                                double                     temperature,
                                                double                     scale_factor)
{
    if (temperature == 0)
    {
        return;
    }
    real T_1 = 1.0/temperature;
    //double hbar  = PLANCK1 / (2 * M_PI);
    for (size_t i = 0; i < frequencies.size(); i++)
    {
        if (frequencies[i] > 0)
        {
            // Scaled frequency
            double sf    = scale_factor*frequencies[i];
            // Vibrational temperature
            // https://en.wikipedia.org/wiki/Vibrational_temperature
            double theta = PLANCK*sf/BOLTZ;
            if (debug)
            {
                fprintf(debug, "scaled freq %10g (1/ps) %10g (1/cm) theta %10g\n", sf, 
                        convertFromGromacs(sf, mpo_unit2(MolPropObservable::FREQUENCY)),
                        theta);
            }
            // Prevent overflow by checking for unreasonably large numbers.
            real tt = theta * T_1;
            if (tt < 100)
            {
                E_[TCComponent::Vibration]  += BOLTZ * theta * (0.5 + 1.0 / (std::expm1(tt)));
                cv_[TCComponent::Vibration] += KILO * BOLTZ * std::exp(tt) * gmx::square(tt / std::expm1(tt));
                S0_[TCComponent::Vibration] += KILO * BOLTZ *  (tt/ std::expm1(tt) - std::log1p(-std::exp(-tt)));
            }
        }
    }
}

double ThermoChemistry::translationalEntropy(double mass,
                                             double temperature,
                                             double pressure)
{
    double ST = 2.5;
    if (temperature > 0)
    {
        double kT = BOLTZ * temperature;
        
        GMX_RELEASE_ASSERT(mass > 0, "Molecular mass should be larger than zero");
        GMX_RELEASE_ASSERT(pressure > 0, "Pressure should be larger than zero");
        // Convert bar to Pascal
        double P   = pressure * 1e5;
        double qT  = (std::pow(2 * M_PI * mass * kT / gmx::square(PLANCK), 1.5) * (kT / P)
                      * (1e30 / AVOGADRO));
        ST        += std::log(qT);
    }
    return universalGasConstant * ST;
}

double ThermoChemistry::rotationalEntropy(double      temperature,
                                          int         natom,
                                          bool        linear,
                                          const rvec &theta,
                                          double      sigma_r)
{
    double SR = 1.5;
    if (linear)
    {
        SR = 1;
    }
    if (temperature > 0)
    {
        GMX_RELEASE_ASSERT(sigma_r > 0, "Symmetry factor should be larger than zero");
        
        if (natom > 1)
        {
            double qR;
            if (linear)
            {
                GMX_RELEASE_ASSERT(theta[0] > 0, "Theta should be larger than zero");
                qR  = temperature / (sigma_r * theta[0]);
            }
            else
            {
                double Q = theta[XX] * theta[YY] * theta[ZZ];
                GMX_RELEASE_ASSERT(Q > 0, "Q should be larger than zero");
                qR  = std::sqrt(M_PI * std::pow(temperature, 3) / Q) / sigma_r;
            }
            SR += std::log(qR);
        }
    }
    return universalGasConstant * SR;
}

static void calcTheta(const MyMol                  *mymol,
                      const std::vector<gmx::RVec> &coords,
                      rvec                          theta)
{
    rvec       xcm;
    const auto atoms = mymol->gmxAtomsConst();
    std::vector<int> index;
    std::vector<gmx::RVec> xx;
    for(int i = 0; i < atoms->nr; i++)
    {
        if (atoms->atom[i].ptype == eptAtom)
        {
            index.push_back(i);
        }
        rvec zzz;
        for(int m = 0; m < DIM; m++)
        {
            zzz[m] = coords[i][m];
        }
        xx.push_back(zzz);
    }
    
    (void) calc_xcm(as_rvec_array(xx.data()),
                    index.size(), index.data(), atoms->atom, xcm, false);
    std::vector<gmx::RVec> x_com;
    x_com.resize(atoms->nr);
    clear_rvec(theta);
    for (int i = 0; i < atoms->nr; i++)
    {
        copy_rvec(coords[i], x_com[i]);
    }
    (void)sub_xcm(as_rvec_array(x_com.data()), index.size(), index.data(), 
                  atoms->atom, xcm, false);

    rvec   inertia;
    matrix trans;
    principal_comp(index.size(), index.data(), atoms->atom, as_rvec_array(x_com.data()),
                   trans, inertia);
    // (kJ/mol ps)^2/(Dalton nm^2 kJ/mol K) =
    // c_kilo kg m^2 ps^2/(s^2 mol g/mol nm^2 K) =
    // c_kilo^2 10^18 / 10^24 K = 1/K
    double rot_const = gmx::square(PLANCK) / (8 * gmx::square(M_PI) * BOLTZ);
    // Rotational temperature (1/K)
    if (mymol->linearMolecule())
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
    if (debug)
    {
        fprintf(debug, "Rotational temperatures (Kelvin)%12.5f%12.5f%12.5f\n", theta[XX], theta[YY], theta[ZZ]);
    }
}

ThermoChemistry::ThermoChemistry(const MyMol                  *mymol,
                                 const std::vector<gmx::RVec> &coords,
                                 const AtomizationEnergy      &atomenergy,
                                 const std::vector<double>    &frequencies,
                                 double                        temperature,
                                 double                        pressure,
                                 double                        scale_factor)
{
    rvec   theta;
    calcTheta(mymol, coords, theta);
    double sigma_r = 1;
    if (mymol->fragments().size() == 1)
    {
        sigma_r = mymol->fragments()[0].symmetryNumber();
    }
    zpe_ = zeroPointEnergy(frequencies, scale_factor);
    cv_[TCComponent::Translation]  = 1.5*RGAS;
    cv_[TCComponent::Rotation]     = mymol->linearMolecule() ? RGAS : 1.5*RGAS;
    E_[TCComponent::Translation]   = 1.5*RGAS*temperature/KILO;
    E_[TCComponent::Rotation]      = (mymol->linearMolecule() ? RGAS : 1.5*RGAS) * temperature/KILO;
    S0_[TCComponent::Translation]  = translationalEntropy(mymol->totalMass(), temperature, pressure);
    S0_[TCComponent::Rotation]     = rotationalEntropy(temperature, mymol->NAtom(), mymol->linearMolecule(),
                                                       theta, sigma_r);
    // Compute componenents of thermochemistry
    calcVibrationalProperties(frequencies, temperature, scale_factor);
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
    // Using the equations from Ochterski2000a
    double sumAtomicEps0 = 0;
    double D0M           = sumAtomicEps0 - mymol->energyTerms()[F_EPOT] - zpe_;
    double sumAtomicH0   = computeAtomizationEnergy(mymol->atomsConst(), atomenergy, 0);
    double dhF0          = sumAtomicH0 - D0M;
    if (temperature == 0)
    {
        dhForm_ = dhF0;
    }
    else
    {
        dhForm_ = dhF0 + E_[TCComponent::Total] - (computeAtomizationEnergy(mymol->atomsConst(), atomenergy, temperature) - sumAtomicH0);
    }
}

}

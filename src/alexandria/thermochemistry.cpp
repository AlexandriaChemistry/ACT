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
    double factor = PLANCK; // (2.0 * M_PI);
    double zpe    = 0;
    for (const auto &omega : frequencies)
    {
        zpe += 0.5 * factor * scale_factor * convertToGromacs(omega, mpo_unit2(MolPropObservable::FREQUENCY));
    }
    return zpe;
}

void ThermoChemistry::calcVibrationalProperties(const std::vector<double> &frequencies,
                                                double                     temperature,
                                                double                     scale_factor)
{
    double hbar  = PLANCK1 / (2 * M_PI);
    for (size_t i = 0; i < frequencies.size(); i++)
    {
        if (frequencies[i] > 0)
        {
            double omega = scale_factor * convertToGromacs(frequencies[i], mpo_unit2(MolPropObservable::FREQUENCY));
            //double omega = scale_factor * frequencies[i];
            // * eigval_to_frequency(eigval[i]);
            double hwkT  = (hbar * omega) / (PICO * BOLTZMANN * temperature);
            if (debug)
            {
                fprintf(debug, "freq %g omega %g hwkT %g\n", frequencies[i], omega, hwkT);
            }
            // Prevent overflow by checking for unreasonably large numbers.
            if (hwkT < 100)
            {
                E_[TCComponent::Vibration]  += temperature * BOLTZ * hwkT * (0.5 + 1.0 / (std::expm1(hwkT)));
                cv_[TCComponent::Vibration] += BOLTZ * std::exp(hwkT) * gmx::square(hwkT / std::expm1(hwkT));
                S0_[TCComponent::Vibration] += BOLTZ *  (hwkT / std::expm1(hwkT) - std::log1p(-std::exp(-hwkT)));
            }
        }
    }
}

double ThermoChemistry::translationalEntropy(double mass,
                                             double temperature,
                                             double pressure)
{
    double kT = BOLTZ * temperature;
    
    GMX_RELEASE_ASSERT(mass > 0, "Molecular mass should be larger than zero");
    GMX_RELEASE_ASSERT(pressure > 0, "Pressure should be larger than zero");
    GMX_RELEASE_ASSERT(temperature > 0, "Temperature should be larger than zero");
    // Convert bar to Pascal
    double P  = pressure * 1e5;
    double qT = (std::pow(2 * M_PI * mass * kT / gmx::square(PLANCK), 1.5) * (kT / P)
                 * (1e30 / AVOGADRO));
    return universalGasConstant * (std::log(qT) + 2.5);
}

double ThermoChemistry::rotationalEntropy(double      temperature,
                                          int         natom,
                                          bool        linear,
                                          const rvec &theta,
                                          double      sigma_r)
{
    GMX_RELEASE_ASSERT(sigma_r > 0, "Symmetry factor should be larger than zero");
    GMX_RELEASE_ASSERT(temperature > 0, "Temperature should be larger than zero");

    double sR = 0;
    if (natom > 1)
    {
        if (linear)
        {
            GMX_RELEASE_ASSERT(theta[0] > 0, "Theta should be larger than zero");
            double qR = temperature / (sigma_r * theta[0]);
            sR        = universalGasConstant * (std::log(qR) + 1);
        }
        else
        {
            double Q = theta[XX] * theta[YY] * theta[ZZ];
            GMX_RELEASE_ASSERT(Q > 0, "Q should be larger than zero");
            double qR = std::sqrt(M_PI * std::pow(temperature, 3) / Q) / sigma_r;
            sR        = universalGasConstant * (std::log(qR) + 1.5);
        }
    }
    return sR;
}

static void calcTheta(const MyMol *mymol,
                      rvec         theta)
{
    rvec   xcm;
    auto  &atoms = mymol->atomsConst();
    std::vector<int> index;
    std::vector<gmx::RVec> xx;
    for(int i = 0; i < atoms.nr; i++)
    {
        if (atoms.atom[i].ptype == eptAtom)
        {
            index.push_back(i);
        }
        rvec zzz;
        for(int m = 0; m < DIM; m++)
        {
            zzz[m] = mymol->x()[i][m];
        }
        xx.push_back(zzz);
    }
    
    (void) calc_xcm(as_rvec_array(xx.data()),
                    index.size(), index.data(), atoms.atom, xcm, false);
    std::vector<gmx::RVec> x_com;
    x_com.resize(atoms.nr);
    clear_rvec(theta);
    for (int i = 0; i < atoms.nr; i++)
    {
        copy_rvec(mymol->x()[i], x_com[i]);
    }
    (void)sub_xcm(as_rvec_array(x_com.data()), index.size(), index.data(), 
                  atoms.atom, xcm, false);

    rvec   inertia;
    matrix trans;
    principal_comp(index.size(), index.data(), atoms.atom, as_rvec_array(x_com.data()),
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
}

ThermoChemistry::ThermoChemistry(const MyMol               *mymol,
                                 const std::vector<double> &frequencies,
                                 double                     temperature,
                                 double                     pressure,
                                 double                     scale_factor)
{
    rvec   theta;
    calcTheta(mymol, theta);
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
    dhForm_ = mymol->energyTerms()[F_EPOT] + E_[TCComponent::Total];

}

}

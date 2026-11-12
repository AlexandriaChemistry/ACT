/*
 * This source file is part of the Alexandria Chemistry Toolkit.
 *
 * Copyright (C) 2026
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
/*! \internal \file
 * \brief
 * The rows of the parameter tables of the GPU kernels, for the processor and the GPU.
 *
 * Each row has 8 numbers, made from the parameters of an interaction in
 * internal units, see gpuPairEntry and gpuBondedEntry in forces_gpu.cpp,
 * which use these functions. They exist only here, such that the GPU can
 * make the rows itself from the parameters, with the same results. They
 * only use arithmetic, in the same order as before, and neither throw nor
 * use containers. The rows of Coulomb, made from the charges, are not here.
 *
 * \author David van der Spoel <david.vanderspoel@icm.uu.se>
 * \ingroup module_gpu
 */
#ifndef ACT_GPU_PAIRROWS_H
#define ACT_GPU_PAIRROWS_H

#include <cstddef>

#include "act/forcefield/forcefield_parametername.h"
#include "act/gpu/gpu.h"

namespace alexandria
{

namespace gpu
{

namespace rows
{

//! \return x to the power 6, as gmx::power6
inline double power6(double x)
{
    double x3 = x*(x*x);
    return x3*x3;
}

/*! \brief The row of an interaction of a pair, except Coulomb
 * \param[in]  kernel  The kernel of the pair
 * \param[in]  params  The parameters of the interaction, in internal units
 * \param[in]  qshell  For Polarization, the charge of the shell, the second particle
 * \param[in]  coulomb For Polarization, the constant of Coulomb's law
 * \param[out] p       The row, 8 numbers, the unused ones zero
 * \return whether the pair interacts, if not its energy is zero
 */
inline bool pairRow(PairKernel kernel, const double *params, double qshell, double coulomb, float *p)
{
    for (int k = 0; k < 8; k++)
    {
        p[k] = 0.0f;
    }
    switch (kernel)
    {
    case PairKernel::Coulomb:
        // From the charges, see gpuPairEntry
        return false;
    case PairKernel::Polarization:
    {
        double alpha = params[polALPHA];
        if (alpha == 0)
        {
            return false;
        }
        p[0] = static_cast<float>(coulomb*qshell*qshell/alpha);
        p[1] = static_cast<float>(params[polRHYPER]);
        p[2] = static_cast<float>(params[polFCHYPER]);
        return true;
    }
    case PairKernel::LJ12_6:
    {
        double sig6 = power6(params[lj12_6SIGMA]);
        double c6   = 4*params[lj12_6EPSILON]*sig6;
        p[0] = static_cast<float>(c6);
        p[1] = static_cast<float>(c6*sig6);
        return true;
    }
    case PairKernel::LJ14_7:
    {
        double sigma   = params[lj14_7SIGMA];
        double epsilon = params[lj14_7EPSILON];
        p[0] = static_cast<float>(1/sigma);
        p[1] = static_cast<float>(epsilon);
        p[2] = static_cast<float>(params[lj14_7GAMMA]);
        p[3] = static_cast<float>(params[lj14_7DELTA]);
        return epsilon > 0 && sigma > 0;
    }
    case PairKernel::WangBuckingham:
    {
        double sigma   = params[wbhSIGMA];
        double epsilon = params[wbhEPSILON];
        double gamma   = params[wbhGAMMA];
        double gamma_3 = 3.0/(gamma + 3.0);
        p[0] = static_cast<float>(power6(sigma));
        p[1] = static_cast<float>(gamma_3);
        p[2] = static_cast<float>(2*epsilon/(1.0 - gamma_3));
        p[3] = static_cast<float>(gamma);
        p[4] = static_cast<float>(1/sigma);
        return epsilon > 0 && gamma > 0 && sigma > 0;
    }
    case PairKernel::DoubleExponential:
        p[0] = static_cast<float>(params[dexpA2] - params[dexpA1]);
        p[1] = static_cast<float>(params[dexpB]);
        return true;
    case PairKernel::HarmonicBond:
        p[0] = static_cast<float>(params[bondLENGTH]);
        p[1] = static_cast<float>(params[bondKB]);
        p[2] = static_cast<float>(params[bondENERGY]);
        return true;
    case PairKernel::CubicBond:
    {
        double b0   = params[cubicLENGTH];
        double rmax = params[cubicRMAX];
        // Not beyond the maximum of the potential, as CubicBond::prepare
        double rhi  = (2*rmax + b0)/3;
        p[0] = static_cast<float>(b0);
        p[1] = static_cast<float>(rmax);
        p[2] = static_cast<float>(params[cubicKB]);
        p[3] = static_cast<float>(params[cubicDE]);
        p[4] = static_cast<float>(rhi*rhi);
        return true;
    }
    case PairKernel::MorseBond:
        p[0] = static_cast<float>(params[morseLENGTH]);
        p[1] = static_cast<float>(params[morseBETA]);
        p[2] = static_cast<float>(params[morseDE]);
        p[3] = static_cast<float>(params[morseD0]);
        return true;
    case PairKernel::HuaBond:
        p[0] = static_cast<float>(params[huaLENGTH]);
        p[1] = static_cast<float>(params[huaB]);
        p[2] = static_cast<float>(params[huaDE]);
        p[3] = static_cast<float>(params[huaC]);
        return true;
    }
    return false;
}

/*! \brief The row of a bonded interaction of three or four particles
 * \param[in]  kernel  The kernel
 * \param[in]  params  The parameters of the interaction, in internal units
 * \param[in]  nparams The number of parameters
 * \param[out] p       The row, 8 numbers, the unused ones zero
 */
inline void bondedRow(BondedKernel kernel, const double *params, size_t nparams, float *p)
{
    for (int k = 0; k < 8; k++)
    {
        p[k] = 0.0f;
    }
    switch (kernel)
    {
    case BondedKernel::HarmonicAngle:
        p[0] = static_cast<float>(params[angleKT]);
        p[1] = static_cast<float>(params[angleANGLE]);
        break;
    case BondedKernel::LinearAngle:
        p[0] = static_cast<float>(params[linangA]);
        p[1] = static_cast<float>(params[linangKLIN]);
        break;
    case BondedKernel::FourierDihedral:
    {
        size_t n = nparams < static_cast<size_t>(fdihNR) ? nparams : static_cast<size_t>(fdihNR);
        for (size_t k = 0; k < n; k++)
        {
            p[k] = static_cast<float>(params[k]);
        }
        p[6] = static_cast<float>(n);
        break;
    }
    case BondedKernel::HarmonicDihedral:
        p[0] = static_cast<float>(params[idihKPHI]);
        break;
    }
}

} // namespace rows

} // namespace gpu

} // namespace alexandria

#endif

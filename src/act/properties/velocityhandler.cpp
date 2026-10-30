/*
 * This source file is part of the Alexandria Chemistry Toolkit.
 *
 * Copyright (C) 2014-2022
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
#include "actpre.h"

#include "velocityhandler.h"

#include <algorithm>
#include <cmath>
#include <map>

#include "act/topology/topology.h"
#include "act/utility/rotation.h"
#include "gromacs/math/units.h"
#include "gromacs/math/vec.h"
#include "gromacs/math/vectypes.h"
#include "gromacs/random/tabulatednormaldistribution.h"
#include "gromacs/random/threefry.h"
#include "gromacs/utility/exceptions.h"
#include "gromacs/utility/fatalerror.h"
#include "gromacs/utility/gmxassert.h"
#include "gromacs/utility/stringutil.h"
#include "gromacs/utility/textwriter.h"

namespace alexandria
{

//! \brief Generate velocities according to Maxwell distribution
static void low_mspeed(real                        tempi,
                       const std::vector<ActAtom> &atoms,
                       std::vector<gmx::RVec>     *v,
                       gmx::ThreeFry2x64<>        *rng)
{
    int                                    nrdf;
    real                                   ekin, temp;
    gmx::TabulatedNormalDistribution<real> normalDist;

    ekin = 0.0;
    nrdf = 0;
    for (size_t i = 0; i < atoms.size(); i++)
    {
        real mass = atoms[i].mass();
        if (mass > 0)
        {
            rng->restart(i, 0);
            real sd = std::sqrt(BOLTZ * tempi / mass);
            for (int m = 0; (m < DIM); m++)
            {
                auto vvv    = sd * normalDist(*rng);
                (*v)[i][m]  = vvv;
                ekin       += 0.5 * mass * vvv * vvv;
            }
            nrdf += DIM;
        }
    }
    temp = (2.0 * ekin) / (nrdf * BOLTZ);
    if (temp > 0)
    {
        real scal = std::sqrt(tempi / temp);
        for (size_t i = 0; (i < atoms.size()); i++)
        {
            for (int m = 0; (m < DIM); m++)
            {
                (*v)[i][m] *= scal;
            }
        }
    }
}

void maxwell_speed(real                        tempi, 
                   unsigned int                seed, 
                   const std::vector<ActAtom> &atoms,
                   std::vector<gmx::RVec>     *v,
                   gmx::TextWriter            *tw)
{
    if (tempi <= 0)
    {
        return;
    }
    if (seed == 0)
    {
        seed = static_cast<int>(gmx::makeRandomSeed());
        if (tw)
        {
            tw->writeLineFormatted("Using random seed %u for generating velocities", seed);
        }
    }
    gmx::ThreeFry2x64<> rng(seed, gmx::RandomDomain::MaxwellVelocities);

    low_mspeed(tempi, atoms, v, &rng);
}

//! \brief Map from center of mass motion removal to string
static const std::map<COMMRemoval, std::string> cr2str = {
    { COMMRemoval::None,    "none"    },
    { COMMRemoval::Linear,  "linear"  },
    { COMMRemoval::Angular, "angular" }
};

const std::string &commRemovalToString(COMMRemoval cr)
{
    return cr2str.find(cr)->second;
}

COMMRemoval stringToCOMMRemoval(const std::string &str)
{
    for(const auto &c : cr2str)
    {
        if (c.second == str)
        {
            return c.first;
        }
    }
    GMX_THROW(gmx::InvalidInputError(gmx::formatString("Unknown center of mass motion removal '%s'", str.c_str())));
}

/*! \brief Compute the center of mass and the inertia tensor around it
 * \param[in]  mass The masses
 * \param[in]  x    The coordinates
 * \param[out] xcm  The center of mass
 * \param[out] I    The inertia tensor
 * \return the total mass
 */
static double centerOfMassInertia(const std::vector<double>    &mass,
                                  const std::vector<gmx::RVec> &x,
                                  dvec                          xcm,
                                  double                        I[DIM][DIM])
{
    double tm = 0;
    clear_dvec(xcm);
    for (size_t i = 0; i < mass.size(); i++)
    {
        tm += mass[i];
        for (int m = 0; m < DIM; m++)
        {
            xcm[m] += mass[i]*x[i][m];
        }
    }
    for (int m = 0; m < DIM; m++)
    {
        xcm[m] /= tm;
        for (int n = 0; n < DIM; n++)
        {
            I[m][n] = 0;
        }
    }
    for (size_t i = 0; i < mass.size(); i++)
    {
        dvec dx;
        for (int m = 0; m < DIM; m++)
        {
            dx[m] = x[i][m] - xcm[m];
        }
        double dx2 = dx[XX]*dx[XX] + dx[YY]*dx[YY] + dx[ZZ]*dx[ZZ];
        for (int m = 0; m < DIM; m++)
        {
            for (int n = 0; n < DIM; n++)
            {
                I[m][n] += mass[i]*((m == n ? dx2 : 0) - dx[m]*dx[n]);
            }
        }
    }
    return tm;
}

COMMRemover::COMMRemover(COMMRemoval                   mode,
                         const std::vector<ActAtom>   &atoms,
                         const std::vector<gmx::RVec> &x) : mode_(mode)
{
    int nMass = 0;
    for (const auto &a : atoms)
    {
        mass_.push_back(a.mass() > 0 ? a.mass() : 0);
        nMass += a.mass() > 0 ? 1 : 0;
    }
    if (nMass <= 1)
    {
        nRotation_ = 0;
    }
    else
    {
        dvec   xcm;
        double I[DIM][DIM], eval[DIM], evec[DIM][DIM];
        centerOfMassInertia(mass_, x, xcm, I);
        diagonalizeSymmetric(I, eval, evec);
        // A linear molecule has a vanishing moment of inertia around its axis
        nRotation_ = (eval[0] < 1e-6*eval[DIM-1]) ? 2 : 3;
    }
}

int COMMRemover::nDOF() const
{
    switch (mode_)
    {
    case COMMRemoval::None:
        return 0;
    case COMMRemoval::Linear:
        return DIM;
    case COMMRemoval::Angular:
        return DIM + nRotation_;
    }
    return 0;
}

void COMMRemover::remove(const std::vector<gmx::RVec> &x,
                         std::vector<gmx::RVec>       *v) const
{
    if (COMMRemoval::None == mode_)
    {
        return;
    }
    GMX_RELEASE_ASSERT(x.size() == mass_.size() && v->size() == mass_.size(),
                       "Number of particles does not match");
    double tm = 0;
    dvec   vcm = { 0, 0, 0 };
    for (size_t i = 0; i < mass_.size(); i++)
    {
        tm += mass_[i];
        for (int m = 0; m < DIM; m++)
        {
            vcm[m] += mass_[i]*(*v)[i][m];
        }
    }
    if (tm == 0)
    {
        return;
    }
    for (int m = 0; m < DIM; m++)
    {
        vcm[m] /= tm;
    }
    dvec omega = { 0, 0, 0 };
    dvec xcm   = { 0, 0, 0 };
    if (COMMRemoval::Angular == mode_ && nRotation_ > 0)
    {
        double I[DIM][DIM], eval[DIM], evec[DIM][DIM];
        centerOfMassInertia(mass_, x, xcm, I);
        // Angular momentum around the center of mass
        dvec L = { 0, 0, 0 };
        for (size_t i = 0; i < mass_.size(); i++)
        {
            dvec dx, dv, a;
            for (int m = 0; m < DIM; m++)
            {
                dx[m] = x[i][m] - xcm[m];
                dv[m] = (*v)[i][m] - vcm[m];
            }
            dcprod(dx, dv, a);
            for (int m = 0; m < DIM; m++)
            {
                L[m] += mass_[i]*a[m];
            }
        }
        // Solve I omega = L for the rotations that are removed, those
        // with the largest moments of inertia.
        diagonalizeSymmetric(I, eval, evec);
        for (int k = DIM - nRotation_; k < DIM; k++)
        {
            if (eval[k] > 0)
            {
                double c = diprod(evec[k], L)/eval[k];
                for (int m = 0; m < DIM; m++)
                {
                    omega[m] += c*evec[k][m];
                }
            }
        }
    }
    for (size_t i = 0; i < mass_.size(); i++)
    {
        if (mass_[i] > 0)
        {
            dvec dx, w;
            for (int m = 0; m < DIM; m++)
            {
                dx[m] = x[i][m] - xcm[m];
            }
            dcprod(omega, dx, w);
            for (int m = 0; m < DIM; m++)
            {
                (*v)[i][m] -= vcm[m] + w[m];
            }
        }
    }
}

} // namespace

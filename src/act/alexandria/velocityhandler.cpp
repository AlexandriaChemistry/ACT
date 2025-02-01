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

#include <cmath>

#include "gromacs/math/units.h"
#include "gromacs/math/vec.h"
#include "gromacs/math/vectypes.h"
#include "gromacs/random/tabulatednormaldistribution.h"
#include "gromacs/random/threefry.h"
#include "gromacs/utility/fatalerror.h"
#include "gromacs/utility/textwriter.h"

namespace alexandria
{

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

static real calc_cm(const std::vector<ActAtom>   &atoms,
                    const std::vector<gmx::RVec> &x,
                    std::vector<gmx::RVec>       *v,
                    rvec                          xcm,
                    rvec                          vcm,
                    rvec                          acm,
                    matrix                        L)
{
    rvec dx, a0;

    clear_rvec(xcm);
    clear_rvec(vcm);
    clear_rvec(acm);
    real tm = 0.0;
    for (size_t i = 0; (i < atoms.size()); i++)
    {
        real m0 = atoms[i].mass();
        tm += m0;
        cprod(x[i], (*v)[i], a0);
        for (int m = 0; (m < DIM); m++)
        {
            xcm[m] += m0 * x[i][m]; /* c.o.m. position */
            vcm[m] += m0 * (*v)[i][m]; /* c.o.m. velocity */
            acm[m] += m0 * a0[m];   /* rotational velocity around c.o.m. */
        }
    }
    cprod(xcm, vcm, a0);
    for (int m = 0; (m < DIM); m++)
    {
        xcm[m] /= tm;
        vcm[m] /= tm;
        acm[m] -= a0[m] / tm;
    }

    clear_mat(L);
    for (size_t i = 0; (i < atoms.size()); i++)
    {
        real m0 = atoms[i].mass();
        for (int m = 0; (m < DIM); m++)
        {
            dx[m] = x[i][m] - xcm[m];
        }
        L[XX][XX] += dx[XX] * dx[XX] * m0;
        L[XX][YY] += dx[XX] * dx[YY] * m0;
        L[XX][ZZ] += dx[XX] * dx[ZZ] * m0;
        L[YY][YY] += dx[YY] * dx[YY] * m0;
        L[YY][ZZ] += dx[YY] * dx[ZZ] * m0;
        L[ZZ][ZZ] += dx[ZZ] * dx[ZZ] * m0;
    }

    return tm;
}

void stop_cm(const std::vector<ActAtom>   &atoms,
             const std::vector<gmx::RVec> &x,
             std::vector<gmx::RVec>       *v)
{
    rvec   xcm, vcm, acm;
    tensor L;

    (void)calc_cm(atoms, x, v, xcm, vcm, acm, L);

    /* Subtract translational center of mass velocity */
    for (size_t i = 0; (i < atoms.size()); i++)
    {
        for (int m = 0; (m < DIM); m++)
        {
            (*v)[i][m] -= vcm[m];
        }
    }
}

} // namespace

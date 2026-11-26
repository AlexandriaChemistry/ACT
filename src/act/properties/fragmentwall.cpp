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

/*! \internal \brief
 * Implements part of the alexandria program.
 * \author David van der Spoel <david.vanderspoel@icm.uu.se>
 */

#include <algorithm>
#include <cmath>
#include <cstddef>
#include <vector>

#include "fragmentwall.h"
#include "gromacs/math/vectypes.h"
#include "gromacs/utility/exceptions.h"
#include "gromacs/utility/stringutil.h"

namespace alexandria
{

FragmentWall::FragmentWall(const std::vector<size_t> &atomStart,
                           const std::vector<double> &mass,
                           double                     radius,
                           double                     forceConstant) :
    atomStart_(atomStart), mass_(mass), radius_(radius), forceConstant_(forceConstant)
{
    if (!std::isfinite(radius) || radius < 0 || !std::isfinite(forceConstant) || forceConstant < 0)
    {
        GMX_THROW(gmx::InvalidInputError(gmx::formatString("The radius %g nm and the force constant %g kJ/mol nm^2 of the wall should be zero or more",
                                                           radius, forceConstant)));
    }
    if (atomStart_.size() < 2 || atomStart_.back() != mass_.size())
    {
        GMX_THROW(gmx::InvalidInputError(gmx::formatString("The fragments end at particle %zu, but there are %zu particles",
                                                           atomStart_.empty() ? 0 : atomStart_.back(), mass_.size())));
    }
    for(size_t f = 0; f + 1 < atomStart_.size(); f++)
    {
        double m = 0;
        for(size_t i = atomStart_[f]; i < atomStart_[f + 1]; i++)
        {
            m += mass_[i];
        }
        if (atomStart_[f + 1] <= atomStart_[f] || m <= 0)
        {
            GMX_THROW(gmx::InvalidInputError(gmx::formatString("Fragment %zu has no particles with a mass", f + 1)));
        }
        fragmentMass_.push_back(m);
    }
}

std::vector<gmx::DVec> FragmentWall::centers(const std::vector<gmx::DVec> &x) const
{
    if (x.size() != mass_.size())
    {
        GMX_THROW(gmx::InternalError(gmx::formatString("There are %zu coordinates for %zu particles", x.size(), mass_.size()).c_str()));
    }
    std::vector<gmx::DVec> com(nfragments(), { 0, 0, 0 });
    for(size_t f = 0; f < nfragments(); f++)
    {
        for(size_t i = atomStart_[f]; i < atomStart_[f + 1]; i++)
        {
            for(int m = 0; m < DIM; m++)
            {
                com[f][m] += mass_[i]*x[i][m]/fragmentMass_[f];
            }
        }
    }
    return com;
}

//! \return the length of a vector
static double length(const gmx::DVec &v)
{
    return std::sqrt(v[XX]*v[XX] + v[YY]*v[YY] + v[ZZ]*v[ZZ]);
}

double FragmentWall::distance(const std::vector<gmx::DVec> &x) const
{
    size_t nfrag = nfragments();
    if (nfrag < 2)
    {
        return 0;
    }
    const auto com = centers(x);
    if (2 == nfrag)
    {
        return length(com[1] - com[0]);
    }
    gmx::DVec all  = { 0, 0, 0 };
    double    mtot = 0;
    for(size_t f = 0; f < nfrag; f++)
    {
        mtot += fragmentMass_[f];
        for(int m = 0; m < DIM; m++)
        {
            all[m] += fragmentMass_[f]*com[f][m];
        }
    }
    double dmax = 0;
    for(size_t f = 0; f < nfrag; f++)
    {
        gmx::DVec d;
        for(int m = 0; m < DIM; m++)
        {
            d[m] = com[f][m] - all[m]/mtot;
        }
        dmax = std::max(dmax, length(d));
    }
    return dmax;
}

double FragmentWall::compute(const std::vector<gmx::DVec> &x,
                             std::vector<gmx::DVec>       *f) const
{
    const gmx::DVec zero = { 0, 0, 0 };
    f->assign(x.size(), zero);
    if (!active())
    {
        return 0;
    }
    size_t                 nfrag = nfragments();
    const auto             com   = centers(x);
    // The force on the center of mass of each fragment
    std::vector<gmx::DVec> fcom(nfrag, zero);
    double                 energy = 0;
    if (2 == nfrag)
    {
        const gmx::DVec d = com[1] - com[0];
        double          r = length(d);
        if (r > radius_)
        {
            energy = 0.5*forceConstant_*(r - radius_)*(r - radius_);
            for(int m = 0; m < DIM; m++)
            {
                double fm   = -forceConstant_*(r - radius_)*d[m]/r;
                fcom[1][m]  = fm;
                fcom[0][m]  = -fm;
            }
        }
    }
    else
    {
        gmx::DVec all  = zero;
        double    mtot = 0;
        for(size_t k = 0; k < nfrag; k++)
        {
            mtot += fragmentMass_[k];
            for(int m = 0; m < DIM; m++)
            {
                all[m] += fragmentMass_[k]*com[k][m];
            }
        }
        // The derivative of the energy to the distance vector of each
        // fragment to the center of mass of all, which moves with the
        // fragments
        gmx::DVec gsum = zero;
        for(size_t k = 0; k < nfrag; k++)
        {
            gmx::DVec d;
            for(int m = 0; m < DIM; m++)
            {
                d[m] = com[k][m] - all[m]/mtot;
            }
            double r = length(d);
            if (r > radius_)
            {
                energy += 0.5*forceConstant_*(r - radius_)*(r - radius_);
                for(int m = 0; m < DIM; m++)
                {
                    double g    = forceConstant_*(r - radius_)*d[m]/r;
                    fcom[k][m] -= g;
                    gsum[m]    += g;
                }
            }
        }
        for(size_t k = 0; k < nfrag; k++)
        {
            for(int m = 0; m < DIM; m++)
            {
                fcom[k][m] += gsum[m]*fragmentMass_[k]/mtot;
            }
        }
    }
    for(size_t k = 0; k < nfrag; k++)
    {
        for(size_t i = atomStart_[k]; i < atomStart_[k + 1]; i++)
        {
            for(int m = 0; m < DIM; m++)
            {
                (*f)[i][m] = fcom[k][m]*mass_[i]/fragmentMass_[k];
            }
        }
    }
    return energy;
}

} // namespace alexandria

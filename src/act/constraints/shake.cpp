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

#include "shake.h"

#include "act/basics/msg_handler.h"
#include "act/forcefield/forcefield_parametername.h"
#include "act/topology/actatom.h"

namespace alexandria
{

int Constrainer::shake(gmx_unused MsgHandler      *msghandler,
                       const TopologyEntryVector  &bonds,
                       const std::vector<ActAtom> &atoms,
                       std::vector<gmx::RVec>     *positions,
                       gmx_unused std::vector<gmx::RVec>*forces)
{
    bool converged = false;
    int  error     = 0;
    // Precompute some stuff
    std::vector<double>    constraint_distance_squared;
    std::vector<double>    half_reduced_mass;
    std::vector<double>    distance_squared_tolerance;
    std::vector<gmx::RVec> initial_displacements;
    std::vector<double>    scaled_lagrange_multiplier(bonds.size(), 0);
    for (const auto &b : bonds)
    {
        auto &indices    = b->atomIndices();
        double hrm       = 0.5*(atoms[indices[0]].invmass() +
                                atoms[indices[1]].invmass());
        half_reduced_mass.push_back(hrm);
        auto &params     = b->params();
        rvec dx;
        rvec_sub((*positions)[b->atomIndices()[0]],
                 (*positions)[b->atomIndices()[1]], dx);
        initial_displacements.push_back(dx);
        auto bondlength2 = params[bondLENGTH]*params[bondLENGTH];
        constraint_distance_squared.push_back(bondlength2);
        distance_squared_tolerance.push_back(toler_*toler_/bondlength2);
    }
    // And some more...
    std::vector<double> invmass;
    for (const auto &a : atoms)
    {
        invmass.push_back(a.invmass());
    }
    for(int iter = 0; iter < maxiter_ && !converged; iter++)
    {
        // Reset at the start of each iteration
        converged = false;
        for (size_t ll = 0; ll < bonds.size() && error == 0; ll++)
        {
            // Get the parameters. We have to know their names to do this.
            // Get the atom indices
            auto &indices   = bonds[ll]->atomIndices();
            auto ai = indices[0];
            auto aj = indices[1];
            // Compute distance before updating
            rvec dx;
            rvec_sub((*positions)[ai], (*positions)[aj], dx);
            const real r_prime_squared                = norm2(dx);
            const real constraint_distance_squared_ll = constraint_distance_squared[ll];
            const real diff = constraint_distance_squared_ll - r_prime_squared;
            /* iconvf is less than 1 when the error is smaller than a bound */
            const real iconvf = std::abs(diff) * distance_squared_tolerance[ll];

            if (iconvf > 1.0_real)
            {
                converged = static_cast<bool>(iconvf);
                const real r_dot_r_prime = iprod(dx, initial_displacements[ll]);
                
                if (r_dot_r_prime < constraint_distance_squared_ll * toler_)
                {
                    error = ll + 1;
                }
                else
                {
                    /* The next line solves equation 5.6 (neglecting
                       the term in g^2), for g */
                    real scaled_lagrange_multiplier_ll =
                        omega_ * diff * half_reduced_mass[ll] / r_dot_r_prime;
                    scaled_lagrange_multiplier[ll] += scaled_lagrange_multiplier_ll;
                    const real xh = dx[XX] * scaled_lagrange_multiplier_ll;
                    const real yh = dx[YY] * scaled_lagrange_multiplier_ll;
                    const real zh = dx[ZZ] * scaled_lagrange_multiplier_ll;
                    const real im = invmass[ai];
                    const real jm = invmass[aj];
                    (*positions)[ai][XX] += xh * im;
                    (*positions)[ai][YY] += yh * im;
                    (*positions)[ai][ZZ] += zh * im;
                    (*positions)[aj][XX] -= xh * jm;
                    (*positions)[aj][YY] -= yh * jm;
                    (*positions)[aj][ZZ] -= zh * jm;
                }
            }
        }
    }
    return error; 
}

}

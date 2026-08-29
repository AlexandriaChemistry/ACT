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

void Constrainer::prepareConstants(MsgHandler                 *msghandler,
                                   Potential                   pot,
                                   const TopologyEntryVector  &bonds,
                                   const std::vector<ActAtom> &atoms)
{
    // Empty arrays
    half_reduced_mass_.clear();
    constraint_distance_squared_.clear();
    distance_squared_tolerance_.clear();
    scaled_lagrange_multiplier_.clear();
    invmass_.clear();
    m2_.clear();
    // Precompute some stuff
    for (const auto &a : atoms)
    {
        invmass_.push_back(a.invmass());
    }
    std::map<Potential, int> p2i = {
        { Potential::HUA_BONDS, huaLENGTH },
        { Potential::MORSE_BONDS, morseLENGTH },
        { Potential::HARMONIC_BONDS, bondLENGTH },
        { Potential::CUBIC_BONDS, cubicLENGTH }
    };
    auto pi = p2i.find(pot);
    if (p2i.end() == pi)
    {
        std::string msg = gmx::formatString("Unsupported bond potential %s in Constrainer",
                                            potentialToString(pot).c_str());
        if (msghandler)
        {
            msghandler->fatal(msg);
        }
        else
        {
            GMX_THROW(gmx::InvalidInputError(msg));
        }
    }
    int blIndex = pi->second;

    // And some more...
    for (const auto &b : bonds)
    {
        auto &indices    = b->atomIndices();
        double hrm       = 0.5*(atoms[indices[0]].invmass() +
                                atoms[indices[1]].invmass());
        half_reduced_mass_.push_back(hrm);
        m2_.push_back(2.0/hrm);
        auto &params     = b->params();
        auto bondlength2 = params[blIndex]*params[blIndex];
        constraint_distance_squared_.push_back(bondlength2);
        distance_squared_tolerance_.push_back(0.5/(toler_*bondlength2));
        scaled_lagrange_multiplier_.push_back(0.0);
    }
}

void Constrainer::computeInitialDisplacements(gmx_unused MsgHandler        *msghandler,
                                              const TopologyEntryVector    &bonds,
                                              const std::vector<gmx::RVec> &positions)
{
    // Empty arrays
    initial_displacements_.clear();

    // Precompute some stuff
    for (const auto &b : bonds)
    {
        auto &indices    = b->atomIndices();
        rvec dx;
        rvec_sub(positions[indices[0]], positions[indices[1]], dx);
        initial_displacements_.push_back(dx);
    }
}                      

int Constrainer::shake(gmx_unused MsgHandler      *msghandler,
                       Potential                   pot,
                       const TopologyEntryVector  &bonds,
                       const std::vector<ActAtom> &atoms,
                       std::vector<gmx::RVec>     *positions)
{
    prepareConstants(msghandler, pot, bonds, atoms);
    computeInitialDisplacements(msghandler, bonds, *positions);
    bool converged = false;
    int  error     = 0;
    for(int iter = 0; iter < maxiter_ && !converged; iter++)
    {
        // Reset at the start of each iteration
        converged = true;
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
            const real r_prime_squared = norm2(dx);
            const real diff            = constraint_distance_squared_[ll] - r_prime_squared;
            /* iconvf is less than 1 when the error is smaller than a bound */
            const real iconvf = std::abs(diff) * distance_squared_tolerance_[ll];
            converged = converged && (iconvf <= 1.0);
            if (iconvf > 1.0_real)
            {
                const real r_dot_r_prime = iprod(dx, initial_displacements_[ll]);
                
                if (r_dot_r_prime < constraint_distance_squared_[ll] * toler_)
                {
                    error = ll + 1;
                }
                else
                {
                    /* The next line solves equation 5.6 (neglecting
                       the term in g^2), for g */
                    real scaled_lagrange_multiplier_ll =
                        omega_ * diff * half_reduced_mass_[ll] / r_dot_r_prime;
                    scaled_lagrange_multiplier_[ll] += scaled_lagrange_multiplier_ll;
                    rvec xyzh;
                    svmul(scaled_lagrange_multiplier_ll, initial_displacements_[ll], xyzh);
                    for(int m = 0; m < DIM; m++)
                    {
                        (*positions)[ai][m] += xyzh[m] * invmass_[ai];
                        (*positions)[aj][m] -= xyzh[m] * invmass_[aj];
                    }
                }
            }
        }
    }
    if (error == 0 && !converged)
    {
        error = -1;
    }
    return error;
}

int Constrainer::rattle(MsgHandler                   *msghandler,
                        Potential                     pot,
                        const TopologyEntryVector    &bonds,
                        const std::vector<ActAtom>   &atoms,
                        const std::vector<gmx::RVec> &positions,
                        std::vector<gmx::RVec>       *velocities)
{
    prepareConstants(msghandler, pot, bonds, atoms);
    computeInitialDisplacements(msghandler, bonds, positions);

    int  iter      = 0;
    int  error     = 0;
    bool converged = false;
    for (iter = 0; (iter < maxiter_) && !converged && (error == 0); iter++)
    {
        converged = true;
        for (size_t ll = 0; (ll < bonds.size()) && (error == 0); ll++)
        {
            auto &indices   = bonds[ll]->atomIndices();
            auto ai = indices[0];
            auto aj = indices[1];
            rvec dv;
            rvec_sub((*velocities)[ai], (*velocities)[aj], dv);

            const real vpijd  = iprod(dv, initial_displacements_[ll]);

            /* iconv is zero when the error is smaller than a bound */
            const real iconvf = std::abs(vpijd) * (distance_squared_tolerance_[ll] / invdt_);
            converged         = converged && (iconvf <= 1.0);

            if (iconvf > 1.0_real)
            {
                const real fac  = omega_ * 2.0_real * half_reduced_mass_[ll] / constraint_distance_squared_[ll];
                const real acor = -fac * vpijd;
                scaled_lagrange_multiplier_[ll] += acor;
                rvec xyzh;
                svmul(acor, initial_displacements_[ll], xyzh);

                for(int m = 0; m < DIM; m++)
                {
                    (*velocities)[ai][m] += xyzh[m] * invmass_[ai];
                    (*velocities)[aj][m] -= xyzh[m] * invmass_[aj];
                }
            }
        }
    }
    if (msghandler && msghandler->debug())
    {
        msghandler->writeDebug(gmx::formatString("Rattled for %d iterations, error code %d",
                                                 iter, error));
    }
    if (error == 0 && !converged)
    {
        error = -1;
    }
    return error;
}

}

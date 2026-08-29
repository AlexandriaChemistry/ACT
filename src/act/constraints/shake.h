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
 * \author Mohammad Mehdi Ghahremanpour <mohammad.ghahremanpour@icm.uu.se>
 */

#ifndef ACT_SHAKE_H
#define ACT_SHAKE_H

#include <vector>

#include "act/molprop/topologyentry.h"
#include "gromacs/math/vec.h"

namespace alexandria
{

class ActAtom;
class ACTMol;
class MsgHandler;

class Constrainer
{
private:
    //! Max number of iterations to try constraining
    int    maxiter_ = 100;
    //! Relative tolerance for bond lengths
    double toler_   = 1e-4;
    //! Overrelaxation
    double omega_   = 1;

public:
    Constrainer(int maxiter, double toler) : maxiter_(maxiter), toler_(toler)
    {}

    /*! \brief Execute the constraint algorithm
     *
     * Original implementation from R.C. van Schaik and W.F. van Gunsteren
     * (ETH Zuerich, June 1992), adapted for GROMACS by David van der
     * Spoel November 1992. Adapted for ACT by DvdS, August 2026.
     *
     * The algorithm here is based section five of Ryckaert, Ciccotti and
     * Berendsen, J Comp Phys, 23, 327, 1977.
     *
     * \param[in]    msghandler  For warnings and errors
     * \param[in]    bonds       The atom identifiers and parameters
     * \param[in]    atoms       The atoms
     * \param[inout] coordinates The coordinates of all particles
     * \param[inout] forces      The forces on all particles
     * \return error code corresponding the (n+1)th constraint or zero if fine
     */
    int shake(MsgHandler                   *msghandler,
              const TopologyEntryVector    &bonds,
              const std::vector<ActAtom>   &atoms,
              std::vector<gmx::RVec>       *coordinates,
              std::vector<gmx::RVec>       *forces);
};

} // namespace alexandria

#endif // ACT_SHAKE_H

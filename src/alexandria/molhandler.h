/*
 * This source file is part of the Alexandria Chemistry Toolkit.
 *
 * Copyright (C) 2014-2022
 *
 * Developers:
 *             Mohammad Mehdi Ghahremanpour,
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
 * \author David van der Spoel <david.vanderspoel@icm.uu.se>
 * \author Julian Marrades <julian.marrades@hotmail.es>
 */

#ifndef ACT_MOLHANDLER_H
#define ACT_MOLHANDLER_H

#include <vector>

#include "act/utility/regression.h"
#include "mymol.h"
#include "gromacs/mdtypes/commrec.h"

namespace alexandria
{

/*! \brief Handles molecules by performing algorithms on them
 * For example, energy minimization and hessian computation.
 * This class can access private members of MyMol as it is a friend class
 */
class MolHandler
{

public:

    /*! \brief Compute the second derivative matrix of the potential energy
     *
     * \param[in]  mol       Molecule to get the hessian for
     * \param[in]  crtmp     Temporary communication record for one core.
     *                       FIXME: another method without this
     * \param[in]  atomIndex Vector containing the indices of the real 
     *                       atoms, not shells or vsites.
     *                       FIXME: Create another method without this argument, then
     *                       get it and call this one
     * \param[out] hessian   MatrixWrapper object that must be pre-
     *                       allocated to NxN where N = 3*atomIndex.size()
     * \param[out] forceZero The forces on the atoms in the input structure,
     *                       that is, not on the shells or vsites. Will be cleared
     *                       and overwritten
     * \return the potential energy of the input structure
     */
    double computeHessian(      MyMol               *mol,
                          const t_commrec           *crtmp,
                          const std::vector<int>    &atomIndex,
                                MatrixWrapper       *hessian,
                                std::vector<double> *forceZero);

};

} // namespace alexandria


#endif // ACT_MOLHANDLER_H
/*
 * This file is part of the GROMACS molecular simulation package.
 *
 * Copyright (c) 1991-2000, University of Groningen, The Netherlands.
 * Copyright (c) 2001-2004, The GROMACS development team.
 * Copyright (c) 2010,2014,2015,2016,2018 by the GROMACS development team.
 * Copyright (c) 2019,2020, by the GROMACS development team, led by
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
#ifndef GMX_GMXANA_PRINC_H
#define GMX_GMXANA_PRINC_H

#include <vector>

#include "act/alexandria/topology.h"
#include "gromacs/math/vectypes.h"
#include "gromacs/utility/basedefinitions.h"

namespace alexandria
{

/*! \brief Rotate all atoms in index using matrix trans 
 * \param[in]    index  List of atom numbers
 * \param[inout] x      The coordinates
 * \param[in]    matrix The 3x3 rotation matrix
 */
void rotate_atoms(const std::vector<int> &index,
                  std::vector<gmx::RVec> *x,
                  const matrix            trans);

/*! \brief Compute principal components
 * Atoms are mass weighted and it is assumed that the center of 
 * mass is in the origin!
 * \param[in]  n       Size of the index array
 * \param[in]  index   List of atom numbers
 * \param[in]  mass    Atomic masses
 * \param[in]  x       The coordinates
 * \param[out] trans   The matrix needed to rotate the molecule to the reference frame
 * \param[out] inertia The moments of inertia
 */
void principal_comp(const std::vector<int>       &index,
                    const std::vector<real>      &mass,
                    const std::vector<gmx::RVec> &x, 
                    matrix                       *trans,
                    gmx::RVec                    *inertia);
                    
/*! \brief Calculate the center of mass of the atoms in index. if bQ then the atoms
 * will be charge weighted rather than mass weighted.
 * \return the total mass/charge.
 */
real calc_xcm(const std::vector<gmx::RVec> &x,
              const std::vector<int>       &index,
              const std::vector<ActAtom>   &atom,
              gmx::RVec                    *xcm,
              bool                          bQ);

/*! \brief Calculate the center of mass and subtract it from all coordinates.
 * Returns the original center of mass in xcm
 * \return the total mass
 */
real sub_xcm(std::vector<gmx::RVec>       *x,
             const std::vector<int>       &index,
             const std::vector<ActAtom>   &atom,
             gmx::RVec                    *xcm,
             bool                          bQ);

void add_xcm(std::vector<gmx::RVec>       *x,
             const std::vector<int>       &index,
             gmx::RVec                    &xcm);
/* Increment all atoms in index with xcm */

} // namespace

#endif

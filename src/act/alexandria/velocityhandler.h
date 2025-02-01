/*
 * This file is part of the GROMACS molecular simulation package.
 *
 * Copyright (c) 1991-2000, University of Groningen, The Netherlands.
 * Copyright (c) 2001-2004, The GROMACS development team.
 * Copyright (c) 2013,2014,2015,2019,2020, by the GROMACS development team, led by
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

#ifndef ACT_VELOCITY_HANDLER
#define ACT_VELOCITY_HANDLER

#include <vector>

#include "topology.h"

namespace gmx
{
class TextWriter;
}

namespace alexandria
{

/*! \brief
 * Generate Maxwellian velocities.
 *
 * \param[in]  tempi  Temperature to generate around. If <= 0 nothing will be done.
 * \param[in]  seed   Random number generator seed
 * \param[in]  atoms  The atoms
 * \param[out] v      Velocities
 * \param[in]  tw     TextWriter
 */
void maxwell_speed(real                        tempi, 
                   unsigned int                seed, 
                   const std::vector<ActAtom> &atoms,
                   std::vector<gmx::RVec>     *v,
                   gmx::TextWriter            *tw);

/*! \brief
 * Remove the center of mass motion in a set of coordinates.
 *
 * \param[in]  atoms  The atoms
 * \param[in]  x      Coordinates
 * \param[out] v      Velocities
 */
void stop_cm(const std::vector<ActAtom>   &atoms,
             const std::vector<gmx::RVec> &x,
             std::vector<gmx::RVec>       *v);

} // namespace alexandria

#endif

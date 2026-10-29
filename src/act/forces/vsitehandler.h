/*
 * This source file is part of the Alexandria Chemistry Toolkit.
 *
 * Copyright (C) 2023,2024
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
#ifndef ACT_VSITEHANDLER_H
#define ACT_VSITEHANDLER_H
    
#include <vector>

#include "gromacs/math/vectypes.h"

namespace alexandria
{

class Topology;

/*! \brief Class to construct vsite positions and distribute vsite forces
 *
 * Periodic boundary conditions are not supported.
 */
class VsiteHandler
{
 public:
    /*! Construct vsite positions
     * \param[in]    top         The molecular topology
     * \param[inout] coordinates The atomic coordinates, the virtual site
     *                           positions are set.
     */
    void constructPositions(const Topology         *top,
                            std::vector<gmx::RVec> *coordinates) const;
                 
    /*! \brief Distribute the forces back to atoms
     * \param[in]    top    The molecular topology
     * \param[in]    coords The atomic coordinates.
     * \param[inout] forces The forces on the atoms
     */
    void distributeForces(const Topology               *top,
                          const std::vector<gmx::RVec> &coords,
                          std::vector<gmx::RVec>       *forces) const;
   
};

} // namespace alexandria

#endif

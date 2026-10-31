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
    
#include <cstddef>
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

    /*! \brief Construct the virtual sites of many conformations at once
     * The coordinates are in the layout of ConformationArrays, component
     * m of particle a in slot c is at (3*a+m)*stride + c. The positions
     * are the same as those of constructPositions.
     * \param[in]    top    The molecular topology
     * \param[inout] x      The coordinates, the virtual site positions are set
     * \param[in]    stride The stride
     * \param[in]    nslot  The number of slots to construct, at most stride
     * \param[in]    nconf  The number of slots that hold conformations, the
     *                      others are padding, of which the floating point
     *                      operations are not counted
     */
    void constructPositions(const Topology *top,
                            double         *x,
                            size_t          stride,
                            size_t          nslot,
                            size_t          nconf) const;

    /*! \brief Distribute the forces of many conformations back to atoms
     * The layout is that of constructPositions with many conformations,
     * and the forces are the same as those of distributeForces.
     * \param[in]    top    The molecular topology
     * \param[in]    x      The coordinates
     * \param[inout] f      The forces, those on the virtual sites are set to zero
     * \param[in]    stride The stride
     * \param[in]    nslot  The number of slots to distribute, at most stride
     * \param[in]    nconf  The number of slots that hold conformations, see
     *                      constructPositions
     */
    void distributeForces(const Topology *top,
                          const double   *x,
                          double         *f,
                          size_t          stride,
                          size_t          nslot,
                          size_t          nconf) const;

};

} // namespace alexandria

#endif

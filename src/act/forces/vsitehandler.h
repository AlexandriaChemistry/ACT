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
    
#include "act/alexandria/topology.h"
#include "gromacs/math/vectypes.h"
#include "gromacs/pbcutil/pbc.h"

namespace alexandria
{

/*! \brief Class to construct vsite positions and distribute vsite forces
 */
class VsiteHandler
{
private:
    //! GROMACS structure for periodic boundary conditions
    t_pbc pbc_     = { 0 };
    //! Integration time step
    real  dt_      = 0.001; // ps
    //! Whether pbc has been initiated
    bool  initPBC_ = false;
 public:
    //! \brief Constructor
    VsiteHandler();
    
    /*! Initialize vars
     * If the box has zero edges (determinant is zero) it will be assumed that no periodic
     * boundary conditions are to be taken into account. If not, the program will crash, since
     * PBC is not supported (yet). https://github.com/dspoel/ACT/issues/127
     * \param[in] box The box edges.
     * \param[in] dt  The integration time step 
     */
    void init(matrix &box, real dt);

    /*! Construct vsite positions
     * \param[in]  top         The molecular topology
     * \param[in]  coordinates The atomic coordinates. Coordinates of
     *                         shell particles may be changed.
     * \param[out] box         The periodic box
     */
    void constructPositions(const Topology         *top,
                            std::vector<gmx::RVec> *coordinates,
                            const matrix           &box) const;
                 
    /*! \brief Distribute the forces back to atoms
     * \param[in]  top    The molecular topology
     * \param[in]  coords The atomic coordinates.
     * \param[out] forces The forces on the atoms
     * \param[out] box    The periodic box
     */
    void distributeForces(const Topology               *top,
                          const std::vector<gmx::RVec> &coords,
                          std::vector<gmx::RVec>       *forces,
                          const matrix                 &box) const;
   
};

} // namespace alexandria

#endif

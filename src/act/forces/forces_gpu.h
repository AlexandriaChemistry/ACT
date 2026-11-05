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
 * Potentials computed for many conformations at once on a GPU.
 *
 * The arrays have the layout of computeExplicitSimd. The parameters of
 * the interactions are prepared here, the GPU code in src/act/gpu only
 * gets plain arrays, see gpu.h. Whether a GPU is used is decided by
 * gpu::useGpu().
 * \author David van der Spoel <david.vanderspoel@icm.uu.se>
 */
#ifndef ACT_FORCES_GPU_H
#define ACT_FORCES_GPU_H

#include <cstddef>
#include <vector>

#include "act/basics/coordinate.h"
#include "act/forcefield/potential.h"
#include "act/molprop/topologyentry.h"
#include "act/topology/actatom.h"

namespace alexandria
{

/*! \brief The least work for the GPU
 * Launching kernels on a GPU and copying the data take tens of
 * microseconds, in which the processor computes some ten thousand
 * interactions. Therefore computeGpu only uses the GPU for at least this
 * number of interactions times conformations. The environment variable
 * ACT_GPU_MIN_WORK sets it.
 * \return the least number of interactions times conformations
 */
size_t gpuMinimumWork();

/*! \brief Set the least work for the GPU, see gpuMinimumWork
 * \param[in] work The least number of interactions times conformations
 * \return the value before
 */
size_t setGpuMinimumWork(size_t work);

/*! \brief Compute a potential for all conformations on the GPU
 * As computeExplicitSimd, for the potentials that have a GPU kernel:
 * COULOMB_POINT and COULOMB_GAUSSIAN. The energy of each interaction is
 * computed in single precision and summed in double precision.
 * \param[in]    mypot    The potential
 * \param[in]    entries  The interactions
 * \param[in]    atoms    The atoms
 * \param[in]    epsilonr Relative dielectric constant
 * \param[in]    stride   Distance between the rows of the arrays
 * \param[in]    nslot    Number of slots to compute, at most stride
 * \param[in]    x        The coordinates
 * \param[inout] f        The forces, added to
 * \param[inout] energies The energy terms, added to
 * \param[in]    forces   Whether to compute forces
 * \return false if this potential is not available on the GPU, no GPU is
 *         used or there is less work than gpuMinimumWork
 */
bool computeGpu(Potential                   mypot,
                TopologyEntryRange          entries,
                const std::vector<ActAtom> &atoms,
                double                      epsilonr,
                size_t                      stride,
                size_t                      nslot,
                const Coordinate           *x,
                float                      *f,
                double                     *energies,
                bool                        forces);

/*! \brief Without arrays in single precision there are no GPU kernels
 * \return false
 */
inline bool computeGpu(Potential, TopologyEntryRange, const std::vector<ActAtom> &, double,
                       size_t, size_t, const Coordinate *, double *, double *, bool)
{
    return false;
}

} // namespace alexandria

#endif

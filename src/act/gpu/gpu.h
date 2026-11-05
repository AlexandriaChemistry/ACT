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
 * Interface to the GPU code, without SYCL types, such that code compiled
 * with any C++ compiler can call it, see cmake/actManageGpu.cmake.
 * \author David van der Spoel <david.vanderspoel@icm.uu.se>
 */
#ifndef ACT_GPU_GPU_H
#define ACT_GPU_GPU_H

#include <cstddef>
#include <cstdint>
#include <string>
#include <vector>

namespace alexandria
{

namespace gpu
{

/*! \brief The programming model of the GPU code
 * \return for instance "SYCL (AdaptiveCpp 25.10.0)", or an empty string if
 *         ACT was built without GPU support
 */
std::string implementation();

/*! \brief The GPUs that the programming model finds
 * The first call initializes the runtime of the programming model.
 * AdaptiveCpp warns about each backend that it cannot load, for instance
 * that of AMD GPUs on a computer with NVIDIA GPUs only, unless the
 * environment variable ACPP_DEBUG_LEVEL is 1 or 0, as ACTRC sets it.
 * \return the names of the GPUs, empty without GPU support or without GPUs
 */
std::vector<std::string> devices();

/*! \brief Whether the GPU computes the potentials that it has kernels for
 * True if ACT was built with GPU support, a GPU was found and neither the
 * environment variable ACT_NO_GPU nor setUseGpu turned it off.
 * \return whether the GPU is used
 */
bool useGpu();

/*! \brief Turn the GPU on or off, for tests and comparisons
 * The GPU is only turned on if useGpu would return true without ACT_NO_GPU.
 * \param[in] use Whether to use the GPU
 * \return whether the GPU was used before
 */
bool setUseGpu(bool use);

/*! \brief A Coulomb interaction between two particles
 * Each particle may have a second charge at the same position with its own
 * width, such that there are up to four terms, see ActAtom::chargeVS. Terms
 * after the first one with a zero charge product are skipped.
 */
struct CoulombPair
{
    //! The first atom
    uint32_t i;
    //! The second atom
    uint32_t j;
    //! The charge products of the terms, including 1/(4 pi eps0 epsilonr)
    float    qq[4];
    //! The widths of the combined Gaussian charges, zero for point charges
    float    zeta[4];
};

/*! \brief Compute Coulomb interactions for all conformations on the GPU
 * The layout is that of computeExplicitSimd: coordinate m of atom a in
 * slot c is x[(3*a+m)*stride + c], in units of coordinateUnit, see
 * coordinate.h, the forces likewise. The forces and the energies of
 * slots 0 to nslot - 1 are added to. The energy of each interaction is
 * computed in single precision and summed in double precision.
 * \param[in]    natom    Number of atoms
 * \param[in]    stride   Distance between the rows of the arrays
 * \param[in]    nslot    Number of slots to compute, at most stride
 * \param[in]    x        The coordinates
 * \param[in]    pairs    The interactions
 * \param[in]    npair    Number of interactions
 * \param[inout] f        The forces, not used if forces is false
 * \param[inout] energies The energy of each slot
 * \param[in]    forces   Whether to compute forces
 * \return false if the GPU is not available
 */
bool computeCoulomb(size_t             natom,
                    size_t             stride,
                    size_t             nslot,
                    const int64_t     *x,
                    const CoulombPair *pairs,
                    size_t             npair,
                    float             *f,
                    double            *energies,
                    bool               forces);

} // namespace gpu

} // namespace alexandria

#endif

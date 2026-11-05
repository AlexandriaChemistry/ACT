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
#include "act/gpu/gpu.h"
#include "act/forcefield/potential.h"
#include "act/molprop/topologyentry.h"
#include "act/topology/actatom.h"

namespace alexandria
{

/*! \brief The least work for the GPU
 * Launching kernels on a GPU and copying the data take tens of
 * microseconds, in which the processor computes some ten thousand
 * interactions. Therefore computeGpu and GpuBatch only use the GPU for at
 * least this number of interactions times conformations, of all molecules
 * of a batch. The environment variable ACT_GPU_MIN_WORK sets it.
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

/*! \brief The GPU kernel of a potential
 * \param[in]  mypot  The potential
 * \param[out] kernel The kernel
 * \return false if the GPU has no kernel for it
 */
bool gpuPairKernel(Potential mypot, gpu::PairKernel *kernel);

/*! \brief The parameters of a pair interaction for the GPU
 * As the prepare functions of the explicit SIMD kernels in
 * forces_explicit_simd.cpp, rounded to single precision.
 * \param[in]  kernel   The kernel, see gpuPairKernel
 * \param[in]  b        The interaction
 * \param[in]  atoms    The atoms
 * \param[in]  epsilonr Relative dielectric constant
 * \param[out] pe       The interaction
 * \return false if there is no interaction, as without polarizability
 */
bool gpuPairEntry(gpu::PairKernel             kernel,
                  const TopologyEntry        &b,
                  const std::vector<ActAtom> &atoms,
                  double                      epsilonr,
                  gpu::PairEntry             *pe);

/*! \brief Potentials of several molecules that are computed on the GPU at once
 * computeEnergies first adds the interactions of all molecules that the GPU
 * can compute, then computes them together, then uses the results where
 * each molecule would compute them, see lowComputeEnergiesT. One GPU
 * calculation for all molecules takes about as long as one for one
 * molecule, since the time to start it dominates for small molecules.
 */
class GpuBatch
{
public:
    /*! \brief Whether the GPU can compute a potential
     * \param[in] mypot The potential
     * \return true if it can, for COULOMB_POINT and COULOMB_GAUSSIAN with a GPU in use
     */
    static bool canCompute(Potential mypot);

    //! Remove all tasks, but keep the memory
    void clear();

    /*! \brief Add the interactions of the conformations of a molecule
     * The arrays have the layout of computeExplicitSimd and have to stay
     * where they are until the results are used.
     * \param[in] mypot    The potential, see canCompute
     * \param[in] entries  The interactions
     * \param[in] atoms    The atoms
     * \param[in] epsilonr Relative dielectric constant
     * \param[in] stride   Distance between the rows of the arrays
     * \param[in] nslot    Number of slots to compute, at most stride
     * \param[in] x        The coordinates
     */
    void add(Potential                   mypot,
             TopologyEntryRange          entries,
             const std::vector<ActAtom> &atoms,
             double                      epsilonr,
             size_t                      stride,
             size_t                      nslot,
             const Coordinate           *x);

    /*! \brief Compute all tasks on the GPU, if there is enough work, see gpuMinimumWork
     * \param[in] forces Whether to compute forces
     * \return whether the results were computed
     */
    bool compute(bool forces);

    //! \return whether the results were computed
    bool computed() const { return computed_; }

    //! \return the number of tasks
    size_t size() const { return ntask_; }

    /*! \brief Add the results of a task
     * \param[in]    k        The task, in the order they were added
     * \param[inout] f        The forces in the layout of the coordinates,
     *                         added to if they were computed
     * \param[inout] energies The energy of each slot, added to
     */
    void addResults(size_t k, float *f, double *energies) const;

private:
    //! The interactions of one molecule
    struct Task
    {
        size_t                         natom  = 0;
        size_t                         stride = 0;
        size_t                         nslot  = 0;
        const Coordinate              *x      = nullptr;
        std::vector<gpu::CoulombPair>  pairs;
        std::vector<float>             f;
        std::vector<double>            energies;
    };
    //! The tasks, of which the first ntask_ are used
    std::vector<Task> tasks_;
    //! Number of tasks in use
    size_t            ntask_    = 0;
    //! Whether the results were computed
    bool              computed_ = false;
    //! Whether forces were computed
    bool              forces_   = false;
};

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

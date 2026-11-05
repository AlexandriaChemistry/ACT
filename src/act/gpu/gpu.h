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

/*! \brief The GPU that this process uses
 * The processes on a node are spread over its GPUs: the number of the
 * process on the node, from the environment variables that Slurm and the
 * MPI libraries set, modulo the number of GPUs.
 * \return the index in devices, -1 without GPU
 */
int gpuIndex();

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

/*! \brief The Coulomb interactions of the conformations of a molecule, see computeCoulombBlocks
 * The layout is that of computeExplicitSimd: coordinate m of atom a in
 * slot c is x[(3*a+m)*stride + c], in units of coordinateUnit, see
 * coordinate.h, the forces likewise.
 */
struct CoulombBlock
{
    //! Number of atoms
    size_t             natom;
    //! Distance between the rows of the arrays
    size_t             stride;
    //! Number of slots to compute, at most stride
    size_t             nslot;
    //! The coordinates
    const int64_t     *x;
    //! The interactions
    const CoulombPair *pairs;
    //! Number of interactions
    size_t             npair;
    //! The forces, added to, not used without forces
    float             *f;
    //! The energy of each slot, added to
    double            *energies;
};

/*! \brief Compute the Coulomb interactions of several molecules at once on the GPU
 * All blocks are computed by the same kernels, such that many small
 * molecules make enough work for a GPU. Blocks with the same coordinate
 * array share its copy on the GPU. The energy of each interaction is
 * computed in single precision and summed in double precision. The
 * results do not depend on the other blocks or slots.
 * \param[in] blocks The molecules
 * \param[in] nblock Number of blocks
 * \param[in] forces Whether to compute forces
 * \return false if the GPU is not available
 */
bool computeCoulombBlocks(const CoulombBlock *blocks,
                          size_t              nblock,
                          bool                forces);

/*! \brief Compute Coulomb interactions for all conformations on the GPU
 * As computeCoulombBlocks for one block.
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

//! The pair potentials of the GPU, see pairInteraction in sycl/pairkernels.h
enum class PairKernel : uint32_t
{
    //! Point and Gaussian charges, with second charges: charge products p[0..3] and widths p[4..7]
    Coulomb,
    //! Shell and core: force constant p[0], distance p[1] and force constant p[2] of the quartic term
    Polarization,
    //! Lennard-Jones 12-6: 4 epsilon sigma^6 p[0] and 4 epsilon sigma^12 p[1]
    LJ12_6,
    //! Buffered 14-7: 1/sigma p[0], epsilon p[1], gamma p[2], delta p[3]
    LJ14_7,
    //! Wang-Buckingham: sigma^6 p[0], 3/(gamma+3) p[1], dispersion prefactor p[2], gamma p[3], 1/sigma p[4]
    WangBuckingham,
    //! MacDaniel-Schmidt, a double exponential: prefactor p[0], exponent p[1]
    DoubleExponential
};

//! Number of energy terms of the pair potentials, at most
constexpr int pairEnergyTerms = 3;

//! A pair interaction with prepared parameters, see PairKernel
struct PairEntry
{
    //! The first particle
    uint32_t i;
    //! The second particle
    uint32_t j;
    //! The parameters
    float    p[8];
};

//! The interactions of one potential of a molecule, see RelaxMolecule
struct PairBlock
{
    //! The potential
    PairKernel       kernel;
    //! The interactions
    const PairEntry *pairs;
    //! Number of interactions
    size_t           npair;
    /*! The energy terms of the last calculation of each slot, term k of
     * slot c at k*stride + c, pairEnergyTerms rows, written for the slots
     * that were computed
     */
    double          *energies;
};

//! A shell that is moved, see ForceComputer::updateShells
struct MovingShell
{
    //! The shell
    uint32_t shell;
    //! Its core
    uint32_t core;
    //! The parent whose distance is limited, if wall
    uint32_t parent;
    //! Whether the distance to the parent is limited
    uint32_t wall;
    //! The harmonic force constant
    double   fcShell;
    //! The force constant of the quartic term
    double   fcHyper;
    //! The distance beyond which the quartic term acts
    double   rHyper;
};

//! A shell whose force is part of the convergence criterion, see shellForces2 in forcecomputer.cpp
struct ForceShell
{
    //! The shell
    uint32_t shell;
    //! The parent whose distance is limited, if wall
    uint32_t parent;
    //! Whether the distance to the parent is limited
    uint32_t wall;
};

//! The shells of a monomer, which converge together
struct ShellFragment
{
    //! The first particle of the monomer
    uint32_t firstAtom;
    //! Number of particles of the monomer
    uint32_t natom;
    //! The first of the shells in RelaxMolecule::forceShells
    uint32_t firstForceShell;
    //! Number of these
    uint32_t nforceShell;
    //! The first of the shells in RelaxMolecule::movingShells
    uint32_t firstMovingShell;
    //! Number of these
    uint32_t nmovingShell;
    //! The number of shells that divides the sum of their squared forces
    double   nShell;
};

/*! \brief The conformations of a molecule whose shells are relaxed, see relaxShells
 * The layout is that of computeExplicitSimd, see CoulombBlock.
 */
struct RelaxMolecule
{
    //! Number of particles
    size_t               natom;
    //! Distance between the rows of the arrays
    size_t               stride;
    //! Number of slots
    size_t               nslot;
    //! The coordinates, the shells are moved, which are the only ones copied back
    int64_t             *x;
    //! The forces of the last calculation, those of the slots computed are replaced
    float               *f;
    //! The forces of the interactions without shells, which do not change
    const float         *fixedForces;
    //! Whether each slot converged, set when it does
    char                *done;
    //! Whether each slot was computed again
    char                *computed;
    //! The shells in the convergence criterion of all fragments
    const ForceShell    *forceShells;
    //! The shells that move of all fragments
    const MovingShell   *movingShells;
    //! The monomers whose shells are relaxed
    const ShellFragment *fragments;
    //! Number of these, one or two
    size_t               nfragment;
    //! The interactions with shells
    const PairBlock     *blocks;
    //! Number of these
    size_t               nblock;
    /*! Whether the forces on all particles are needed, otherwise only
     * those on the shells are copied to f
     */
    bool                 allForces;
};

//! Settings of relaxShells, see ForceComputer
struct RelaxSettings
{
    //! The mean square force on the shells below which they converged
    double msForceToler;
    //! Square of the largest distance of a shell from its parent, 0 for none
    double maxShellDistance2;
    //! The smallest fraction of the Newton step
    double minStep;
    //! The largest fraction of the Newton step
    double maxStep;
    //! Iterations, including the first calculation, as ForceComputer::maxiter
    int    maxIter;
};

/*! \brief Relax the shells of the conformations of several molecules on the GPU
 * As ForceComputer::relaxMonomer and relaxDimer: each iteration moves the
 * shells of the slots that did not converge, see ForceComputer::updateShells,
 * then computes the interactions with shells for these, which with the
 * forces of the interactions without shells give the forces. All data
 * stay on the GPU until the shells of all slots converged or maxIter is
 * reached. The step of each fragment starts at one.
 * \param[inout] mols       The molecules
 * \param[in]    nmol       Number of molecules
 * \param[in]    settings   The settings
 * \param[out]   nComputed  The number of slots computed in all iterations
 * \param[out]   notFinite  Whether a shell got a displacement that is not
 *                          finite, then it was not moved
 * \return false if the GPU is not available
 */
bool relaxShells(const RelaxMolecule *mols,
                 size_t               nmol,
                 const RelaxSettings &settings,
                 size_t              *nComputed,
                 bool                *notFinite);

} // namespace gpu

} // namespace alexandria

#endif

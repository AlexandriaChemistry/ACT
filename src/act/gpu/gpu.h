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
#include <new>
#include <string>
#include <utility>
#include <vector>

namespace alexandria
{

//! See act/forcefield/parameterrecipes.h
struct SetRecipe;

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
 * environment variable ACPP_DEBUG_LEVEL is 1 or 0, which the GPU library
 * of ACT sets when it is loaded, unless it is set already.
 * \return the names of the GPUs, empty without GPU support or without GPUs
 */
std::vector<std::string> devices();

/*! \brief Why no GPU was found, for the log file
 * Looks at the backends of the programming model, the GPU device files of
 * the node and the GPUs that the queue system gave the job.
 * \return lines with the likely causes and what to do, empty with GPUs
 */
std::vector<std::string> noGpuReasons();

//! An estimate of the peak performance of a GPU
struct DevicePeak
{
    //! Single precision GFlop/s, 0 if not known
    double      gflops = 0;
    //! How it was estimated
    std::string basis;
};

/*! \brief Estimate the peak single precision performance of the GPUs
 * From what the programming model reports, see estimatePeak.
 * \return the estimates, in the order of devices
 */
std::vector<DevicePeak> devicePeaks();

/*! \brief Estimate the peak single precision performance of a GPU
 * The number of compute units times the number of single precision lanes
 * of a compute unit, which depends on the vendor and the type of GPU,
 * times two for fused multiply-add instructions, times the clock frequency.
 * NVIDIA streaming multiprocessors have 128 lanes, those of P100, V100,
 * A100 and A30 64. AMD compute units have 64 lanes, without packed
 * instructions.
 * \param[in] vendor The vendor of the GPU
 * \param[in] name   The name of the GPU
 * \param[in] units  The number of compute units
 * \param[in] ghz    The maximum clock frequency in GHz
 * \return the estimate, with gflops 0 for other vendors or without
 *         compute units or frequency
 */
DevicePeak estimatePeak(const std::string &vendor, const std::string &name,
                        unsigned units, double ghz);

//! What the message of OutOfMemory says about the node, see outOfMemoryFacts
struct OutOfMemoryFacts
{
    //! The number of processes on this node, 0 if not known
    int  nodeProcesses = 0;
    //! The number of GPUs on this node, 0 if not known
    int  nodeGpus      = 0;
    //! Whether the GPU is an NVIDIA GPU, for which the multi-process service exists
    bool nvidia        = false;
    //! Whether the NVIDIA multi-process service runs
    bool mpsRuns       = false;
};

/*! \brief The facts about this node for the message of OutOfMemory
 * The number of processes on the node comes from the environment variables
 * of the MPI libraries and Slurm, that of the GPUs from what
 * selectGpuOfProcess found or from the NVIDIA device files, and the
 * multi-process service runs if the process in its file
 * nvidia-cuda-mps-control.pid exists, in CUDA_MPS_PIPE_DIRECTORY or else
 * in /tmp/nvidia-mps.
 * \param[in] device The name of the GPU
 * \return the facts
 */
OutOfMemoryFacts outOfMemoryFacts(const std::string &device);

/*! \brief The message of OutOfMemory
 * \param[in] device         The name of the GPU
 * \param[in] totalBytes     Its memory
 * \param[in] allocatedBytes What this process has allocated on it
 * \param[in] requestedBytes What it asked for
 * \param[in] facts          What is known about the node
 * \return the message, in lines
 */
std::string outOfMemoryMessage(const std::string &device, double totalBytes,
                               size_t allocatedBytes, size_t requestedBytes,
                               const OutOfMemoryFacts &facts);

/*! \brief The GPU has no memory left for what a process puts on it
 * Thrown by the GPU code instead of using an allocation that failed. The
 * message tells how much memory the process has on the GPU, how many
 * processes share it, whether the NVIDIA multi-process service runs, and
 * what to do, see outOfMemoryMessage. As a std::bad_alloc, the program
 * stops with "Memory allocation failed" and this message.
 */
class OutOfMemory : public std::bad_alloc
{
public:
    /*! \brief Constructor
     * \param[in] device         The name of the GPU
     * \param[in] totalBytes     Its memory
     * \param[in] allocatedBytes What this process has allocated on it for the molecules
     * \param[in] requestedBytes What it asked for
     */
    OutOfMemory(const std::string &device, double totalBytes,
                size_t allocatedBytes, size_t requestedBytes);

    //! \return the message
    const char *what() const noexcept override { return message_.c_str(); }

private:
    //! The message, see the constructor
    std::string message_;
};

//! What noGpuReasons looks at on a node without usable GPU, see explainNoGpu
struct NoGpuFacts
{
    //! CUDA_VISIBLE_DEVICES, ROCR_VISIBLE_DEVICES and HIP_VISIBLE_DEVICES, those that are set, with their values
    std::vector<std::pair<std::string, std::string>> visibleDevices;
    //! Whether the node has the device files of NVIDIA GPUs
    bool                                             nvidiaDevices = false;
    //! Whether the node has the device file of AMD GPUs
    bool                                             amdDevices    = false;
    //! The directory of the AdaptiveCpp runtime library, empty if not known
    std::string                                      runtimeDirectory;
    //! The backends of AdaptiveCpp in that directory, such as cuda, hip, ocl and omp
    std::vector<std::string>                         backends;
};

/*! \brief Why no GPU was found, from what the node and the job look like
 * \param[in] facts What noGpuReasons found
 * \return lines with the likely causes and what to do
 */
std::vector<std::string> explainNoGpu(const NoGpuFacts &facts);

/*! \brief Let this process see only the GPU of the node that it uses
 * The processes on a node are spread over its GPUs by their number on the
 * node, from the environment variables of the MPI libraries and Slurm.
 * If CUDA_VISIBLE_DEVICES, ROCR_VISIBLE_DEVICES or HIP_VISIBLE_DEVICES
 * lists several GPUs, as Slurm sets them, only that of this process is
 * kept, and without such a list on a node with several NVIDIA GPUs,
 * CUDA_VISIBLE_DEVICES is set to it. Otherwise, the GPU runtime makes a
 * context on each GPU in each process, and the NVIDIA multi-process
 * service takes only 48 or 60 processes per GPU. Must be called before anything
 * else uses the GPU, at the start of main.
 */
void selectGpuOfProcess();

//! \return the GPU of the node that selectGpuOfProcess chose, empty if it did not choose
std::string gpuOfProcess();

/*! \brief The GPU that this process uses
 * The processes on a node are spread over its GPUs: the number of the
 * process on the node, from the environment variables that the MPI
 * libraries and Slurm set, modulo the number of GPUs. After
 * selectGpuOfProcess, each process sees only its own.
 * \return the index in devices, -1 without GPU
 */
int gpuIndex();

/*! \brief Whether the GPU computes the potentials that it has kernels for
 * True if ACT was built with GPU support, a GPU was found and neither
 * turnOffGpu nor setUseGpu turned it off.
 * \return whether the GPU is used
 */
bool useGpu();

/*! \brief Never use the GPU in this process, see the option -nogpu of MsgHandler
 * Called before the GPU is used, the GPU is not started at all, such that
 * it takes no place among the processes a GPU may have, and none of its memory.
 * \param[in] reason Why, for the log file, completing "turned off ..."
 */
void turnOffGpu(const std::string &reason = "with -nogpu");

//! \return whether turnOffGpu was called
bool gpuTurnedOff();

//! \return the reason given to turnOffGpu
std::string gpuOffReason();

/*! \brief Turn the GPU on or off, for tests and comparisons
 * The GPU is only turned on if useGpu would return true without turnOffGpu.
 * \param[in] use Whether to use the GPU
 * \return whether the GPU was used before
 */
bool setUseGpu(bool use);

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
    DoubleExponential,
    //! Harmonic bond: length p[0], force constant p[1], energy p[2]
    HarmonicBond,
    //! Cubic bond: length p[0], largest length p[1], force constant p[2], energy p[3], square of the largest distance p[4]
    CubicBond,
    //! Morse bond: length p[0], beta p[1], De p[2], D0 p[3]
    MorseBond,
    //! Hua bond: length p[0], b p[1], De p[2], c p[3]
    HuaBond
};

//! Number of energy terms of the pair potentials, at most
constexpr int pairEnergyTerms = 3;

/*! \brief A pair interaction with prepared parameters, see PairKernel
 * For the potentials other than Coulomb, a nonzero p[7] marks an
 * interaction without parameters, which has no energy and force.
 */
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
 * The layout is that of computeExplicitSimd, see PairBatchBlock.
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
    //! Steps without halving the mean square force after which the shells stop, 0 for never, see ForceComputer::setStallIterations
    int    stallIterations;
    //! Largest ratio of the force across the core-shell vector to the force for a shell at the wall, see shellWallAngleTolerance
    double wallAngleToler;
    /*! Whether relaxBlocks asks the GPU every few iterations whether all
     * slots are done, and returns the statistics. Otherwise it submits all
     * iterations without waiting, and adds the statistics to those that
     * relaxStatistics returns.
     */
    bool   poll = true;
};

/*! \brief The interactions of one potential in the conformations of a molecule, see computePairBlocks
 * The layout is that of computeExplicitSimd: coordinate m of atom a in
 * slot c is x[(3*a+m)*stride + c], in units of coordinateUnit, see
 * coordinate.h, the forces likewise.
 */
struct PairBatchBlock
{
    //! Number of atoms
    size_t           natom;
    //! Distance between the rows of the arrays
    size_t           stride;
    //! Number of slots to compute, at most stride
    size_t           nslot;
    //! The coordinates
    const int64_t   *x;
    //! The potential
    PairKernel       kernel;
    //! The interactions
    const PairEntry *pairs;
    //! Number of interactions
    size_t           npair;
    //! The forces, added to, not used without forces
    float           *f;
    //! The energy terms of each slot, term k of slot c at k*stride + c, pairEnergyTerms rows, added to
    double          *energies;
};

/*! \brief Compute pair interactions of several molecules and potentials at once on the GPU
 * All blocks are computed by the same kernels, such that many small
 * molecules make enough work for a GPU. Blocks with the same coordinate
 * array share its copy on the GPU. The energy of each interaction is
 * computed in single precision and summed in double precision. The
 * results do not depend on the other blocks or slots.
 * \param[in] blocks The interactions
 * \param[in] nblock Number of blocks
 * \param[in] forces Whether to compute forces
 * \return false if the GPU is not available
 */
bool computePairBlocks(const PairBatchBlock *blocks,
                       size_t                nblock,
                       bool                  forces);

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

//! The bonded potentials of the GPU, see bondedInteraction in sycl/bondedkernels.h
enum class BondedKernel : uint32_t
{
    //! Harmonic angle i-j-k: force constant p[0], angle p[1] in radians
    HarmonicAngle,
    //! Linear angle i-j-k: a p[0], force constant p[1]
    LinearAngle,
    //! Fourier dihedral i-j-k-l: coefficients p[0..5], their number p[6]
    FourierDihedral,
    //! Harmonic improper dihedral i-j-k-l: force constant p[0]
    HarmonicDihedral
};

//! A bonded interaction with prepared parameters, see BondedKernel
struct BondedEntry
{
    //! The particles, three for angles
    uint32_t a[4];
    //! The parameters
    float    p[8];
};

//! How a virtual site is constructed from its parents, see vsitePosition in sycl/bondedkernels.h
enum class VsiteKind : uint32_t
{
    //! On the first parent
    Copy,
    //! On the line of two parents: a p[0]
    Linear2,
    //! At a fixed distance p[0] from the first parent towards the second
    Fd2,
    //! In the plane of three parents: a p[0], b p[1]
    Linear3,
    //! In the plane of three parents at a fixed distance: a p[0], b p[1]
    Fd3,
    //! Out of the plane of three parents: a p[0], b p[1], c p[2] including the sign
    Out3,
    //! A combination of four parents: a p[0], b p[1], c p[2]
    Linear4
};

/*! \brief A virtual site with prepared parameters
 * The parents are ordered as the construction routines of vsitehandler.cpp
 * get them, such that the site is constructed relative to the first.
 */
struct VsiteEntry
{
    //! The construction
    VsiteKind kind;
    //! Number of parents
    uint32_t  nparent;
    //! The parents
    uint32_t  parent[4];
    //! The site
    uint32_t  site;
    //! The parameters
    float     p[3];
};

/*! \brief A list of interactions of one potential in a block, see SystemBlock
 * The interactions without shells are computed once for each calculation,
 * those with shells also in each iteration of relaxing the shells.
 */
struct InteractionList
{
    //! Whether the interactions are bonded, otherwise pairs
    bool               bonded;
    //! The potential, a PairKernel or a BondedKernel
    uint32_t           kernel;
    //! Whether these are interactions with shells
    bool               shells;
    //! The energy array of the block that the energies are added to, one per interaction type
    uint32_t           energyIndex;
    //! Number of energy terms of the potential, the others are zero
    uint32_t           nterm = pairEnergyTerms;
    //! Number of interactions
    size_t             n;
    //! The pair interactions, if not bonded
    const PairEntry   *pairs;
    //! The bonded interactions, if bonded
    const BondedEntry *bonds;
    /*! \brief The row of each interaction in the parameter table of the list, see updateParameters
     * Not used for Coulomb, of which the parameters follow from the charges.
     * May be nullptr if the parameters are only given with the interactions.
     */
    const uint32_t    *index  = nullptr;
    //! Number of rows of the parameter table
    size_t             ntable = 0;
};

/*! \brief A topology and its conformations, resident on the GPU, see createSystem
 * The layout is that of computeExplicitSimd, see PairBatchBlock.
 */
struct SystemBlock
{
    //! Number of particles
    size_t                 natom;
    //! Distance between the rows of the coordinates
    size_t                 stride;
    //! Number of slots, the conformations
    size_t                 nslot;
    //! The coordinates
    const int64_t         *x;
    //! The interactions
    const InteractionList *lists;
    //! Number of these
    size_t                 nlist;
    //! Number of energy arrays, the interaction types with energies
    size_t                 nenergy;
    //! The virtual sites, in the order to construct them
    const VsiteEntry      *vsites;
    //! Number of these
    size_t                 nvsite;
    //! The shells in the convergence criterion of all fragments, see ShellFragment
    const ForceShell      *forceShells;
    //! The shells that move of all fragments
    const MovingShell     *movingShells;
    //! The fragments whose shells can be relaxed, their shells index the arrays above
    const ShellFragment   *fragments;
    //! Number of these
    size_t                 nfragment;
    //! The shells with one parent and the parent, to put them on their parents, two per shell
    const uint32_t        *shellParents;
    //! Number of these shells
    size_t                 nshellParent;
    /*! \brief An earlier block with the same particles and interactions, or -1
     * Such a block, a copy for another parameter set, shares the
     * interactions of that block on the GPU, and has its own parameters,
     * coordinates and forces. Its lists have to match those of that block.
     */
    int64_t                topologyOf = -1;
};

//! The data of a system resident on the GPU, see createSystem
class System;

/*! \brief Copy topologies and their conformations to the GPU
 * They stay there until destroySystem. The parameters of the interactions
 * can be changed with updateParameters, the coordinates of the slots are
 * changed on the GPU only, by virtual sites and shells. The interactions
 * of a block with SystemBlock::topologyOf are stored only once.
 * \param[in] blocks The topologies
 * \param[in] nblock Number of these
 * \return the system, or nullptr if the GPU is not available
 */
System *createSystem(const SystemBlock *blocks, size_t nblock);

//! Free the memory of a system on the GPU
void destroySystem(System *system);

/*! \brief Replace the moving shells of a block, with new force constants
 * \param[in] system The system
 * \param[in] block  The block
 * \param[in] shells The moving shells of all fragments, as for createSystem
 * \param[in] n      Their number, which has to be that of createSystem
 */
void updateShells(System *system, size_t block, const MovingShell *shells, size_t n);

/*! \brief Set the parameters of all interactions of a system, with one upload
 *
 * The parameters of each interaction are a row of the table of its list,
 * see InteractionList::index, except for Coulomb, of which they are
 * computed on the device from the charges and widths of the particles,
 * as gpuPairEntry does.
 * \param[in] sys      The system
 * \param[in] tables   The rows of 8 parameters of the tables of all lists
 *                     of all blocks, in the order of createSystem
 * \param[in] charges  For all particles of all blocks in order: the charge,
 *                     the second charge and the widths of both, 4 per
 *                     particle, or nullptr to keep them, see createChargeRecipes
 * \param[in] coulomb  The Coulomb constant, ONE_4PI_EPS0
 * \param[in] epsilonr The relative dielectric constant
 */
void updateParameters(System *sys, const float *tables, const double *charges,
                      double coulomb, double epsilonr);

/*! \brief How the device makes a row of a parameter table from the parameters of a recipe
 * See createRowRecipes and SetRecipe in act/forcefield/parameterrecipes.h.
 */
struct RowRecipe
{
    //! The row in the tables of all blocks of all replicas, see updateParameters
    uint64_t row     = 0;
    //! The set of parameters, see SetRecipe
    int32_t  set     = 0;
    //! The replica, of which the values are used
    uint32_t replica = 0;
    //! The kernel, a PairKernel, or a BondedKernel if bonded
    uint32_t kernel  = 0;
    //! Whether the interaction is bonded, see BondedKernel
    uint32_t bonded  = 0;
    //! For PairKernel::Polarization, the position of the charge of the shell in the charges, see updateParameters, otherwise -1
    int64_t  charge  = -1;
};

/*! \brief Let the device make the rows of the parameter tables from the values of the force field
 * After this, updateParameters makes the rows of the recipes from the values
 * that updateRecipeValues gave last, after the charges, such that the rows of
 * the shells use the charges of the charge generation too, and the tables
 * passed to it may be nullptr. The rows without recipe are those uploaded last.
 * \param[in] sys       The system
 * \param[in] recipes   The recipes of the sets of parameters, see ParameterRecipes::recipes
 * \param[in] nset      Their number
 * \param[in] rows      The rows made from the sets
 * \param[in] nrow      Their number
 * \param[in] factors   The factor of the unit of each value, see ParameterRecipes::factors
 * \param[in] nvalue    The number of values of each replica
 * \param[in] nreplica  The number of replicas
 */
void createRowRecipes(System *sys, const SetRecipe *recipes, size_t nset, const RowRecipe *rows, size_t nrow,
                      const double *factors, size_t nvalue, size_t nreplica);

/*! \brief How the device makes the force constants of a moving shell, see MovingShell
 * As the rows of the polarization, from the parameters of a set of the
 * polarization and the charge of the shell, see createShellVsiteRecipes.
 */
struct ShellRecipe
{
    //! The moving shell, in the order of updateMovingShells
    uint64_t shell   = 0;
    //! The set of parameters, see SetRecipe
    int32_t  set     = 0;
    //! The replica, of which the values are used
    uint32_t replica = 0;
    //! The position of the charge of the shell in the charges, see updateParameters
    int64_t  charge  = 0;
};

/*! \brief How the device makes the parameters of a virtual site, see VsiteEntry
 * Parameter k is factor[k] times parameter index[k] of the set, or unchanged
 * if index[k] is negative, see gpuVsiteParameters in act/forces/forces_gpu.h.
 */
struct VsiteRecipe
{
    //! The virtual site, in the order of updateVsites
    uint64_t vsite    = 0;
    //! The set of parameters, see SetRecipe
    int32_t  set      = 0;
    //! The replica, of which the values are used
    uint32_t replica  = 0;
    //! The parameter of the set of each parameter of the site
    int32_t  index[3] = { -1, -1, -1 };
    //! The factor of each
    double   factor[3] = { 1, 1, 1 };
};

/*! \brief Let the device make the moving shells and the virtual sites from the values too
 * After createRowRecipes, updateParameters makes the force constants of
 * these shells and the parameters of these sites as well, after which
 * updateMovingShells and updateVsites are only needed for the others.
 * \param[in] sys     The system
 * \param[in] shells  The recipes of the moving shells
 * \param[in] nshell  Their number
 * \param[in] vsites  The recipes of the virtual sites
 * \param[in] nvsite  Their number
 */
void createShellVsiteRecipes(System *sys, const ShellRecipe *shells, size_t nshell,
                             const VsiteRecipe *vsites, size_t nvsite);

/*! \brief How the device makes a fixed part of the charge of a particle
 * One of the four numbers of a particle in the charges, see
 * updateParameters, is a value of the recipes, or its value in the units of
 * the program, as FragmentHandler::updateFixedCharges gives them on the
 * processor. See createChargeRecipes.
 */
struct ChargeRecipe
{
    //! The position in the charges of all blocks of all replicas
    uint64_t position = 0;
    //! The replica, of which the values are used
    uint32_t replica  = 0;
    //! The value, see createRowRecipes
    int32_t  slot     = 0;
    //! Whether the value is converted to the units of the program, see ForceFieldParameter::internalValue
    uint32_t internal = 0;
};

/*! \brief Let the device make the fixed parts of the charges from the values too
 * After createRowRecipes, updateParameters writes these numbers of the
 * charges when the values changed, after which the charges passed to it
 * may be nullptr. The charges are generated again only if the values of
 * the charge generation changed, see updateAcm, or acmChargesChanged.
 * \param[in] sys     The system
 * \param[in] charges The recipes
 * \param[in] n       Their number
 */
void createChargeRecipes(System *sys, const ChargeRecipe *charges, size_t n);

/*! \brief Generate the charges again at the next updateParameters
 * For instance when fixed charges of the recipes changed, since the
 * charges that downloadAcmCharges gives include them.
 * \param[in] sys The system
 */
void acmChargesChanged(System *sys);

/*! \brief The values of the force field of all replicas for the recipes, see createRowRecipes
 * \param[in] sys    The system
 * \param[in] values The values, those of each replica after those of the previous one
 */
void updateRecipeValues(System *sys, const double *values);

/*! \brief Replace the moving shells of all blocks, in the order of createSystem
 * \param[in] sys    The system
 * \param[in] shells The moving shells
 * \param[in] n      Their number, which must not change
 */
void updateMovingShells(System *sys, const MovingShell *shells, size_t n);

/*! \brief Replace the virtual sites of all blocks, in the order of createSystem
 * \param[in] sys    The system
 * \param[in] vsites The virtual sites
 * \param[in] n      Their number, which must not change
 */
void updateVsites(System *sys, const VsiteEntry *vsites, size_t n);

//! No block, see AcmMolecule
constexpr uint32_t noBlock = 0xffffffff;

//! No parameter, zero, see AcmSystem
constexpr uint32_t noValue = 0xffffffff;

//! A fragment whose charges are generated on the GPU, see QgenAcm::solveSQE
struct AcmFragment
{
    //! Number of particles
    uint32_t natom;
    //! Its first particle in the molecule
    uint32_t atomStart;
    //! The first particle of its molecule, AcmMolecule::firstAtom
    uint32_t moleculeAtom;
    //! The block of the fragment alone
    uint32_t block;
    //! The particles whose charges are generated, in AcmSystem::acmParticles
    uint32_t firstAcm, nacm;
    //! The bonds, in AcmSystem::bonds
    uint32_t firstBond, nbond;
    //! The fixed charges, in AcmSystem::fixedValues
    uint32_t firstFixed, nfixed;
    //! The total charge
    int32_t  qtotal;
};

//! A molecule whose charges are generated on the GPU
struct AcmMolecule
{
    //! Number of particles
    uint32_t natom;
    //! Its first particle in AcmSystem::symmetry
    uint32_t firstAtom;
    //! Its fragments, in AcmSystem::fragments
    uint32_t firstFragment, nfragment;
    //! The block with all fragments, or noBlock
    uint32_t dimerBlock;
};

//! Charge generation with SQE on the GPU, see createAcm
struct AcmSystem
{
    const AcmMolecule *molecules;
    size_t             nmolecule;
    const AcmFragment *fragments;
    size_t             nfragment;
    //! The particles in their fragment whose charges are generated
    const uint32_t    *acmParticles;
    size_t             nacm;
    //! The coordinates in nm of the ACM particles in the structure the charges are for
    const double      *x;
    //! For each ACM particle its electronegativity, hardness and width in the values, see updateAcm
    const uint32_t    *particleValues;
    //! The two particles of each bond, as indices in the ACM particles of their fragment
    const uint32_t    *bonds;
    size_t             nbond;
    //! For each bond its electronegativity difference and hardness in the values
    const uint32_t    *bondValues;
    //! For each bond the factor of its electronegativity difference
    const double      *bondFactors;
    //! The fixed charges of the fragments in the values
    const uint32_t    *fixedValues;
    size_t             nfixed;
    //! Number of values, see updateAcm
    size_t             nvalue;
    //! For each particle of each molecule which charges are equal, see apply_symmetrized_charges
    const int32_t     *symmetry;
    //! Gaussian charges, otherwise point charges
    bool               gaussian;
    //! The Coulomb constant ONE_4PI_EPS0 and the electron volt ELECTRONVOLT, as in QgenAcm::calcJ
    double             coulomb;
    double             electronVolt;
};

/*! \brief Generate the charges with SQE on the GPU from now on
 *
 * In updateParameters, after uploading the charges, the charges of the
 * ACM particles are generated, made symmetric and stored in the blocks
 * of their fragments and molecule, before the Coulomb parameters follow.
 * \param[in] sys The system
 * \param[in] acm The molecules and their fragments
 */
void createAcm(System *sys, const AcmSystem &acm);

/*! \brief Set the values of the parameters of the charge generation, used by the next updateParameters
 * \param[in] sys    The system
 * \param[in] values The force field parameters that AcmSystem refers to
 */
void updateAcm(System *sys, const double *values);

/*! \brief Set the coordinates of the ACM particles, see AcmSystem::x
 * \param[in] sys The system
 * \param[in] x   Three per ACM particle
 */
void updateAcmCoordinates(System *sys, const double *x);

/*! \brief A virtual site of the structure of a molecule for the charge generation
 * Its parameters follow from a set of the recipes as those of VsiteRecipe,
 * see createAcmVsites.
 */
struct AcmVsite
{
    //! The construction, with the particles of the molecule, see AcmMolecule
    VsiteEntry entry;
    //! The replica, of which the values are used
    uint32_t   replica  = 0;
    //! The set of parameters, see SetRecipe
    int32_t    set      = 0;
    //! The parameter of the set of each parameter of the site, -1 if unused
    int32_t    index[3] = { -1, -1, -1 };
    //! The factor of each
    double     factor[3] = { 1, 1, 1 };
};

/*! \brief Let the device construct the virtual sites of the structures of the charge generation
 * After createRowRecipes, the device makes the coordinates of the ACM
 * particles, see updateAcmCoordinates, from the structures of the molecules
 * and their virtual sites, with the parameters of the recipes, when
 * acmVsitesChanged asks for it. As ActMol::GenerateAcmCharges, the virtual
 * sites are constructed in fixed point and converted to nm in single
 * precision.
 * \param[in] sys        The system
 * \param[in] x          The structures of the molecules in fixed point, three
 *                       per particle, in the order of AcmMolecule::firstAtom
 * \param[in] vsiteStart The first virtual site of each molecule, and the number
 *                       of sites after the last, AcmSystem::nmolecule + 1 numbers
 * \param[in] vsites     The virtual sites
 * \param[in] nvsite     Their number
 */
void createAcmVsites(System *sys, const int64_t *x, const uint32_t *vsiteStart,
                     const AcmVsite *vsites, size_t nvsite);

/*! \brief Construct the coordinates of the charge generation again, see createAcmVsites
 * At the next updateParameters, after the parameters of the recipes.
 * \param[in] sys The system
 */
void acmVsitesChanged(System *sys);

/*! \brief Download the charges and widths, as updateParameters takes them
 * \param[in]  sys     The system
 * \param[out] charges Four per particle of all blocks
 */
void downloadCharges(System *sys, double *charges);

/*! \brief Download the generated charges of the molecules of createAcm
 * \param[in]  sys     The system
 * \param[out] charges For each molecule from AcmMolecule::firstAtom on, the charge of each particle
 */
void downloadAcmCharges(System *sys, double *charges);

/*! \brief Construct the virtual sites of the slots of the blocks
 * \param[in] system The system
 */
void constructVsites(System *system);

/*! \brief Compute all interactions of all slots of all blocks
 * The forces stay on the GPU. The energies of each energy array of each
 * block, term t of slot c at (e*pairEnergyTerms + t)*stride + c for energy
 * array e, are written to energies[b] for block b.
 * \param[in]  system   The system
 * \param[out] energies The energies of each block
 */
void computeAll(System *system, double *const *energies);

/*! \brief Compute all interactions of the slots of some blocks
 * As computeAll, for the blocks with use set, the others are not changed.
 * \param[in]  system   The system
 * \param[in]  use      Whether to compute each block
 * \param[out] energies The energies of each block computed, or nullptr to
 *                      compute only the forces, without waiting for the GPU
 * \param[in]  forces   Whether to compute the forces too, which relaxBlocks needs
 */
void computeBlocks(System *system, const char *use, double *const *energies, bool forces = true);

/*! \brief Compute all interactions of the slots of some blocks, with the energies where they arrive
 * As computeBlocks with energies, without copying them: energies[b] points
 * to those of block b, if it was computed, as they arrived on the host.
 * They stay there until the next calculation of energies.
 * \param[in]  system   The system
 * \param[in]  use      Whether to compute each block
 * \param[out] energies Where the energies of each block computed are, one pointer for each block
 * \param[in]  forces   Whether to compute the forces too, which relaxBlocks needs
 */
void computeBlocksInPlace(System *system, const char *use, const double **energies, bool forces);

/*! \brief Compute all interactions of the slots of some blocks, keeping the energies on the GPU
 * As computeBlocksInPlace followed by keepEnergies, without copying the
 * energies to the host and without waiting for the GPU, for computeDeviations.
 * \param[in] system The system
 * \param[in] use    Whether to compute each block
 * \param[in] forces Whether to compute the forces too, which relaxBlocks needs
 * \param[in] before Whether these are the energies before relaxing, see keepEnergies
 */
void computeBlocksKeep(System *system, const char *use, bool forces, bool before);

/*! \brief Put the shells on their parents, see SystemBlock::shellParents
 * \param[in] system The system
 * \param[in] use    Whether to do this for each block
 */
void resetShells(System *system, const char *use);

/*! \brief Keep the shell coordinates of some blocks, or put the kept ones back
 * To start a relaxation from the shells of the previous one. For the slots
 * whose shells did not converge in their last relaxation, see relaxBlocks,
 * the positions of the parents are kept instead.
 * \param[in] system  The system
 * \param[in] use     Whether to do this for each block
 * \param[in] restore Whether to put the kept shells back, otherwise keep them
 */
void keepShells(System *system, const char *use, bool restore);

/*! \brief Relax the shells of some blocks, as relaxShells
 * The forces and coordinates are those on the GPU, the energies are not
 * computed, see computeBlocks.
 * \param[in]  system    The system
 * \param[in]  use       Whether to relax each block
 * \param[in]  relax     For each block, whether to relax the shells of each
 *                        fragment, two per block
 * \param[in]  settings  The settings
 * \param[out] nComputed The number of slots computed in all iterations
 * \param[out] notFinite Whether a displacement was not finite
 * \param[out] nFailed   The number of slots whose shells did not converge,
 *                       those that stopped without progress and those left
 *                       after the last iteration
 * Without RelaxSettings::poll the outputs are zero, see relaxStatistics.
 */
void relaxBlocks(System *system, const char *use, const char *relax,
                 const RelaxSettings &settings, size_t *nComputed, bool *notFinite,
                 size_t *nFailed);

//! The statistics of the relaxations submitted without waiting, see RelaxSettings::poll
struct RelaxStatistics
{
    //! The number of slots computed in all iterations
    uint64_t computed  = 0;
    //! The number of slots whose shells did not converge
    uint64_t failed    = 0;
    //! Whether a displacement was not finite
    bool     notFinite = false;
};

/*! \brief The statistics of the relaxations without polling since the last call
 * Waits for the GPU, the counters start again.
 * \param[in]  system The system
 * \param[out] stats  The statistics
 */
void relaxStatistics(System *system, RelaxStatistics *stats);

/*! \brief Whether a displacement was not finite in the last relaxation, on the GPU
 * For work that depends on it without waiting, such as keepShells.
 * \param[in] system The system
 * \return The flag on the device, nonzero if not finite
 */
const int *deviceNotFinite(System *system);


/*! \brief Copy the coordinates of particles from one block to another, in all slots
 * As copyRelaxedMonomerShells in forcecomputer.cpp, which copies the shells
 * of a monomer to the dimer.
 * \param[in] system The system
 * \param[in] from   The block to copy from
 * \param[in] to     The block to copy to
 * \param[in] pairs  The particle in from and in to of each coordinate to copy, two each
 * \param[in] n      Number of particles
 * \param[in] slots  For each slot of to, the slot of from with the same conformation
 */
void copyCoordinates(System *system, size_t from, size_t to, const uint32_t *pairs, size_t n,
                     const uint32_t *slots);

/*! \brief Keep the coordinates of a block, to restore them with restoreCoordinates
 * \param[in] system The system
 * \param[in] block  The block
 */
void saveCoordinates(System *system, size_t block);

/*! \brief Restore the coordinates of a block that saveCoordinates kept
 * \param[in] system The system
 * \param[in] block  The block
 */
void restoreCoordinates(System *system, size_t block);

/*! \brief Copy the coordinates and forces of a block to the host
 * For the slots, in the layout of SystemBlock::x.
 * \param[in]  system The system
 * \param[in]  block  The block
 * \param[out] x      The coordinates
 * \param[out] f      The forces, if not nullptr
 */
void downloadBlock(System *system, size_t block, int64_t *x, float *f);

//! How the energy of a type of a block is made from its energy arrays, see DeviationType
enum DeviationKind : uint32_t
{
    //! One term of one energy array, DeviationType::source
    deviationTerm     = 0,
    //! The potential energy, the sum of all terms of all energy arrays in order
    deviationEpot     = 1,
    //! The exchange plus the van der Waals correction, DeviationType::source2
    deviationExchange = 2,
    //! Zero, the van der Waals correction once added to the exchange
    deviationZero     = 3
};

//! What an energy type is for the deviations, flags of DeviationType::role
enum DeviationRole : uint32_t
{
    //! An electrostatic term, of which the change by relaxing is induction
    deviationElec                = 1,
    //! The exchange
    deviationExchangeRole        = 2,
    //! An induction from a potential
    deviationInduction           = 4,
    //! An induction correction from a potential
    deviationInductionCorrection = 8
};

//! An energy type of a block, see DeviationSystem
struct DeviationType
{
    //! The interaction type, as an integer, increasing within a block
    int32_t  type    = 0;
    //! How it is made, DeviationKind
    uint32_t kind    = deviationTerm;
    //! Energy array times pairEnergyTerms plus term, -1 if none
    int32_t  source  = -1;
    //! For deviationExchange, the source of the van der Waals correction
    int32_t  source2 = -1;
    //! DeviationRole flags
    uint32_t role    = 0;
};

//! A dimer of one replica with its reference data, see DeviationSystem
struct DeviationDimer
{
    //! The blocks of the dimer and of its monomers, noBlock if none
    uint32_t dimer      = 0;
    uint32_t monomer[2] = { 0xffffffff, 0xffffffff };
    //! The replica
    uint32_t replica    = 0;
    //! The conformations with a loss
    uint32_t nconf      = 0;
    //! For each type of the dimer the index of that type in each monomer, -1 if none, two per type
    uint64_t matchOffset = 0;
    //! The first reference conformation, see DeviationSystem::weights
    uint64_t refOffset   = 0;
    //! The first conformation of this dimer among those of all dimers, in their order
    uint64_t confOffset  = 0;
};

//! A topology of one replica with the reference potential energies, see DeviationSystem
struct DeviationMonomer
{
    //! The block
    uint32_t block     = 0;
    //! The replica
    uint32_t replica   = 0;
    //! The conformations with a loss
    uint32_t nconf     = 0;
    //! The first reference conformation, see DeviationSystem::epotWeights
    uint64_t refOffset = 0;
    //! The first conformation of this topology among those of all, in their order
    uint64_t confOffset = 0;
};

/*! \brief What computeDeviations needs, the energy types and the reference data
 * As ForceEnergyDevComputer::calcDeviationCompact. The dimers and the
 * monomers of each replica are in the order of the molecules.
 */
struct DeviationSystem
{
    //! The energy types of all blocks
    const DeviationType    *types          = nullptr;
    //! The first type of each block, and the number of types at the end
    const uint32_t         *blockFirstType = nullptr;
    //! The slot of each conformation of each block
    const uint32_t         *slotOfConf     = nullptr;
    //! The first conformation of each block in slotOfConf, and the number at the end
    const uint64_t         *blockFirstConf = nullptr;
    //! The dimers
    const DeviationDimer   *dimers         = nullptr;
    size_t                  ndimer         = 0;
    //! The type matches of the dimers, see DeviationDimer::matchOffset
    const int32_t          *matches        = nullptr;
    size_t                  nmatch         = 0;
    //! The targets, interaction types as integers
    const int32_t          *targets        = nullptr;
    size_t                  ntarget        = 0;
    //! The interaction types of the sums and of the induction
    int32_t                 allElecType    = 0;
    int32_t                 exchIndType    = 0;
    int32_t                 inductionType  = 0;
    //! The weight of each reference conformation
    const double           *weights        = nullptr;
    size_t                  nweight        = 0;
    //! The reference energy of each target of each conformation, at refOffset*ntarget + j*nconf + c
    const double           *eqm            = nullptr;
    //! Whether there is a reference energy, as eqm
    const char             *have           = nullptr;
    //! The topologies with the potential energy target, none if there is no such target
    const DeviationMonomer *monomers       = nullptr;
    size_t                  nmonomer       = 0;
    //! Their weights, reference energies and whether there are those, per conformation
    const double           *epotWeights    = nullptr;
    const double           *epotEqm        = nullptr;
    const char             *epotHave       = nullptr;
    size_t                  nepot          = 0;
    //! The number of replicas
    size_t                  nreplica       = 1;
    //! The loss function: 0 mean squared, 1 mean absolute, 2 Huber, 3 asinh
    int                     lossFunction   = 0;
};

/*! \brief Copy what computeDeviations needs to the GPU, until the next call or destroySystem
 * \param[in] system     The system
 * \param[in] deviations The energy types and the reference data
 */
void createDeviations(System *system, const DeviationSystem &deviations);

/*! \brief Keep the energies computed last on the GPU, for computeDeviations
 * \param[in] system The system
 * \param[in] use    The blocks computed last, see computeBlocksInPlace
 * \param[in] before Whether these are the energies before relaxing, otherwise after
 */
void keepEnergies(System *system, const char *use, bool before);

/*! \brief The deviations of the energies kept by keepEnergies, as ForceEnergyDevComputer::calcDeviationCompact
 * Only of the dimers and topologies whose block is used. The sums over
 * the conformations of a dimer and the terms of the conformations of a
 * topology are added in their order by the caller.
 * \param[in]  system    The system
 * \param[in]  use       Whether to compute each block
 * \param[out] dimerSums For each dimer and target, the sum of the weights and that
 *                       of the losses over its conformations, in order, 2*ntarget per dimer
 * \param[out] epotTerms For each conformation of each topology of DeviationSystem::monomers,
 *                       the weight and the loss of the potential energy, 2 per conformation
 * Without dimerSums the results stay on the GPU, for chainComputeChi2, and
 * nothing waits for the GPU.
 */
void computeDeviations(System *system, const char *use, double *dimerSums, double *epotTerms);

/*! \brief The start of the particles of a block in the arrays of all blocks
 * The charges of the blocks have four values per particle there, the
 * charge, the second charge and their widths.
 * \param[in] system The system
 * \param[in] block  The block
 * \return The index of its first particle
 */
uint64_t blockAtomOffset(const System *system, size_t block);

//! A molecule of the Monte Carlo chain on the GPU, see ChainSystem
struct ChainMolecule
{
    //! The dimer among DeviationSystem::dimers, -1 if none
    int32_t  dimer       = -1;
    //! The topology among DeviationSystem::monomers, -1 if none
    int32_t  monomer     = -1;
    //! The first term of its potential energy among those of all topologies, and their number
    uint64_t epotOffset  = 0;
    uint32_t nepot       = 0;
    //! The parameters its blocks depend on, in ChainSystem::deps
    uint64_t depOffset   = 0;
    uint32_t ndep        = 0;
    //! The bounds of its charges, in ChainSystem::bounds
    uint64_t boundOffset = 0;
    uint32_t nbound      = 0;
    //! Whether it may be skipped when its parameters did not change
    uint32_t skippable   = 0;
};

//! The bound of the charge of a particle and those that are added to it
struct ChainBound
{
    //! The charge, among those of generateAcmCharges
    uint64_t charge       = 0;
    //! The second charge, among the charges of the blocks, -1 for none
    int64_t  secondCharge = -1;
    //! The bounds
    double   minimum      = 0;
    double   maximum      = 0;
    //! The charges added, of shells and virtual sites, in ChainSystem::extras
    uint64_t extraOffset  = 0;
    uint32_t nextra       = 0;
};

//! What a term of the total chi squared is, see ChainTerm
enum ChainTermKind : int32_t
{
    //! A target of the dimers, ChainTerm::index among DeviationSystem::targets
    chainTermDimer  = 0,
    //! The potential energy of the monomers
    chainTermEpot   = 1,
    //! The bounds of the charges
    chainTermBounds = 2
};

//! A term of the total chi squared, in the order of the fitting targets
struct ChainTerm
{
    //! ChainTermKind
    int32_t kind   = chainTermDimer;
    int32_t index  = 0;
    //! The weight of the target
    double  weight = 0;
};

/*! \brief What the Monte Carlo chain on the GPU needs to compute the chi squared
 * As ACMFitnessComputer::prepareDeviation and finishDeviation for the
 * deviations of the GPU, of replica 0: which molecules are skipped since
 * the values of the parameters they depend on are those of one of the
 * two sets kept, the bounds of the charges, the sums of the targets and
 * the total, and the sets kept of each molecule.
 */
struct ChainSystem
{
    const ChainMolecule *mols      = nullptr;
    size_t               nmol      = 0;
    //! The indices of the parameters
    const uint32_t      *deps      = nullptr;
    size_t               ndep      = 0;
    const ChainBound    *bounds    = nullptr;
    size_t               nbound    = 0;
    const uint64_t      *extras    = nullptr;
    size_t               nextra    = 0;
    const ChainTerm     *terms     = nullptr;
    size_t               nterm     = 0;
    //! For each block its molecule, -1 if not of replica 0
    const int32_t       *blockMol  = nullptr;
    //! For each block what it is, 0 or 1 a monomer, 2 a dimer
    const int32_t       *blockWhat = nullptr;
    //! The number of parameters
    size_t               nparam    = 0;
    //! For each parameter its position among the values of the recipes, see updateRecipeValues, -1 if none
    const int32_t       *recipeSlot = nullptr;
    //! For each parameter its position among the values of updateAcm, -1 if none
    const int32_t       *acmIndex   = nullptr;
    //! The loss function, as DeviationSystem::lossFunction
    int                  lossFunction = 0;
};

/*! \brief Copy what the chain needs to the GPU, the sets kept are empty
 * Needs createDeviations first.
 * \param[in] system The system
 * \param[in] chain  The molecules, bounds and terms
 */
void createChain(System *system, const ChainSystem &chain);

/*! \brief The number of values of the sets kept of a molecule
 * In this order: the values of its parameters, the sums of its dimer
 * (two per target), the terms of its potential energy (two per
 * conformation) and its bound.
 * \param[in] system The system
 * \param[in] mol    The molecule
 * \return The number of values of one set
 */
size_t chainEntrySize(const System *system, size_t mol);

/*! \brief Set the sets kept of all molecules
 * \param[in] system The system
 * \param[in] values For each molecule two sets of chainEntrySize, one after the other
 * \param[in] valid  For each molecule whether each set is valid, two per molecule
 */
void chainSetCache(System *system, const double *values, const char *valid);

/*! \brief Get the sets kept of all molecules, as chainSetCache
 * \param[in]  system The system
 * \param[out] values The values
 * \param[out] valid  Whether each is valid
 */
void chainGetCache(System *system, double *values, char *valid);

/*! \brief Set the values of the parameters on the GPU
 * \param[in] system The system
 * \param[in] params The values, ChainSystem::nparam
 */
void chainSetParameters(System *system, const double *params);

/*! \brief Copy values of the parameters on the GPU to the values of the recipes and of the charge generation
 * Those of replica 0, as the processor would with updateRecipeValues and
 * updateAcm. The recipes are made again by updateParameters, and so are
 * the charges and the virtual sites of the charge generation if asked,
 * until updateParameters made them.
 * \param[in] system  The system
 * \param[in] param   The parameter, -1 for all
 * \param[in] charges Whether the charges must be generated again
 * \param[in] vsites  Whether the virtual sites of the charge generation must be made again
 */
void chainApply(System *system, int param, bool charges = true, bool vsites = true);

/*! \brief Make the recipes again in updateParameters, and the charges and virtual sites if asked, see chainApply
 * \param[in] system  The system
 * \param[in] charges Whether the charges must be generated again
 * \param[in] vsites  Whether the virtual sites of the charge generation must be made again
 */
void chainMarkChanged(System *system, bool charges, bool vsites);

/*! \brief The values of the recipes and of the charge generation of replica 0 on the GPU
 * \param[in]  system The system
 * \param[out] recipe The values of the recipes, as many as there are of a replica
 * \param[out] acm    The values of the charge generation, as many as there are of a replica
 */
void downloadValues(System *system, std::vector<double> *recipe, std::vector<double> *acm);

//! A step of the chain, with the random numbers the processor drew, see MCMCMutator::proposeStep
struct ChainStep
{
    //! The parameter changed
    uint32_t param  = 0;
    //! The random number of the change, see MCMCMutator::changeParam
    double   change = 0;
    //! The random number of the Metropolis criterion
    double   accept = 0;
    //! The inverse temperature
    double   beta   = 0;
};

/*! \brief Set up the steps of the chain
 * \param[in] system      The system
 * \param[in] steps       The steps
 * \param[in] nstep       Their number
 * \param[in] stepSize    For each parameter, the largest change
 * \param[in] lower       For each parameter, the lower bound
 * \param[in] upper       For each parameter, the upper bound
 * \param[in] bounded     For each parameter, whether it has bounds
 * \param[in] temperature For each parameter, its weight of the temperature
 * \param[in] chi2        The chi squared of the parameters on the GPU now
 * \param[in] haveBest    Whether there is a best chi squared
 * \param[in] bestChi2    The best chi squared
 */
void chainSetSteps(System *system, const ChainStep *steps, size_t nstep,
                   const double *stepSize, const double *lower, const double *upper,
                   const char *bounded, const double *temperature,
                   double chi2, bool haveBest, double bestChi2);

/*! \brief Before the evaluation of a step: change its parameter, as MCMCMutator::changeParam, and its values, as chainApply
 * The processor tells updateParameters what to make again, see chainMarkChanged.
 * \param[in] system The system
 * \param[in] step   The step
 */
void chainStepBefore(System *system, size_t step);

/*! \brief After the evaluation of a step: the chi squared and the sets kept, as chainComputeChi2, and the Metropolis criterion, as MCMCMutator::finishStep, which undoes a rejected change
 * \param[in] system The system
 * \param[in] step   The step
 */
void chainStepAfter(System *system, size_t step);

//! The results of the steps of the chain, see chainResults
struct ChainResults
{
    //! For each step whether it was accepted and whether it gave a new lowest chi squared
    std::vector<char>   accepted, newBest;
    //! For each step the chi squared and the value of the parameter tried
    std::vector<double> chi2, trial;
    //! The parameters at the end and those with the lowest chi squared
    std::vector<double> params, best;
    //! The chi squared at the end, and the lowest
    double              finalChi2 = 0, bestChi2 = 0;
    bool                haveBest  = false;
};

/*! \brief The results of the steps, waits for the GPU
 * \param[in]  system  The system
 * \param[out] results The results
 */
void chainResults(System *system, ChainResults *results);

/*! \brief Which molecules to compute, from the values of the parameters on the GPU
 * As ACMFitnessComputer::prepareDeviation: a molecule is skipped if
 * skipping is allowed, it may be skipped, and the values of its
 * parameters are those of a set kept. The masks of the blocks, see
 * chainMask, follow.
 * \param[in] system  The system
 * \param[in] canSkip Whether molecules may be skipped
 */
void chainSkip(System *system, bool canSkip);

/*! \brief The mask of the blocks to compute after chainSkip, on the GPU
 * It can be passed as use to the functions of the system.
 * \param[in] system The system
 * \param[in] what   -1 for all blocks, -2 for the monomers, 2 for the dimers
 * \return The mask on the GPU
 */
const char *chainMask(System *system, int what);

/*! \brief Which molecules are computed and which set kept each has, see chainSkip
 * \param[in]  system The system
 * \param[out] active For each molecule whether it is computed
 * \param[out] hit    For each molecule the set kept with its values, -1 if none
 */
void chainActive(System *system, char *active, int *hit);

/*! \brief The total chi squared from the deviations of computeDeviations and the sets kept
 * As ACMFitnessComputer::finishDeviation and StaticIndividualInfo::sumChiSquared.
 * With keep, the deviations of the molecules computed become the first
 * set kept of each, the one before the second, those of the molecules
 * skipped are swapped to the set with their values. Nothing waits for the
 * GPU, see chainChi2.
 * \param[in] system The system
 * \param[in] keep   Whether to keep the deviations
 */
void chainComputeChi2(System *system, bool keep);

/*! \brief The chi squared of the last chainComputeChi2
 * \param[in] system The system
 * \return The chi squared
 */
double chainChi2(System *system);

} // namespace gpu

} // namespace alexandria

#endif

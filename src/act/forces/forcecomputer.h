/*
 * This source file is part of the Alexandria Chemistry Toolkit.
 *
 * Copyright (C) 2021-2026
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
#ifndef ACT_FORCECOMPUTER_H
#define ACT_FORCECOMPUTER_H

#include <algorithm>
#include <functional>
#include <memory>
#include <set>
#include <vector>

#include "act/basics/msg_handler.h"
#include "act/forces/gpusystem.h"
#include "act/forces/vsitehandler.h"
#include "act/forcefield/forcefield.h"
#include "act/forcefield/potential.h"
#include "act/topology/actmol.h"
#include "act/topology/coords_energy.h"
#include "act/topology/topology.h"
#include "gromacs/math/vectypes.h"

namespace alexandria
{

class MsgHandler;
class QtypeProps;
class Topology;

/*! \brief Largest ratio of the force across the core-shell vector to the force on a shell at the wall
 * A shell held at the largest distance from its core by a force F pointing
 * away from it can only move across the core-shell vector. With the
 * component of the force in that direction at most this fraction of F, the
 * energy it can still gain is about F d_max theta^2/2, some 1e-6 kJ/mol for
 * F = 1e4 kJ/(mol nm) and d_max = 0.02 nm, and its force counts as zero in
 * the convergence criterion. The force across the core-shell vector is the
 * difference of two large numbers, so it often cannot be reduced below the
 * tolerance on the mean square force, see shellForce2 in forcecomputer.cpp.
 */
constexpr double shellWallAngleTolerance = 1e-4;

//! \brief How to treat the shells when computing total energies, see ForceComputer::computeTotal
struct ShellOptions
{
    //! Whether to relax the shells, otherwise they stay where they are
    bool          relax       = true;
    //! Whether to put the shells on top of their core first
    bool          resetShells = true;
    //! With resetShells, whether to start the relaxation from the input
    //! positions of the shells. The energies before relaxation, from
    //! which the induction follows, are still those with the shells on
    //! top of their core. Useful if the input positions are close to the
    //! relaxed ones, for instance predicted from the previous step in MD.
    bool          shellsFromInput = false;
    //! Indices of the shells to relax, if empty all shells are relaxed
    std::set<int> relaxShells;
    //! Electric field applied to the charges
    gmx::RVec     field       = { 0.0, 0.0, 0.0 };
    //! Whether to add INDUCTION, ALLELEC and EXCHIND as compute, otherwise
    //! the energies are those of the potentials, as in computeOnce
    bool          energyTerms = true;
    //! Whether to construct the virtual sites and spread their forces
    bool          virtualSites = true;
};

/*! \brief Conformations of one topology, with the data for computing them
 * The data of the conformations is kept between calculations, which is
 * faster if only the coordinates of some particles change. Use
 * setCoordinates or invalidate after changing the coordinates.
 */
struct ConformationSet
{
    //! The topology
    const Topology   *top = nullptr;
    //! The conformations, with the data for computing them
    ConformationData  data;

    /*! \brief Replace the conformations
     * \param[in] x The coordinates of all particles of each conformation
     */
    void setCoordinates(const std::vector<std::vector<gmx::RVec>> &x)
    {
        data.confs.resize(x.size());
        for(size_t c = 0; c < x.size(); c++)
        {
            data.setX(c, x[c]);
        }
        data.invalidate();
    }

    //! \copydoc setCoordinates(const std::vector<std::vector<gmx::RVec>> &)
    void setCoordinates(const std::vector<Coordinates> &x)
    {
        data.confs.resize(x.size());
        for(size_t c = 0; c < x.size(); c++)
        {
            data.setX(c, x[c].view());
        }
        data.invalidate();
    }

    //! \brief Mark the coordinates as changed
    void invalidate() { data.invalidate(); }
};

/*! \brief Interaction energies and forces of the conformations of a dimer
 * The coordinates and forces are stored as in ConformationArrays: component
 * m of particle a of conformation c is at (DIM*a+m)*stride + c.
 */
struct DimerInteractions
{
    //! Number of conformations
    size_t              nconf  = 0;
    //! Distance between the rows of x and f, at least nconf
    size_t              stride = 0;
    //! Coordinates of all particles, with the relaxed shells and the virtual sites
    std::vector<double> x;
    //! Interaction forces, the forces in the dimer minus those in the monomers
    std::vector<double> f;
    //! Interaction energy of each conformation
    std::vector<double> einter;
    //! Potential energy of each conformation, the interaction energy plus the energies of the monomers
    std::vector<double> epot;
};

//! Name of the option of the polarization for the maximum distance between a shell and its core
const char shellMaxDistanceName[] = "shellmaxdistance";

/*! \brief Return the maximum distance between a shell and its core
 * This is part of the force field, since the parameters are trained with
 * it. It is the option shellMaxDistanceName of the polarization, and
 * 0.02 nm if the force field does not have it. Zero means that there is
 * no maximum. A shell held at this distance by a force pointing away from
 * its core counts as relaxed when the force perpendicular to the
 * core-shell vector is small enough.
 * \param[in] pd The force field
 * \return The distance in nm
 */
double shellMaxDistance(const ForceField *pd);

/*! \brief Class to compute all the forces in a molecule or complex.
 * It is imperative that there is only one force computer in each program,
 * since results may depend on the flags given when initiating.
 */
class ForceComputer
{
private:
    //! Convergence criterium for minimizing shells: mean square force
    double         msForceToler_     = 1e-6;
    //! Maximum number of iterations to spend on minimizing shells
    int            maxiter_          = 25;
    //! Whether to start relaxing from the previous shells, see setKeepShells
    bool           keepShells_       = false;
    //! Whether to compute on the GPU, see setGpuResident
    bool           gpuResident_      = false;
    //! The number of parameter sets on the GPU, see setGpuReplicas
    size_t         gpuReplicas_      = 1;
    //! Whether to keep the energies as the GPU gives them, see setKeepGpuEnergies
    bool           keepGpuEnergies_  = false;
    //! Whether the GPU could not compute the molecules, see prepareGpuResident
    mutable bool   gpuResidentUnavailable_ = false;
    //! Why the GPU does not compute the molecules, see gpuResidentWhyNot
    mutable std::string gpuResidentWhyNot_;
    //! \return whether to try computing on the GPU, see setGpuResident
    bool useGpuResident() const;
    //! Steps without progress after which the shells stop, 0 for never, see setStallIterations
    int            stallIterations_  = 0;
    //! Relative dielectric constant
    double         epsilonr_         = 1.0;
    //! Virtual site handler
    VsiteHandler   vsiteHandler_;
    //! Total number of energy/force evaluations performed (calls to computeOnce)
    mutable size_t nEval_      = 0;
    //! Total number of shell-minimization iterations performed across all compute() calls
    mutable size_t nShellIter_ = 0;
    //! Total number of failed shell-minimizations across all compute() calls
    mutable size_t nShellConvergeFailed_ = 0;
    //! The molecules on the GPU, see computeRelaxedEnergiesGpu, shared by copies
    mutable std::shared_ptr<GpuTrainingSystem> gpuSystem_;

    /*! \brief Compute the energies of relaxed shells on the GPU
     * As computeRelaxedEnergies, with all topologies and conformations on
     * the GPU, see GpuTrainingSystem: the parameters are copied to the GPU,
     * the energies of each calculation back. Without forces.
     * \param[in]    msghandler  For messages
     * \param[in]    pd          The force field
     * \param[inout] allmols     The molecules, prepared on the GPU
     * \param[in]    resetShells Whether to put the shells on their parents first
     * \param[in]    separateInductionCorrection Whether the induction correction is separate
     * \param[in]    checkCharges Whether to compare the charges of the GPU with those of the molecules
     * \param[in]    stage       With several replicas, puts the parameters of a replica in
     *                           the force field and the molecules, see computeRelaxedEnergiesReplicas,
     *                           otherwise nullptr
     * \param[in]    deviationsOnly Whether only the deviations are needed, see computeGpuDeviations
     */
    void computeRelaxedEnergiesGpu(MsgHandler                          *msghandler,
                                   const ForceField                    *pd,
                                   const std::vector<ACTMol *>         &allmols,
                                   bool                                 resetShells,
                                   bool                                 separateInductionCorrection,
                                   bool                                 checkCharges,
                                   const std::function<void(size_t)> *stage = nullptr,
                                   bool                                 deviationsOnly = false) const;
    /*! \brief Update shell positions if needed.
     * This will work on a whole ACT mol or just one of its monomers.
     * The force constants of the shells are computed from the current
     * polarizabilities and shell charges in the topology.
     * \param[in]    top      The topology
     * \param[in]    offset   Index of the first atom of the topology in the coordinates
     * \param[in]    nShell   The number of shells
     * \param[in]    isShell  Whether each atom is a shell
     * \param[in]    fragment The fragment the shells belong to, for keeping
     *                         the step of the shells of each fragment
     * \param[in]    maxShellDistance The maximum distance between a shell and
     *                         its core, see shellMaxDistance
     * \param[inout] cdata    The conformations
     * \param[in]    active   If not nullptr, the conformations to update, the
     *                         others are not visited and not counted
     * \return The number of conformations that are done, of those in active
     *         if it is given
     */
    size_t updateShells(const Topology            *top,
                        size_t                     offset,
                        size_t                     nShell,
                        const std::vector<bool>   &isShell,
                        size_t                     fragment,
                        double                     maxShellDistance,
                        ConformationData          *cdata,
                        const std::vector<size_t> *active = nullptr) const;
    /*! \brief Relax the shells of one monomer in all conformations of all molecules
     * Iterates shell updates and energy computations until all shells
     * converged or the max number of iterations is reached.
     * \param[in] msg_handler For debugging, may be nullptr
     * \param[in] pd          The force field structure
     * \param[in] im          Index of the monomer (0 or 1)
     * \param[in] allmols     The molecules to relax
     * \param[in] allForces   Whether the forces on all particles are needed,
     *                        otherwise only those on the shells are computed
     */
    void relaxMonomer(MsgHandler                  *msg_handler,
                      const ForceField            *pd,
                      size_t                       im,
                      const std::vector<ACTMol *> &allmols,
                      bool                         allForces) const;
    /*! \brief Relax the shells in all dimer conformations of all molecules
     * Iterates shell updates and energy computations until all shells
     * converged or the max number of iterations is reached.
     * \param[in] msg_handler For debugging, may be nullptr
     * \param[in] pd          The force field structure
     * \param[in] fixMonomer  For each monomer whether its shells are kept fixed
     * \param[in] allmols     The molecules to relax
     * \param[in] allForces   Whether the forces on all particles are needed,
     *                        otherwise only those on the shells are computed
     */
    void relaxDimer(MsgHandler                  *msg_handler,
                    const ForceField            *pd,
                    bool                         fixMonomer[2],
                    const std::vector<ACTMol *> &allmols,
                    bool                         allForces) const;
    /*! \brief Replace the conformations of a dimer and compute their energies and forces
     * See computeInteraction.
     * \param[in]    msghandler                  For debugging, may be nullptr
     * \param[in]    pd                          The force field
     * \param[inout] mol                         The dimer, the conformations are stored in it
     * \param[in]    coords                      Coordinates of all particles of each conformation
     * \param[in]    separateInductionCorrection Whether to store InductionCorrection separately
     */
    void computeDimerConformations(MsgHandler                                *msghandler,
                                   const ForceField                          *pd,
                                   ACTMol                                    *mol,
                                   const std::vector<std::vector<gmx::RVec>> &coords,
                                   bool                                       separateInductionCorrection) const;

 public:
    //! \brief Default constructor
    ForceComputer() {}
    /*! \brief Constructor with initialization of variables
     * \param[in] msForce          The tolerance for the mean square force on shells
     * \param[in] maxiter          The maximum number of iterations for shell minimization
     * \param[in] epsilonr         Relative dielectric constant
     */
    ForceComputer(double   msForce,
                  int      maxiter,
                  double   epsilonr)
    {
        init(msForce, maxiter, epsilonr);
    }
    /*! \brief Initialization of variables
     * \param[in] msForce          The tolerance for the mean square force on shells
     * \param[in] maxiter          The maximum number of iterations for shell minimization
     * \param[in] epsilonr         Relative dielectric constant
     */
    void init(double   msForce,
              int      maxiter,
              double   epsilonr);

    /*! \brief Generate coordinates for virtual sites
     * \param[in]    top         Molecular topology containing vsite information
     * \param[inout] coordinates The atomic coordinates, will be updated only for vsites
     */
    void constructVsiteCoordinates(const Topology         *top,
                                   std::vector<gmx::RVec> *coordinates) const;
    /*! \brief Distribute forces on virtual sites to atoms
     * \param[in]    top         Molecular topology containing vsite information
     * \param[inout] coordinates The atomic coordinates
     * \param[inout] forces      The atomic forces will be updated and vsite forces set to zero
     */
    void spreadVsiteForces(const Topology               *top,
                           const std::vector<gmx::RVec> *coordinates,
                           std::vector<gmx::RVec>       *forces) const;

    /*! \brief Construct virtual site coordinates from integer coordinates
     * As the other constructVsiteCoordinates, from the coordinates relative
     * to the parents, which are exact, see coordinate.h.
     * \param[in]    top         Molecular topology containing vsite information
     * \param[inout] coordinates The coordinates, the virtual sites are set
     */
    void constructVsiteCoordinates(const Topology *top,
                                   Coordinates    *coordinates) const;

    /*! \brief Distribute forces on virtual sites to atoms with integer coordinates
     * \param[in]    top         Molecular topology containing vsite information
     * \param[in]    coordinates The coordinates
     * \param[inout] forces      The forces will be updated and vsite forces set to zero
     */
    void spreadVsiteForces(const Topology         *top,
                           const Coordinates      *coordinates,
                           std::vector<gmx::RVec> *forces) const;

    /*! Do one actual computations.
     * Will do one force/energy computation.
     * \param[in]  msg_handler For debugging, may be nullptr
     * \param[in]  pd          The force field structure
     * \param[in]  top         The molecular topology
     * \param[in]  coordinates The atomic coordinates. Coordinates of
     *                         shell particles may be changed.
     * \param[out] forces      The atomic forces
     * \param[out] energies    The energy components
     * \param[in]  field       Optional electric field to be applied
     */
    void computeOnce(MsgHandler                        *msg_handler,
                     const ForceField                  *pd,
                     const Topology                    *top,
                     std::vector<gmx::RVec>            *coordinates,
                     std::vector<gmx::RVec>            *forces,
                     std::map<InteractionType, double> *energies,
                     const gmx::RVec                   &field) const;

    /*! \brief Do one computation with integer coordinates
     * As the other computeOnce, without rounding the coordinates to real.
     * \param[in]  msg_handler For debugging, may be nullptr
     * \param[in]  pd          The force field structure
     * \param[in]  top         The molecular topology
     * \param[in]  coordinates The coordinates
     * \param[out] forces      The atomic forces
     * \param[out] energies    The energy components
     * \param[in]  field       Optional electric field to be applied
     */
    void computeOnce(MsgHandler                        *msg_handler,
                     const ForceField                  *pd,
                     const Topology                    *top,
                     const Coordinates                 *coordinates,
                     std::vector<gmx::RVec>            *forces,
                     std::map<InteractionType, double> *energies,
                     const gmx::RVec                   &field) const;
    /*! \brief Compute the total energies and forces of many conformations
     *
     * The conformations are computed together with the SIMD routines. The
     * virtual sites are constructed, the shells relaxed as in the options
     * and the forces on the virtual sites spread over their parents. The
     * energy terms are those of compute, including INDUCTION, ALLELEC and
     * EXCHIND.
     * \param[in]    msghandler For debugging, may be nullptr
     * \param[in]    pd         The force field
     * \param[inout] set        The conformations, on return with the forces
     *                          and the relaxed shells
     * \param[in]    options    How to treat the shells, and the electric field
     * \param[out]   energies   The energy terms of each conformation
     */
    void computeTotal(MsgHandler                                     *msghandler,
                      const ForceField                               *pd,
                      ConformationSet                                *set,
                      const ShellOptions                             &options,
                      std::vector<std::map<InteractionType, double>> *energies) const;

    /*! \brief Compute interaction energies of many dimer conformations
     *
     * The conformations are computed together with the SIMD routines, see
     * computeRelaxedEnergies. This replaces the conformations stored in the
     * molecule, their energies are in ConformationData::energies.
     * The output vectors may hold the results of a previous call, their
     * memory is used again.
     * \param[in]    msghandler                  For debugging, may be nullptr
     * \param[in]    pd                          The force field
     * \param[inout] mol                         The dimer, the conformations are stored in it
     * \param[inout] coords                      Coordinates of all particles of each
     *                                           conformation, on return with the relaxed shells
     * \param[in]    separateInductionCorrection Whether to store InductionCorrection separately
     *                                           or add it to Induction
     * \param[out]   einter                      The interaction energy terms of each conformation
     * \param[out]   epot                        The potential energy of each dimer conformation
     * \param[out]   forces                      The interaction forces of each conformation, the
     *                                           forces in the dimer minus those in the monomers
     */
    void computeInteraction(MsgHandler                                     *msghandler,
                            const ForceField                               *pd,
                            ACTMol                                         *mol,
                            std::vector<std::vector<gmx::RVec>>            *coords,
                            bool                                            separateInductionCorrection,
                            std::vector<std::map<InteractionType, double>> *einter,
                            std::vector<double>                            *epot,
                            std::vector<std::vector<gmx::RVec>>            *forces) const;

    /*! \brief Compute interaction energies and forces of many dimer conformations
     *
     * As the other computeInteraction, but only the total interaction energy
     * is returned, and the coordinates and forces are returned with the
     * conformations in the inner loop, taken from the arrays of the
     * conformations stored in the molecule.
     * \param[in]    msghandler                  For debugging, may be nullptr
     * \param[in]    pd                          The force field
     * \param[inout] mol                         The dimer, the conformations are stored in it
     * \param[in]    coords                      Coordinates of all particles of each conformation
     * \param[in]    separateInductionCorrection Whether to store InductionCorrection separately
     *                                           or add it to Induction
     * \param[out]   result                      The energies, coordinates and forces, the memory
     *                                           of a previous call is used again
     */
    void computeInteraction(MsgHandler                                *msghandler,
                            const ForceField                          *pd,
                            ACTMol                                    *mol,
                            const std::vector<std::vector<gmx::RVec>> &coords,
                            bool                                       separateInductionCorrection,
                            DimerInteractions                         *result) const;

    /*! Do complete energy/force computation.
     * If shells are present their positions will be minimized.
     * If the shell minimization does not converge a warning will
     * be written.
     * \param[in]  msg_handler For debugging, may be nullptr
     * \param[in]  pd          Pointer to force field structure
     * \param[in]  top         The molecular topology
     * \param[in]  coordinates The atomic coordinates. Coordinates of
     *                         shell particles may be changed.
     * \param[out] forces      The atomic forces
     * \param[out] energies    The energy components
     * \param[in]  field       Optional electric field to be applied
     * \param[in]  resetShells Set the position of the shells to that of the connecting atoms
     * \param[in]  relax       Indices of the shells to relax, if empty all shells are relaxed
     * \param[in]  shellsFromInput With resetShells, start the relaxation from the
     *                         input positions of the shells, see ShellOptions
     */
    void compute(MsgHandler                        *msg_handler,
                 const ForceField                  *pd,
                 const Topology                    *top,
                 std::vector<gmx::RVec>            *coordinates,
                 std::vector<gmx::RVec>            *forces,
                 std::map<InteractionType, double> *energies,
                 const gmx::RVec                   &field = { 0.0, 0.0, 0.0 },
                 bool                               resetShells = true,
                 std::set<int>                      relax = {},
                 bool                               shellsFromInput = false) const;

    /*! \brief Do complete energy/force computation with integer coordinates
     * As the other compute, without rounding the coordinates to real.
     * \param[in]  msg_handler For debugging, may be nullptr
     * \param[in]  pd          Pointer to force field structure
     * \param[in]  top         The molecular topology
     * \param[in]  coordinates The coordinates. Coordinates of shell
     *                         particles may be changed.
     * \param[out] forces      The atomic forces
     * \param[out] energies    The energy components
     * \param[in]  field       Optional electric field to be applied
     * \param[in]  resetShells Set the position of the shells to that of the connecting atoms
     * \param[in]  relax       Indices of the shells to relax, if empty all shells are relaxed
     * \param[in]  shellsFromInput With resetShells, start the relaxation from the
     *                         input positions of the shells, see ShellOptions
     */
    void compute(MsgHandler                        *msg_handler,
                 const ForceField                  *pd,
                 const Topology                    *top,
                 Coordinates                       *coordinates,
                 std::vector<gmx::RVec>            *forces,
                 std::map<InteractionType, double> *energies,
                 const gmx::RVec                   &field = { 0.0, 0.0, 0.0 },
                 bool                               resetShells = true,
                 std::set<int>                      relax = {},
                 bool                               shellsFromInput = false) const;
                 
    /*! \brief Whether to compute the relaxed energies of training on the GPU
     * With a GPU, a single precision build and a polarizable force field
     * of which the GPU supports the potentials, computeRelaxedEnergies
     * then keeps the molecules on the GPU, generates their charges and
     * relaxes their shells there, see GpuTrainingSystem. Otherwise the
     * processor computes them, as without this.
     * \param[in] resident Whether to use the GPU like this
     */
    void setGpuResident(bool resident) { gpuResident_ = resident; }

    /*! \brief Set the number of parameter sets to compute at once on the GPU
     * The molecules are on the GPU once for each, see GpuTrainingSystem and
     * computeRelaxedEnergiesReplicas. Other evaluations compute all
     * replicas with the same parameters.
     * \param[in] nreplica The number of replicas, at least 1
     */
    void setGpuReplicas(size_t nreplica) { gpuReplicas_ = std::max<size_t>(1, nreplica); }

    //! \return the number of replicas on the GPU, see setGpuReplicas
    size_t gpuReplicas() const { return gpuReplicas_; }

    /*! \brief Keep the energies of the calculations on the GPU, for gpuTypeEnergies
     * \param[in] keep Whether to keep them
     */
    void setKeepGpuEnergies(bool keep) { keepGpuEnergies_ = keep; }

    /*! \brief The energies of each interaction type of a topology from the last GPU evaluation
     * See GpuTrainingSystem::typeEnergies, needs setKeepGpuEnergies.
     * \param[in]  mol      The molecule
     * \param[in]  what     0 or 1 for a monomer, 2 for the dimer
     * \param[out] energies The interaction types, in increasing order
     * \return false if there are none
     */
    bool gpuTypeEnergies(const ACTMol *mol, int what, std::vector<GpuTypeEnergies> *energies) const
    {
        return gpuSystem_ && gpuSystem_->typeEnergies(mol, what, energies);
    }

    /*! \brief Compute the deviations of these energy targets on the GPU, in the evaluations from now on
     * See GpuTrainingSystem::setDeviations. The molecules must be on the GPU,
     * see prepareGpuResident.
     * \param[in] setup The targets
     * \return false if the GPU cannot compute them
     */
    bool setGpuDeviations(const GpuDeviationSetup &setup) const
    {
        return gpuSystem_ && gpuSystem_->setDeviations(setup);
    }

    /*! \brief The deviations of the last evaluation on the GPU, see GpuTrainingSystem::deviationSums
     * \param[out] sums The sums of the weights and of the losses of each target
     * \return false if there are none
     */
    bool gpuDeviationSums(std::vector<std::pair<double, double>> *sums) const
    {
        return gpuSystem_ && gpuSystem_->deviationSums(sums);
    }

    /*! \brief The deviations of one molecule in the last evaluation on the GPU
     * See GpuTrainingSystem::moleculeDeviations.
     * \param[in]  mol The index of the molecule in those on the GPU
     * \param[out] dev The deviations
     * \return false if there are none
     */
    bool gpuMoleculeDeviations(size_t mol, GpuMoleculeDeviations *dev) const
    {
        return gpuSystem_ && gpuSystem_->moleculeDeviations(mol, dev);
    }

    /*! \brief Compute only some molecules on the GPU in the next evaluation
     * See GpuTrainingSystem::setActiveMolecules.
     * \param[in] replica The replica, or -1 for all
     * \param[in] active  Whether each molecule on the GPU is computed
     */
    void setGpuActiveMolecules(int replica, const std::vector<char> &active) const
    {
        if (gpuSystem_)
        {
            gpuSystem_->setActiveMolecules(replica, active);
        }
    }

    /*! \brief A signature of the parameters of each molecule on the GPU
     * See GpuTrainingSystem::moleculeSignatures.
     * \param[in] pd The force field
     * \return the signatures, empty without molecules on the GPU
     */
    std::vector<uint64_t> gpuMoleculeSignatures(const ForceField *pd) const
    {
        return gpuSystem_ ? gpuSystem_->moleculeSignatures(pd, this) : std::vector<uint64_t>();
    }

    //! \return a number that changes when the molecules on the GPU are made again
    size_t gpuVersion() const
    {
        return gpuSystem_ ? gpuSystem_->version() : 0;
    }

    //! \return whether the GPU makes the rows of its parameter tables from recipes, see GpuTrainingSystem::recipeRows
    bool gpuRecipeRows() const
    {
        return gpuSystem_ && !gpuResidentUnavailable_ ? gpuSystem_->recipeRows() : false;
    }

    /*! \return why the GPU does not compute the molecules of the training
     * although it was asked to, see prepareGpuResident, or an empty string
     */
    const std::string &gpuResidentWhyNot() const { return gpuResidentWhyNot_; }

    //! \return whether the GPU makes the fixed parts of the charges, see GpuTrainingSystem::chargeRecipes
    bool gpuChargeRecipes() const
    {
        return gpuSystem_ && !gpuResidentUnavailable_ ? gpuSystem_->chargeRecipes() : false;
    }

    //! \return a number that changes when the molecules or the deviations on the GPU are made again
    size_t gpuDeviationsVersion() const
    {
        return gpuSystem_ ? gpuSystem_->deviationsVersion() : 0;
    }

    //! \return whether the GPU computes the deviations of the energies, see setGpuDeviations
    bool gpuDeviationsReady() const
    {
        return gpuSystem_ && gpuSystem_->deviationsReady();
    }

    /*! \brief Make the system of these molecules on the GPU, see setGpuResident
     * If the GPU cannot compute them, the processor computes everything
     * from then on.
     * \param[in] pd       The force field
     * \param[in] allmols  The molecules
     * \return whether the GPU computes these molecules
     */
    bool prepareGpuResident(const ForceField *pd, const std::vector<ACTMol *> &allmols) const;

    /*! \brief Compute the relaxed energies for several sets of parameters at once on the GPU
     * As computeRelaxedEnergies with the shells reset and without forces,
     * for each replica of the molecules on the GPU, see setGpuReplicas.
     * Afterwards, selectGpuReplica puts the energies and the charges of a
     * replica in the molecules.
     * \param[in]    msghandler  For messages
     * \param[in]    pd          The force field
     * \param[inout] allmols     The molecules
     * \param[in]    stage       Puts the parameters of a replica in the force field and
     *                           the topologies of the molecules, and the fixed charges
     *                           in their atoms, called for each replica in order
     * \param[in]    separateInductionCorrection Whether the induction correction is separate
     * \param[in]    deviationsOnly Whether only the deviations are needed, see computeGpuDeviations
     * \return false if the GPU cannot compute these molecules, then nothing was computed
     */
    bool computeRelaxedEnergiesReplicas(MsgHandler                          *msghandler,
                                        const ForceField                    *pd,
                                        const std::vector<ACTMol *>         &allmols,
                                        const std::function<void(size_t)> &stage,
                                        bool                                 separateInductionCorrection,
                                        bool                                 deviationsOnly = false) const;

    /*! \brief Compute only the deviations of the energies on the GPU, see setGpuDeviations
     * As computeRelaxedEnergies with the shells reset and without forces, but
     * the energies stay on the GPU, which computes the deviations of the
     * energy targets, see gpuDeviationSums. The energies of the molecules do
     * not change, their charges do.
     * \param[in]    msghandler For messages
     * \param[in]    pd         The force field
     * \param[inout] allmols    The molecules, with the deviations set up
     * \return false if this is not possible, then nothing was computed
     */
    bool computeGpuDeviations(MsgHandler                  *msghandler,
                              const ForceField            *pd,
                              const std::vector<ACTMol *> &allmols) const;

    /*! \brief Put the energies and charges of a replica in the molecules
     * See computeRelaxedEnergiesReplicas.
     * \param[in] replica The replica
     */
    void selectGpuReplica(size_t replica) const;

    /*! \brief Whether the GPU generates the charges in computeRelaxedEnergies
     * With setGpuResident or ACT_GPU_RESIDENT_CHECK, a GPU, a polarizable
     * force field and SQE with Gaussian or point charges, the GPU generates
     * the charges of the molecules that FragmentHandler::sqeOnGpu allows,
     * when the energies are computed without forces.
     * \param[in] pd The force field
     */
    bool gpuCharges(const ForceField *pd) const;

    /*! \brief Whether only the GPU generates the charges, see gpuCharges
     * Then the processor need not generate them before computeRelaxedEnergies,
     * except for ACT_GPU_RESIDENT_CHECK, which compares the charges.
     * \param[in] pd The force field
     */
    bool gpuChargesOnly(const ForceField *pd) const;

    /*! \brief Compute energies and forces for monomers and/or dimers.
     * Despite the name of the function it is possible to select RelaxShells::None
     * which results in a single energy calculation. This will happen in case of
     * a non-polarizable force field.
     * Different molecules and conformation may take different number of 
     * iterations, but the intention of the code is to keep running just until
     * the position of the selected shells have been minimized, or maxiter has been
     * reached.
     * The ACT energies are stored in the conformation data of the molecules,
     * see ACTMol::conformationData, the energy maps of the conformations hold
     * the QM energies.
     * \param[in] msghandler     For debugging and messages
     * \param[in] pd             An ACT force field
     * \param[in] allmols        The molecules to compute
     * \param[in] resetShells    Whether or not to reset the shell positions before relaxation.
     *                           This is an optimization that may or may not work. 
     * \param[in] separateInductionCorrection Whether or not to treat induction correction separately
     * \param[in] computeForces  Whether the forces are needed. If not, they are not
     *                           computed and the forces of the conformations are not
     *                           meaningful afterwards. Polarizable force fields always
     *                           compute forces, since the shells are relaxed with them.
     */
    void computeRelaxedEnergies(MsgHandler                  *msghandler,
                                const ForceField            *pd,
                                const std::vector<ACTMol *> &allmols,
                                bool                         resetShells,
                                bool                         separateInductionCorrection,
                                bool                         computeForces = true) const;

    /*! \brief Return the ACT potential used
     * In practice this converts the InteractionType to the ftype
     * used within the force computer.
     * \param[in] pd      Pointer to force field structure
     * \param[in] itype The interaction type
     * \return the potential type
     */
    Potential ftype(const ForceField *pd,
                    InteractionType   itype) const;
    
    //! \return the force tolerance
    double forceTolerance() const { return msForceToler_; }

    /*! \brief Set the shell minimization RMS tolerance
     * \param[in] toler The new tolerance
     */
    void setForceTolerance(double toler) { msForceToler_ = toler; }

    //! \return total number of force/energy evaluations performed so far
    size_t numEvaluations() const { return nEval_; }

    //! \return total number of shell minimization iterations performed so far
    size_t numShellIterations() const { return nShellIter_; }

    //! \return Number of failed shell convergence trials
    size_t numShellConvergenceFailed() const { return nShellConvergeFailed_; }

    //! \brief Reset both statistics counters to zero
    void resetStatistics() const { nEval_ = 0; nShellIter_ = 0; nShellConvergeFailed_ = 0; }
    
    //! \return the max number of shell iterations
    int maxIter() const { return maxiter_; }

    /*! \brief Set the max number of shell iterations
     * \param[in] maxiter The number!
     */
    void setMaxIter(int maxiter) { maxiter_ = maxiter; }

    //! \return the steps without progress after which the shells stop, 0 for never
    int stallIterations() const { return stallIterations_; }

    /*! \brief Stop the shells of a fragment that make no progress
     * Progress is a mean square force on the shells of at most half the
     * last one that was. Shells that a strong field holds at the largest
     * distance from their core can keep a force across the core-shell vector
     * well above the tolerance, in a cycle or stuck, for as many steps as
     * there are. A conformation that stops counts as not converged. In
     * training, stopping after 3 steps gave clearly worse force fields, since
     * many conformations that converge later stopped too, therefore no
     * program uses this by default.
     * \param[in] n The number of steps without halving the mean square
     *              force on the shells, 0 (the default) to relax until
     *              convergence or the maximum number of iterations
     */
    void setStallIterations(int n) { stallIterations_ = n; }

    /*! \brief Whether to start relaxing the shells from those of the previous relaxation
     * In computeRelaxedEnergies the shells are still put on their parents
     * for the energies before relaxation, from which the induction follows.
     * The relaxation then starts from where the shells were relaxed the
     * previous time, which needs far fewer iterations when the parameters
     * changed little. The converged positions differ within the tolerance
     * on the force, so a smaller tolerance is needed for the same precision.
     * \param[in] keep Whether to keep the shells
     */
    void setKeepShells(bool keep) { keepShells_ = keep; }

    //! \return whether the relaxation starts from the shells of the previous one
    bool keepShells() const { return keepShells_; }

    /*! \brief Plot the potential functions
     * This plots the potential functions corresponding to
     * InteractionType.
     * \param[in] msghandler Message handler
     * \param[in] pd    Pointer to force field structure
     * \param[in] itype The interaction type
     */
    void plot(MsgHandler       *msghandler,
              const ForceField *pd,
              InteractionType   itype) const;

    /*! \brief Generate vsites
     * \param[in] top         Topology
     * \param[in] coordinates The coordinates to update
     */
    void generateVsites(const Topology         *top,
                        std::vector<gmx::RVec> *coordinates)
    {
        vsiteHandler_.constructPositions(top, coordinates);
    }
};

/*! \brief Return the force scaling when hyper polarizability is used. 
 * 
 * The code in this function was derived using Mathematica:
 * V[x_] := (1/2) k1 x^2 + k2 x^4
 * F[x_] := -D[V[x], x]
 * MyX = Assuming[{k1 > 0, k2 > 0, F0 >= 0}, Solve[F0 == F[x], x]]
 * Simplify[First[MyX]]
 *
 * \param[in] k1 The harmonic force constant
 * \param[in] k2 The quartic force constant
 * \param[in] F  The norm of the force
 * \return Force scaling constant when hyper polarizability is used
 */
double calc_fscale(double k1, double k2, double F);

/*! \brief Return the ratio of the polarization force and the core-shell distance
 *
 * The polarization potential is 1/2 k1 r^2 + k2 (r - rh)^4, where the
 * quartic term acts for distances r larger than rh only.
 *
 * \param[in] k1 The harmonic force constant
 * \param[in] k2 The quartic force constant
 * \param[in] rh The distance beyond which the quartic term acts
 * \param[in] r  The core-shell distance, larger than zero
 * \return The force divided by the distance
 */
double shellForceConstant(double k1, double k2, double rh, double r);

} // namespace alexandria

#endif

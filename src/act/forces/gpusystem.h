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
 * The topologies and conformations of all molecules of a process resident
 * on the GPU, for computing the energies of relaxed shells, see
 * ForceComputer::computeRelaxedEnergies, with little communication: per
 * evaluation the parameters go to the GPU and the energies come back.
 * \author David van der Spoel <david.vanderspoel@icm.uu.se>
 */
#ifndef ACT_FORCES_GPUSYSTEM_H
#define ACT_FORCES_GPUSYSTEM_H

#include <memory>
#include <string>
#include <vector>

#include "act/basics/interactiontype.h"
#include "act/gpu/gpu.h"

namespace alexandria
{

class ACTMol;
class ForceComputer;
class ForceField;
struct DeviationReference;

/*! \brief The deviations of one molecule of one replica, see GpuTrainingSystem::moleculeDeviations
 * The caller adds those of all molecules in their order, as
 * ForceEnergyDevComputer::calcDeviationGpu.
 */
struct GpuMoleculeDeviations
{
    //! Whether it is a dimer with interaction energy targets
    bool                dimer = false;
    //! For each target, the sum of the weights and that of the losses over the conformations
    std::vector<double> targets;
    //! Whether its first topology has the potential energy target
    bool                epot  = false;
    //! For each conformation, the weight and the loss of the potential energy
    std::vector<double> epotTerms;
};

/*! \brief The energy targets that GpuTrainingSystem::computeDeviations computes
 * As ForceEnergyDevComputer::calcDeviationCompact, with the induction
 * correction not separate and without forces.
 */
struct GpuDeviationSetup
{
    //! The interaction energy targets of the dimers, in the order of the fitting targets
    std::vector<InteractionType>            targets;
    //! The molecules
    std::vector<const ACTMol *>             mols;
    //! The reference interaction energies of each molecule, nullptr if it is no dimer
    std::vector<const DeviationReference *> interaction;
    //! The number of conformations with a loss of each molecule
    std::vector<size_t>                     interactionConfs;
    //! Whether there is a potential energy target, of the first topology of each molecule
    bool                                    epot = false;
    //! The reference potential energies of each molecule, and the number of conformations with a loss
    std::vector<const DeviationReference *> epotReference;
    std::vector<size_t>                     epotConfs;
    //! The loss function, as LossFunction
    int                                     lossFunction = 0;
};

//! The energies of one interaction type of a topology, see GpuTrainingSystem::typeEnergies
struct GpuTypeEnergies
{
    //! The interaction type
    InteractionType     type = InteractionType::EPOT;
    //! The energy of each conformation after relaxing the shells
    std::vector<double> act;
    //! The energy of each conformation before relaxing, with the shells on their parents
    std::vector<double> before;
};

/*! \brief All topologies of the molecules of a process on the GPU
 * There is a block for each topology of each molecule: the monomers and,
 * for dimers, the dimer, as computeEnergies computes them.
 *
 * With several replicas, all blocks are there once for each replica, each
 * with its own parameters, charges and shells, to compute the molecules
 * for several sets of force field parameters at once: stageParameters
 * takes the parameters of each replica, uploadParameters sends them all,
 * and after computing, selectReplica puts the energies and charges of one
 * replica in the molecules.
 */
class GpuTrainingSystem
{
public:
    GpuTrainingSystem();
    ~GpuTrainingSystem();
    GpuTrainingSystem(const GpuTrainingSystem &) = delete;
    GpuTrainingSystem &operator=(const GpuTrainingSystem &) = delete;

    /*! \brief Make the system for these molecules, unless it was made already
     * \param[in] pd               The force field
     * \param[in] allmols          The molecules
     * \param[in] maxShellDistance Largest distance of a shell from its parent, 0 for none
     * \param[in] charges          Whether to generate the charges with SQE on the GPU,
     *                             for the molecules that FragmentHandler::sqeOnGpu allows
     * \param[in] fc               To construct the virtual sites of the structures the
     *                             charges are generated for
     * \param[in] nreplica         The number of replicas, at least 1
     * \return false if the GPU cannot compute these molecules, then nothing
     *         else may be called
     */
    bool prepare(const ForceField *pd, const std::vector<ACTMol *> &allmols, double maxShellDistance,
                 bool charges, const ForceComputer *fc, size_t nreplica = 1);

    //! \return the number of replicas, see prepare
    size_t nreplica() const;

    /*! \brief Copy the parameters of the topologies to the GPU, and generate the charges there
     * All replicas get the same parameters.
     * \param[in] pd The force field
     * \param[in] fc To construct the virtual sites of the structures the charges are generated for
     */
    void updateParameters(const ForceField *pd, const ForceComputer *fc);

    /*! \brief Collect the parameters of one replica, for uploadParameters
     * The parameters are those of the force field and of the topologies of
     * the molecules, and the charges those of their atoms, now.
     * \param[in] replica The replica
     * \param[in] pd      The force field
     * \param[in] fc      To construct the virtual sites of the structures the charges are generated for
     */
    void stageParameters(size_t replica, const ForceField *pd, const ForceComputer *fc);

    /*! \brief Copy the parameters that stageParameters collected to the GPU,
     * and generate the charges there
     * \param[in] pd The force field, for the dielectric constant
     */
    void uploadParameters(const ForceField *pd);

    /*! \brief Put the energies and the charges of a replica in the molecules
     * The energies of the other replicas are kept here, compute stores
     * the energies of each replica where they are.
     * \param[in] replica The replica
     */
    void selectReplica(size_t replica);

    /*! \brief The largest difference between the charges on the GPU and those of the molecules
     * \return the difference and the charge, for ACT_GPU_RESIDENT_CHECK
     */
    std::pair<double, double> chargeDifference() const;

    /*! \brief Copy the charges generated on the GPU to the molecules and their fragments,
     * for the processor, for instance for BoundsDevComputer, those of the
     * replica of selectReplica
     */
    void copyChargesToMolecules();

    //! Put the shells on their parents, as ForceComputer::computeRelaxedEnergies
    void resetShells();

    /*! \brief Keep the shells of all topologies, for restoreShells
     * Nothing is kept when a displacement was not finite since the last call.
     */
    void keepShells();

    /*! \brief Put back the shells that keepShells kept, to start the relaxation there
     * \param[in] what As for compute
     * \return whether there were shells kept
     */
    bool restoreShells(int what);

    //! Construct the virtual sites
    void constructVsites();

    /*! \brief Keep the energies of the calculations as the GPU gives them, for typeEnergies
     * \param[in] keep Whether to keep them
     */
    void setKeepEnergies(bool keep);

    /*! \brief The energies of each interaction type of a topology of a molecule
     * As compute stores them in the conformation data, for the replica of
     * selectReplica: those after relaxing are from the last calculation, those
     * before relaxing from the last one with backup. Needs setKeepEnergies.
     * \param[in]  mol      The molecule
     * \param[in]  what     0 or 1 for a monomer, 2 for the dimer
     * \param[out] energies The interaction types, in increasing order
     * \return false if there are no energies of this topology
     */
    bool typeEnergies(const ACTMol *mol, int what, std::vector<GpuTypeEnergies> *energies) const;

    /*! \brief Set what computeDeviations computes, the energy targets and their reference data
     * The reference data are copied to the GPU, unless they are those of the
     * last call. The energies are kept on the GPU from then on.
     * \param[in] setup The targets
     * \return false if the GPU cannot compute these, then computeDeviations does nothing
     */
    bool setDeviations(const GpuDeviationSetup &setup);

    //! \return whether computeDeviations computes deviations, see setDeviations
    bool deviationsReady() const;

    /*! \brief Keep the energies of compute on the GPU only, for computeDeviations
     * Then compute does not copy the energies to the molecules, nor wait for the GPU.
     * \param[in] only Whether to do so, needs setDeviations
     */
    void setDeviationsOnly(bool only);

    /*! \brief Compute the deviations of the energies of all replicas on the GPU
     * From the energies of the last calculations, as typeEnergies.
     */
    void computeDeviations();

    /*! \brief The deviations of computeDeviations of the replica of selectReplica
     * Of all molecules, which must all have been active, see setActiveMolecules.
     * \param[out] sums For each target of setDeviations and then for the potential
     *                  energy: the sum of the weights and that of the losses
     * \return false if there are none
     */
    bool deviationSums(std::vector<std::pair<double, double>> *sums) const;

    /*! \brief The deviations of one molecule of the replica of selectReplica
     * From the last computeDeviations, in which the molecule was active.
     * \param[in]  mol The molecule, its index in the molecules of prepare
     * \param[out] dev The deviations
     * \return false if there are none
     */
    bool moleculeDeviations(size_t mol, GpuMoleculeDeviations *dev) const;

    /*! \brief Compute only some molecules of a replica in the next evaluation
     * Until the end of computeRelaxedEnergies, after which all are active
     * again. The others keep their shells and are not staged, see stageParameters.
     * \param[in] replica The replica, or -1 for all
     * \param[in] active  Whether each molecule of prepare is computed
     */
    void setActiveMolecules(int replica, const std::vector<char> &active);

    /*! \brief A signature of the parameters of each molecule as they go to the GPU
     * Stages the parameters of all molecules of the first replica, from the
     * force field and the molecules now. Molecules with the same signature
     * get the same results.
     * \param[in] pd The force field
     * \param[in] fc To construct the virtual sites of the structures the charges are generated for
     * \return a hash of the parameters of each molecule
     */
    std::vector<uint64_t> moleculeSignatures(const ForceField *pd, const ForceComputer *fc);

    //! \return a number that changes when the system is made again, see prepare
    size_t version() const;

    /*! \return whether the GPU makes the rows of the parameter tables from
     * recipes, see gpu::createRowRecipes, and the moving shells, the virtual
     * sites and the structures of the charge generation as well
     */
    bool recipeRows() const;

    /*! \return whether the GPU makes the fixed parts of the charges as
     * well, see gpu::createChargeRecipes, such that the processor need not
     * update them, see FragmentHandler::updateFixedCharges
     */
    bool chargeRecipes() const;

    //! \return why prepare failed, empty if it did not
    const std::string &whyNot() const;

    //! \return a number that changes when the system or the deviations are made again, see setDeviations
    size_t deviationsVersion() const;

    /*! \brief Compute all interactions, the energies go to the conformation data
     * \param[in] what   0 or 1 for a monomer, -2 for both, 2 for the dimers, -1 for all
     * \param[in] forces Whether to compute the forces too, needed before relax
     * \param[in] backup Whether the energies of all molecules of all replicas
     *                   are backed up as well, as the energies before relaxing
     */
    void compute(int what, bool forces = true, bool backup = false);

    /*! \brief Compute only the forces, for relax, without the energies
     * Does not wait for the GPU, the energies of the molecules do not change.
     * \param[in] what As for compute
     */
    void computeForces(int what);

    //! \return whether keepShells kept shells that restoreShells can put back
    bool haveKeptShells() const;

    /*! \brief Relax the shells
     * \param[in]  what         -1 for both monomers, 2 for the dimers
     * \param[in]  relaxMonomer For dimers, whether to relax the shells of each monomer
     * \param[in]  settings     The settings
     * \param[out] nFailed      The number of conformations whose shells did not converge
     * \return the number of conformations computed in all iterations
     */
    size_t relax(int what, const bool relaxMonomer[2], const gpu::RelaxSettings &settings,
                 size_t *nFailed);

    //! Copy the shells of the relaxed monomers to the dimers, as copyRelaxedMonomerShells
    void copyRelaxedMonomerShells();

    //! Keep the coordinates of the dimers
    void saveDimers();

    //! Restore the coordinates of the dimers kept by saveDimers
    void restoreDimers();

private:
    class Impl;
    std::unique_ptr<Impl> impl_;
};

} // namespace alexandria

#endif

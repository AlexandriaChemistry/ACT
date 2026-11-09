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
#include <vector>

#include "act/gpu/gpu.h"

namespace alexandria
{

class ACTMol;
class ForceComputer;
class ForceField;

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

    /*! \brief Compute all interactions, the energies go to the conformation data
     * \param[in] what   0 or 1 for a monomer, -2 for both, 2 for the dimers, -1 for all
     * \param[in] forces Whether to compute the forces too, needed before relax
     */
    void compute(int what, bool forces = true);

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

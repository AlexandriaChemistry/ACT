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
 * Implements part of the alexandria program.
 * \author David van der Spoel <david.vanderspoel@icm.uu.se>
 */

#ifndef ACT_RINGPOLYMER_H
#define ACT_RINGPOLYMER_H

#include <cstddef>
#include <functional>
#include <map>
#include <random>
#include <utility>
#include <vector>

#include "act/basics/interactiontype.h"
#include "act/basics/msg_handler.h"
#include "act/forcefield/forcefield.h"
#include "act/forces/forcecomputer.h"
#include "act/topology/coords_energy.h"
#include "act/topology/topology.h"
#include "gromacs/math/vectypes.h"

namespace alexandria
{

/*! \brief Settings of the path integral molecular dynamics of a ring polymer
 * See RingPolymerSampler.
 */
struct RingPolymerSettings
{
    //! The number of beads, copies of the molecule. One gives classical statistics.
    int    nbeads      = 1;
    //! The temperature in K
    double temperature = 298.15;
    //! The time step in ps
    double timeStep    = 0.00025;
    //! The time constant in ps of the thermostat of the centroid of the beads
    double centroidTau = 0.05;
    //! Seed of the random numbers of the momenta and the thermostat, 0 for a random one
    int    seed        = 1993;
    //! Whether to use the thermostat. Without it the Hamiltonian of the
    //! ring polymer is conserved, for testing.
    bool   thermostat  = true;

    /*! \brief Check the settings
     * \throws gmx::InvalidInputError for settings that are out of range
     */
    void check() const;
};

/*! \brief One configuration of a ring polymer
 * The copies of a molecule that represent its quantum mechanical
 * distribution at a temperature, with their energies.
 */
struct RingPolymer
{
    /*! The coordinates of all particles of each bead in nm, with the relaxed
     * shells and the virtual sites. They are relative to the center of mass
     * of the average of the beads.
     */
    std::vector<std::vector<gmx::RVec>> x;
    //! The potential energy of each bead in kJ/mol, with relaxed shells
    std::vector<double>                 epot;
};

//! Averages over the configurations of a ring polymer simulation
struct RingPolymerAverages
{
    //! The number of configurations
    double n            = 0;
    //! Sum of the potential energy, the average over the beads (kJ/mol)
    double epot         = 0;
    //! Sum of the squared potential energy
    double epot2        = 0;
    //! Sum of the kinetic energy, see RingPolymerSampler::kineticEnergy (kJ/mol)
    double ekin         = 0;
    //! Sum of the squared kinetic energy
    double ekin2        = 0;
    //! Sum of the squared radius of gyration of the beads, see RingPolymerSampler::gyration2
    double gyration2    = 0;

    //! \return the average potential energy
    double potentialEnergy() const { return n > 0 ? epot/n : 0.0; }
    //! \return the average kinetic energy
    double kineticEnergy() const { return n > 0 ? ekin/n : 0.0; }
    //! \return the average squared radius of gyration of the beads
    double gyrationSquared() const { return n > 0 ? gyration2/n : 0.0; }
};

/*! \brief A potential of the centroids of the atoms, in addition to the force field
 *
 * It gets the centroids, the averages over the beads, of all particles in
 * nm, of which only those of particles with a mass are set, and returns
 * the energy in kJ/mol and the forces on the centroids in kJ/mol nm. See
 * RingPolymerSampler::setCentroidPotential.
 */
using CentroidPotential = std::function<double(const std::vector<gmx::DVec> &, std::vector<gmx::DVec> *)>;

/*! \brief Path integral molecular dynamics of one molecule
 *
 * The quantum mechanical distribution of the atoms of a molecule at a
 * temperature T is sampled as the classical distribution of a ring
 * polymer: P copies (beads) of the molecule, in which each atom is
 * connected to itself in the next bead by a spring with the frequency
 * \f$\omega_P = P k_B T/\hbar\f$. The ring polymer is simulated at the
 * temperature PT with the Hamiltonian
 * \f[ H_P = \sum_{k=1}^P \left[ \sum_i \frac{{\bf p}_{i,k}^2}{2m_i} +
 * \frac{m_i\omega_P^2}{2}({\bf q}_{i,k}-{\bf q}_{i,k+1})^2 + V({\bf q}_k) \right] \f]
 * With one bead this is classical molecular dynamics.
 *
 * Only particles with a mass are part of the ring polymer. Shells are
 * relaxed for each bead, and virtual sites are constructed for each bead,
 * such that V is the energy of a bead with relaxed shells. All beads are
 * computed together by the force computer.
 *
 * The integrator is BAOAB, with the motion in the springs integrated
 * exactly in the normal modes of the free ring polymer (A) and a Langevin
 * thermostat on each normal mode (O), as in the path integral Langevin
 * equation of Ceriotti et al., J. Chem. Phys. 133 (2010) 124104. The
 * thermostat of the internal modes is critically damped, that of the
 * centroid has the time constant of the settings. The center of mass of
 * the centroid does not move.
 *
 * A potential of the centroids can be added, for instance a wall that
 * keeps the compounds of a dimer together. With the energy U of the
 * centroids \f$\bar{\bf q}\f$, the Hamiltonian of the ring polymer has the
 * term \f$P U(\bar{\bf q})\f$ in addition, such that each bead of an atom gets
 * the force on its centroid. This does not change the distribution of the
 * beads around their centroid.
 */
class RingPolymerSampler
{
public:
    /*! \brief Constructor
     * \param[in] pd        The force field
     * \param[in] top       The topology of the molecule
     * \param[in] forceComp The force computer
     * \param[in] settings  The settings
     * \throws gmx::InvalidInputError for settings that are out of range, or
     * a time step that is too large for the springs of the ring polymer
     */
    RingPolymerSampler(const ForceField          *pd,
                       const Topology            *top,
                       const ForceComputer       *forceComp,
                       const RingPolymerSettings &settings);

    /*! \brief Add a potential of the centroids to the force field
     * Call it before start.
     * \param[in] potential The potential, an empty one for none
     */
    void setCentroidPotential(CentroidPotential potential);

    /*! \brief Start with all beads at the same coordinates
     * The momenta are drawn from the Maxwell distribution.
     * \param[in] msghandler  For warnings and debugging, may be nullptr
     * \param[in] x           The coordinates of all particles of the molecule
     * \param[in] energyTerms Whether to compute the energy terms, see energyTerms()
     */
    void start(MsgHandler                   *msghandler,
               const std::vector<gmx::RVec> &x,
               bool                          energyTerms = false);

    /*! \brief Do a number of time steps
     * \param[in] msghandler  For warnings and debugging, may be nullptr
     * \param[in] nsteps      The number of steps
     * \param[in] energyTerms Whether to compute the energy terms in the last
     *                        step, see energyTerms(). With shells this takes
     *                        one more energy calculation.
     */
    void run(MsgHandler *msghandler,
             int         nsteps,
             bool        energyTerms = false);

    /*! \brief Sample configurations of the ring polymer
     * \param[in]  msghandler For warnings and debugging, may be nullptr
     * \param[in]  nequil     The number of steps before the first configuration
     * \param[in]  stride     The number of steps between configurations. With
     *                        zero, all configurations are the same.
     * \param[in]  npolymers  The number of configurations
     * \param[out] pool       The configurations are added to it
     * \return the averages over the configurations
     */
    RingPolymerAverages sample(MsgHandler               *msghandler,
                               int                       nequil,
                               int                       stride,
                               size_t                    npolymers,
                               std::vector<RingPolymer> *pool);

    //! \return the present configuration
    RingPolymer polymer() const;

    //! \return the number of beads
    size_t nbeads() const { return nbeads_; }

    /*! \brief The potential energy of the force field
     * \return the energy in kJ/mol, the average over the beads, without
     *         the potential of the centroids
     */
    double potentialEnergy() const;

    //! \return the energy in kJ/mol of the potential of the centroids, see setCentroidPotential
    double centroidEnergy() const { return ecentroid_; }

    /*! \brief The energy terms of the present configuration
     * As ForceComputer::compute gives them, with the induction energy
     * with respect to the shells on their parent particles.
     * \return the terms in kJ/mol, the averages over the beads. Empty,
     *         unless they were asked for in start or in the last step of run.
     */
    const std::map<InteractionType, double> &energyTerms() const { return terms_; }

    /*! \brief The centroids of all particles
     * \param[out] x The averages over the beads of the coordinates in nm,
     *               of the shells and virtual sites as well
     */
    void centroids(std::vector<gmx::DVec> *x) const;

    /*! \brief The centroid virial estimator of the kinetic energy
     * \f[ K = \frac{3N}{2}k_B T - \frac{1}{2P}\sum_{k=1}^P\sum_i
     * ({\bf q}_{i,k}-\bar{\bf q}_i)\cdot{\bf F}_{i,k} \f]
     * for N atoms with the centroids \f$\bar{\bf q}_i\f$, including the
     * translation of the molecule.
     * \return the kinetic energy in kJ/mol
     */
    double kineticEnergy() const;

    /*! \brief The spread of the beads
     * \return the squared distance of the atoms to their centroid in nm^2,
     *         averaged over the beads and the atoms with a mass
     */
    double gyration2() const;

    /*! \brief The Hamiltonian of the ring polymer
     * \return \f$H_P\f$ in kJ/mol, which is conserved without thermostat,
     *         with P times the potential of the centroids
     */
    double hamiltonian() const;

    /*! \brief The kinetic temperature of the ring polymer
     * \return the temperature in K that corresponds to the momenta of the
     *         beads, which is P times the temperature on average
     */
    double beadTemperature() const;

private:
    //! The force field
    const ForceField                      *pd_        = nullptr;
    //! The topology of the molecule
    const Topology                        *top_       = nullptr;
    //! The force computer
    const ForceComputer                   *forceComp_ = nullptr;
    //! The settings
    RingPolymerSettings                    settings_;
    //! Number of beads
    size_t                                 nbeads_    = 1;
    //! The particles with a mass
    std::vector<size_t>                    massive_;
    //! Their masses, one for each component
    std::vector<double>                    mass_;
    //! The shells with one parent, and that parent
    std::vector<std::pair<size_t, size_t>> shells_;
    //! k_B times the temperature of the ring polymer, P times that of the molecule
    double                                 kTP_       = 0;
    //! The frequency of the springs
    double                                 omegaP_    = 0;
    //! The frequencies of the normal modes of the free ring polymer
    std::vector<double>                    omega_;
    //! The transformation from beads to normal modes, bead j and mode k at P*j+k
    std::vector<double>                    cmat_;
    //! Positions, momenta and forces: component d of bead j at n*j+d, with n components
    std::vector<double>                    q_, p_, f_;
    //! Positions and momenta in normal modes, component d of mode k at n*k+d
    std::vector<double>                    qn_, pn_;
    //! Potential energy of each bead
    std::vector<double>                    epot_;
    //! The potential of the centroids, may be empty
    CentroidPotential                      centroidPotential_;
    //! Its energy
    double                                 ecentroid_ = 0;
    //! The centroids of all particles and the forces on them, for the potential of the centroids
    std::vector<gmx::DVec>                 xcentroid_, fcentroid_;
    //! The energy terms, the averages over the beads, if they were computed
    std::map<InteractionType, double>      terms_;
    //! The beads as conformations of the molecule, with all particles
    ConformationSet                        set_;
    //! Coordinates of one bead, kept to use the memory again
    Coordinates                            scratch_;
    //! Whether the beads have coordinates
    bool                                   started_   = false;
    //! Random number generator
    std::mt19937_64                        gen_;
    //! Normal distribution
    std::normal_distribution<double>       normal_;

    /*! \brief Transform from beads to normal modes
     * \param[in]  beads Component d of bead j at n*j+d
     * \param[out] modes Component d of mode k at n*k+d
     */
    void toNormalModes(const std::vector<double> &beads,
                       std::vector<double>       *modes) const;

    /*! \brief Transform from normal modes to beads
     * \param[in]  modes Component d of mode k at n*k+d
     * \param[out] beads Component d of bead j at n*j+d
     */
    void fromNormalModes(const std::vector<double> &modes,
                         std::vector<double>       *beads) const;

    /*! \brief Move the free ring polymer exactly, in normal modes
     * \param[in] dt The time
     */
    void freeMotion(double dt);

    /*! \brief Apply the thermostat to the momenta of the normal modes
     * \param[in] dt The time
     */
    void thermostat(double dt);

    //! \brief Remove the motion of the center of mass of the centroid, in normal modes
    void stopCenterOfMass();

    /*! \brief Compute the forces and energies of all beads
     * The shells move with their parents and are relaxed, the virtual
     * sites are constructed.
     * \param[in] msghandler  For warnings and debugging, may be nullptr
     * \param[in] first       Whether these are the starting coordinates,
     *                        then the shells start from their parents
     * \param[in] energyTerms Whether to compute the energy terms
     */
    void computeForces(MsgHandler *msghandler,
                       bool        first,
                       bool        energyTerms);

    /*! \brief Do one time step
     * \param[in] msghandler  For warnings and debugging, may be nullptr
     * \param[in] energyTerms Whether to compute the energy terms
     */
    void step(MsgHandler *msghandler,
              bool        energyTerms);
};

} // namespace alexandria

#endif

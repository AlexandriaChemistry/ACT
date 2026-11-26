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

#ifndef ACT_PATHINTEGRALMD_H
#define ACT_PATHINTEGRALMD_H

#include <cstddef>
#include <map>
#include <string>
#include <vector>

#include "act/basics/interactiontype.h"
#include "gromacs/math/vectypes.h"

struct gmx_output_env_t;

namespace alexandria
{

class ForceComputer;
class ForceField;
class JsonTree;
class MsgHandler;
class SimulationConfigHandler;
class Topology;
struct SimulationObserver;

//! An average with its standard error
struct AverageError
{
    //! The average
    double average = 0;
    //! The standard error from blocks, zero if there are fewer than two complete blocks
    double error   = 0;
    //! The number of values
    size_t n       = 0;
};

/*! \brief Average of a series of values with a standard error from blocks
 *
 * The values of a simulation are correlated in time. The series is divided
 * into consecutive blocks of equal length, and the standard error of the
 * average follows from the averages of the blocks, which are independent
 * if the blocks are longer than the correlation time. The values are not
 * kept.
 */
class BlockAverager
{
public:
    /*! \brief Constructor
     * \param[in] nexpected The number of values that will be added
     * \param[in] nblocks   The number of blocks
     */
    explicit BlockAverager(size_t nexpected = 0, size_t nblocks = 8);

    /*! \brief Add a value
     * \param[in] value The value
     */
    void add(double value);

    /*! \brief The average of the values added so far
     * \return the average of all values, and the standard error from the
     *         blocks that are complete
     */
    AverageError result() const;

private:
    //! The number of values that will be added
    size_t              nexpected_ = 0;
    //! The number of values added
    size_t              n_         = 0;
    //! The sum of all values
    double              total_     = 0;
    //! The sum of the values of each block
    std::vector<double> sum_;
    //! The number of values of each block
    std::vector<size_t> count_;
};

/*! \brief Covariance of two series of values with a standard error from blocks
 *
 * As BlockAverager, for the covariance of two quantities of a simulation.
 * The covariance is that of all values, the standard error follows from
 * the covariances within the blocks.
 */
class CovarianceAverager
{
public:
    /*! \brief Constructor
     * \param[in] nexpected The number of pairs of values that will be added
     * \param[in] nblocks   The number of blocks
     */
    explicit CovarianceAverager(size_t nexpected = 0, size_t nblocks = 8);

    /*! \brief Add a pair of values
     * \param[in] a The value of the first quantity
     * \param[in] b The value of the second quantity
     */
    void add(double a, double b);

    /*! \brief The covariance of the values added so far
     * \return the average of the product minus the product of the averages
     *         of all values, and the standard error from the blocks that
     *         are complete
     */
    AverageError result() const;

private:
    //! Sums of a block
    struct Sums
    {
        //! The number of pairs
        size_t n  = 0;
        //! The sum of the first quantity
        double a  = 0;
        //! The sum of the second quantity
        double b  = 0;
        //! The sum of the product
        double ab = 0;
        //! \return the covariance, zero without values
        double covariance() const { return n > 0 ? ab/n - (a/n)*(b/n) : 0.0; }
    };
    //! The number of pairs that will be added
    size_t            nexpected_ = 0;
    //! The sums of all values
    Sums              total_;
    //! The sums of each block
    std::vector<Sums> block_;
};

//! The results of a simulation with ring polymers, see pathIntegralMD
struct PathIntegralResults
{
    //! The number of beads
    int                                     nbeads       = 0;
    //! The temperature in K
    double                                  temperature  = 0;
    //! The number of steps done
    int                                     nsteps       = 0;
    //! The number of steps at the start that are not in the averages
    int                                     nequil       = 0;
    //! Whether all steps were done, false if an observer ended the simulation
    bool                                    finished     = true;
    //! The number of fragments (compounds)
    size_t                                  nfragments   = 1;
    //! The potential energy of the force field in kJ/mol, the average over the beads
    AverageError                            epot;
    //! The kinetic energy in kJ/mol, from the centroid virial estimator, with the translation
    AverageError                            ekin;
    //! Their sum
    AverageError                            etot;
    //! The kinetic energy in kJ/mol from the primitive estimator, a check of the number of beads and the time step
    AverageError                            ekinPrimitive;
    /*! \brief The heat capacity at constant volume in kJ/mol K
     *
     * From the fluctuations of the energy,
     * \f$C_V = k_B\beta^2(\langle\epsilon_V\epsilon_T\rangle -
     * \langle\epsilon_V\rangle\langle\epsilon_T\rangle) + \frac{3N}{2}k_B\f$
     * with the total energy from the centroid virial estimator
     * \f$\epsilon_V\f$ and from the primitive estimator \f$\epsilon_T\f$,
     * including the translation of the molecule. It has no meaning with a
     * wall, and it converges slowly with the number of beads, see
     * computeThermochemistry for the correction.
     */
    AverageError                            heatCapacity;
    //! The spread of the beads in nm, see RingPolymerSampler::gyration2
    AverageError                            gyration;
    //! The distance between the fragments in nm, see FragmentWall::distance
    AverageError                            distance;
    //! The energy terms in kJ/mol, averaged over the steps in which the energies were written
    std::map<InteractionType, AverageError> terms;
    //! Whether there was a wall
    bool                                    wall         = false;
    //! The fraction of the steps in which the fragments were at the wall
    double                                  wallFraction = 0;
    //! Whether the fragments moved apart, without a wall
    bool                                    movedApart   = false;
    //! The width in nm of the bins of the histogram of the distance
    double                                  binWidth     = 0.002;
    //! The number of steps with the distance in each bin, for two fragments
    std::vector<double>                     histogram;
};

//! The files of a simulation with ring polymers, each may be nullptr for no file
struct PathIntegralFiles
{
    //! The trajectory of the centroids (pdb)
    const char             *trajectory = nullptr;
    //! The energies (xvg)
    const char             *energy     = nullptr;
    //! The progress, a line with a JSON object each time the energies are written
    const char             *progress   = nullptr;
    //! The histogram of the distance between two fragments (xvg)
    const char             *histogram  = nullptr;
    //! The output environment for the xvg files
    const gmx_output_env_t *oenv       = nullptr;
};

/*! \brief Simulate a compound, dimer or cluster at constant temperature with ring polymers
 *
 * Path integral molecular dynamics with RingPolymerSampler, with the
 * number of beads, the temperature, the time step, the number of steps
 * and the output intervals of the configuration. With one bead this is
 * classical dynamics with a Langevin thermostat. The center of mass does
 * not move. With a wall, see FragmentWall, the fragments are kept together.
 *
 * The energy file has the energy terms as averages over the beads, the
 * energy of the wall if there is one (WALL), the kinetic energy from the
 * centroid virial estimator (EKIN), which includes the translation of the
 * center of mass, the sum of the potential and the kinetic energy (ETOT),
 * the temperature of the momenta of the beads (TEMPERATURE), the spread of
 * the beads (RGYR, nm) and, with more than one fragment, the distance
 * between them (DISTANCE, nm). The trajectory has the centroids.
 *
 * \param[in] msghandler For output and warnings
 * \param[in] pd         The force field
 * \param[in] top        The topology
 * \param[in] atomStart  The first particle of each fragment, and one more
 *                       than the last particle of the last fragment
 * \param[in] x          The starting coordinates of all particles
 * \param[in] forceComp  The force computer
 * \param[in] simConfig  The settings
 * \param[in] seed       The seed of the momenta and the thermostat, 0 for a random one
 * \param[in] files      The files to write
 * \param[in] observer   Gets the energies and the centroids when they are
 *                       due, and can end the simulation, may be nullptr
 * \param[in] planckScale The factor for the constant of Planck, see
 *                       RingPolymerSettings::planckScale
 * \return the averages over the steps after the equilibration
 * \throws gmx::InvalidInputError for settings that are out of range, see
 * RingPolymerSettings::check
 */
PathIntegralResults pathIntegralMD(MsgHandler                    *msghandler,
                                   const ForceField              *pd,
                                   const Topology                *top,
                                   const std::vector<size_t>     &atomStart,
                                   const std::vector<gmx::RVec>  &x,
                                   const ForceComputer           *forceComp,
                                   const SimulationConfigHandler &simConfig,
                                   int                            seed,
                                   const PathIntegralFiles       &files,
                                   const SimulationObserver      *observer = nullptr,
                                   double                         planckScale = 1.0);

/*! \brief Write the averages of a simulation with ring polymers
 * \param[in]    msghandler For the log file
 * \param[in]    title      What was simulated
 * \param[in]    results    The results
 * \param[inout] jtree      The averages are added to it, may be nullptr
 */
void reportPathIntegral(MsgHandler                *msghandler,
                        const std::string         &title,
                        const PathIntegralResults &results,
                        JsonTree                  *jtree);

/*! \brief The energy of binding from simulations of a cluster and its fragments
 * \param[in]  cluster   The results for all fragments together
 * \param[in]  fragments The results for each fragment alone
 * \param[out] deltaH    The enthalpy of binding of two fragments in kJ/mol,
 *                       the energy minus k_B T for the change of pV of ideal gases
 * \return the energy of the cluster minus those of the fragments in kJ/mol,
 *         with the standard error from those of the simulations
 */
AverageError bindingEnergy(const PathIntegralResults              &cluster,
                           const std::vector<PathIntegralResults> &fragments,
                           AverageError                           *deltaH);

} // namespace alexandria

#endif

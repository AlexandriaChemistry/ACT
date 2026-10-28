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

#include <set>
#include <vector>

#include "act/basics/msg_handler.h"
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

/*! \brief Class to compute all the forces in a molecule or complex.
 * It is imperative that there is only one force computer in each program,
 * since results may depend on the flags given when initiating. This goes
 * in particule for the maxShellDistance.
 */
class ForceComputer
{
private:
    //! Convergence criterium for minimizing shells: mean square force
    double         msForceToler_     = 1e-6;
    //! Maximum number of iterations to spend on minimizing shells
    int            maxiter_          = 25;
    //! Maximum allowed distance for shells to be away from their core
    double         maxShellDistance_ = 0.04;
    //! Relative dielectric constant
    double         epsilonr_         = 1.0;
    //! Box for periodic boundaries
    matrix         box_   = { { 0 } };
    //! Virtual site handler
    VsiteHandler   vsiteHandler_;
    //! Total number of energy/force evaluations performed (calls to computeOnce)
    mutable size_t nEval_      = 0;
    //! Total number of shell-minimization iterations performed across all compute() calls
    mutable size_t nShellIter_ = 0;
    //! Total number of failed shell-minimizations across all compute() calls
    mutable size_t nShellConvergeFailed_ = 0;
    /*! \brief Update shell positions if needed.
     * This will work on a whole ACT mol or just one of its monomers.
     * \param[in]    mol      The ACT molecule 
     * \param[in]    top      The topology
     * \param[inout] mycoords The coordinates and energy array
     * \return The number of shells that converged
     */
    size_t updateShells(const Topology            *top,
                        size_t                     offset,
                        size_t                     nShell,
                        const std::vector<bool>    isShell,
                        const std::vector<double>  fcShell1,
                        const std::vector<double>  fcHyper,
                        std::vector<CoordsEnergy> *myCoords) const;
    void relaxMonomer(MsgHandler                  *msg_handler,
                      const ForceField            *pd,
                      size_t                       im,
                      const std::vector<ACTMol *> &allmols) const;
    void relaxDimer(MsgHandler                  *msg_handler,
                    const ForceField            *pd,
                    bool                         fixMonomer[2],
                    const std::vector<ACTMol *> &allmols) const;

 public:
    //! \brief Default constructor
    ForceComputer() {}
    /*! \brief Constructor with initialization of variables
     * \param[in] msForce          The tolerance for the mean square force on shells
     * \param[in] maxiter          The maximum number of iterations for shell minimization
     * \param[in] maxShellDistance Maximum allowed distance (nm) for shells to be away from their core
     * \param[in] epsilonr         Relative dielectric constant
     */
    ForceComputer(double   msForce,
                  int      maxiter,
                  double   maxShellDistance,
                  double   epsilonr)
    {
        init(msForce, maxiter, maxShellDistance, epsilonr);
    }
    /*! \brief Initialization of variables
     * \param[in] msForce          The tolerance for the mean square force on shells
     * \param[in] maxiter          The maximum number of iterations for shell minimization
     * \param[in] maxShellDistance Maximum allowed distance (nm) for shells to be away from their core
     */
    void init(double   msForce,
              int      maxiter,
              double   maxShellDistance,
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
    /*! Do one actual computations.
     * Will do one force/energy computation for many coordinate sets
     * \param[in]  msg_handler For debugging, may be nullptr
     * \param[in]  pd          The force field structure
     * \param[in]  mols        The molecules and/or molecule dimers
     * \param[in]  field       Optional electric field to be applied
     */
    void computeManyOnce(MsgHandler          *msg_handler,
                         const ForceField    *pd,
                         std::vector<ACTMol> *mols,
                         const gmx::RVec     &field) const;
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
     */
    void compute(MsgHandler                        *msg_handler,
                 const ForceField                  *pd,
                 const Topology                    *top,
                 std::vector<gmx::RVec>            *coordinates,
                 std::vector<gmx::RVec>            *forces,
                 std::map<InteractionType, double> *energies,
                 const gmx::RVec                   &field = { 0.0, 0.0, 0.0 },
                 bool                               resetShells = true,
                 std::set<int>                      relax = {}) const;
                 
    /*! \brief Compute energies and forces for monomers and/or dimers.
     * Despite the name of the function it is possible to select RelaxShells::None
     * which results in a single energy calculation. This will happen in case of
     * a non-polarizable force field.
     * Different molecules and conformation may take different number of 
     * iterations, but the intention of the code is to keep running just until
     * the position of the selected shells have been minimized, or maxiter has been
     * reached.
     * The ACT energies are stored in the conformation data of the molecules,
     * see ACTMol::conformationData, not in the energy maps of the conformations.
     * ACTMol::storeEnergies copies them to the energy maps.
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
                                                                
    /*! Do many complete energy/force computation.
     * If shells are present their positions will be minimized.
     * If the shell minimization does not converge a warning will
     * be written.
     * \param[in]  msg_handler For debugging, may be nullptr
     * \param[in]  pd          Pointer to force field structure
     * \param[in]  mols        The molecules and/or molecule dimers
     * \param[in]  field       Optional electric field to be applied
     * \param[in]  resetShells Set the position of the shells to that of the connecting atoms
     * \param[in]  relax       Specify the shell indices that should be relaxed. If empty, all shells will be relaxed.
     */
    void computeMany(MsgHandler          *msg_handler,
                     const ForceField    *pd,
                     std::vector<ACTMol> *mols,
                     const gmx::RVec     &field = { 0.0, 0.0, 0.0 },
                     bool                 resetShells = true,
                     std::set<int>        relax = {}) const;

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

    //! \return the maximum allowed distance for a shell to be displaced from its core
    double maxShellDistance() const { return maxShellDistance_; }

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
        vsiteHandler_.constructPositions(top, coordinates, box_);
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
real calc_fscale(real k1, real k2, real F);

} // namespace alexandria

#endif

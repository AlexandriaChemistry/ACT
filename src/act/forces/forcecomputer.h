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
    const Topology           *top = nullptr;
    //! The conformations
    std::vector<CoordsEnergy> confs;
    //! The data for computing them
    ConformationData          data;

    /*! \brief Replace the conformations
     * \param[in] x The coordinates of all particles of each conformation
     */
    void setCoordinates(const std::vector<std::vector<gmx::RVec>> &x)
    {
        confs.resize(x.size());
        for(size_t c = 0; c < x.size(); c++)
        {
            confs[c].x = x[c];
        }
        data.invalidate();
    }

    //! \brief Mark the coordinates as changed
    void invalidate() { data.invalidate(); }
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
     * \param[inout] myCoords The coordinates and energy array
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
                        std::vector<CoordsEnergy> *myCoords,
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
     * molecule, their energy maps are left empty, see ACTMol::storeEnergies.
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
real calc_fscale(real k1, real k2, real F);

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

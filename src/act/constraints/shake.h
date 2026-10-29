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
 * \author Mohammad Mehdi Ghahremanpour <mohammad.ghahremanpour@icm.uu.se>
 */

#ifndef ACT_SHAKE_H
#define ACT_SHAKE_H

#include <array>
#include <cmath>
#include <string>
#include <vector>

#include "gromacs/math/vec.h"

namespace alexandria
{

class ActAtom;
class ForceField;
class MsgHandler;
class Topology;

//! \brief What to constrain in a simulation
enum class ConstraintType {
    //! No constraints
    No,
    //! Bond lengths
    Bonds,
    //! Bond lengths and the distances between the outer atoms of angles
    Angles
};

/*! \brief Map constraint type to string
 * \param[in] ct The constraint type
 * \return string corresponding to the constraint type
 */
const std::string &constraintTypeToString(ConstraintType ct);

/*! \brief Map string to constraint type
 * \param[in] str The string
 * \return the constraint type
 * \throws if the string does not correspond to a constraint type
 */
ConstraintType stringToConstraintType(const std::string &str);

/*! \brief Distance constraints using SHAKE and RATTLE
 *
 * Constrains bond lengths and, for angle constraints, the 1-3 distances
 * of angles. All constraints are solved together using Newton iterations
 * (matrix SHAKE).
 */
class Constrainer
{
private:
    //! Max number of iterations to try constraining
    int    maxiter_ = 100;
    //! Relative tolerance for bond lengths
    double toler_   = 1e-4;
    //! Inverse time step (1/0.001 ps -> 1000)
    double invdt_   = 1000;
    //! The atom pairs with a constrained distance
    std::vector<std::array<int, 2>> pairs_;
    //! Number of bond constraints, these come first in pairs_
    size_t                 nbonds_            = 0;
    //! Number of angles that were not constrained since they are redundant
    size_t                 nredundant_        = 0;
    //! Squared constrained distances
    std::vector<double>    constraint_distance_squared_;
    //! Inverse tolerances on the squared distances
    std::vector<double>    distance_squared_tolerance_;
    //! Distance vectors of the constrained pairs in the reference coordinates
    std::vector<gmx::RVec> initial_displacements_;
    //! Inverse masses of all particles
    std::vector<double>    invmass_;
    //! Two constraints that share an atom
    struct CoupledPair
    {
        //! The first constraint
        int    a;
        //! The second constraint
        int    b;
        //! How moving the atoms of b changes the distance vector of a
        double coupling;
    };
    //! All pairs of coupled constraints, including each constraint with itself
    std::vector<CoupledPair> coupled_;

    /*! \brief Make the list of constraints, see makeConstraints
     * \param[in] msghandler For warnings and errors
     * \param[in] ctype      What to constrain
     * \param[in] pd         The force field
     * \param[in] top        The topology
     * \param[in] positions  The starting coordinates
     */
    void makeConstraintList(MsgHandler                   *msghandler,
                            ConstraintType                ctype,
                            const ForceField             *pd,
                            const Topology               *top,
                            const std::vector<gmx::RVec> &positions);

    //! \brief Find the coupled constraints
    void makeCoupledPairs();

    /*! \brief Prepare for shaking or rattling by computing the initial displacements
     * \param[in] positions  The coordinates of all particles
     */
    void computeInitialDisplacements(const std::vector<gmx::RVec> &positions);

    /*! \brief Add a constraint
     * \param[in] ai     The first atom
     * \param[in] aj     The second atom
     * \param[in] length The distance between the atoms
     */
    void addConstraint(int ai, int aj, double length);

public:
    /*! \brief Constructor
     * \param[in] maxiter Max number of iterations to try constraining
     * \param[in] toler   Relative tolerance for the constrained distances
     * \param[in] invdt   Inverse time step, used to turn SHAKE displacements into velocity corrections
     */
    Constrainer(int maxiter, double toler, double invdt) : maxiter_(maxiter), toler_(toler), invdt_(invdt)
    {}

    /*! \brief Make the list of constraints
     *
     * With bond constraints, all bonds are constrained to the reference
     * bond length. With angle constraints, the distance between the outer
     * atoms i and k of each angle i-j-k is constrained as well, to the
     * distance given by the reference bond lengths and angle. Around an
     * atom with n bonded neighbours at most 2n-3 angles are constrained,
     * since more would be redundant, and two around a planar atom with
     * three neighbours. Angles of more than 175 degrees are not
     * constrained either. Finally, an angle constraint is left out if it
     * is (nearly) linearly dependent on the other constraints in the
     * starting coordinates, for instance in rings.
     *
     * \param[in] msghandler For warnings and errors
     * \param[in] ctype      What to constrain
     * \param[in] pd         The force field
     * \param[in] top        The topology
     * \param[in] positions  The starting coordinates
     */
    void makeConstraints(MsgHandler                   *msghandler,
                         ConstraintType                ctype,
                         const ForceField             *pd,
                         const Topology               *top,
                         const std::vector<gmx::RVec> &positions);

    /*! \brief Keep atoms at their positions
     *
     * Frozen atoms are not moved by the constraints. Constraints between
     * two frozen atoms are removed. Call this after makeConstraints.
     * \param[in] frozen The indices of the frozen atoms
     */
    void setFrozen(const std::vector<int> &frozen);

    //! \return the number of constraints
    size_t numConstraints() const { return pairs_.size(); }

    //! \return the number of bond constraints
    size_t numBondConstraints() const { return nbonds_; }

    //! \return the number of angle constraints
    size_t numAngleConstraints() const { return pairs_.size() - nbonds_; }

    //! \return the number of angles that were left out since they are redundant
    size_t numRedundantAngles() const { return nredundant_; }

    //! \return the atom pairs with a constrained distance
    const std::vector<std::array<int, 2>> &pairs() const { return pairs_; }

    /*! \brief Return the constrained distance of a pair
     * \param[in] index The index of the constraint
     * \return the constrained distance
     */
    double length(size_t index) const { return std::sqrt(constraint_distance_squared_[index]); }

    /*! \brief Execute the Shake algorithm
     *
     * The constraint equations of Ryckaert, Ciccotti and Berendsen,
     * J Comp Phys, 23, 327, 1977 are solved for all constraints
     * together by Newton iterations, as in the matrix SHAKE algorithm of
     * Kraeutler, van Gunsteren and Huenenberger, J Comput Chem 22, 501,
     * 2001. Unlike iterating over the constraints one at a time, this
     * converges quickly for coupled constraints, as with angle constraints.
     *
     * \param[in]    msghandler For warnings and errors
     * \param[inout] positions  The coordinates of all particles
     * \return error code corresponding the (n+1)th constraint or zero if fine
     */
    int shake(MsgHandler             *msghandler,
              std::vector<gmx::RVec> *positions)
    {
        return shake(msghandler, *positions, positions, nullptr);
    }

    /*! \brief Execute the Shake algorithm as the first part of rattle
     *
     * The constraint forces act along the constrained pairs in the reference
     * coordinates, that should satisfy the constraints. The displacements
     * due to the constraints divided by the time step are added to the
     * velocities, as in the first part of the rattle algorithm.
     * See Andersen, J Comp Phys, 52, 24, 1983.
     *
     * \param[in]    msghandler For warnings and errors
     * \param[in]    reference  The coordinates before the integration step
     * \param[inout] positions  The coordinates of all particles
     * \param[inout] velocities The velocities for all particles, may be nullptr
     * \return error code corresponding the (n+1)th constraint or zero if fine
     */
    int shake(MsgHandler                   *msghandler,
              const std::vector<gmx::RVec> &reference,
              std::vector<gmx::RVec>       *positions,
              std::vector<gmx::RVec>       *velocities);

    /*! \brief Execute the second part of the rattle algorithm
     *
     * The velocity components along the constraints are removed by
     * solving the linear equations for all constraints together.
     * \param[in]    msghandler For warnings and errors
     * \param[in]    positions  The coordinates of all particles
     * \param[inout] velocities The velocities for all particles
     * \return zero if fine, 1 if the equations could not be solved
     */
    int rattle(MsgHandler                   *msghandler,
               const std::vector<gmx::RVec> &positions,
               std::vector<gmx::RVec>       *velocities);

    /*! \brief Remove the force components along the constraints
     *
     * Forces of constraints are added such that the accelerations of the
     * atoms do not change the constrained distances. The result is zero
     * at a minimum of the energy under the constraints, for instance for
     * the convergence check in energy minimization.
     * \param[in]    msghandler For warnings and errors
     * \param[in]    positions  The coordinates of all particles
     * \param[inout] forces     The forces on all particles
     * \return zero if fine
     */
    int projectForces(MsgHandler                   *msghandler,
                      const std::vector<gmx::RVec> &positions,
                      std::vector<gmx::RVec>       *forces);
};

} // namespace alexandria

#endif // ACT_SHAKE_H

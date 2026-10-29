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

#include "shake.h"

#include <map>
#include <set>

#pragma GCC diagnostic push
#pragma GCC diagnostic ignored "-Wpragmas"
#pragma GCC diagnostic ignored "-Wunknown-warning-option"
#pragma GCC diagnostic ignored "-Wmaybe-uninitialized"
#pragma GCC diagnostic ignored "-Wdeprecated-redundant-constexpr-static-def"
#pragma GCC diagnostic ignored "-Wdeprecated-copy-with-dtor"
#pragma GCC diagnostic ignored "-Wdeprecated-enum-enum-conversion"
#include "Eigen/Dense"
#pragma GCC diagnostic pop

#include "act/basics/msg_handler.h"
#include "act/forcefield/forcefield.h"
#include "act/forcefield/forcefield_parametername.h"
#include "act/topology/actatom.h"
#include "act/topology/topology.h"
#include "gromacs/math/units.h"
#include "gromacs/utility/exceptions.h"
#include "gromacs/utility/stringutil.h"

namespace alexandria
{

//! \brief Map from constraint type to string
static const std::map<ConstraintType, std::string> ct2str = {
    { ConstraintType::No,     "no"     },
    { ConstraintType::Bonds,  "bonds"  },
    { ConstraintType::Angles, "angles" }
};

const std::string &constraintTypeToString(ConstraintType ct)
{
    return ct2str.find(ct)->second;
}

ConstraintType stringToConstraintType(const std::string &str)
{
    for(const auto &c : ct2str)
    {
        if (c.second == str)
        {
            return c.first;
        }
    }
    GMX_THROW(gmx::InvalidInputError(gmx::formatString("Unknown constraint type '%s'", str.c_str())));
}

//! \brief Throw or report a fatal error
static void constraintError(MsgHandler *msghandler, const std::string &msg)
{
    if (msghandler)
    {
        msghandler->fatal(msg);
    }
    else
    {
        GMX_THROW(gmx::InvalidInputError(msg));
    }
}

/*! \brief Find the index of a parameter in a potential
 * \param[in] msghandler For warnings and errors
 * \param[in] pot        The potential
 * \param[in] p2i        Map from supported potentials to the parameter index
 * \return the parameter index
 */
static int parameterIndex(MsgHandler                     *msghandler,
                          Potential                       pot,
                          const std::map<Potential, int> &p2i)
{
    auto pi = p2i.find(pot);
    if (p2i.end() == pi)
    {
        constraintError(msghandler,
                        gmx::formatString("Unsupported potential %s in Constrainer",
                                          potentialToString(pot).c_str()));
        return -1;
    }
    return pi->second;
}

void Constrainer::addConstraint(int ai, int aj, double length)
{
    pairs_.push_back({ ai, aj });
    auto length2 = length*length;
    constraint_distance_squared_.push_back(length2);
    distance_squared_tolerance_.push_back(0.5/(toler_*length2));
}

/*! \brief Add the gradient of a distance to an orthonormal basis if it is independent
 *
 * The gradient of the distance between atoms ai and aj is orthogonalized
 * with respect to the basis using Gram-Schmidt. If the norm of what is
 * left, relative to the norm of the gradient, is larger than a tolerance,
 * the normalized remainder is added to the basis.
 * \param[in]    positions The coordinates
 * \param[in]    ai        The first atom
 * \param[in]    aj        The second atom
 * \param[inout] basis     The orthonormal basis of the gradients so far
 * \return whether the gradient was independent and added to the basis
 */
static bool addIndependent(const std::vector<gmx::RVec>     &positions,
                           int                               ai,
                           int                               aj,
                           std::vector<std::vector<double>> *basis)
{
    // Relative norm below which a gradient counts as dependent
    const double        indepToler = 0.05;
    std::vector<double> grad(DIM*positions.size(), 0.0);
    rvec                dx;
    rvec_sub(positions[ai], positions[aj], dx);
    double              norm = std::sqrt(2*norm2(dx));
    if (norm == 0)
    {
        return false;
    }
    for(int m = 0; m < DIM; m++)
    {
        grad[DIM*ai+m] =  dx[m]/norm;
        grad[DIM*aj+m] = -dx[m]/norm;
    }
    for(const auto &b : *basis)
    {
        double proj = 0;
        for(size_t i = 0; i < grad.size(); i++)
        {
            proj += grad[i]*b[i];
        }
        for(size_t i = 0; i < grad.size(); i++)
        {
            grad[i] -= proj*b[i];
        }
    }
    double rest = 0;
    for(auto g : grad)
    {
        rest += g*g;
    }
    rest = std::sqrt(rest);
    if (rest < indepToler)
    {
        return false;
    }
    for(auto &g : grad)
    {
        g /= rest;
    }
    basis->push_back(std::move(grad));
    return true;
}

void Constrainer::makeConstraints(MsgHandler                   *msghandler,
                                  ConstraintType                ctype,
                                  const ForceField             *pd,
                                  const Topology               *top,
                                  const std::vector<gmx::RVec> &positions)
{
    makeConstraintList(msghandler, ctype, pd, top, positions);
    makeCoupledPairs();
}

void Constrainer::makeConstraintList(MsgHandler                   *msghandler,
                                     ConstraintType                ctype,
                                     const ForceField             *pd,
                                     const Topology               *top,
                                     const std::vector<gmx::RVec> &positions)
{
    pairs_.clear();
    constraint_distance_squared_.clear();
    distance_squared_tolerance_.clear();
    invmass_.clear();
    nbonds_     = 0;
    nredundant_ = 0;
    for (const auto &a : top->atoms())
    {
        invmass_.push_back(a.invmass());
    }
    if (ConstraintType::No == ctype || !top->hasEntry(InteractionType::BONDS))
    {
        return;
    }
    // Basis of the constraint gradients, to detect redundant constraints
    std::vector<std::vector<double>> basis;
    // Bond lengths by atom pair, the lowest index first
    std::map<std::pair<int, int>, double> bondLength;
    auto blIndex = parameterIndex(msghandler,
                                  pd->findForcesConst(InteractionType::BONDS).potential(),
                                  { { Potential::HUA_BONDS, huaLENGTH },
                                    { Potential::MORSE_BONDS, morseLENGTH },
                                    { Potential::HARMONIC_BONDS, bondLENGTH },
                                    { Potential::CUBIC_BONDS, cubicLENGTH } });
    for (const auto &b : top->entry(InteractionType::BONDS))
    {
        auto &indices = b->atomIndices();
        auto  length  = b->params()[blIndex];
        bondLength.insert({ std::minmax(indices[0], indices[1]), length });
        addIndependent(positions, indices[0], indices[1], &basis);
        addConstraint(indices[0], indices[1], length);
    }
    nbonds_ = pairs_.size();
    if (ConstraintType::Angles != ctype || !top->hasEntry(InteractionType::ANGLES))
    {
        return;
    }
    auto thIndex = parameterIndex(msghandler,
                                  pd->findForcesConst(InteractionType::ANGLES).potential(),
                                  { { Potential::HARMONIC_ANGLES, angleANGLE },
                                    { Potential::UREY_BRADLEY_ANGLES, ubANGLE } });
    auto &angles = top->entry(InteractionType::ANGLES);
    // Number of bonded neighbours and sum of the reference angles per central atom
    std::vector<int>    nbonded(invmass_.size(), 0);
    std::vector<double> angleSum(invmass_.size(), 0.0);
    for (const auto &b : bondLength)
    {
        nbonded[b.first.first]  += 1;
        nbonded[b.first.second] += 1;
    }
    for (const auto &a : angles)
    {
        angleSum[a->atomIndices()[1]] += a->params()[thIndex];
    }
    // Number of angle constraints per central atom so far
    std::vector<int> nconstrained(invmass_.size(), 0);
    for (const auto &a : angles)
    {
        auto &indices = a->atomIndices();
        auto  aj      = indices[1];
        // With n neighbours, 2n-3 angles fix the positions of the
        // neighbours with respect to each other, except for a planar
        // atom with three neighbours, where the third angle follows
        // from the other two. Linear angles do not fix anything.
        int  nmax     = 2*nbonded[aj] - 3;
        if (3 == nbonded[aj] && angleSum[aj] > 355*DEG2RAD)
        {
            nmax = 2;
        }
        if (nconstrained[aj] >= nmax || a->params()[thIndex] > 175*DEG2RAD)
        {
            nredundant_ += 1;
            continue;
        }
        auto  bij     = bondLength.find(std::minmax(indices[0], indices[1]));
        auto  bjk     = bondLength.find(std::minmax(indices[1], indices[2]));
        if (bondLength.end() == bij || bondLength.end() == bjk)
        {
            constraintError(msghandler,
                            gmx::formatString("Angle %d-%d-%d without bonds in Constrainer",
                                              indices[0]+1, indices[1]+1, indices[2]+1));
            continue;
        }
        if (!addIndependent(positions, indices[0], indices[2], &basis))
        {
            nredundant_ += 1;
            continue;
        }
        nconstrained[aj] += 1;
        // Law of cosines, the angle is in radians
        double rij   = bij->second;
        double rjk   = bjk->second;
        double theta = a->params()[thIndex];
        addConstraint(indices[0], indices[2],
                      std::sqrt(rij*rij + rjk*rjk - 2*rij*rjk*std::cos(theta)));
    }
}

void Constrainer::computeInitialDisplacements(const std::vector<gmx::RVec> &positions)
{
    initial_displacements_.clear();
    for (const auto &p : pairs_)
    {
        rvec dx;
        rvec_sub(positions[p[0]], positions[p[1]], dx);
        initial_displacements_.push_back(dx);
    }
}

/*! \brief Compute the coupling coefficient of two constraints
 *
 * When constraint b moves its atoms along a vector s, as
 * x_ib += invmass_ib g s and x_jb -= invmass_jb g s, the distance
 * vector of constraint a changes by coupling(a, b) g s.
 * \param[in] pa      The atoms of constraint a
 * \param[in] pb      The atoms of constraint b
 * \param[in] invmass The inverse masses
 * \return the coupling coefficient
 */
static double coupling(const std::array<int, 2> &pa,
                       const std::array<int, 2> &pb,
                       const std::vector<double> &invmass)
{
    double c = 0;
    if (pa[0] == pb[0]) { c += invmass[pa[0]]; }
    if (pa[0] == pb[1]) { c -= invmass[pa[0]]; }
    if (pa[1] == pb[0]) { c -= invmass[pa[1]]; }
    if (pa[1] == pb[1]) { c += invmass[pa[1]]; }
    return c;
}

void Constrainer::makeCoupledPairs()
{
    // Constraints sharing an atom are coupled
    std::vector<std::vector<int>> atomConstraints(invmass_.size());
    for (size_t a = 0; a < pairs_.size(); a++)
    {
        atomConstraints[pairs_[a][0]].push_back(a);
        atomConstraints[pairs_[a][1]].push_back(a);
    }
    coupled_.clear();
    for (size_t a = 0; a < pairs_.size(); a++)
    {
        std::set<int> bs;
        for (auto ai : pairs_[a])
        {
            bs.insert(atomConstraints[ai].begin(), atomConstraints[ai].end());
        }
        for (auto b : bs)
        {
            coupled_.push_back({ static_cast<int>(a), b,
                                 coupling(pairs_[a], pairs_[b], invmass_) });
        }
    }
}

/*! \brief Move the atoms of the constraints along vectors
 * \param[in]    pairs   The constrained atom pairs
 * \param[in]    invmass The inverse masses
 * \param[in]    vecs    The vector per constraint
 * \param[in]    factor  The factor per constraint
 * \param[inout] x       The coordinates or velocities to change
 */
static void displace(const std::vector<std::array<int, 2>> &pairs,
                     const std::vector<double>             &invmass,
                     const std::vector<gmx::RVec>          &vecs,
                     const Eigen::VectorXd                 &factor,
                     std::vector<gmx::RVec>                *x)
{
    for (size_t b = 0; b < pairs.size(); b++)
    {
        auto ai = pairs[b][0];
        auto aj = pairs[b][1];
        for (int m = 0; m < DIM; m++)
        {
            (*x)[ai][m] += invmass[ai]*factor[b]*vecs[b][m];
            (*x)[aj][m] -= invmass[aj]*factor[b]*vecs[b][m];
        }
    }
}

int Constrainer::shake(gmx_unused MsgHandler        *msghandler,
                       const std::vector<gmx::RVec> &reference,
                       std::vector<gmx::RVec>       *positions,
                       std::vector<gmx::RVec>       *velocities)
{
    if (pairs_.empty())
    {
        return 0;
    }
    computeInitialDisplacements(reference);
    // Keep the unconstrained positions to correct the velocities
    std::vector<gmx::RVec> unconstrained;
    if (velocities)
    {
        unconstrained = *positions;
    }
    // Newton iterations on all constraints together, see
    // Kraeutler, van Gunsteren and Huenenberger, J Comput Chem 22, 501 (2001).
    // The constraint forces act along the reference displacements s_b.
    // The step is halved as long as it increases the deviations, since
    // far from the solution a full step may overshoot.
    auto                   nc = pairs_.size();
    Eigen::MatrixXd        A(nc, nc);
    Eigen::VectorXd        diff(nc);
    std::vector<gmx::RVec> dx(nc);
    std::vector<gmx::RVec> start;
    // Compute the deviations, return the sum of the squared relative deviations
    auto deviations = [&]()
    {
        double sum = 0;
        for (size_t a = 0; a < nc; a++)
        {
            rvec_sub((*positions)[pairs_[a][0]], (*positions)[pairs_[a][1]], dx[a]);
            diff[a]    = constraint_distance_squared_[a] - norm2(dx[a]);
            double rel = diff[a] * distance_squared_tolerance_[a];
            sum       += rel*rel;
        }
        return sum;
    };
    double dev   = deviations();
    int    error = -1;
    for (int iter = 0; iter < maxiter_; iter++)
    {
        /* The error is smaller than the bound if this is at most 1 */
        bool converged = true;
        for (size_t a = 0; a < nc && converged; a++)
        {
            converged = (std::abs(diff[a]) * distance_squared_tolerance_[a] <= 1.0);
        }
        if (converged)
        {
            error = 0;
            break;
        }
        // Linearized change of the squared distances
        A.setZero();
        for (const auto &c : coupled_)
        {
            A(c.a, c.b) = 2*c.coupling*iprod(dx[c.a], initial_displacements_[c.b]);
        }
        Eigen::VectorXd g = A.partialPivLu().solve(diff);
        if (!g.allFinite())
        {
            error = 1;
            break;
        }
        start = *positions;
        for (int halve = 0; halve < 20; halve++)
        {
            displace(pairs_, invmass_, initial_displacements_, g, positions);
            double newDev = deviations();
            if (newDev < dev)
            {
                dev = newDev;
                break;
            }
            *positions = start;
            g         *= 0.5;
        }
        if (*positions == start)
        {
            // No step reduces the deviations
            deviations();
            error = 1;
            break;
        }
    }
    if (velocities)
    {
        for(size_t i = 0; i < positions->size(); i++)
        {
            for(int m = 0; m < DIM; m++)
            {
                (*velocities)[i][m] += ((*positions)[i][m] - unconstrained[i][m])*invdt_;
            }
        }
    }
    return error;
}

void Constrainer::setFrozen(const std::vector<int> &frozen)
{
    for (auto f : frozen)
    {
        invmass_[f] = 0;
    }
    std::vector<std::array<int, 2>> pairs;
    std::vector<double>             d2, tol;
    size_t                          nbonds = 0;
    for (size_t a = 0; a < pairs_.size(); a++)
    {
        if (invmass_[pairs_[a][0]] > 0 || invmass_[pairs_[a][1]] > 0)
        {
            pairs.push_back(pairs_[a]);
            d2.push_back(constraint_distance_squared_[a]);
            tol.push_back(distance_squared_tolerance_[a]);
            if (a < nbonds_)
            {
                nbonds += 1;
            }
        }
    }
    pairs_                       = pairs;
    constraint_distance_squared_ = d2;
    distance_squared_tolerance_  = tol;
    nbonds_                      = nbonds;
    makeCoupledPairs();
}

/*! \brief Remove the components along the constraints from vectors
 *
 * Solves the linear equations for the multipliers h such that
 * r_a . (u_a1 + invmass_a1 Delta_a1 - u_a2 - invmass_a2 Delta_a2) = 0
 * for all constraints a, where Delta_i = sum_b sigma_ib h_b r_b.
 * \param[in]    pairs   The constrained atom pairs
 * \param[in]    coupled The coupled constraints
 * \param[in]    invmass The inverse masses
 * \param[in]    r       The distance vectors of the constraints
 * \param[in]    u       The vectors (velocities) per atom
 * \param[out]   h       The multipliers
 * \return whether the solution is finite
 */
template <typename CoupledPair>
static bool solveProjection(const std::vector<std::array<int, 2>> &pairs,
                            const std::vector<CoupledPair>        &coupled,
                            const std::vector<gmx::RVec>          &r,
                            const std::vector<gmx::RVec>          &u,
                            Eigen::VectorXd                       *h)
{
    auto            nc = pairs.size();
    Eigen::MatrixXd A  = Eigen::MatrixXd::Zero(nc, nc);
    Eigen::VectorXd uproj(nc);
    for (size_t a = 0; a < nc; a++)
    {
        rvec du;
        rvec_sub(u[pairs[a][0]], u[pairs[a][1]], du);
        uproj[a] = -iprod(du, r[a]);
    }
    for (const auto &c : coupled)
    {
        A(c.a, c.b) = c.coupling*iprod(r[c.a], r[c.b]);
    }
    *h = A.partialPivLu().solve(uproj);
    return h->allFinite();
}

int Constrainer::rattle(MsgHandler                   *msghandler,
                        const std::vector<gmx::RVec> &positions,
                        std::vector<gmx::RVec>       *velocities)
{
    if (pairs_.empty())
    {
        return 0;
    }
    computeInitialDisplacements(positions);
    // The velocity constraints are linear, such that one solution of the
    // linear equations removes all velocity components along the constraints.
    Eigen::VectorXd h;
    int error = 0;
    if (solveProjection(pairs_, coupled_, initial_displacements_, *velocities, &h))
    {
        displace(pairs_, invmass_, initial_displacements_, h, velocities);
    }
    else
    {
        error = 1;
    }
    if (msghandler && msghandler->debug())
    {
        msghandler->writeDebug(gmx::formatString("Rattled with error code %d", error));
    }
    return error;
}

int Constrainer::projectForces(gmx_unused MsgHandler        *msghandler,
                               const std::vector<gmx::RVec> &positions,
                               std::vector<gmx::RVec>       *forces)
{
    if (pairs_.empty())
    {
        return 0;
    }
    computeInitialDisplacements(positions);
    // The accelerations, the forces divided by the masses
    std::vector<gmx::RVec> acc(forces->size());
    for (size_t i = 0; i < forces->size(); i++)
    {
        svmul(invmass_[i], (*forces)[i], acc[i]);
    }
    Eigen::VectorXd h;
    if (!solveProjection(pairs_, coupled_, initial_displacements_, acc, &h))
    {
        return 1;
    }
    // Add the constraint forces, without dividing by the masses
    std::vector<double> ones(invmass_.size(), 1.0);
    displace(pairs_, ones, initial_displacements_, h, forces);
    return 0;
}

} // namespace alexandria

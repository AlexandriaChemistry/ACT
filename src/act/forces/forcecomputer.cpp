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
#include "forcecomputer.h"

#include <algorithm>
#include <limits>
#include <set>
#include <utility>
#include <vector>

#include <cstdlib>

#include "act/basics/chargemodel.h"
#include "act/basics/msg_handler.h"
#include "act/forcefield/forcefield_parametername.h"
#include "act/forces/conformationdata.h"
#include "act/forces/forcecomputerimpl.h"
#include "act/forces/forces_simd.h"
#include "act/qgen/qtype.h"
#include "act/topology/actmol.h"
#include "act/topology/fragmenthandler.h"
#include "act/topology/topology.h"
#include "act/utility/timings.h"
#include "gromacs/gmxpreprocess/grompp-impl.h"
#include "gromacs/math/vec.h"
#include "gromacs/math/vectypes.h"
#include "gromacs/utility/futil.h"

namespace alexandria
{

/*! \brief Compute the sum of the squared forces on the shells
 * Used to compute the mean square force on shell particles when
 * minimizing the energy with respect to the shell positions.
 * A shell at the maximum distance from its core with a force pointing
 * away from the core is held there, see update_shell_positions. For such
 * a shell only the force perpendicular to the core-shell vector counts,
 * since that is the only force that can still be relaxed.
 * \param[in] atoms             The particles
 * \param[in] isShell           Whether or not the particles are shells
 * \param[in] x                 The coordinates
 * \param[in] f                 The forces
 * \param[in] offset            Index of the first particle in x and f
 * \param[in] maxShellDistance2 Square of the maximum core-shell distance, 0 if there is none
 * \return The sum of the squared shell forces
 */
static double shellForce2(const std::vector<ActAtom>   &atoms,
                          const std::vector<bool>      &isShell,
                          const std::vector<gmx::RVec> &x,
                          const std::vector<gmx::RVec> &f,
                          size_t                        offset,
                          double                        maxShellDistance2)
{
    // Shells are put on the wall by scaling, allow for rounding
    const double atWall = 0.9999*maxShellDistance2;
    double       f2sum  = 0;
    for(size_t i = 0; i < isShell.size(); i++)
    {
        if (!isShell[i])
        {
            continue;
        }
        const auto &fi = f[i+offset];
        double      f2 = iprod(fi, fi);
        if (maxShellDistance2 > 0 && atoms[i].parents().size() == 1)
        {
            rvec dx;
            rvec_sub(x[i+offset], x[atoms[i].parents()[0]+offset], dx);
            double dx2 = iprod(dx, dx);
            double fdx = iprod(fi, dx);
            if (dx2 >= atWall && fdx > 0)
            {
                f2 -= fdx*fdx/dx2;
            }
        }
        f2sum += f2;
    }
    return f2sum;
}

void ForceComputer::init(double  msForce,
                         int     maxiter,
                         double  maxShellDistance,
                         double  epsilonr)
{
    msForceToler_     = msForce;
    maxiter_          = maxiter;
    maxShellDistance_ = maxShellDistance;
    epsilonr_         = epsilonr;
    clear_mat(box_);
    real dt = 0.001;
    vsiteHandler_.init(box_, dt);
}

void ForceComputer::constructVsiteCoordinates(const Topology         *top,
                                              std::vector<gmx::RVec> *coordinates) const
{
    // Construct virtual site coordinates
    vsiteHandler_.constructPositions(top, coordinates, box_);
}

void ForceComputer::spreadVsiteForces(const Topology               *top,
                                      const std::vector<gmx::RVec> *coordinates,
                                      std::vector<gmx::RVec>       *forces) const
{
    // Spread virtual site forces
    vsiteHandler_.distributeForces(top, *coordinates, forces, box_);
}

real calc_fscale(real k1, real k2, real F0)
{
    if (F0 == 0)
    {
        return 0;
    }
    if (k2 == 0)
    {
        // F0 = -k1 Delta x => Delta x = -F0/k1
        // However, F0 is taken out because we multiply by the
        // force vector down below.
        return -1/k1;
    }
    else
    {
        // Solution from Mathematica:
        // V[x_] := (1/2) k1 x^2 + k2 x^4
        // F[x_] := -D[V[x], x]
        // MyX = Assuming[{k1 > 0, k2 > 0}, Solve[F0 == F[x], x]]
        // three solutions are obtained, but only the first one is real.
        // Note that this solution assumes the force is independent of shell
        // position, which will never be the case.
        // Also, this assumes a 1D problem where the force on the shell is
        // parallel to the core-shell distance vector.
        static const real third = 1.0/3.0;
        static const real three = std::pow(3.0, third);
        static const real sqrt3 = std::sqrt(3.0);
        real k12    = k1*k2;
        real F0k22  = F0*k2*k2;
        real throot = std::pow(-9*F0k22 +
                               sqrt3*std::sqrt(k12*k12*k12 + 27*F0k22*F0k22), third);
        auto k = (-three*three*k12 + three*throot*throot)/(6*k2*throot);
        return k/F0;
    }
}

double shellForceConstant(double k1, double k2, double rh, double r)
{
    double fc = k1;
    if (k2 > 0 && r > rh)
    {
        double drh = r - rh;
        fc        += 4*k2*drh*drh*drh/r;
    }
    return fc;
}

static inline void update_shell_pos(int                           shell,
                                    std::vector<gmx::RVec>       *coordinates,
                                    const std::vector<gmx::RVec> &forces,
                                    double                        fcShell,
                                    double                        fcHyper,
                                    double                        step)
{
    auto Fnorm  = norm(forces[shell]);
    auto fscale = step*calc_fscale(fcShell, fcHyper, Fnorm);
    for(int m = 0; m < DIM; m++)
    {
        (*coordinates)[shell][m] -= forces[shell][m] * fscale;
    }
}

/*! \brief Move a shell towards the minimum of its potential
 * \param[in]    core        The core of the shell
 * \param[in]    shell       The shell
 * \param[inout] coordinates The coordinates
 * \param[in]    forces      The forces
 * \param[in]    fcShell     The harmonic force constant of the shell
 * \param[in]    fcHyper     The quartic force constant of the shell
 * \param[in]    rHyper      The distance beyond which the quartic term acts
 * \param[in]    step        Fraction of the displacement that would bring
 *                            an isolated shell to its minimum
 */
static void update_shell_positions(int                           core,
                                   int                           shell,
                                   std::vector<gmx::RVec>       *coordinates,
                                   const std::vector<gmx::RVec> &forces,
                                   double                        fcShell,
                                   double                        fcHyper,
                                   double                        rHyper,
                                   double                        step)
{
    // Core shell vector
    gmx::RVec r_core_shell;
    rvec_sub((*coordinates)[shell], (*coordinates)[core], r_core_shell);
    double rcs    = norm(r_core_shell);
    //double toler  = 1e-2;
    // Check for co-linearity of shell displacement and the force
    //double cosang = cos_angle(r_core_shell, forces[shell]);
    if (rcs == 0)
    {
        // On top of the core the quartic term only contributes if
        // it starts at zero distance.
        update_shell_pos(shell, coordinates, forces, fcShell,
                         rHyper > 0 ? 0 : fcHyper, step);
    }
    else
    {
        // Now we have the following geometry:
        //
        //                 Ftotal
        //                   ^
        //                    ^
        //          Fspring    ^
        //  Core    <-------- Shell ---> r_core_shell
        //
        // We want to determine what distance x along the F_total vector the shell
        // has to be displaced. To do so, we assume that
        // r_core_shell / Fspring = x / Ftotal
        // Ratio between distance and polarization force
        double fac  = step/shellForceConstant(fcShell, fcHyper, rHyper, rcs);
        for(int m = 0; m < DIM; m++)
        {
            (*coordinates)[shell][m] += fac*forces[shell][m];
        }
    }
}

void ForceComputer::compute(MsgHandler                        *msg_handler,
                            const ForceField                  *pd,
                            const Topology                    *top,
                            std::vector<gmx::RVec>            *coordinates,
                            std::vector<gmx::RVec>            *forces,
                            std::map<InteractionType, double> *energies,
                            const gmx::RVec                   &field,
                            bool                               resetShells,
                            std::set<int>                      relax) const
{
    nEval_ += 1;
    constructVsiteCoordinates(top, coordinates);
    // Short-cut
    auto &atoms = top->atoms();
    // Reset shells if needed and if possible
    if (resetShells)
    {
        for(size_t i = 0; i < top->nAtoms(); i++)
        {
            if (atoms[i].pType() == ActParticle::Shell &&
                atoms[i].parents().size() == 1)
            {
                auto parent = atoms[i].parents()[0];
                copy_rvec((*coordinates)[parent], (*coordinates)[i]);
            }
        }
    }
    // Do first calculation every time.
    computeOnce(msg_handler, pd, top, coordinates, forces, energies, field);
    // Store total electrostatics energy in independent copy.
    std::map<InteractionType, double> eBefore = *energies;

    // Now let's have a look whether we are polarizable
    auto itype = InteractionType::POLARIZATION;
    if (pd->polarizable() && top->hasEntry(itype))
    {
        // mean square shell force
        double msForce = 0;
        int    iter    = 1;
        // Is this particle a shell?
        std::vector<bool>   isShell;
        // One over force constant for this particle
        std::vector<double> fcShell_1;
        std::vector<double> fcHyper;
        std::vector<double> rHyper;
        auto &ffpl  = pd->findForcesConst(itype);
        auto pol_name = potentialToParameterName(ffpl.potential());
        int  nshell = 0;
        for(auto &aa : atoms)
        {
            bool bIS = aa.pType() == ActParticle::Shell;
            isShell.push_back(bIS);
            double fc_1    = 0;
            double fchyper = 0;
            double rhyper  = 0;
            if (bIS)
            {
                Identifier atID(aa.ffType());
                auto params = ffpl.findParameterMapConst(atID);
                auto alpha  = params.find(pol_name[polALPHA])->second.internalValue();
                auto ffp = params.find(pol_name[polFCHYPER]);
                if (ffp != params.end())
                {
                    fchyper = ffp->second.internalValue();
                }
                auto ffr = params.find(pol_name[polRHYPER]);
                if (ffr != params.end())
                {
                    rhyper = ffr->second.internalValue();
                }
                auto q     = aa.charge();
                if (alpha > 0 && q != 0)
                {
                    fc_1 = (q*q*ONE_4PI_EPS0)/alpha;
                }
                nshell += 1;
            }
            fcHyper.push_back(fchyper);
            rHyper.push_back(rhyper);
            fcShell_1.push_back(fc_1);
        }
        auto &pols = top->entry(itype);
        // Prevent shells from flying off.
        real maxShellDistance2 = gmx::square(maxShellDistance_);
        msForce    = shellForce2(atoms, isShell, *coordinates, *forces, 0, maxShellDistance2)/nshell;
        // Fraction of the step to the minimum of an isolated shell
        double step = 1;
        while (msForce > msForceToler_ && iter < maxiter_)
        {
            // Loop over polarizabilities
            for(const auto &p : pols)
            {
                // Displace the shells according to the force
                // Since the potential is harmonic we use Hooke's law
                // F = k dx -> dx = F / k
                // Note that this only works for purely harmonic potentials,
                // see function calc_fscale for quartic potentials.
                int shell = p->atomIndex(1);
                if (relax.empty() || relax.end() != relax.find(shell))
                {
                    update_shell_positions(p->atomIndex(0), p->atomIndex(1),
                                           coordinates, *forces,
                                           fcShell_1[shell], fcHyper[shell],
                                           rHyper[shell], step);
                    // Check distance from parent, if there is only one
                    if (maxShellDistance2 > 0 && atoms[shell].parents().size() == 1)
                    {
                        int core = atoms[shell].parents()[0];
                        rvec dx;
                        rvec_sub((*coordinates)[shell], (*coordinates)[core], dx);
                        real dx2 = iprod(dx, dx);
                        if (dx2 > maxShellDistance2)
                        {
                            // Move back the shell/drude to the wall distance
                            real scale = std::sqrt(maxShellDistance2/dx2);
                            for (int m = 0; m < DIM; m++)
                            {
                                (*coordinates)[shell][m] = (*coordinates)[core][m] + scale*dx[m];
                            } 
                        }
                    }
                }
            }
            // Do next calculation
            computeOnce(msg_handler, pd, top, coordinates, forces, energies, field);
            double msForceNew = shellForce2(atoms, isShell, *coordinates, *forces, 0, maxShellDistance2)/nshell;
            // Coupled shells can overshoot, then the step is reduced.
            // Smaller steps can freeze shells far from the minimum.
            if (msForceNew > msForce)
            {
                step = 0.5;
            }
            msForce  = msForceNew;
            iter    += 1;
            nShellIter_ += 1;
        }
        if (msForce > msForceToler_)
        {
            nShellConvergeFailed_ += 1;
            if (msg_handler && msg_handler->debug())
            {
                msg_handler->msg(ACTStatus::Debug,
                                 gmx::formatString("Shell optimization did not converge. RMS force is %g, iter %d/%d, total evals %zu, total shell iters %zu", std::sqrt(msForce), iter, maxiter_, nEval_, nShellIter_));
            }
        }
    }
    {
        // Induction energy
        bool   haveInduction = false;
        double eInduction    = 0;
        // Extract electrostatics once more
        // Note that the INDUCTIONCORRECTION is treated in the calling routine
        std::set<InteractionType> eTerms = {
            InteractionType::ELECTROSTATICS,
            InteractionType::POLARIZATION,
            InteractionType::CHARGETRANSFER
        };
        for(const auto et : eTerms)
        {
            auto tt = energies->find(et);
            if (energies->end() != tt &&
                eBefore.end() != eBefore.find(et))
            {
                eInduction   += tt->second - eBefore[et];
                tt->second    = eBefore[et];
                haveInduction = true;
            }
        }
        if (haveInduction)
        {
            energies->insert_or_assign(InteractionType::INDUCTION, eInduction);
        }
    }
    {
        // Sum of all electrostatic terms
        double allelec = 0;
        for(const auto &itype : { InteractionType::ELECTROSTATICS,
                                  InteractionType::POLARIZATION,
                                  InteractionType::INDUCTION,
                                  InteractionType::INDUCTIONCORRECTION })
        {
            auto ee = energies->find(itype);
            if (energies->end() != ee)
            {
                allelec += ee->second;
            }
        }
        energies->insert_or_assign(InteractionType::ALLELEC, allelec);
    }
    {
        // Sum of exchange and induction terms
        double exchind = 0;
        for(const auto &itype : { InteractionType::EXCHANGE, InteractionType::INDUCTION, InteractionType::INDUCTIONCORRECTION })
        {
            auto ee = energies->find(itype);
            if (energies->end() != ee)
            {
                exchind += ee->second;
            }
        }
        energies->insert_or_assign(InteractionType::EXCHIND, exchind);
    }
    // Spread forces to atoms
    spreadVsiteForces(top, coordinates, forces);
}

/*! \brief The shells with one parent particle, with that parent
 * \param[in] atoms The particles
 * \return The shells and their parents
 */
static std::vector<std::pair<size_t, size_t>> shellParents(const std::vector<ActAtom> &atoms)
{
    std::vector<std::pair<size_t, size_t>> shells;
    for(size_t i = 0; i < atoms.size(); i++)
    {
        if (atoms[i].pType() == ActParticle::Shell &&
            atoms[i].parents().size() == 1)
        {
            shells.push_back({ i, atoms[i].parents()[0] });
        }
    }
    return shells;
}

/*! \brief Put shells on top of their parent particle
 * \param[in]    shells The shells and their parents, see shellParents
 * \param[inout] coords The coordinates
 */
static void resetShells(const std::vector<std::pair<size_t, size_t>> &shells,
                        std::vector<gmx::RVec>                       *coords)
{
    for(const auto &sp : shells)
    {
        copy_rvec((*coords)[sp.second], (*coords)[sp.first]);
    }
}

size_t ForceComputer::updateShells(const Topology            *top,
                                   size_t                     offset,
                                   size_t                     nShell,
                                   const std::vector<bool>    isShell,
                                   size_t                     fragment,
                                   std::vector<CoordsEnergy> *myCoords) const
{
    auto   itype      = InteractionType::POLARIZATION;
    size_t numberDone = 0;
    double maxShellDistance2 = maxShellDistance_*maxShellDistance_;
    for(auto &dim : *myCoords)
    {
        if (dim.done)
        {
            numberDone += 1;
            continue;
        }
        // Mean square shell force
        double msForce = shellForce2(top->atoms(), isShell, dim.x, dim.fcalc,
                                     offset, maxShellDistance2);
        if (nShell > 0)
        {
            msForce /= nShell;
        }
        if (msForce <= msForceToler_)
        {
            dim.done    = true;
            numberDone += 1;
        }
        else
        {
            // Coupled shells can overshoot, then the step is reduced, see compute
            if (msForce > dim.shellForce[fragment])
            {
                dim.shellStep[fragment] = 0.5;
            }
            dim.shellForce[fragment] = msForce;
            auto &pols  = top->entry(itype);
            auto &atoms = top->atoms();
            // Loop over polarizabilities
            for(const auto &p : pols)
            {
                // Displace the shells according to the force
                // Since the potential is harmonic we use Hooke's law
                // F = k dx -> dx = F / k
                // Note that this only works for purely harmonic potentials,
                // see function calc_fscale for quartic potentials.
                // The atom indices stored in the fragment topology are local,
                // whereas the coordinates/forces of the dimer are global, so
                // we add the offset to map local -> global indices.
                int shellLocal   = p->atomIndex(1);
                int coreLocal    = p->atomIndex(0);
                int shellGlobal  = shellLocal + static_cast<int>(offset);
                int coreGlobal   = coreLocal + static_cast<int>(offset);
                // The force constant from the current polarizability and
                // shell charge, since both may be trained, as in compute.
                auto   &params = p->params();
                double  alpha  = params[polALPHA];
                double  qshell = atoms[shellLocal].charge();
                double  fc1    = 0;
                if (alpha > 0 && qshell != 0)
                {
                    fc1 = (qshell*qshell*ONE_4PI_EPS0)/alpha;
                }
                update_shell_positions(coreGlobal, shellGlobal,
                                       &dim.x, dim.fcalc,
                                       fc1, params[polFCHYPER], params[polRHYPER],
                                       dim.shellStep[fragment]);
                // Check distance from parent, if there is only one
                if (maxShellDistance2 > 0 && atoms[shellLocal].parents().size() == 1)
                {
                    int coreParentLocal  = atoms[shellLocal].parents()[0];
                    int coreParentGlobal = coreParentLocal + static_cast<int>(offset);
                    rvec dx;
                    rvec_sub(dim.x[shellGlobal], dim.x[coreParentGlobal], dx);
                    real dx2 = iprod(dx, dx);
                    if (dx2 > maxShellDistance2)
                    {
                        // Move back the shell/drude to the wall distance
                        real scale = std::sqrt(maxShellDistance2/dx2);
                        for (int m = 0; m < DIM; m++)
                        {
                            dim.x[shellGlobal][m] = dim.x[coreParentGlobal][m] + scale*dx[m];
                        }
                    }
                }
            }
        }
    }
    return numberDone;
}

//! \return Whether an interaction type changes upon shell relaxation
static inline bool isElecTerm(InteractionType itype)
{
    return (InteractionType::ELECTROSTATICS == itype ||
            InteractionType::POLARIZATION   == itype ||
            InteractionType::CHARGETRANSFER == itype);
}

/*! \brief Compute the induction, all electrostatics and exchange plus induction
 * The induction is the change of the electrostatic terms since their backup,
 * and these terms are restored to their backup.
 * \param[inout] cener The energies of all conformations
 */
static void calcInduction(ConformationEnergies *cener)
{
    // The terms present, in the same order as in an energy map.
    // Note that the INDUCTIONCORRECTION is treated in the calling routine.
    std::vector<std::pair<InteractionType, ConformationEnergies::Term>> terms;
    for(auto t : cener->types())
    {
        if (isElecTerm(t) || InteractionType::EXCHANGE == t)
        {
            terms.push_back({ t, cener->find(t) });
        }
    }
    auto induction = cener->find(InteractionType::INDUCTION);
    auto eAllElec  = cener->add(InteractionType::ALLELEC);
    auto eExchInd  = cener->add(InteractionType::EXCHIND);
    for(size_t c = 0; c < cener->nconf(); c++)
    {
        bool   haveInduction = false;
        double eInduction    = 0;
        double allElec       = 0;
        double exchInd       = 0;
        for(auto &tt : terms)
        {
            auto &et = tt.second;
            if (!et.has(c))
            {
                continue;
            }
            if (isElecTerm(tt.first))
            {
                // Without backups nothing was relaxed
                if (cener->haveBackups() && et.before[c] != et.act[c])
                {
                    eInduction += et.act[c] - et.before[c];
                    et.restoreBackup(c);
                    haveInduction = true;
                }
                // Copy new value
                allElec += et.act[c];
            }
            else
            {
                exchInd += et.act[c];
            }
        }
        if (haveInduction)
        {
            // This is the total induction including the induction correction,
            // which may be split off later. The induction and induction
            // correction are not up to date yet.
            if (!induction.have)
            {
                induction = cener->add(InteractionType::INDUCTION);
            }
            induction.set(c, eInduction);
            allElec += eInduction;
            exchInd += eInduction;
        }
        else if (induction.has(c))
        {
            // Induction from a potential function rather than from shells
            exchInd += induction.act[c];
        }
        eAllElec.set(c, allElec);
        eExchInd.set(c, exchInd);
    }
}

/*! \brief Sum the ACT energy terms that change upon shell relaxation
 * \param[inout] cener The energies of all conformations
 * \return The sum of electrostatics, polarization and charge transfer for each conformation
 */
static std::vector<double> sumElecTerms(ConformationEnergies *cener)
{
    std::vector<double> esum(cener->nconf(), 0.0);
    for(const auto itype : { InteractionType::ELECTROSTATICS,
                             InteractionType::POLARIZATION,
                             InteractionType::CHARGETRANSFER })
    {
        auto ee = cener->find(itype);
        for(size_t c = 0; c < cener->nconf(); c++)
        {
            if (ee.has(c))
            {
                esum[c] += ee.act[c];
            }
        }
    }
    return esum;
}

/*! \brief Back up all ACT energies
 * \param[inout] cener The energies of all conformations
 */
static void backupEnergies(ConformationEnergies *cener)
{
    // The backups are only used for conformations with an energy, so
    // copying all of them is the same as copying only those.
    for(auto t : cener->types())
    {
        auto et = cener->find(t);
        std::copy(et.act, et.act + cener->nconf(), et.before);
    }
}

/*! \brief Subtract the energies of the monomers from those of the dimer
 * The backups are subtracted as well.
 * \param[inout] dimer    The energies of the dimer conformations
 * \param[inout] monomers The energies of the monomer conformations
 */
static void subtractMonomers(ConformationEnergies *dimer, ConformationEnergies *monomers[2])
{
    // Without backups, they are the same as the energies
    bool backups = dimer->haveBackups();
    for(auto t : dimer->types())
    {
        auto ed = dimer->find(t);
        ConformationEnergies::Term em[2] = { monomers[0]->find(t), monomers[1]->find(t) };
        for(size_t c = 0; c < dimer->nconf(); c++)
        {
            if (!ed.has(c))
            {
                continue;
            }
            double edim = ed.act[c];
            for(size_t j = 0; j < 2; j++)
            {
                if (c < monomers[j]->nconf() && em[j].has(c))
                {
                    edim -= em[j].act[c];
                }
            }
            ed.act[c] = edim;
            if (backups)
            {
                double edb4 = ed.before[c];
                for(size_t j = 0; j < 2; j++)
                {
                    if (c < monomers[j]->nconf() && em[j].has(c))
                    {
                        edb4 -= em[j].before[c];
                    }
                }
                ed.before[c] = edb4;
            }
        }
    }
}

/*! \brief Check whether the virtual sites have to be constructed
 * The atoms do not move, therefore the virtual sites only have to be
 * constructed the first time and when their parameters changed, unless
 * shells are used to construct them.
 * \param[in]    top   The topology
 * \param[inout] cdata The conformation data, updated if the virtual sites
 *                     have to be constructed
 * \return Whether the virtual sites have to be constructed
 */
static bool vsitesOutdated(const Topology *top, ConformationData *cdata)
{
    bool outdated = (!cdata->haveVsites || cdata->vsiteVersion != top->vsiteVersion());
    if (!outdated)
    {
        auto &atoms = top->atoms();
        for(const auto &entry : top->entries())
        {
            if (!isVsite(entry.first))
            {
                continue;
            }
            for(const auto &vs : entry.second)
            {
                for(auto a : vs->atomIndices())
                {
                    outdated = outdated || atoms[a].pType() == ActParticle::Shell;
                }
            }
        }
    }
    if (outdated)
    {
        cdata->haveVsites   = true;
        cdata->vsiteVersion = top->vsiteVersion();
        cdata->arrays.invalidateVsites();
        cdata->coordinates.update();
    }
    return outdated;
}

/*! \brief Copy the relaxed monomer shell positions into the dimer
 * \param[inout] allmols The molecules
 */
static void copyRelaxedMonomerShells(const std::vector<ACTMol *> &allmols)
{
    for(auto molptr : allmols)
    {
        auto &mol = *molptr;
        if (!mol.isDimer())
        {
            continue;
        }
        auto  fh     = mol.fragmentHandler();
        auto &tops   = fh->topologies();
        auto &astart = fh->atomStart();
        auto &dimers = *mol.dimerCoords();
        for(size_t im = 0; im < tops.size(); im++)
        {
            auto &atoms = tops[im].atoms();
            auto &mons  = *mol.monomerCoords(im);
            for(size_t i = 0; i < dimers.size(); i++)
            {
                for(size_t j = 0; j < atoms.size(); j++)
                {
                    if (atoms[j].pType() == ActParticle::Shell)
                    {
                        copy_rvec(mons[i].xrelax[j], dimers[i].x[astart[im]+j]);
                    }
                }
            }
        }
    }
}

/*! \brief Start the shell relaxation with full steps
 * \param[inout] coords The conformations
 */
static void resetShellSteps(std::vector<CoordsEnergy> *coords)
{
    for(auto &ce : *coords)
    {
        for(size_t i = 0; i < 2; i++)
        {
            ce.shellStep[i]  = 1;
            ce.shellForce[i] = std::numeric_limits<double>::max();
        }
    }
}

void ForceComputer::relaxMonomer(MsgHandler                  *msghandler,
                                 const ForceField            *pd,
                                 size_t                       im,
                                 const std::vector<ACTMol *> &allmols) const
{
    static const size_t timer = Timings::instance().index("Shell relaxation");
    ScopedTimer         scoped(timer);
    for(auto molptr : allmols)
    {
        resetShellSteps(molptr->monomerCoords(im));
    }
    // We have already done one energy calculation
    int iter     = 1;
    // Prevent shells from flying off.
    bool allDone = false;
    while (iter < maxiter_ && !allDone)
    {
        // Set loop variable to true
        allDone = true;
        for(auto molptr : allmols)
        {
            auto &mol = *molptr;
            // Shortcuts...
            auto  fh    = mol.fragmentHandler();
            auto &tops  = fh->topologies();
            auto &pols  = fh->polInfo();
            // Relax monomers first.
            if (mol.monomerCoords(im)->size() !=
                updateShells(&tops[im], 0, pols[im].nShell(), pols[im].isShell(),
                             0, mol.monomerCoords(im)))
            {
                allDone = false;
                // Not all done, recompute
            }
        }
        if (!allDone)
        {
            // Only for the conformations that did not converge
            computeEnergies(msghandler, pd, allmols, im, true, true);
        }
        iter++;
    }
}

void ForceComputer::relaxDimer(MsgHandler                  *msghandler,
                               const ForceField            *pd,
                               bool                         fixMonomer[2],
                               const std::vector<ACTMol *> &allmols) const
{
    static const size_t timer = Timings::instance().index("Shell relaxation");
    ScopedTimer         scoped(timer);
    for(auto molptr : allmols)
    {
        resetShellSteps(molptr->dimerCoords());
    }
    // We have already done one energy calculation
    int iter     = 1;
    // Prevent shells from flying off.
    bool allDone = false;
    while (iter < maxiter_ && !allDone)
    {
        // Set loop variable to true
        allDone = true;
        for(auto molptr : allmols)
        {
            auto &mol = *molptr;
            // Shortcuts...
            auto  fh    = mol.fragmentHandler();
            auto &tops  = fh->topologies();
            auto &pols  = fh->polInfo();
            // The done flag is per dimer conformation, but convergence
            // has to be checked for each monomer separately. A conformation
            // is done when the shells of all monomers that are relaxed converged.
            auto              &dims = *mol.dimerCoords();
            std::vector<char>  wasDone(dims.size()), done(dims.size());
            for(size_t c = 0; c < dims.size(); c++)
            {
                wasDone[c] = dims[c].done;
                done[c]    = 1;
            }
            for(size_t im = 0; im < 2; im++)
            {
                if (!fixMonomer[im])
                {
                    for(size_t c = 0; c < dims.size(); c++)
                    {
                        dims[c].done = wasDone[c];
                    }
                    updateShells(&tops[im], fh->atomStart()[im],
                                 pols[im].nShell(), pols[im].isShell(),
                                 im, &dims);
                    for(size_t c = 0; c < dims.size(); c++)
                    {
                        done[c] = done[c] && dims[c].done;
                    }
                }
            }
            for(size_t c = 0; c < dims.size(); c++)
            {
                dims[c].done = done[c];
                allDone      = allDone && done[c];
            }
        }
        if (!allDone)
        {
            // Not all done, recompute the conformations that did not converge
            computeEnergies(msghandler, pd, allmols, 2, true, true);
        }
        iter++;
    }
}

void ForceComputer::computeRelaxedEnergies(MsgHandler                  *msghandler,
                                           const ForceField            *pd,
                                           const std::vector<ACTMol *> &allmols,
                                           bool                         resetTheShells,
                                           bool                         separateInductionCorrection,
                                           bool                         computeForces) const
{
    static const size_t timer = Timings::instance().index("Energy calculation overhead");
    ScopedTimer         scoped(timer);
    // Shells are relaxed using the forces
    bool forces = computeForces || pd->polarizable();
    // Without polarization, the interaction energies of dimers are those of
    // the intermolecular interactions, the monomers are not needed then.
    // The forces would lack the intramolecular contributions.
    bool intermolecular = !forces;
    // Loop over molecules
    for(auto molptr : allmols)
    {
        auto &mol = *molptr;
        bool  inter = intermolecular && mol.isDimer() && mol.topology()->hasIntermolecular();
        if (mol.isDimer())
        {
            // We have a dimer!
            auto top    = mol.topology();
            bool vsites = vsitesOutdated(top, mol.conformationData(2));
            // Only the shells, if any, are reset
            std::vector<std::pair<size_t, size_t>> shells;
            if (resetTheShells)
            {
                shells = shellParents(top->atoms());
            }
            // Virtual sites for dimer
            for(auto &dim : *mol.dimerCoords())
            {
                nEval_ += 1;
                if (vsites)
                {
                    constructVsiteCoordinates(top, &dim.x);
                }
                // Reset shells if needed and if possible
                resetShells(shells, &dim.x);
            }
        }
        // Shortcuts...
        auto &tops  = mol.fragmentHandler()->topologies();
        // Virtual sites for monomers
        for(size_t index = 0; index < tops.size() && !inter; index++)
        {
            auto &top    = tops[index];
            bool  vsites = vsitesOutdated(&top, mol.conformationData(index));
            std::vector<std::pair<size_t, size_t>> shells;
            if (resetTheShells)
            {
                shells = shellParents(top.atoms());
            }
            for(auto &mon : *mol.monomerCoords(index))
            {
                nEval_ += 1;
                if (vsites)
                {
                    constructVsiteCoordinates(&top, &mon.x);
                }
                // Reset shells if needed and if possible
                resetShells(shells, &mon.x);
            }
        }
    }
    // Do first calculation every time.
    computeEnergies(msghandler, pd, allmols, 0, forces, false, intermolecular); // Monomer 0
    computeEnergies(msghandler, pd, allmols, 1, forces, false, intermolecular); // Monomer 1
    computeEnergies(msghandler, pd, allmols, 2, forces, false, intermolecular); // Complete topology

    // Store the original (unrelaxed) energies so that we can extract the
    // induction later. Only terms that actually have an ACT energy are backed up.
    // Only if shells are relaxed the energies change.
    bool relax = pd->polarizable();
    for(auto molptr : allmols)
    {
        auto &mol = *molptr;
        for(size_t what = 0; what < 3; what++)
        {
            auto &energies = mol.conformationData(what)->energies;
            energies.setBackups(relax);
            if (relax)
            {
                backupEnergies(&energies);
            }
        }
    }
    // Now let's have a look whether we are polarizable
    if (pd->polarizable())
    {
        for(size_t im = 0; im < 2; im++)
        {
            relaxMonomer(msghandler, pd, im, allmols);
        }
        // Backup the now-relaxed monomer coordinates so that we can return to
        // them later when computing the second order induction.
        for(auto molptr : allmols)
        {
            auto &mol = *molptr;
            // Shortcuts...
            auto &tops  = mol.fragmentHandler()->topologies();
            for(size_t index = 0; index < tops.size(); index++)
            {
                for(auto &mon : *mol.monomerCoords(index))
                {
                    mon.xrelax = mon.x;
                }
            }
        }
    }
    // Second order induction energy for each dimer conformation of each molecule
    std::vector<std::vector<double>> eInduc2(allmols.size());
    bool doInduc2 = pd->polarizable() && separateInductionCorrection;
    // From now on, only polarizable models relax the dimer further, but the
    // post-processing below runs in all cases.
    if (pd->polarizable())
    {
        if (doInduc2)
        {
            // Store the original dimer coordinates, such that the final
            // relaxation below does not depend on what we do here.
            std::vector<std::vector<std::vector<gmx::RVec>>> xdimer(allmols.size());
            for(size_t k = 0; k < allmols.size(); k++)
            {
                auto &mol = (*allmols[k]);
                if (mol.isDimer())
                {
                    for(auto &dim : *mol.dimerCoords())
                    {
                        xdimer[k].push_back(dim.x);
                    }
                    eInduc2[k].resize(mol.dimerCoords()->size(), 0.0);
                }
            }
            // We have to compute second order induction first, starting
            // from the shells in the relaxed monomer positions.
            copyRelaxedMonomerShells(allmols);
            computeEnergies(msghandler, pd, allmols, 2);
            // Reference energy without mutual induction
            std::vector<std::vector<double>> eRef(allmols.size());
            for(size_t k = 0; k < allmols.size(); k++)
            {
                eRef[k] = sumElecTerms(&(*allmols[k]).conformationData(2)->energies);
            }
            for(size_t im = 0; im < 2; im++)
            {
                if (im > 0)
                {
                    // Reset coordinates to the relaxed monomer shells
                    copyRelaxedMonomerShells(allmols);
                    computeEnergies(msghandler, pd, allmols, 2);
                }
                // Relax monomer im in the presence of the other monomer
                bool fixMonomer[2] = { im != 0, im != 1 };
                relaxDimer(msghandler, pd, fixMonomer, allmols);
                // Extract energy
                for(size_t k = 0; k < allmols.size(); k++)
                {
                    auto esum = sumElecTerms(&(*allmols[k]).conformationData(2)->energies);
                    for(size_t i = 0; i < eInduc2[k].size() && i < esum.size(); i++)
                    {
                        eInduc2[k][i] += esum[i] - eRef[k][i];
                    }
                }
            }
            // Restore the original dimer coordinates and energies
            for(size_t k = 0; k < allmols.size(); k++)
            {
                auto &dimers = *(*allmols[k]).dimerCoords();
                for(size_t i = 0; i < dimers.size(); i++)
                {
                    dimers[i].x = xdimer[k][i];
                }
            }
            computeEnergies(msghandler, pd, allmols, 2);
        }
        // Final, total energy
        bool ff[2] = { false, false };
        relaxDimer(msghandler, pd, ff, allmols);
    }
    
    // Now post-process the energies to extract induction and allelec
    for(const auto &molptr : allmols)
    {
        auto &mol = *molptr;
        auto *dimer = &mol.conformationData(2)->energies;
        ConformationEnergies *monomers[2] = { &mol.conformationData(0)->energies,
                                              &mol.conformationData(1)->energies };
        bool inter = intermolecular && mol.isDimer() && mol.topology()->hasIntermolecular();
        if (mol.isDimer())
        {
            // First, subtract monomer energies from dimer energies,
            // unless only the intermolecular interactions were computed
            if (!inter)
            {
                subtractMonomers(dimer, monomers);
            }
            calcInduction(dimer);
            if (doInduc2)
            {
                // Split the total induction into second order induction and
                // the remaining higher order induction correction.
                size_t k     = &molptr - allmols.data();
                auto   eit   = dimer->find(InteractionType::INDUCTION);
                std::vector<double> eInd(dimer->nconf(), 0.0);
                for(size_t i = 0; i < dimer->nconf(); i++)
                {
                    if (eit.has(i))
                    {
                        eInd[i] = eit.act[i];
                    }
                }
                auto induction  = dimer->add(InteractionType::INDUCTION);
                auto correction = dimer->add(InteractionType::INDUCTIONCORRECTION);
                for(size_t i = 0; i < dimer->nconf() && i < eInduc2[k].size(); i++)
                {
                    induction.set(i, eInduc2[k][i]);
                    correction.set(i, eInd[i] - eInduc2[k][i]);
                }
            }
        }
        auto &tops = mol.fragmentHandler()->topologies();
        for(size_t i = 0; i < tops.size() && !inter; i++)
        {
            calcInduction(monomers[i]);
        }
    }
    if (!forces)
    {
        return;
    }
    // Spread forces to atoms
    // Loop over molecules
    for(auto molptr : allmols)
    {
        auto &mol = *molptr;
        // Shortcuts...
        auto &tops  = mol.fragmentHandler()->topologies();
        if (tops.size() == 2)
        {
            // We have a dimer!
            auto top    = mol.topology();
            // Virtual sites for dimer
            for(auto &dim : *mol.dimerCoords())
            {
                spreadVsiteForces(top, &dim.x, &dim.fcalc);
            }
        }
        // Virtual sites for monomers
        for(size_t index = 0; index < tops.size(); index++)
        {
            auto &top   = tops[index];
            for(auto &mon : *mol.monomerCoords(index))
            {
                spreadVsiteForces(&top, &mon.x, &mon.fcalc);
            }
        }
    }
}

void ForceComputer::computeOnce(MsgHandler                        *msg_handler,
                                const ForceField                  *pd,
                                const Topology                    *top,
                                std::vector<gmx::RVec>            *coordinates,
                                std::vector<gmx::RVec>            *forces,
                                std::map<InteractionType, double> *energies,
                                const gmx::RVec                   &field) const
{
    // Clear energies
    energies->clear();
    // Clear forces
    auto &atoms = top->atoms();
    for(size_t ff = 0; ff < forces->size(); ++ff)
    {
        real fac = FIELDFAC*atoms[ff].charge();
        svmul(fac, field, (*forces)[ff]);
    }
    double epot = 0;
    for(const auto &entry : top->entries())
    {
        if (entry.second.empty())
        {
            continue;
        }
        // Not all potentials are available in the force field,
        // in particular restraint potentials.
        Potential mypot;
        if (pd->interactionPresent(entry.first))
        {
            // Force field parameter list
            auto &ffpl = pd->findForcesConst(entry.first);
            // The function we need to do the math
            mypot = ffpl.potential();
        }
        else
        {
            mypot = defaultPotential(entry.first);
        }
        auto bfc   = getBondForceComputer(mypot);
        if (bfc)
        {
            // Now do the calculations and store the energy
            std::map<InteractionType, double> my_energy;
            auto ener = bfc(msg_handler, entry.second, atoms, coordinates, forces, &my_energy, epsilonr_);
            if (my_energy.size() > 1)
            {
                for(const auto &me : my_energy)
                {
                    if (energies->find(me.first) != energies->end())
                    {
                        GMX_THROW(gmx::InternalError(gmx::formatString("Energy term %s occurs twice",
                                                                       interactionTypeToString(me.first).c_str()).c_str()));
                    }
                    energies->insert_or_assign( me.first, me.second );
                    epot += me.second;
                }
            }
            else
            {
                energies->insert_or_assign( entry.first, ener );

                epot += ener;
            }
        }
        else if (debug && !isVsite(entry.first))
        {
            fprintf(debug, "Please implement a force function for type %s\n",
                    potentialToString(mypot).c_str());
        }
    }
    auto ivdwcorr = energies->find(InteractionType::VDWCORRECTION);
    if (energies->end() != ivdwcorr)
    {
        energies->find(InteractionType::EXCHANGE)->second += ivdwcorr->second;
        ivdwcorr->second = 0;
    }
    energies->insert_or_assign( InteractionType::EPOT, epot );
}

Potential ForceComputer::ftype(const ForceField *pd,
                               InteractionType   itype) const
{
    Potential ftype = Potential::NONE;
    if (pd->interactionPresent(itype))
    {
        ftype = pd->findForcesConst(itype).potential();
    }
    return ftype;
}

void ForceComputer::plot(MsgHandler        *msghandler,
                         const ForceField  *pd,
                         InteractionType    itype) const
{
    if (!pd->interactionPresent(itype))
    {
        fprintf(stderr, "No such interaction %s in the force field.\n",
                interactionTypeToString(itype).c_str());
        return;
    }
    std::string btype("bondtype");

    std::map<InteractionType, const std::string> i2s = {
        { InteractionType::BONDS,              btype },
        { InteractionType::ANGLES,             btype },
        { InteractionType::LINEAR_ANGLES,      btype },
        { InteractionType::PROPER_DIHEDRALS,   btype },
        { InteractionType::IMPROPER_DIHEDRALS, btype },
        { InteractionType::VDW,                "vdwtype" },
        { InteractionType::ELECTROSTATICS,     "acmtype" }
    };
    auto &fs   = pd->findForcesConst(itype);
    // The function we need to do the math
    auto bfc   = getBondForceComputer(fs.potential());
    if (nullptr == bfc)
    {
        fprintf(stderr, "Please implement a force function for type %s\n",
                potentialToString(fs.potential()).c_str());
    }
    else
    {
        std::vector<std::pair<int, int>> linear_bonds = {
            { 0, 1}, {1, 2}, {2, 3}
        };
        std::vector<std::pair<int, int>> idih_bonds = {
            { 0, 1}, {0, 2}, {0, 3}
        };
        for(const auto &f : fs.parametersConst())
        {
            unsigned int ntrain  = 0;
            for(const auto &pp : f.second)
            {
                ntrain = std::max(ntrain, pp.second.ntrain());
            }
            if (ntrain == 0)
            {
                continue;
            }
            auto bos   = f.first.bondOrders();
            std::vector<Bond> bbb;
            for(size_t i = 0; i < bos.size(); i++)
            {
                if (InteractionType::IMPROPER_DIHEDRALS == itype)
                {
                    bbb.push_back(Bond(idih_bonds[i].first, idih_bonds[i].second, bos[i]));
                }
                else if (InteractionType::VDW == itype ||
                         InteractionType::ELECTROSTATICS == itype)
                {
                    ;
                }
                else
                {
                    bbb.push_back(Bond(linear_bonds[i].first, linear_bonds[i].second, bos[i]));
                }
            }
            Topology top;
            top.init(bbb);

            auto subtype = i2s.find(itype);
            if (i2s.end() != subtype)
            {
                auto ptypes = pd->particleTypesConst();
                // Loop over subtypes (depending on the InteractionType)
                for(const auto &atomsub : f.first.atoms())
                {
                    for (const auto &pt : ptypes)
                    {
                        if (pt.second.hasOption(subtype->second) &&
                            pt.second.optionValue(subtype->second) == atomsub)
                        {
                            ActAtom a(pt.second);
                            a.setCharge(1);
                            top.addAtom(a);
                            break;
                        }
                    }
                }
            }
            if (top.nAtoms() < 2)
            {
                continue;
            }
            // Open outfile
            std::string filename = gmx::formatString("%s_%s.xvg",
                                                     interactionTypeToString(itype).c_str(),
                                                     f.first.id().c_str());
            std::vector<gmx::RVec> forces;
            gmx::RVec rvnul = { 0, 0, 0 };
            FILE *fp = gmx_ffopen(filename.c_str(), "w");
            std::map<InteractionType, double> energies;
            switch (itype)
            {
            case InteractionType::BONDS:
            case InteractionType::VDW:
            case InteractionType::ELECTROSTATICS:
                {
                    std::vector<gmx::RVec> coordinates = { { 0, 0, 0 }, { 1, 0, 0 } };
                    // top.build will resize coordinates array.
                    top.build(msghandler,
                              pd, &coordinates, 175.0, 5.0, missingParameters::Error);
                    forces.resize(top.nAtoms(), rvnul);
                    // First atom is zero, second must be the other particle
                    size_t jatom = 1+top.atoms()[0].shells().size()+top.atoms()[0].vsites().size();
                    if (jatom == 0)
                    {
                        GMX_THROW(gmx::InternalError(gmx::formatString("Could not find a second atom to make a plot, there are %zu atoms, interactionType %s, id %s", top.nAtoms(), interactionTypeToString(itype).c_str(), f.first.id().c_str()).c_str()));
                    }
                    {
                        // Now remove everything except what we want to compute.
                        std::map<InteractionType, TopologyEntryVector> *entries = top.entriesPtr();
                        for (const auto &ik : i2s)
                        {
                            if (itype != ik.first)
                            {
                                entries->erase(ik.first);
                            }
                        }
                    }
                    std::vector<double> rr, vv, ff;
                    // Now do the calculations and store the energy
                    double r0 = 0.05, r1 = 1.0, delta = 0.001;
                    int    nsteps = (r1-r0)/delta+1;
                    for(int i = 0; i < nsteps; i++)
                    {
                        double x = r0+i*delta;
                        coordinates[jatom][0] = x;
                        // Move shell position
                        for(int j : top.atoms()[jatom].shells())
                        {
                            coordinates[j][0] = x;
                        }
                        // Move vsite position
                        for(int j : top.atoms()[jatom].vsites())
                        {
                            coordinates[j][0] = x;
                        }
                        rr.push_back(x);
                        energies.clear();
                        for(size_t k = 0; k < top.nAtoms(); k++)
                        {
                            copy_rvec(rvnul, forces[k]);
                        }
                        // We need to compute the sum of all relevant forces and run shell minimization.
                        compute(msghandler, pd, &top, &coordinates, &forces, &energies);

                        auto ener = energies[InteractionType::EPOT];

                        if (ener == 0 && InteractionType::VDW == itype)
                        {
                            auto irep = InteractionType::EXCHANGE;
                            auto idsp = InteractionType::DISPERSION;
                            if (energies.find(irep) != energies.end() &&
                                energies.find(idsp) != energies.end())
                            {
                                ener = energies[irep] + energies[idsp];
                            }
                            auto iec = InteractionType::VDWCORRECTION;
                            if (energies.find(iec) != energies.end())
                            {
                                ener += energies[iec];
                            }
                        }
                        fprintf(fp, "%10g  %10g\n", x, ener);
                        vv.push_back(ener);
                        ff.push_back(forces[0][0]);
                    }
                    // Check whether force is derivative of energy
                    for(size_t i = 1; i < vv.size()-1; i++)
                    {
                        if (std::abs(ff[i]) > 1e-6)
                        {
                            double fnumeric = (vv[i+1]-vv[i-1])/(2*delta);
                            double relerror = (fnumeric-ff[i])/ff[i];
                            if (std::abs(relerror) > 1e-1)
                            {
                                printf("%s: Force %g, expected %g. Relative error %g, r = %g v+ %g v- %g delta %g\n",
                                       filename.c_str(),
                                       ff[i], fnumeric, relerror, rr[i],
                                       vv[i+1], vv[i-1], 2*delta);
                            }
                        }
                    }
                }
                break;
            case InteractionType::ANGLES:
                {
                    std::vector<gmx::RVec> acoords     = { { 0, 0, 0 }, { 1, 0, 0 }, { 1, 1, 0 } };
                    std::vector<gmx::RVec> coordinates = acoords;
                    // top.build will resize coordinates array.
                    top.build(msghandler,
                              pd, &coordinates, 175.0, 5.0, missingParameters::Error);
                    forces.resize(top.nAtoms(), rvnul);
                    if (!top.hasEntry(itype))
                    {
                        break;
                    }
                    size_t ak = 0;
                    for(const auto &entry : top.entry(itype))
                    {
                        const auto &ind = entry->atomIndices();
                        GMX_RELEASE_ASSERT(ind.size() == acoords.size(), "Number of atoms in angle");
                        for(size_t m = 0; m < ind.size(); m++)
                        {
                            coordinates[ind[m]] = acoords[m];
                            ak = ind[m];
                        }
                        break;
                    }
                    double th0 = 45, th1 = 145, delta = 1;
                    int    nsteps = (th1-th0)/delta+1;
                    for(int i = 0; i < nsteps; i++)
                    {
                        double theta = (th0+i*delta);
                        coordinates[ak][0] = 1+std::cos(theta*DEG2RAD);
                        coordinates[ak][1] = std::sin(theta*DEG2RAD);
                        energies.clear();
                        bfc(msghandler, top.entry(itype), top.atoms(), &coordinates, &forces, &energies, epsilonr_);
                        fprintf(fp, "%10g  %10g\n", theta, energies[itype]);
                    }
                }
                break;
            case InteractionType::LINEAR_ANGLES:
                {
                    //! \todo take a into account
                    std::vector<gmx::RVec> coordinates = { { 0, 0, 0 }, { 1, 0, 0 }, { 2, 0, 0 } };
                    top.build(msghandler,
                              pd, &coordinates, 175.0, 5.0, missingParameters::Error);
                    forces.resize(top.nAtoms(), rvnul);
                    double r0 = 0.0, r1 = 0.1, delta = 0.001;
                    int    nsteps = (r1-r0)/delta+1;
                    for(int i = 0; i < nsteps; i++)
                    {
                        double xx = (r0+i*delta);
                        coordinates[1][1] = xx;
                        energies.clear();
                        bfc(msghandler, top.entry(itype), top.atoms(), &coordinates, &forces, &energies, epsilonr_);
                        fprintf(fp, "%10g  %10g\n", xx, energies[itype]);
                    }

                }
                break;
            case InteractionType::PROPER_DIHEDRALS:
                {
                    std::vector<gmx::RVec> coordinates = { { 0, 0, 0 }, { 1, 0, 0 }, { 1, 1, 0 }, { 1, 1, 1 } };
                    top.build(msghandler,
                              pd, &coordinates, 175.0, 5.0, missingParameters::Error);
                    forces.resize(top.nAtoms(), rvnul);
                    double th0 = 0, th1 = 360, delta = 2;
                    int    nsteps = (th1-th0)/delta+1;
                    for(int i = 0; i < nsteps; i++)
                    {
                        double theta = (th0+i*delta);
                        coordinates[3][0] = 1+std::cos(theta*DEG2RAD);
                        coordinates[3][1] = 1+std::sin(theta*DEG2RAD);
                        energies.clear();
                        bfc(msghandler, top.entry(itype), top.atoms(), &coordinates, &forces, &energies, epsilonr_);
                        fprintf(fp, "%10g  %10g\n", theta, energies[itype]);
                    }
                }
                break;
            case InteractionType::IMPROPER_DIHEDRALS:
                {
                    std::vector<gmx::RVec> coordinates = { { 1, 0.5, 0 }, { 0, 0, 0 }, { 2, 0, 0 }, { 1, 1.5, 0 } };
                    top.build(msghandler,
                              pd, &coordinates, 175.0, 5.0, missingParameters::Error);
                    forces.resize(top.nAtoms(), rvnul);
                    double th0 = -0.02, th1 = 0.02, delta = 0.001;
                    int    nsteps = (th1-th0)/delta+1;
                    for(int i = 0; i < nsteps; i++)
                    {
                        double theta = (th0+i*delta);
                        coordinates[3][2] = theta;
                        energies.clear();
                        bfc(msghandler, top.entry(itype), top.atoms(), &coordinates, &forces, &energies, epsilonr_);
                        fprintf(fp, "%10g  %10g\n", theta, energies[itype]);
                    }
                }
                break;
            default: // will not plot this interaction
                break;
            }
            gmx_ffclose(fp);
        }
    }
}

} // namespace alexandria

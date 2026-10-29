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
}

void ForceComputer::constructVsiteCoordinates(const Topology         *top,
                                              std::vector<gmx::RVec> *coordinates) const
{
    // Construct virtual site coordinates
    vsiteHandler_.constructPositions(top, coordinates);
}

void ForceComputer::spreadVsiteForces(const Topology               *top,
                                      const std::vector<gmx::RVec> *coordinates,
                                      std::vector<gmx::RVec>       *forces) const
{
    // Spread virtual site forces
    vsiteHandler_.distributeForces(top, *coordinates, forces);
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
        // We want to determine the displacement of the shell that cancels
        // the total force, assuming the other forces do not change, that is
        // a Newton step for the spring. For a spring potential V(r) the
        // stiffness along the core-shell vector is V''(r), perpendicular
        // to it V'(r)/r. These differ for a quartic spring.
        double kperp = shellForceConstant(fcShell, fcHyper, rHyper, rcs);
        double kpar  = fcShell;
        if (fcHyper > 0 && rcs > rHyper)
        {
            double drh = rcs - rHyper;
            kpar      += 12*fcHyper*drh*drh;
        }
        double fpar = iprod(forces[shell], r_core_shell)/rcs;
        for(int m = 0; m < DIM; m++)
        {
            double fparm = fpar*r_core_shell[m]/rcs;
            (*coordinates)[shell][m] += step*(fparm/kpar + (forces[shell][m] - fparm)/kperp);
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
    // One conformation, computed with the SIMD routines
    ConformationSet set;
    set.top = top;
    set.setCoordinates({ *coordinates });
    ShellOptions options;
    options.field       = field;
    options.resetShells = resetShells;
    options.relaxShells = relax;
    std::vector<std::map<InteractionType, double>> emaps;
    computeTotal(msg_handler, pd, &set, options, &emaps);
    *coordinates = set.confs[0].x;
    *forces      = set.confs[0].fcalc;
    *energies    = emaps[0];
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
                // Only the shells in the mask are relaxed
                if (!isShell[shellLocal])
                {
                    continue;
                }
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

/*! \brief Take the induction correction of a potential into account
 * A potential for the INDUCTIONCORRECTION, e.g. MACDANIEL_SCHMIDT, is part
 * of all electrostatics and of exchange plus induction, as in
 * ForceComputer::compute. If it is not treated separately, it is moved to
 * the induction.
 * \param[inout] cener      The energies of all conformations
 * \param[in]    correction The induction correction of the potential per conformation
 * \param[in]    toInduction Whether to move the correction to the induction
 */
static void addCorrection(ConformationEnergies      *cener,
                          const std::vector<double> &correction,
                          bool                       toInduction)
{
    auto eAllElec = cener->add(InteractionType::ALLELEC);
    auto eExchInd = cener->add(InteractionType::EXCHIND);
    for(size_t c = 0; c < cener->nconf() && c < correction.size(); c++)
    {
        if (eAllElec.has(c))
        {
            eAllElec.act[c] += correction[c];
        }
        if (eExchInd.has(c))
        {
            eExchInd.act[c] += correction[c];
        }
    }
    if (toInduction)
    {
        auto induction = cener->add(InteractionType::INDUCTION);
        auto icorr     = cener->find(InteractionType::INDUCTIONCORRECTION);
        for(size_t c = 0; c < cener->nconf() && c < correction.size(); c++)
        {
            if (icorr.has(c))
            {
                induction.set(c, (induction.has(c) ? induction.act[c] : 0.0) + correction[c]);
                icorr.have[c] = 0;
            }
        }
    }
}

/*! \brief The induction correction of a potential in all conformations
 * \param[in] cener The energies of all conformations
 * \return The correction of each conformation, empty if there is none
 */
static std::vector<double> potentialCorrection(ConformationEnergies *cener)
{
    std::vector<double> correction;
    auto icorr = cener->find(InteractionType::INDUCTIONCORRECTION);
    if (icorr.have)
    {
        correction.resize(cener->nconf(), 0.0);
        for(size_t c = 0; c < cener->nconf(); c++)
        {
            if (icorr.has(c))
            {
                correction[c] = icorr.act[c];
            }
        }
    }
    return correction;
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
                                 const std::vector<ACTMol *> &allmols,
                                 bool                         allForces) const
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
            for(auto molptr : allmols)
            {
                for(const auto &mon : *molptr->monomerCoords(im))
                {
                    nShellIter_ += mon.done ? 0 : 1;
                }
            }
            computeEnergies(msghandler, pd, allmols, im, true, true, false, allForces);
        }
        iter++;
    }
    // Check the conformations computed last
    double maxShellDistance2 = maxShellDistance_*maxShellDistance_;
    for(auto molptr : allmols)
    {
        auto  fh   = molptr->fragmentHandler();
        auto &tops = fh->topologies();
        auto &pols = fh->polInfo();
        if (im >= tops.size() || pols[im].nShell() == 0)
        {
            continue;
        }
        for(const auto &mon : *molptr->monomerCoords(im))
        {
            if (!mon.done &&
                shellForce2(tops[im].atoms(), pols[im].isShell(), mon.x, mon.fcalc,
                            0, maxShellDistance2)/pols[im].nShell() > msForceToler_)
            {
                nShellConvergeFailed_ += 1;
            }
        }
    }
}

void ForceComputer::relaxDimer(MsgHandler                  *msghandler,
                               const ForceField            *pd,
                               bool                         fixMonomer[2],
                               const std::vector<ACTMol *> &allmols,
                               bool                         allForces) const
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
            for(auto molptr : allmols)
            {
                for(const auto &dim : *molptr->dimerCoords())
                {
                    nShellIter_ += dim.done ? 0 : 1;
                }
            }
            computeEnergies(msghandler, pd, allmols, 2, true, true, false, allForces);
        }
        iter++;
    }
    // Check the conformations computed last, the shells of the monomers that
    // are relaxed have to be converged
    double maxShellDistance2 = maxShellDistance_*maxShellDistance_;
    for(auto molptr : allmols)
    {
        auto  fh   = molptr->fragmentHandler();
        auto &tops = fh->topologies();
        auto &pols = fh->polInfo();
        for(const auto &dim : *molptr->dimerCoords())
        {
            if (dim.done)
            {
                continue;
            }
            bool failed = false;
            for(size_t im = 0; im < tops.size() && im < 2; im++)
            {
                if (!fixMonomer[im] && pols[im].nShell() > 0 &&
                    shellForce2(tops[im].atoms(), pols[im].isShell(), dim.x, dim.fcalc,
                                fh->atomStart()[im], maxShellDistance2)/pols[im].nShell() > msForceToler_)
                {
                    failed = true;
                }
            }
            nShellConvergeFailed_ += failed ? 1 : 0;
        }
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
    // Without forces in the result, only those on the shells are needed
    // for relaxing them
    bool allForces = computeForces;
    // Do first calculation every time.
    computeEnergies(msghandler, pd, allmols, 0, forces, false, intermolecular, allForces); // Monomer 0
    computeEnergies(msghandler, pd, allmols, 1, forces, false, intermolecular, allForces); // Monomer 1
    computeEnergies(msghandler, pd, allmols, 2, forces, false, intermolecular, allForces); // Complete topology

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
            relaxMonomer(msghandler, pd, im, allmols, allForces);
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
            computeEnergies(msghandler, pd, allmols, 2, true, false, false, allForces);
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
                    computeEnergies(msghandler, pd, allmols, 2, true, false, false, allForces);
                }
                // Relax monomer im in the presence of the other monomer
                bool fixMonomer[2] = { im != 0, im != 1 };
                relaxDimer(msghandler, pd, fixMonomer, allmols, allForces);
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
            computeEnergies(msghandler, pd, allmols, 2, true, false, false, allForces);
        }
        // Final, total energy
        bool ff[2] = { false, false };
        relaxDimer(msghandler, pd, ff, allmols, allForces);
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
            // The induction correction from a potential, before the
            // induction correction from the shells is added
            auto correction = potentialCorrection(dimer);
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
                auto induction   = dimer->add(InteractionType::INDUCTION);
                auto icorrection = dimer->add(InteractionType::INDUCTIONCORRECTION);
                for(size_t i = 0; i < dimer->nconf() && i < eInduc2[k].size(); i++)
                {
                    double ecorr = i < correction.size() ? correction[i] : 0.0;
                    induction.set(i, eInduc2[k][i]);
                    icorrection.set(i, ecorr + eInd[i] - eInduc2[k][i]);
                }
            }
            if (!correction.empty())
            {
                addCorrection(dimer, correction, !separateInductionCorrection);
            }
        }
        auto &tops = mol.fragmentHandler()->topologies();
        for(size_t i = 0; i < tops.size() && !inter; i++)
        {
            calcInduction(monomers[i]);
            auto correction = potentialCorrection(monomers[i]);
            if (!correction.empty())
            {
                addCorrection(monomers[i], correction, false);
            }
        }
    }
    if (!allForces)
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

/*! \brief Add the forces of an electric field on the charges
 * \param[in]    atoms  The particles
 * \param[in]    field  The electric field
 * \param[inout] coords The conformations, only those not done are changed
 */
static void addFieldForces(const std::vector<ActAtom>  &atoms,
                           const gmx::RVec             &field,
                           std::vector<CoordsEnergy>   *coords)
{
    if (0 == field[XX] && 0 == field[YY] && 0 == field[ZZ])
    {
        return;
    }
    for(auto &ce : *coords)
    {
        if (ce.done)
        {
            continue;
        }
        for(size_t i = 0; i < atoms.size() && i < ce.fcalc.size(); i++)
        {
            real fac = FIELDFAC*atoms[i].charge();
            for(int m = 0; m < DIM; m++)
            {
                ce.fcalc[i][m] += fac*field[m];
            }
        }
    }
}

void ForceComputer::computeTotal(MsgHandler                                     *msghandler,
                                 const ForceField                               *pd,
                                 ConformationSet                                *set,
                                 const ShellOptions                             &options,
                                 std::vector<std::map<InteractionType, double>> *energies) const
{
    auto  top   = set->top;
    auto &confs = set->confs;
    auto  cdata = &set->data;
    auto &atoms = top->atoms();
    // Construct the virtual sites and reset the shells
    bool vsites = options.virtualSites && vsitesOutdated(top, cdata);
    std::vector<std::pair<size_t, size_t>> shells;
    if (options.resetShells)
    {
        shells = shellParents(atoms);
    }
    for(auto &conf : confs)
    {
        nEval_ += 1;
        if (vsites)
        {
            constructVsiteCoordinates(top, &conf.x);
        }
        resetShells(shells, &conf.x);
    }
    computeConformations(msghandler, pd, top, &confs, cdata, true, false);
    addFieldForces(atoms, options.field, &confs);
    // The energies before relaxing the shells give the induction energy
    auto &cener = cdata->energies;
    auto  itype = InteractionType::POLARIZATION;
    bool  relax = options.relax && pd->polarizable() && top->hasEntry(itype);
    cener.setBackups(relax);
    if (relax)
    {
        backupEnergies(&cener);
        // The shells to relax
        std::vector<bool> isShell(atoms.size(), false);
        size_t            nShell = 0;
        for(size_t i = 0; i < atoms.size(); i++)
        {
            if (atoms[i].pType() == ActParticle::Shell &&
                (options.relaxShells.empty() ||
                 options.relaxShells.end() != options.relaxShells.find(i)))
            {
                isShell[i] = true;
                nShell    += 1;
            }
        }
        resetShellSteps(&confs);
        // We have already done one energy calculation
        int  iter    = 1;
        bool allDone = false;
        while (iter < maxiter_ && !allDone)
        {
            allDone = (updateShells(top, 0, nShell, isShell, 0, &confs) == confs.size());
            if (!allDone)
            {
                // Only for the conformations that did not converge
                for(const auto &conf : confs)
                {
                    nShellIter_ += conf.done ? 0 : 1;
                }
                computeConformations(msghandler, pd, top, &confs, cdata, true, true);
                addFieldForces(atoms, options.field, &confs);
            }
            iter++;
        }
        // Check the conformations computed last
        double maxShellDistance2 = maxShellDistance_*maxShellDistance_;
        for(auto &conf : confs)
        {
            if (conf.done)
            {
                continue;
            }
            double msForce = shellForce2(atoms, isShell, conf.x, conf.fcalc, 0, maxShellDistance2);
            if (nShell > 0)
            {
                msForce /= nShell;
            }
            if (msForce > msForceToler_)
            {
                nShellConvergeFailed_ += 1;
            }
            conf.done = true;
        }
    }
    // The energy terms, as in compute
    energies->clear();
    for(size_t c = 0; c < confs.size(); c++)
    {
        std::map<InteractionType, double> emap;
        for(auto t : cener.types())
        {
            auto et = cener.find(t);
            if (et.has(c))
            {
                emap.insert({ t, et.act[c] });
            }
        }
        if (!options.energyTerms)
        {
            energies->push_back(emap);
            continue;
        }
        // The change of the electrostatic terms due to the relaxation
        // of the shells is the induction energy.
        bool   haveInduction = false;
        double eInduction    = 0;
        for(auto et : { InteractionType::ELECTROSTATICS,
                        InteractionType::POLARIZATION,
                        InteractionType::CHARGETRANSFER })
        {
            auto tt = emap.find(et);
            if (emap.end() != tt)
            {
                double before = relax ? cener.find(et).before[c] : tt->second;
                eInduction   += tt->second - before;
                tt->second    = before;
                haveInduction = true;
            }
        }
        if (haveInduction)
        {
            emap.insert_or_assign(InteractionType::INDUCTION, eInduction);
        }
        auto sum = [&emap](std::initializer_list<InteractionType> itypes)
        {
            double esum = 0;
            for(auto it : itypes)
            {
                auto ee = emap.find(it);
                if (emap.end() != ee)
                {
                    esum += ee->second;
                }
            }
            return esum;
        };
        emap.insert_or_assign(InteractionType::ALLELEC,
                              sum({ InteractionType::ELECTROSTATICS,
                                    InteractionType::POLARIZATION,
                                    InteractionType::INDUCTION,
                                    InteractionType::INDUCTIONCORRECTION }));
        emap.insert_or_assign(InteractionType::EXCHIND,
                              sum({ InteractionType::EXCHANGE,
                                    InteractionType::INDUCTION,
                                    InteractionType::INDUCTIONCORRECTION }));
        energies->push_back(emap);
    }
    // Spread forces to atoms
    for(auto &conf : confs)
    {
        if (options.virtualSites)
        {
            spreadVsiteForces(top, &conf.x, &conf.fcalc);
        }
    }
}

void ForceComputer::computeInteraction(MsgHandler                                     *msghandler,
                                       const ForceField                               *pd,
                                       ACTMol                                         *mol,
                                       std::vector<std::vector<gmx::RVec>>            *coords,
                                       bool                                            separateInductionCorrection,
                                       std::vector<std::map<InteractionType, double>> *einter,
                                       std::vector<double>                            *epot,
                                       std::vector<std::vector<gmx::RVec>>            *forces) const
{
    einter->clear();
    epot->clear();
    forces->clear();
    if (!mol->isDimer())
    {
        return;
    }
    // Replace the conformations of the dimer and its monomers
    for(size_t what = 0; what < 3; what++)
    {
        mol->conformationData(what)->invalidate();
    }
    auto  &dimers = *mol->dimerCoords();
    auto  &mons0  = *mol->monomerCoords(0);
    auto  &mons1  = *mol->monomerCoords(1);
    size_t natom1 = mol->fragmentHandler()->atomStart()[1];
    dimers.clear();
    mons0.clear();
    mons1.clear();
    for(const auto &x : *coords)
    {
        CoordsEnergy ce = {};
        ce.x = x;
        dimers.push_back(ce);
        CoordsEnergy m0 = {}, m1 = {};
        m0.x.assign(x.begin(), x.begin()+natom1);
        m1.x.assign(x.begin()+natom1, x.end());
        m0.xrelax = m0.x;
        m1.xrelax = m1.x;
        mons0.push_back(m0);
        mons1.push_back(m1);
    }
    computeRelaxedEnergies(msghandler, pd, { mol }, true, separateInductionCorrection, true);
    mol->storeEnergies();
    for(size_t c = 0; c < dimers.size(); c++)
    {
        std::map<InteractionType, double> emap;
        for(const auto &e : dimers[c].eMap)
        {
            // The polarization is part of the induction
            if (e.second.haveACT() && InteractionType::POLARIZATION != e.first)
            {
                emap.insert({ e.first, e.second.eact() });
            }
        }
        // The dimer energy is the interaction energy plus the energies of
        // the relaxed monomers, and the interaction forces are the dimer
        // forces minus the forces in the monomers.
        auto   ep = emap.find(InteractionType::EPOT);
        double e  = (emap.end() != ep) ? ep->second : 0;
        auto   f  = dimers[c].fcalc;
        for(size_t k = 0; k < 2; k++)
        {
            auto &mon = (*mol->monomerCoords(k))[c];
            auto  mep = mon.eMap.find(InteractionType::EPOT);
            if (mon.eMap.end() != mep && mep->second.haveACT())
            {
                e += mep->second.eact();
            }
            size_t offset = (k == 0) ? 0 : natom1;
            for(size_t i = 0; i < mon.fcalc.size(); i++)
            {
                rvec_dec(f[offset+i], mon.fcalc[i]);
            }
        }
        (*coords)[c] = dimers[c].x;
        einter->push_back(emap);
        epot->push_back(e);
        forces->push_back(f);
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
    // One conformation, computed with the SIMD routines, the shells
    // and the virtual sites stay where they are
    ConformationSet set;
    set.top = top;
    set.setCoordinates({ *coordinates });
    ShellOptions options;
    options.relax        = false;
    options.resetShells  = false;
    options.field        = field;
    options.energyTerms  = false;
    options.virtualSites = false;
    std::vector<std::map<InteractionType, double>> emaps;
    computeTotal(msg_handler, pd, &set, options, &emaps);
    // The evaluations are counted by the callers of computeOnce
    nEval_      -= 1;
    *forces      = set.confs[0].fcalc;
    *energies    = emaps[0];
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

/*! \brief Compute the energy of one potential for one structure
 * \param[in] msghandler For warnings and errors
 * \param[in] pot        The potential
 * \param[in] entries    The interactions
 * \param[in] atoms      The atoms
 * \param[in] x          The coordinates
 * \param[in] epsilonr   Relative dielectric constant
 * \return The energy
 */
static double potentialEnergy(MsgHandler                   *msghandler,
                              Potential                     pot,
                              const TopologyEntryVector    &entries,
                              const std::vector<ActAtom>   &atoms,
                              const std::vector<gmx::RVec> &x,
                              double                        epsilonr)
{
    std::vector<CoordsEnergy> ce(1);
    ce[0].x = x;
    computePotential(msghandler, pot, entries, atoms, &ce, epsilonr, false);
    auto ep = ce[0].eMap.find(InteractionType::EPOT);
    return (ce[0].eMap.end() != ep && ep->second.haveACT()) ? ep->second.eact() : 0;
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
                        double ener = potentialEnergy(msghandler, fs.potential(), top.entry(itype),
                                                      top.atoms(), coordinates, epsilonr_);
                        fprintf(fp, "%10g  %10g\n", theta, ener);
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
                        double ener = potentialEnergy(msghandler, fs.potential(), top.entry(itype),
                                                      top.atoms(), coordinates, epsilonr_);
                        fprintf(fp, "%10g  %10g\n", xx, ener);
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
                        double ener = potentialEnergy(msghandler, fs.potential(), top.entry(itype),
                                                      top.atoms(), coordinates, epsilonr_);
                        fprintf(fp, "%10g  %10g\n", theta, ener);
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
                        double ener = potentialEnergy(msghandler, fs.potential(), top.entry(itype),
                                                      top.atoms(), coordinates, epsilonr_);
                        fprintf(fp, "%10g  %10g\n", theta, ener);
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

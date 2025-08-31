/*
 * This source file is part of the Alexandria Chemistry Toolkit.
 *
 * Copyright (C) 2021-2025
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

#include <set>

#include <cstdlib>

#include "act/basics/chargemodel.h"
#include "act/basics/msg_handler.h"
#include "act/forcefield/forcefield_parametername.h"
#include "act/forces/forcecomputerimpl.h"
#include "act/qgen/qtype.h"
#include "gromacs/gmxpreprocess/grompp-impl.h"
#include "gromacs/math/vec.h"
#include "gromacs/utility/futil.h"

namespace alexandria
{

static double dotProdRvec(const std::vector<bool>      &isShell,
                          const std::vector<gmx::RVec> &rv)
{
    double dpr = 0;
    int    i   = 0;
    for(const auto &rr : rv)
    {
        if (isShell[i++])
        {
            dpr += iprod(rr, rr);
        }
    }
    return dpr;
}

ForceComputer::ForceComputer(double  msForce,
                             int     maxiter,
                             double  maxShellDistance)
{
    msForceToler_     = msForce;
    maxiter_          = maxiter;
    maxShellDistance_ = maxShellDistance;
    clear_mat(box_);
    real dt = 0.001;
    vsiteHandler_ = new VsiteHandler(box_, dt);
}

ForceComputer::~ForceComputer()
{
    if (vsiteHandler_)
    {
        delete vsiteHandler_;
    }
}

void ForceComputer::constructVsiteCoordinates(const Topology         *top,
                                              std::vector<gmx::RVec> *coordinates) const
{
    // Construct virtual site coordinates
    vsiteHandler_->constructPositions(top, coordinates, box_);
}

void ForceComputer::spreadVsiteForces(const Topology         *top,
                                      std::vector<gmx::RVec> *coordinates,
                                      std::vector<gmx::RVec> *forces) const
{
    // Spread virtual site forces
    vsiteHandler_->distributeForces(top, *coordinates, forces, box_);
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
    constructVsiteCoordinates(top, coordinates);
    // Short-cut
    auto &atoms = top->atoms();
    // Reset shells if needed
    if (resetShells)
    {
        for(size_t i = 0; i < top->nAtoms(); i++)
        {
            for(int sh : atoms[i].shells())
            {
                copy_rvec((*coordinates)[i], (*coordinates)[sh]);
            }
        }
    }
    // Do first calculation every time.
    computeOnce(pd, top, coordinates, forces, energies, field);
    // Store total electrostatics energy in independent copy.
    std::map<InteractionType, double> eBefore = *energies;

    // Now let's have a look whether we are polarizable
    auto itype = InteractionType::POLARIZATION;
    // mean square shell force
    double msForce = 0;
    if (pd->polarizable() && top->hasEntry(itype))
    {
        // Is this particle a shell?
        std::vector<bool>   isShell;
        // One over force constant for this particle
        std::vector<double> fcShell_1;
        auto &ffpl  = pd->findForcesConst(itype);
        auto pol_name = potentialToParameterName(ffpl.potential());
        int  nshell = 0;
        for(auto &aa : atoms)
        {
            bool bIS = aa.pType() == ActParticle::Shell;
            isShell.push_back(bIS);
            double fc_1 = 0;
            if (bIS)
            {
                Identifier atID(aa.ffType());
                auto alpha = ffpl.findParameterTypeConst(atID, pol_name[polALPHA]).internalValue();
                auto q     = aa.charge();
                if (alpha > 0 && q != 0)
                {
                    fc_1 = alpha/(q*q*ONE_4PI_EPS0);
                }
                nshell += 1;
            }
            fcShell_1.push_back(fc_1);
        }
        msForce    = dotProdRvec(isShell, *forces)/nshell;
        auto &pols = top->entry(itype);
        int   iter = 1;
        // Golden ratio, may be used for overrelaxation
        // double gold     = 0.5*(1+std::sqrt(5.0));
        // Prevent shells from flying off.
        real maxShellDistance2 = gmx::square(maxShellDistance_);
        while (msForce > msForceToler_ && iter < maxiter_)
        {
            // Loop over polarizabilities
            for(const auto &p : pols)
            {
                // Displace the shells according to the force
                // Since the potential is harmonic we use Hooke's law
                // F = k dx -> dx = F / k
                // TODO Optimize this protocol using overrelaxation
                int shell = p->atomIndex(1);
                if (relax.empty() || relax.end() != relax.find(shell))
                {
                    for(int m = 0; m < DIM; m++)
                    {
                        (*coordinates)[shell][m] += (*forces)[shell][m] * fcShell_1[shell];
                    }
                    // Check distance from core, if there is only one
                    if (atoms[shell].cores().size() == 1)
                    {
                        int core = atoms[shell].cores()[0];
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
            computeOnce(pd, top, coordinates, forces, energies, field);
            msForce  = dotProdRvec(isShell, *forces)/nshell;
            iter    += 1;
        }
        if (msg_handler && msg_handler->debug())
        {
            if (msForce > msForceToler_)
            {
                msg_handler->msg(ACTStatus::Warning,
                                 gmx::formatString("Shell optimization did not converge. RMS force is %g", std::sqrt(msForce/pols.size())));
            }
            else
            {
                msg_handler->msg(ACTStatus::Debug, "Shell optimization converged");
            }
        }
    }
    {
        // Induction energy
        double eInduction = 0;
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
                eInduction += tt->second - eBefore[et];
                tt->second = eBefore[et];
            }
        }
        if (eInduction != 0)
        {
            energies->insert_or_assign(InteractionType::INDUCTION, eInduction);
        }
    }
    {
        // Sum of all electrostatic terms
        double allelec = 0;
        for(const auto &itype : { InteractionType::ELECTROSTATICS, InteractionType::POLARIZATION, InteractionType::INDUCTION, InteractionType::INDUCTIONCORRECTION })
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

void ForceComputer::computeOnce(const ForceField                  *pd,
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
        // Force field parameter list
        auto &ffpl = pd->findForcesConst(entry.first);
        // The function we need to do the math
        auto bfc   = getBondForceComputer(ffpl.potential());
        if (bfc)
        {
            // Now do the calculations and store the energy
            std::map<InteractionType, double> my_energy;
            auto ener = bfc(entry.second, atoms, coordinates, forces, &my_energy);
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
                    potentialToString(ffpl.potential()).c_str());
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
                else if (InteractionType::VDW == itype)
                {
                    ;
                }
                else
                {
                    bbb.push_back(Bond(linear_bonds[i].first, linear_bonds[i].second, bos[i]));
                }
            }
            Topology               top(bbb);

            auto subtype = i2s.find(itype);
            if (i2s.end() != subtype)
            {
                for(const auto &atomname : f.first.atoms())
                {
                    if (pd->hasParticleType(atomname))
                    {
                        auto p = pd->findParticleType(atomname);
                        if (p->hasOption(subtype->second))
                        {
                            ActAtom a(*p);
                            a.setCharge(1);
                            top.addAtom(a);
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
                    top.build(msghandler,
                              pd, &coordinates, 175.0, 5.0, missingParameters::Error);
                    forces.resize(top.nAtoms(), rvnul);
                    // First atom is zero, second must be the other particle
                    size_t jatom = 1+top.atoms()[0].shells().size()+top.atoms()[0].vsites().size();
                    if (jatom == 0)
                    {
                        GMX_THROW(gmx::InternalError(gmx::formatString("Could not find a second atom to make a plot, there are %zu atoms, interactionType %s, id %s", top.nAtoms(), interactionTypeToString(itype).c_str(), f.first.id().c_str()).c_str()));
                    }
                    std::vector<double> rr, vv, ff;
                    // Now do the calculations and store the energy
                    double r0 = 0.05, r1 = 1.0, delta = 0.001;
                    int    nsteps = (r1-r0)/delta+1;
                    for(int i = 0; i < nsteps; i++)
                    {
                        double x = r0+i*delta;
                        coordinates[jatom][0] = x;
                        rr.push_back(x);
                        energies.clear();
                        for(size_t k = 0; k < top.nAtoms(); k++)
                        {
                            copy_rvec(rvnul, forces[k]);
                        }
                        bfc(top.entry(itype), top.atoms(), &coordinates, &forces, &energies);
                        auto ener = energies[itype];
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
                    std::vector<gmx::RVec> coordinates = { { 0, 0, 0 }, { 1, 0, 0 }, { 1, 1, 0 } };
                    top.build(msghandler,
                              pd, &coordinates, 175.0, 5.0, missingParameters::Error);
                    forces.resize(top.nAtoms(), rvnul);
                    double th0 = 0, th1 = 180, delta = 1;
                    int    nsteps = (th1-th0)/delta+1;
                    for(int i = 0; i < nsteps; i++)
                    {
                        double theta = (th0+i*delta);
                        coordinates[2][0] = 1+std::cos(theta*DEG2RAD);
                        coordinates[2][1] = std::sin(theta*DEG2RAD);
                        energies.clear();
                        bfc(top.entry(itype), top.atoms(), &coordinates, &forces, &energies);
                        fprintf(fp, "%10g  %10g\n", theta, energies[itype]);
                    }
                }
                break;
            case InteractionType::LINEAR_ANGLES:
                {
                    // TODO take a into account
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
                        bfc(top.entry(itype), top.atoms(), &coordinates, &forces, &energies);
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
                        bfc(top.entry(itype), top.atoms(), &coordinates, &forces, &energies);
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
                        bfc(top.entry(itype), top.atoms(), &coordinates, &forces, &energies);
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

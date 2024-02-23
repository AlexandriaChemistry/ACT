/*
 * This source file is part of the Alexandria Chemistry Toolkit.
 *
 * Copyright (C) 2021-2024
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

#include "act/forces/forcecomputerimpl.h"
#include "gromacs/gmxpreprocess/grompp-impl.h"
#include "gromacs/math/vec.h"
#include "gromacs/utility/futil.h"
#include "act/qgen/qtype.h"

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

ForceComputer::ForceComputer(double   msForce,
                             int      maxiter)
{
    msForceToler_ = msForce;
    maxiter_      = maxiter;
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

double ForceComputer::compute(const ForceField                  *pd,
                              const Topology                    *top,
                              std::vector<gmx::RVec>            *coordinates,
                              std::vector<gmx::RVec>            *forces,
                              std::map<InteractionType, double> *energies,
                              const gmx::RVec                   &field) const
{
    // Spread virtual sites
    vsiteHandler_->constructPositions(top, coordinates, box_);
    // Do first calculation every time.
    computeOnce(pd, top, coordinates, forces, energies, field);
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
        int  nshell = 0;
        for(auto &aa : top->atoms())
        {
            bool bIS = aa.pType() == eptShell;
            isShell.push_back(bIS);
            double fc_1 = 0;
            if (bIS)
            {
                Identifier atID(aa.ffType());
                auto fc = ffpl.findParameterTypeConst(atID, "kshell").internalValue();
                if (fc != 0)
                {
                    fc_1 = 1.0/fc;
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
                for(int m = 0; m < DIM; m++)
                {
                    (*coordinates)[shell][m] += (*forces)[shell][m] * fcShell_1[shell];
                }
            }
            // Do next calculation
            computeOnce(pd, top, coordinates, forces, energies, field);
            msForce  = dotProdRvec(isShell, *forces)/nshell;
            iter    += 1;
        }
    }
    // Spread forces to atoms
    vsiteHandler_->distributeForces(top, *coordinates, forces, box_);
    return msForce;
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
    auto atoms = top->atoms();
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
        auto bfc   = getBondForceComputer(ffpl.gromacsType());
        if (bfc)
        {
            // Now do the calculations and store the energy
            std::map<InteractionType, double> my_energy;
            bfc(entry.second, top->atoms(), coordinates, forces, &my_energy);
            for(const auto &me : my_energy)
            {
                if (energies->find(me.first) != energies->end())
                {
                    GMX_THROW(gmx::InternalError(gmx::formatString("Energy term %s occurs twice",
                                                                   interactionTypeToString(me.first).c_str()).c_str()));
                }
                energies->insert({ me.first, me.second });
                epot += me.second;
            }
        }
        else if (!isVsite(entry.first))
        {
            fprintf(stderr, "Please implement a force function for type %s\n",
                    interaction_function[ffpl.gromacsType()].name);
        }
    }
    energies->insert({ InteractionType::EPOT, epot });
}

int ForceComputer::ftype(const ForceField   *pd,
                         InteractionType  itype) const
{
    int ftype = F_EPOT;
    if (pd->interactionPresent(itype))
    {
        ftype = pd->findForcesConst(itype).gromacsType();
    }
    return ftype;
}

void ForceComputer::plot(const ForceField  *pd,
                         InteractionType    itype) const
{
    if (!pd->interactionPresent(itype))
    {
        fprintf(stderr, "No such interaction %s in the force field.\n",
                interactionTypeToString(itype).c_str());
        return;
    }
    std::string btype("bondtype");
    std::string vdwtype("vdwtype");

    std::map<InteractionType, const std::string> i2s = {
        { InteractionType::BONDS,              btype },
        { InteractionType::ANGLES,             btype },
        { InteractionType::LINEAR_ANGLES,      btype },
        { InteractionType::PROPER_DIHEDRALS,   btype },
        { InteractionType::IMPROPER_DIHEDRALS, btype },
        { InteractionType::VDW,                vdwtype },
        { InteractionType::COULOMB,            vdwtype }
    };
    auto &fs = pd->findForcesConst(itype);
    // The function we need to do the math
    auto bfc = getBondForceComputer(fs.gromacsType());
    if (nullptr == bfc)
    {
        fprintf(stderr, "Please implement a force function for type %s\n",
                interaction_function[fs.gromacsType()].name);
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
            int  ntrain  = 0;
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
            case InteractionType::COULOMB:
                {
                    std::vector<gmx::RVec> coordinates = { { 0, 0, 0 }, { 1, 0, 0 } };
                    top.build(pd, &coordinates, 175.0, 5.0, missingParameters::Error);
                    forces.resize(top.nAtoms(), rvnul);
                    size_t jatom = top.nAtoms()/2;
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
                            auto irep = InteractionType::REPULSION;
                            auto idsp = InteractionType::DISPERSION;
                            if (energies.find(irep) != energies.end() &&
                                energies.find(idsp) != energies.end())
                            {
                                ener = energies[irep] + energies[idsp];
                            }
                        }
                        fprintf(fp, "%10g  %10g\n", x, ener);
                        vv.push_back(ener);
                        ff.push_back(forces[0][0]);
                    }
                    // Check whether force is derivative of energy
                    if (debug)
                    {
                        for(size_t i = 1; i < vv.size()-1; i++)
                        {
                            if (std::abs(ff[i]) > 1e-6)
                            {
                                double fnumeric = (vv[i+1]-vv[i-1])/(2*delta);
                                double relerror = (fnumeric-ff[i])/ff[i];
                                if (std::abs(relerror) > 1e-1)
                                {
                                    fprintf(debug, "%s: Force %g, expected %g. Relative error %g, r = %g v+ %g v- %g delta %g\n",
                                            filename.c_str(),
                                            ff[i], fnumeric, relerror, rr[i],
                                            vv[i+1], vv[i-1], 2*delta);
                                }
                            }
                        }
                    }
                }
                break;
            case InteractionType::ANGLES:
                {
                    std::vector<gmx::RVec> coordinates = { { 0, 0, 0 }, { 1, 0, 0 }, { 1, 1, 0 } };
                    top.build(pd, &coordinates, 175.0, 5.0, missingParameters::Error);
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
                    top.build(pd, &coordinates, 175.0, 5.0, missingParameters::Error);
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
                    top.build(pd, &coordinates, 175.0, 5.0, missingParameters::Error);
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
                    top.build(pd, &coordinates, 175.0, 5.0, missingParameters::Error);
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
            default:
                break;
            }
            gmx_ffclose(fp);
        }
    }
}

} // namespace alexandria

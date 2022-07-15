/*
 * This source file is part of the Alexandria Chemistry Toolkit.
 *
 * Copyright (C) 2014-2022
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

#include <cstdlib>

#include "alexandria/topology.h"
#include "act/forces/forcecomputerimpl.h"
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

double ForceComputer::compute(const Topology                    *top,
                              std::vector<gmx::RVec>            *coordinates,
                              std::vector<gmx::RVec>            *forces,
                              std::map<InteractionType, double> *energies) const
{
    gmx::RVec field = { 0, 0, 0 };
    // Do first calculation every time.
    computeOnce(top, coordinates, forces, energies, field);
    // Now let's have a look whether we are polarizable
    auto itype = InteractionType::POLARIZATION;
    if (!pd_->polarizable() || !top->hasEntry(itype))
    {
        return 0;
    }
    // Is this particle a shell?
    std::vector<bool>   isShell;
    // One over force constant for this particle
    std::vector<double> fcShell_1;
    auto &ffpl  = pd_->findForcesConst(itype);
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
    double msForceMax = (rmsForce_*rmsForce_)*nshell;
    double msForce    = dotProdRvec(isShell, *forces);
    auto pols         = top->entry(itype);
    int    iter       = 1;
    // Golden ratio, may be used for overrelaxation
    // double gold     = 0.5*(1+std::sqrt(5.0));
    while (msForce > msForceMax && iter < maxiter_)
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
        computeOnce(top, coordinates, forces, energies, field);
        msForce  = dotProdRvec(isShell, *forces);
        iter    += 1;
    }
    return std::sqrt(msForce/nshell);
}

void ForceComputer::computeOnce(const Topology                    *top,
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
        svmul(atoms[ff].charge(), field, (*forces)[ff]);
    }
    double epot = 0;
    for(const auto &entry : top->entries())
    {
        if (entry.second.empty())
        {
            continue;
        }
        // Force field parameter list
        auto &ffpl = pd_->findForcesConst(entry.first);
        // The function we need to do the math
        auto bfc   = getBondForceComputer(ffpl.fType());
        if (nullptr == bfc)
        {
            fprintf(stderr, "Please implement a force function for type %s\n", interaction_function[ffpl.fType()].name);
        }
        else
        {
            // Now do the calculations and store the energy
            double eee = bfc(entry.second, top->atoms(), coordinates, forces);
            energies->insert({ entry.first, eee });
            epot += eee;
        }
    }
    energies->insert({ InteractionType::EPOT, epot });
}

void ForceComputer::calcPolarizability(const Topology         *top,
                                       std::vector<gmx::RVec> *coordinates,
                                       QtypeProps             *qtp) const
{
    std::vector<gmx::RVec>            forces(coordinates->size());
    std::map<InteractionType, double> energies;
    gmx::RVec  field = { 0, 0, 0 };
    
    computeOnce(top, coordinates, &forces, &energies, field);
    std::vector<double> q;
    for (auto &at : top->atoms())
    {
        q.push_back(at.charge());
    }
    qtp->setQ(q);
    qtp->setX(*coordinates);
    qtp->calcMoments();
    auto mpo = MolPropObservable::DIPOLE;
    if (!qtp->hasMultipole(mpo))
    {
        GMX_THROW(gmx::InternalError("No dipole to compute."));
    }
    auto mu_ref = qtp->getMultipole(mpo);
    // Convert from e nm2/V to cubic nm
    double enm2_V = E_CHARGE*1e6*1e-18/(4*M_PI*EPSILON0_SI)*1e21;
    tensor alpha  = { { 0 } };
    double efield = 0.1;
    for (auto m = 0; m < DIM; m++)
    {
        field[m] = efield;
        computeOnce(top, coordinates, &forces, &energies, field);
        qtp->setX(*coordinates);
        field[m] = 0;
        qtp->calcMoments();
        auto qmu = qtp->getMultipole(mpo);
        for (auto n = 0; n < DIM; n++)
        {
            alpha[n][m] = enm2_V*((qmu[n]-mu_ref[n])/efield);
        }
    }
    // Store the tensor
    qtp->setPolarizabilityTensor(alpha);
    // Reset energies etc.
    computeOnce(top, coordinates, &forces, &energies, field);
}

int ForceComputer::ftype(InteractionType itype) const
{
    int ftype = F_EPOT;
    if (pd_->interactionPresent(itype))
    {
        ftype = pd_->findForcesConst(itype).fType();
    }
    return ftype;
}

void ForceComputer::plot(InteractionType itype) const
{
    if (!pd_->interactionPresent(itype))
    {
        fprintf(stderr, "No such interaction %s in the force field.\n",
                interactionTypeToString(itype).c_str());
        return;
    }
    auto &fs = pd_->findForcesConst(itype);
    // The function we need to do the math
    auto bfc = getBondForceComputer(fs.fType());
    if (nullptr == bfc)
    {
        fprintf(stderr, "Please implement a force function for type %s\n",
                interaction_function[fs.fType()].name);
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
            std::vector<gmx::RVec> forces;
            gmx::RVec rvnul = { 0, 0, 0 };
            for(const auto &atomname : f.first.atoms())
            {
                auto p = pd_->findParticleType(itype, atomname);
                if (p != pd_->particleTypesConst().end())
                {
                    top.addAtom(ActAtom(*p));
                }
            }
            if (top.nAtoms() < 2)
            {
                continue;
            }
            forces.resize(top.nAtoms(), rvnul);
            // Open outfile
            std::string filename = gmx::formatString("%s_%s.xvg", 
                                                     interactionTypeToString(itype).c_str(),
                                                     f.first.id().c_str());
            FILE *fp = gmx_ffopen(filename.c_str(), "w");
            switch (itype)
            {
            case InteractionType::BONDS:
            case InteractionType::VDW:
                {
                    std::vector<gmx::RVec> coordinates = { { 0, 0, 0 }, { 1, 0, 0 } };
                    top.build(pd_, coordinates, 175.0, 5.0, missingParameters::Error);
                    top.setIdentifiers(pd_);
                    top.fillParameters(pd_);
                    
                    // Now do the calculations and store the energy
                    double r0 = 0.05, r1 = 1.0, delta = 0.01;
                    int    nsteps = (r1-r0)/delta+1;
                    for(int i = 0; i < nsteps; i++)
                    {
                        double x = r0+i*delta;
                        coordinates[1][0] = x;
                        double eee = bfc(top.entry(itype), top.atoms(), &coordinates, &forces);
                        fprintf(fp, "%10g  %10g\n", x, eee);
                    }
                }
                break;
            case InteractionType::ANGLES:
                {
                    std::vector<gmx::RVec> coordinates = { { 0, 0, 0 }, { 1, 0, 0 }, { 1, 1, 0 } };
                    top.build(pd_, coordinates, 175.0, 5.0, missingParameters::Error);
                    top.setIdentifiers(pd_);
                    top.fillParameters(pd_);
                    double th0 = 0, th1 = 180, delta = 1;
                    int    nsteps = (th1-th0)/delta+1;
                    for(int i = 0; i < nsteps; i++)
                    {
                        double theta = (th0+i*delta);
                        coordinates[2][0] = 1+std::cos(theta*DEG2RAD);
                        coordinates[2][1] = std::sin(theta*DEG2RAD);
                        double eee = bfc(top.entry(itype), top.atoms(), &coordinates, &forces);
                        fprintf(fp, "%10g  %10g\n", theta, eee);
                    }
                }
                break;
            case InteractionType::LINEAR_ANGLES:
                {
                    // TODO take a into account
                    std::vector<gmx::RVec> coordinates = { { 0, 0, 0 }, { 1, 0, 0 }, { 2, 0, 0 } };
                    top.build(pd_, coordinates, 175.0, 5.0, missingParameters::Error);
                    top.setIdentifiers(pd_);
                    top.fillParameters(pd_);
                    double r0 = 0.0, r1 = 0.1, delta = 0.001;
                    int    nsteps = (r1-r0)/delta+1;
                    for(int i = 0; i < nsteps; i++)
                    {
                        double xx = (r0+i*delta);
                        coordinates[1][1] = xx;
                        double eee = bfc(top.entry(itype), top.atoms(), &coordinates, &forces);
                        fprintf(fp, "%10g  %10g\n", xx, eee);
                    }
                    
                }
                break;
            case InteractionType::PROPER_DIHEDRALS:
                {
                    std::vector<gmx::RVec> coordinates = { { 0, 0, 0 }, { 1, 0, 0 }, { 1, 1, 0 }, { 1, 1, 1 } };
                    top.build(pd_, coordinates, 175.0, 5.0, missingParameters::Error);
                    top.setIdentifiers(pd_);
                    top.fillParameters(pd_);
                    double th0 = 0, th1 = 360, delta = 2;
                    int    nsteps = (th1-th0)/delta+1;
                    for(int i = 0; i < nsteps; i++)
                    {
                        double theta = (th0+i*delta);
                        coordinates[3][0] = 1+std::cos(theta*DEG2RAD);
                        coordinates[3][1] = 1+std::sin(theta*DEG2RAD);
                        double eee = bfc(top.entry(itype), top.atoms(), &coordinates, &forces);
                        fprintf(fp, "%10g  %10g\n", theta, eee);
                    }
                }
                break;
            case InteractionType::IMPROPER_DIHEDRALS:
                {
                    std::vector<gmx::RVec> coordinates = { { 1, 0.5, 0 }, { 0, 0, 0 }, { 2, 0, 0 }, { 1, 1.5, 0 } };
                    top.build(pd_, coordinates, 175.0, 5.0, missingParameters::Error);
                    top.setIdentifiers(pd_);
                    top.fillParameters(pd_);
                    double th0 = -0.02, th1 = 0.02, delta = 0.001;
                    int    nsteps = (th1-th0)/delta+1;
                    for(int i = 0; i < nsteps; i++)
                    {
                        double theta = (th0+i*delta);
                        coordinates[3][2] = theta;
                        double eee = bfc(top.entry(itype), top.atoms(), &coordinates, &forces);
                        fprintf(fp, "%10g  %10g\n", theta, eee);
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

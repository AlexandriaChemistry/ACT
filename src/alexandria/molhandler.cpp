/*
 * This source file is part of the Alexandria Chemistry Toolkit.
 *
 * Copyright (C) 2014-2022
 *
 * Developers:
 *             Mohammad Mehdi Ghahremanpour,
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
 * \author David van der Spoel <david.vanderspoel@icm.uu.se>
 * \author Julian Marrades <julian.marrades@hotmail.es>
 */

#include "molhandler.h"

#include "gromacs/gmxlib/network.h"
#include "gromacs/math/do_fit.h"

#include "gromacs/mdtypes/inputrec.h"
#include "gromacs/mdtypes/enerdata.h"
#include "gromacs/topology/topology.h"

namespace alexandria
{

double MolHandler::computeHessian(      MyMol               *mol,
                                  const t_commrec           *crtmp,
                                  const std::vector<int>    &atomIndex,
                                        MatrixWrapper       *hessian,
                                        std::vector<double> *forceZero)
{
    std::vector<double> fneg;
    fneg.resize(DIM*atomIndex.size(), 0.0);
    
    double shellForceRMS = 0;
    (void) mol->calculateEnergy(crtmp, &shellForceRMS);
    double epot0  = mol->enerd_->term[F_EPOT];

    // Store central force vector
    forceZero->clear();
    for(auto &atom : atomIndex)
    {
        for(int m = 0; m < DIM; m++)
        {
            forceZero->push_back(mol->f_[atom][m]);
        }
    }
    double    stepSize     = 1e-6; // 1 pm
    for(size_t ai = 0; ai < atomIndex.size(); ai++)
    {
        auto atomI = atomIndex[ai];
        for(int atomXYZ = 0; atomXYZ < DIM; atomXYZ++)
        {
            double xyzRef = mol->state_->x[atomI][atomXYZ];
            for(int delta = 0; delta <= 1; delta++)
            {
                mol->state_->x[atomI][atomXYZ] = xyzRef + (2*delta-1)*stepSize;
                (void) mol->calculateEnergy(crtmp, &shellForceRMS);
                if (delta == 0)
                {
                    for(size_t aj = 0; aj < atomIndex.size(); aj++)
                    {
                        int atomJ = atomIndex[aj];
                        for (int d = 0; d < DIM; d++)
                        {
                            fneg[aj*DIM+d] = mol->f_[atomJ][d];
                        }
                    }
                }
            }
            mol->state_->x[atomI][atomXYZ] = xyzRef;  // Restore positions
            {
                int col = ai*DIM+atomXYZ;
                for(size_t aj = 0; aj < atomIndex.size(); aj++)
                {
                    int    atomJ = atomIndex[aj];
                    for(int d = 0; d < DIM; d++)
                    {
                        int    row   = aj*DIM+d;
                        double value = -(mol->f_[atomJ][d]-fneg[row])/(2*stepSize);
                        if (false && debug)
                        {
                            fprintf(debug, "Setting H[%2d][%2d] = %10g\n", row, col, value); 
                        }
                        hessian->set(row, col, value);
                    }
                }
            }
        }
    }
    return epot0;
}

immStatus MolHandler::minimizeCoordinates(MyMol  *mol,
                                          double *rmsd)
{
    mol->updateMDAtoms();
    // We need to use a single core system for minimizing shells
    t_commrec *crtmp = init_commrec();
    crtmp->nnodes = 1;
    crtmp->nodeid = 0;

    // TODO: check if this is really necessary
    mol->restoreCoordinates();
    immStatus imm          = immStatus::OK;
    auto      mdatoms      = mol->MDatoms_->get()->mdatoms();
    // Below is a Newton-Rhapson algorithm
    bool      converged    = false;
    double    myForceToler = mol->inputrec_->em_tol; // kJ/mol nm
    int       myIter       = 0;
    int       maxIter      = 100;
    // List of atoms (not shells) and weighting factors
    std::vector<int> theAtoms;
    std::vector<real> w_rls;
    for(int atom = 0; atom < mdatoms->nr; atom++)
    {
        if (mdatoms->ptype[atom] == eptAtom)
        {
            theAtoms.push_back(atom);
            w_rls.push_back(1);
        }
        else
        {
            w_rls.push_back(0);
        }
    }
    MatrixWrapper       Hessian(DIM*theAtoms.size(), DIM*theAtoms.size());
    std::vector<double> f0, f00;
    double              epot0     = 0;
    bool                firstStep = true;
    // Now start the minimization loop.
    do
    {
        double newEpot = computeHessian(mol, crtmp, theAtoms, &Hessian, &f0);
        if (firstStep)
        {
            epot0     = newEpot;
            f00       = f0;
            firstStep = false;
        }
        std::vector<double> deltaX(DIM*theAtoms.size(), 0.0);
        // Solve H delta X = -grad (E) = force(E)
        int result = Hessian.solve(f0, &deltaX);
        if (0 == result)
        {
            // Set a maximum displacement to prevent exploding molecules
            double deltaXTolerance = 0.01; // nm
            double normDeltaX = 0;
            for(auto &dx : deltaX)
            {
                normDeltaX += dx*dx;
            }
            double rmsDeltaX = std::sqrt(normDeltaX/deltaX.size());
            double scaleDeltaX = 1;
            if (rmsDeltaX > deltaXTolerance)
            {
                scaleDeltaX = deltaXTolerance/rmsDeltaX;
            }
            int  i    = 0;
            for (auto &atomI : theAtoms)
            {
                for (int m = 0; m < DIM; m++)
                {
                    mol->state_->x[atomI][m] += scaleDeltaX*deltaX[i++];
                }
            }
        }
        else
        {
            break;
        }
                
        // One more energy and force calculation with the new coordinates
        double shellForceRMS = 0;
        imm = mol->calculateEnergy(crtmp, &shellForceRMS);
        
        double msAtomForce  = 0;
        for(size_t kk = 0; kk < theAtoms.size(); kk++)
        {
            int atomI = theAtoms[kk];
            msAtomForce  += iprod(mol->f_[atomI], mol->f_[atomI]);
        }
        msAtomForce /= theAtoms.size();
        converged = msAtomForce <= myForceToler;
        
        if (debug)
        {
            fprintf(debug, "%s rmsForce %10g Epot before %10g now %10g\n", mol->getMolname().c_str(),
                    std::sqrt(msAtomForce), epot0, mol->enerd_->term[F_EPOT]);
        }
        myIter += 1;
    }
    while (!converged && myIter < maxIter);
    std::vector<gmx::RVec> xmin(mol->mtop_->natoms);
    for (auto &kk : theAtoms)
    {
        copy_rvec(mol->state_->x[kk], xmin[kk]);
        if (debug)
        {
            fprintf(debug, "force[%2d] = %10g %10g %10g\n", kk, mol->f_[kk][XX], mol->f_[kk][YY], mol->f_[kk][ZZ]);
        }
    }
    // Fetch back the input structure.
    mol->restoreCoordinates();
    // Restore the forces
    for(size_t kk = 0; kk < theAtoms.size(); kk++)
    {
        for(int m = 0; m < DIM; m++)
        {
            mol->f_[theAtoms[kk]][m] = f00[DIM*kk+m];
        }
    }
    std::vector<gmx::RVec> xp(mol->mtop_->natoms);
    for (auto &kk : theAtoms)
    {
        copy_rvec(mol->state_->x[kk], xp[kk]);
    }
    // Compute RMSD
    do_fit(w_rls.size(), w_rls.data(),
           as_rvec_array(xp.data()), as_rvec_array(xmin.data()));
    double msd  = 0;
    for (auto &kk : theAtoms)
    {
        rvec dx;
        rvec_sub(xmin[kk], mol->state_->x[kk], dx);
        msd += iprod(dx, dx);
    }
    *rmsd = std::sqrt(msd/theAtoms.size());
    done_commrec(crtmp);
    return imm;
}

} // namespace alexandria

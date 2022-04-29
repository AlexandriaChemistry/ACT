/*
 * This source file is part of the Alexandria Chemistry Toolkit.
 *
 * Copyright (C) 2022
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

#include "act/molprop/molpropobservable.h"
#include "act/utility/units.h"
#include "gromacs/math/units.h"
#include "gromacs/gmxlib/network.h"
#include "gromacs/linearalgebra/eigensolver.h"
#include "gromacs/math/do_fit.h"
#include "gromacs/mdtypes/enerdata.h"
#include "gromacs/mdtypes/inputrec.h"
#include "gromacs/topology/topology.h"

namespace alexandria
{

double MolHandler::computeHessian(      MyMol               *mol,
                                  const t_commrec           *crtmp,
                                  const std::vector<int>    &atomIndex,
                                        MatrixWrapper       *hessian,
                                        std::vector<double> *forceZero) const
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
                            int row   = aj*DIM+d;
                            fneg[row] = mol->f_[atomJ][d];
                        }
                    }
                }
            }
            mol->state_->x[atomI][atomXYZ] = xyzRef;  // Restore positions
            
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
    (void) mol->calculateEnergy(crtmp, &shellForceRMS);
    return epot0;
}

void MolHandler::nma(MyMol               *mol,
                     std::vector<double> *frequencies,
                     std::vector<double> *intensities,
                     FILE                *fp) const
{

    // Init single use commrec
    t_commrec *crtmp = init_commrec();
    crtmp->nnodes = 1;
    crtmp->nodeid = 0;

    // Get the indices of the real atoms of the molecule (not shells and such)
    auto atoms = mol->atoms();
    std::vector<int> atomIndex;
    for(int atom = 0; atom < atoms->nr; atom++)
    {
        if (atoms->atom[atom].ptype == eptAtom)
        {
            atomIndex.push_back(atom);
        }
    }

    // Compute and average hessian
    const int matrixSide = DIM*atomIndex.size();
    MatrixWrapper       hessian(matrixSide, matrixSide);
    std::vector<double> f0;
    computeHessian(mol, crtmp, atomIndex, &hessian, &f0);
    hessian.averageTriangle();

    // Dispose single use commrec
    done_commrec(crtmp);

    // divide elements hessian[i][j] by sqrt(mass[i])*sqrt(mass[j])
    double massFac;
    for (size_t i = 0; (i < atomIndex.size()); i++)
    {
        size_t ai = atomIndex[i];
        for (size_t j = 0; (j < DIM); j++)
        {
            for (size_t k = 0; (k < atomIndex.size()); k++)
            {
                size_t ak = atomIndex[k];
                massFac   = gmx::invsqrt(atoms->atom[ai].m * atoms->atom[ak].m);
                for (size_t l = 0; (l < DIM); l++)
                {
                    hessian.mult(i*DIM+j, k*DIM+l, massFac);
                }
            }
        }
    }

    // Get the vibrational frequencies
    int rot_trans = 6;
    if (mol->linearMolecule())
    {
        rot_trans = 5;
    }
    // Call diagonalization routine
    // fprintf(stderr, "\nDiagonalizing to find vectors...\n");
    auto hessianFlat = hessian.flatten();
    std::vector<double> eigenvalues(matrixSide);
    std::vector<double> eigenvectors(matrixSide*matrixSide);
    eigensolver(hessianFlat.data(), matrixSide, 0, matrixSide - 1,
                eigenvalues.data(), eigenvectors.data());
    // Scale the output eigenvectors
    for (int i = 0; i < matrixSide; i++)
    {
        for (size_t j = 0; j < atomIndex.size(); j++)
        {
            size_t aj = atomIndex[j];
            massFac   = gmx::invsqrt(atoms->atom[aj].m);
            for (size_t k = 0; (k < DIM); k++)
            {
                eigenvectors[i * matrixSide + j * DIM + k] *= massFac;
            }
        }
    }

    frequencies->clear();
    intensities->clear();
    for (size_t i = rot_trans; i < eigenvalues.size(); i++)
    {
        auto val = eigenvalues[i];
        if (val < 0)
        {
            if (fp)
            {
                fprintf(fp, "Warning: negative eigenvalue %zu = %g\n", i, val);
            }
            // We need to store something to be able to compare
            frequencies->push_back(0);
        }
        else
        {
            #define SOL_CM_PER_S 29979245800
            frequencies->push_back(
                1.0/(2.0*M_PI*SOL_CM_PER_S) * 1E12 * std::sqrt(val)
            );
        }
    }

    // Get the eigenvectors into a MatrixWrapper
    MatrixWrapper eigenvecMat(eigenvectors, matrixSide);

    if (fp)
    {
        std::vector<GenericProperty *> harm;
        auto mpo = MolPropObservable::FREQUENCY;
        for (auto &ee : mol->experimentConst())
        {
            if (ee.hasMolPropObservable(mpo))
            {
                harm = ee.propertyConst(mpo);
                break;
            }
        }
        const char *unit = mpo_unit2(mpo);
        double delta = 0;
        if (!harm.empty())
        {
            fprintf(fp, "Electronic vibrational frequencies: (%s):\n", unit);
            size_t k = 0;
            for(auto &ff : harm[0]->getVector())
            {
                fprintf(fp, "  %8.3f", convertFromGromacs(ff, unit));
                delta += gmx::square(ff - (*frequencies)[k]);
                k++;
            }
            fprintf(fp, "\n");
        }
        fprintf(fp, "Alexandria vibrational frequencies: (%s):\n", unit);
        for (const auto &freq : *frequencies)
        {
            fprintf(fp, "  %8.3f", convertFromGromacs(freq, unit));
        }
        if (delta > 0)
        {
            fprintf(fp, "\nFrequency RMSD %g\n",
                    convertFromGromacs(std::sqrt(delta/frequencies->size()),
                                       unit));
        }
        else
        {
            fprintf(fp, "\n\n");
        }
    }
    // Print eigenvalues and eigenvectors to debug files
    if (debug)
    {
        // Print vibrational frequencies to file
        const int FLOAT_SIZE  = 13;
        fprintf(debug, "\nMolecule: %s\nHessian eigenvalues:\n[ ", mol->getMolname().c_str());
        for (const auto &val : eigenvalues)
        {
            fprintf(debug, "%-*g ", FLOAT_SIZE, val);
        }
        fprintf(debug, "]\nHessian eigenvectors:\n%s\n", eigenvecMat.toString().c_str());
    }

}

immStatus MolHandler::minimizeCoordinates(MyMol  *mol,
                                          double *rmsd) const
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
    // The original coordinates, including shell positions
    std::vector<gmx::RVec> xold(mol->mtop_->natoms);
    // Create Hessian, just containings atoms, not shells
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
            // Store energy
            epot0     = newEpot;
            // Store forces
            f00       = f0;
            // Store coordinates, including shell positions
            for (int kk = 0; kk < mol->mtop_->natoms; kk++)
            {
                copy_rvec(mol->state_->x[kk], xold[kk]);
            }
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
    // Re-compute the energy one last time.
    double shellForceRMS = 0;
    (void) mol->calculateEnergy(crtmp, &shellForceRMS);
    // Compute RMSD
    std::vector<gmx::RVec> xmin(mol->mtop_->natoms);
    for (auto &kk : theAtoms)
    {
        copy_rvec(mol->state_->x[kk], xmin[kk]);
    }
    do_fit(w_rls.size(), w_rls.data(),
           as_rvec_array(xold.data()), as_rvec_array(xmin.data()));
    double msd  = 0;
    for (auto &kk : theAtoms)
    {
        rvec dx;
        rvec_sub(mol->state_->x[kk], xold[kk], dx);
        msd += iprod(dx, dx);
    }
    *rmsd = std::sqrt(msd/theAtoms.size());
    done_commrec(crtmp);
    return imm;
}

} // namespace alexandria

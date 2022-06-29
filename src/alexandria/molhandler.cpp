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
#include "alexandria/velocityhandler.h"
#include "gromacs/fileio/pdbio.h"
#include "gromacs/utility/futil.h"
#include "gromacs/math/units.h"
#include "gromacs/gmxlib/network.h"
#include "gromacs/linearalgebra/eigensolver.h"
#include "gromacs/fileio/xvgr.h"
#include "gromacs/math/do_fit.h"
#include "gromacs/math/vectypes.h"
#include "gromacs/mdtypes/enerdata.h"
#include "gromacs/mdtypes/inputrec.h"
#include "gromacs/topology/topology.h"

namespace alexandria
{

double MolHandler::computeHessian(      MyMol               *mol,
                                  const ForceComputer       *forceComp,
                                  const std::vector<int>    &atomIndex,
                                        MatrixWrapper       *hessian,
                                        std::vector<double> *forceZero) const
{
    std::vector<double> fneg;
    fneg.resize(DIM*atomIndex.size(), 0.0);
    
    (void) mol->calculateEnergy(forceComp);
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
                (void) mol->calculateEnergy(forceComp);
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
    (void) mol->calculateEnergy(forceComp);

    return epot0;
}

void MolHandler::nma(MyMol               *mol,
                     const ForceComputer *forceComp,
                     std::vector<double> *frequencies,
                     std::vector<double> *intensities,
                     FILE                *fp) const
{
    // Get the indices of the real atoms of the molecule (not shells and such)
    auto atoms = mol->topology()->atoms();
    std::vector<int> atomIndex;
    for(size_t atom = 0; atom < atoms.size(); atom++)
    {
        if (atoms[atom].pType() == eptAtom)
        {
            atomIndex.push_back(atom);
        }
    }

    // Compute and average hessian
    const int matrixSide = DIM*atomIndex.size();
    MatrixWrapper       hessian(matrixSide, matrixSide);
    std::vector<double> f0;
    computeHessian(mol, forceComp, atomIndex, &hessian, &f0);
    hessian.averageTriangle();

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
                massFac   = gmx::invsqrt(atoms[ai].mass() * atoms[ak].mass());
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
            massFac   = gmx::invsqrt(atoms[aj].mass());
            for (size_t k = 0; (k < DIM); k++)
            {
                eigenvectors[i * matrixSide + j * DIM + k] *= massFac;
            }
        }
    }

    frequencies->clear();
    intensities->clear();
    auto mpo = MolPropObservable::FREQUENCY;
    const char *unit = mpo_unit2(mpo);
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
            double f = 1.0/(2.0*M_PI*SOL_CM_PER_S) * 1E12 * std::sqrt(val);
            frequencies->push_back(convertToGromacs(f, unit));
        }
    }

    // Get the eigenvectors into a MatrixWrapper
    MatrixWrapper eigenvecMat(eigenvectors, matrixSide);

    if (fp)
    {
        std::vector<GenericProperty *> harm;
        for (auto &ee : mol->experimentConst())
        {
            if (ee.hasMolPropObservable(mpo))
            {
                harm = ee.propertyConst(mpo);
                break;
            }
        }
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

immStatus MolHandler::minimizeCoordinates(MyMol               *mol,
                                          const ForceComputer *forceComp) const
{
    // TODO: check if this is really necessary
    mol->restoreCoordinates(coordSet::Minimized);  // Is minimized defined??? I guess it makes sense: if not defined, nothing changes
    immStatus imm          = immStatus::OK;
    // Below is a Newton-Rhapson algorithm
    bool      converged    = false;
    double    myForceToler = forceComp->convergenceTolerance();
    int       myIter       = 0;
    int       maxIter      = 100;
    // List of atoms (not shells) and weighting factors
    auto      myatoms      = mol->atomsConst();
    std::vector<int> theAtoms;
    for(size_t atom = 0; atom < myatoms.size(); atom++)
    {
        if (myatoms[atom].pType() == eptAtom)
        {
            theAtoms.push_back(atom);
        }
    }
    // Create Hessian, just containings atoms, not shells
    MatrixWrapper       Hessian(DIM*theAtoms.size(), DIM*theAtoms.size());
    std::vector<double> f0, f00;
    double              epot0     = 0;
    bool                firstStep = true;
    // Now start the minimization loop.
    do
    {
        double newEpot = computeHessian(mol, forceComp, theAtoms, &Hessian, &f0);
        if (firstStep)
        {
            // Store energy
            epot0     = newEpot;
            // Store forces
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
        (void) mol->calculateEnergy(forceComp);
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
    (void) mol->calculateEnergy(forceComp);
    mol->backupCoordinates(coordSet::Minimized);
    return imm;
}

double MolHandler::coordinateRmsd(MyMol                                       *mol,
                                  std::map<coordSet, std::vector<gmx::RVec> > *x) const
{
    // Compute RMSD
    auto   myatoms = mol->atomsConst();
    double tmass   = mol->totalMass();
    std::vector<real> myMass;
    x->clear();
    for (const auto &cs : { coordSet::Original, coordSet::Minimized })
    {
        if (!mol->hasCoordinateSet(cs))
        {
            return 0.0;
        }
        std::vector<gmx::RVec> xcs = mol->coordinateSet(cs);
        rvec                   com = { 0, 0, 0 };
        for (size_t i = 0; i < xcs.size(); ++i)
        {
            myMass.push_back(myatoms[i].mass());
            double relativeMass = myatoms[i].mass()/tmass;
            for(int m = 0; m < DIM; m++)
            {
                com[m] += relativeMass*xcs[i][m];
            }
        }
        for (size_t i = 0; i < xcs.size(); ++i)
        {
            for(int m = 0; m < DIM; m++)
            {
                xcs[i][m] -= com[m];
            }
        }
        x->insert({ cs, xcs });
    }
    do_fit(myatoms.size(), myMass.data(),
           as_rvec_array((*x)[coordSet::Original].data()),
           as_rvec_array((*x)[coordSet::Minimized].data()));
    
    double msd  = 0;
    for (size_t i = 0; i < myatoms.size(); i++)
    {
        if (myMass[i] > 0)
        {
            rvec dx;
            rvec_sub((*x)[coordSet::Original][i], (*x)[coordSet::Minimized][i], dx);
            msd += iprod(dx, dx);
        }
    }
    return std::sqrt(msd/mol->nRealAtoms());
}

static void energyLegend(FILE                                    *exvg, 
                         const std::map<InteractionType, double> &energies,
                         const gmx_output_env_t                  *oenv)
{
    std::vector<const char *> enames;
    for(auto &e : energies)
    {
        enames.push_back(interactionTypeToString(e.first).c_str());
    }
    enames.push_back("EKIN");
    enames.push_back("ETOT");
    enames.push_back("TEMPERATURE");
    xvgr_legend(exvg, enames.size(), enames.data(), oenv);
}

static double rvecNormSquared(const gmx::RVec &v)
{
    return v[XX]*v[XX]+v[YY]*v[YY]+v[ZZ]*v[ZZ];
}

static void printEnergies(FILE                                    *exvg, 
                          int                                      step, 
                          double                                   deltat, 
                          const std::map<InteractionType, double> &energies,
                          double                                   ekin,
                          int                                      nDOF,
                          const gmx_output_env_t                  *oenv)
{
    // For energy writing first time around
    if (step == 0)
    {
        energyLegend(exvg, energies, oenv);
    }
    fprintf(exvg, "%10g", step*deltat);
    for(auto &e : energies)
    {
        fprintf(exvg, "  %10g", e.second);
    }
    auto epot = energies.find(InteractionType::EPOT)->second;
    auto temp = ekin*2/(3*nDOF*BOLTZ);
    fprintf(exvg, " %10g %10g %10g\n", ekin, ekin+epot, temp);
}

static void initArrays(const MyMol            *mol,
                       std::vector<gmx::RVec> *coordinates, 
                       std::vector<gmx::RVec> *velocities, 
                       std::vector<gmx::RVec> *forcesCur,
                       std::vector<gmx::RVec> *forcesPrev,
                       std::vector<double>    *inv2mass,
                       std::vector<double>    *mass_2)
{
    auto &xmol = mol->x();
    coordinates->resize(xmol.size());
    velocities->resize(xmol.size());
    forcesPrev->resize(xmol.size());
    forcesCur->resize(xmol.size());
    gmx::RVec zero = { 0, 0, 0 };
    auto &atoms = mol->atomsConst();
    for (long i = 0; i < xmol.size(); i++)
    {
        copy_rvec(xmol[i], (*coordinates)[i]);
        copy_rvec(zero, (*velocities)[i]);
        copy_rvec(zero, (*forcesPrev)[i]);
        copy_rvec(zero, (*forcesCur)[i]);
        mass_2->push_back(atoms[i].mass()*0.5);
        if (atoms[i].mass() > 0)
        {
            inv2mass->push_back(0.5/atoms[i].mass());
        }
        else
        {
            inv2mass->push_back(0);
        }
    }
}
                       
void MolHandler::simulate(MyMol                         *mol,
                          const ForceComputer           *forceComp,
                          const SimulationConfigHandler &simConfig,
                          FILE                          *logFile,
                          const char                    *trajectoryFile,
                          const char                    *energyFile,
                          const gmx_output_env_t        *oenv) const
{
    FILE                              *traj = gmx_ffopen(trajectoryFile, "w");
    FILE                              *exvg = xvgropen(energyFile, "ACT Energies",
                                                       "Time (ps)", "Energy (kJ/mol)", oenv);
    std::vector<gmx::RVec>             coordinates, velocities, forces[2];
    std::map<InteractionType, double>  energies;
    std::vector<double>                inv2mass;
    std::vector<double>                mass_2;

    // To swap between forces efficiently
    int cur = 0;
#define prev (1-cur)

    // Initiate coordinates, etc.
    initArrays(mol, &coordinates, &velocities,
               &forces[cur], &forces[prev], &inv2mass, &mass_2);
               
    // Generate velocities
    if (simConfig.temperature() > 0)
    {
        maxwell_speed(simConfig.temperature(), simConfig.seed(),
                      mol->atomsConst(), &velocities, logFile);
        stop_cm(mol->atomsConst(), coordinates, &velocities);
    }
    auto xmol    = mol->x();
    auto deltat  = simConfig.deltat();
    auto deltat2 = deltat*deltat;
    // For trajectory writing
    matrix      box   = { { 4, 0, 0 }, { 0, 4, 0 }, { 0, 0, 4 } };
    char        chain = ' ';
    const char *title = "ACT trajectory";
    int         nDOF  = DIM*mol->nRealAtoms();
    fprintf(logFile, "\nWill start simulation.\n");
    for (int step = 0; step < simConfig.nsteps(); step++)
    {
        // Check whether we need to print stuff
        bool   doPrintEner = simConfig.nstener() > 0 && step % simConfig.nstener() == 0;
        double ekin        = 0;
        // Integrate positions using velocity Verlet
        for(long i = 0; i < xmol.size(); i++)
        {
            if (inv2mass[i] > 0)
            {
                double dt2m = deltat2*inv2mass[i];
                for(int m = 0; m < DIM; m++)
                {
                    coordinates[i][m] += (deltat*velocities[i][m] + 
                                          dt2m*forces[prev][i][m]);
                }
                if (doPrintEner)
                {
                    ekin += mass_2[i]*rvecNormSquared(velocities[i]);
                }
            }
        }
        // Compute forces and energies
        forceComp->compute(mol->topology(), &coordinates,
                           &forces[cur], &energies);
        // Write energies if requested
        if (doPrintEner)
        {
            printEnergies(exvg, step, deltat, energies, ekin, nDOF, oenv);
        }
        // Integrate velocities using velocity Verlet
        for(long i = 0; i < xmol.size(); i++)
        {
            if (inv2mass[i] > 0)
            {
                for(int m = 0; m < DIM; m++)
                {
                    velocities[i][m] += deltat*inv2mass[i]*(forces[cur][i][m] +
                                                            forces[prev][i][m]);
                }
            }
        }
        // Write output if needed
        if (simConfig.nstxout() > 0 && step % simConfig.nstxout() == 0)
        {
            write_pdbfile(traj, title, mol->gmxAtoms(),
                          as_rvec_array(coordinates.data()), epbcNONE,
                          box, chain, step+1, nullptr, false);
        }
        // Swap force arrays
        cur = prev;
    }
    
    gmx_ffclose(traj);
    xvgrclose(exvg);
    fprintf(logFile, "Simulation finished correctly.\n");
}

} // namespace alexandria

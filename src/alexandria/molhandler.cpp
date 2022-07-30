/*
 * This source file is part of the Alexandria Chemistry Toolkit.
 *
 * Copyright (C) 2022
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
#include "gromacs/math/vec.h"
#include "gromacs/gmxlib/network.h"
#include "gromacs/linearalgebra/eigensolver.h"
#include "gromacs/fileio/xvgr.h"
#include "gromacs/math/do_fit.h"
#include "gromacs/math/vectypes.h"
#include "gromacs/mdtypes/enerdata.h"
#include "gromacs/mdtypes/inputrec.h"
#include "gromacs/topology/topology.h"

#include "Eigen/Eigenvalues"

namespace alexandria
{

std::map<eMinimizeStatus, std::string> eMinStat2String = {
    { eMinimizeStatus::OK,           "OK"             },
    { eMinimizeStatus::TooManySteps, "Too many steps" },
    { eMinimizeStatus::Solver,       "Solver failed"  }
};

const std::string &eMinimizeStatusToString(eMinimizeStatus e)
{
    return eMinStat2String[e];
}

void MolHandler::computeHessian(const MyMol                       *mol,
                                const ForceComputer               *forceComp,
                                std::vector<gmx::RVec>            *coords,
                                const std::vector<int>            &atomIndex,
                                MatrixWrapper                     *hessian,
                                std::vector<double>               *forceZero,
                                std::map<InteractionType, double> *energyZero,
                                std::vector<gmx::RVec>            *dpdq) const
{
    std::vector<double> fneg;
    fneg.resize(DIM*atomIndex.size(), 0.0);
    
    const auto &atoms    = mol->atomsConst();
    std::vector<gmx::RVec> fzero(atoms.size());
    (void) mol->calculateEnergy(forceComp, coords, &fzero, energyZero);
    forceZero->clear();
    for(auto &atom : atomIndex)
    {
        for(int m = 0; m < DIM; m++)
        {
            forceZero->push_back(fzero[atom][m]);
        }
    }
    
    double      stepSize = 1e-8; // 0.001 pm
#define GMX_DOUBLE_EPS 2.2204460492503131e-16
    stepSize = 10.0 * std::sqrt(GMX_DOUBLE_EPS);

    if (dpdq)
    {
        gmx::RVec zero = { 0, 0, 0 };
        dpdq->resize(atomIndex.size(), zero);
    }
    std::vector<gmx::RVec> forces(atoms.size());
    for(size_t ai = 0; ai < atomIndex.size(); ai++)
    {
        gmx::RVec mu[2] = { { 0, 0, 0 }, { 0, 0, 0 } };
        auto atomI = atomIndex[ai];
        for(int atomXYZ = 0; atomXYZ < DIM; atomXYZ++)
        {
            double xyzRef = (*coords)[atomI][atomXYZ];
            for(int delta = 0; delta <= 1; delta++)
            {
                (*coords)[atomI][atomXYZ] = xyzRef + (2*delta-1)*stepSize;
                std::map<InteractionType, double> energies;
                (void) mol->calculateEnergy(forceComp, coords, &forces, &energies);
                if (dpdq)
                {
                    // To compute dipole we need to take shells into account as well!
                    // The same goes for virtual sites.
                    for(size_t ii = 0; ii < atoms.size(); ii++)
                    {
                        for(int m = 0; m < DIM; m++)
                        {
                            mu[delta][m] += atoms[ii].charge()*(*coords)[ii][m];
                        }
                    }
                }
                if (delta == 0)
                {
                    for(size_t aj = 0; aj < atomIndex.size(); aj++)
                    {
                        int atomJ = atomIndex[aj];
                        for (int d = 0; d < DIM; d++)
                        {
                            int row   = aj*DIM+d;
                            fneg[row] = forces[atomJ][d];
                        }
                    }
                }
            }
            (*coords)[atomI][atomXYZ] = xyzRef;  // Restore positions
            
            int col = ai*DIM+atomXYZ;
            for(size_t aj = 0; aj < atomIndex.size(); aj++)
            {
                int    atomJ = atomIndex[aj];
                for(int d = 0; d < DIM; d++)
                {
                    int    row   = aj*DIM+d;
                    double value = -(forces[atomJ][d]-fneg[row])/(2*stepSize);
                    if (false && debug)
                    {
                        fprintf(debug, "Setting H[%2d][%2d] = %10g\n", row, col, value); 
                    }
                    hessian->set(row, col, value);
                }
            }
            if (dpdq)
            {
                for(int d = 0; d < DIM; d++)
                {
                    (*dpdq)[ai][d] += (mu[1][d]-mu[0][d])/(2*stepSize);
                }
            }
        }
    }
    (void) mol->calculateEnergy(forceComp, coords, &fzero, energyZero);
}

static void computeFrequencies(const std::vector<double> &eigenvalues,
                               size_t                     rot_trans,
                               std::vector<double>       *frequencies,
                               std::vector<std::string>  *output)
{
    frequencies->clear();
    auto        mpo  = MolPropObservable::FREQUENCY;
    const char *unit = mpo_unit2(mpo);
    for (size_t i = rot_trans; i < eigenvalues.size(); i++)
    {
        auto val = eigenvalues[i];
        if (val < 0)
        {
            if (output)
            {
                output->push_back(gmx::formatString("Warning: negative eigenvalue %zu = %g\n", i+1, val));
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
}

static void computeIntensities(const std::vector<double>    &eigenvalues,
                               const std::vector<double>    &eigenvectors,
                               const std::vector<int>       &atomIndex,
                               const std::vector<ActAtom>   &atoms,
                               const std::vector<gmx::RVec> &dpdq,
                               std::vector<double>          *intensities)
{
    intensities->clear();
    auto matrixSide = eigenvalues.size();
    for (size_t i = 0; i < eigenvalues.size(); i++)
    {
        double In = 0;
        for (size_t j = 0; j < atomIndex.size(); j++)
        {
            size_t aj      = atomIndex[j];
            double massFac = gmx::invsqrt(atoms[aj].mass());
            for(int d = 0; d < DIM; d++)
            {
                auto evindex = i * matrixSide + j * DIM + d;
                GMX_RELEASE_ASSERT(evindex < eigenvectors.size(), "Range check");
                auto myev = eigenvectors[evindex];
                if (myev != 0)
                {
                    In += gmx::square(dpdq[j][d]*massFac*myev);
                }
                else if (debug)
                {
                    fprintf(debug, "Eigenvector[%zu][%zu][%d] = %g\n",
                            i, j, d, eigenvectors[evindex]);
                }
            } 
        }
        intensities->push_back(In);
    }
}

static void outputFreqInten(const MyMol               *mol,
                            const std::vector<double> &frequencies,
                            const std::vector<double> &intensities,
                            std::vector<std::string>  *output)
{
    std::vector<GenericProperty *> harm;
    auto mpo = MolPropObservable::FREQUENCY;
    const char *unit = mpo_unit2(mpo);
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
        output->push_back(gmx::formatString("Electronic vibrational frequencies: (%s):", unit));
        size_t k = 0;
        std::string freq;
        for(auto &ff : harm[0]->getVector())
        {
            freq += gmx::formatString("  %10g", convertFromGromacs(ff, unit));
            delta += gmx::square(ff - frequencies[k]);
            k++;
        }
        output->push_back(freq);
    }
    output->push_back(gmx::formatString("Alexandria vibrational frequencies: (%s):", unit));
    std::string actfreq;
    for (const auto &freq : frequencies)
    {
        actfreq += gmx::formatString("  %10g", convertFromGromacs(freq, unit));
    }
    output->push_back(actfreq);
    // Now print intensities
    auto mpoi = MolPropObservable::INTENSITY;
    const char *uniti = mpo_unit2(mpoi);
    output->push_back(gmx::formatString("Alexandria IR intensities: (%s):", uniti));
    std::string actinten;
    for (const auto &inten : intensities)
    {
        actinten += gmx::formatString("  %10g", convertFromGromacs(inten, uniti));
    }
    output->push_back(actinten);
    if (delta > 0)
    {
        output->push_back(gmx::formatString("Frequency RMSD %g",
                                            convertFromGromacs(std::sqrt(delta/frequencies.size()),
                                                               unit)));
    }
}

static void solveEigen(const MatrixWrapper          &hessian,
                       const std::vector<int>       &atomIndex,
                       const std::vector<ActAtom>   &atoms,
                       const std::vector<gmx::RVec> &dpdq,
                       int                           rot_trans,
                       std::vector<double>          *frequencies,
                       std::vector<double>          *intensities,
                       std::vector<std::string>     *output,
                       bool                          debugNMA)
{
    if (output)
    {
        output->push_back("Using Eigen to solve eigenvalue problem.");
    }
    int matrixSide = hessian.nColumn();
    Eigen::MatrixXd mat(matrixSide, matrixSide);
    for(int col = 0; col < matrixSide; col++)
    {
        for(int row = 0; row < matrixSide; row++)
        {
            mat(col, row) = hessian.get(col, row);
        }
    }
    bool computeEigenVectors = true;
    Eigen::EigenSolver<Eigen::MatrixXd> solver(mat, computeEigenVectors);
    auto evecs = solver.eigenvectors();
    if (output)
    {
        for(Eigen::Index i =  0; i < matrixSide; i++)
        {
            double sumReal = 0, sumComplex = 0;
            auto evv = evecs.col(i);
            std::string evi;
            for(Eigen::Index j = 0; j < matrixSide; j++)
            {
                auto evec   = evv[j];
                sumReal    += gmx::square(evec.real());
                sumComplex += gmx::square(std::abs(evec));
                evi += gmx::formatString("  %10g", evec.real());
            }
            output->push_back(gmx::formatString("|v[%2zu]|^2 = %8g (real) %8g (complex) Eigenval = %g + i%g.",
                                                i, sumReal, sumComplex,
                                                solver.eigenvalues()[i].real(),
                                                solver.eigenvalues()[i].imag()));
        }
        if (debugNMA)
        {
            for(Eigen::Index i =  0; i < matrixSide; i++)
            {
                auto evv = evecs.col(i);
                std::string evi;
                for(Eigen::Index j = 0; j < matrixSide; j++)
                {
                    evi += gmx::formatString("  %10g", evv[j].real());
                }
                output->push_back(gmx::formatString("EV%ld %s", i, evi.c_str()));
            }
        }
    }

    // Eigen does not sort the output according to increasing eigenvalue
    // so we have to do the sorting ourselves.
    typedef struct
    {
        Eigen::Index index;
        double       value;
    } iv_t;
    std::vector<iv_t> iv;
    auto &evals    = solver.eigenvalues();
    int nnegative = 0;
    for(Eigen::Index i = 0; i < matrixSide; i++)
    {
        // Both eigenvalues and eigenvectors can be complex numbers
        // If we take the real part only, we reproduce the LAPACK
        // result, but is this correct?
        auto eval = evals[i].real();
        //auto eval = std::abs(evals[i]);
        iv.push_back({ i, eval });
        if (eval < 0)
        {
            nnegative += 1;
        }
    }
    if (output && nnegative > rot_trans)
    {
        output->push_back(gmx::formatString("There are %d negative eigenvalues.",
                                            nnegative-rot_trans));
    }
    std::sort(iv.begin(), iv.end(), [](const iv_t &a, const iv_t &b)
              { return a.value < b.value; });

    std::vector<double> In(matrixSide, 0);
    for(Eigen::Index col = 0; col < matrixSide; col++)
    {
        auto   icol = iv[col].index;
        auto   evec = evecs.col(icol);
        for (int row = 0; row < matrixSide; row++)
        {
            auto irow = iv[row].index;
            // Eigenvectors are complex numbers too. The eigenvectors
            // should be normalized (length 1). This is indeed the case
            // if we take the imaginary part into account, but whether
            // this gives the correct physical result is unclear.
            // auto myev = std::abs(evec[irow]);
            auto myev  = evec[irow].real();
            auto atomI = irow / DIM;
            int k      = irow % DIM;
            size_t ai  = atomIndex[atomI];
            // The components of the root-mass weighted Hessian
            // have unit kJ/(Da mol nm^2). Even though the eigenvectors
            // are normalized, they have the same unit. Since kJ =
            // 1000 kg m^2/s^2 the unit reduces to
            // 1000 * 1000 * 1e18/s^2 = 1e24/s^2
            // dpdq has unit of e
            // so the output unit of the intensities here is e^2 =
            // nm^2 g/mol. To convert this to km/mol we have to assume
            // a line width for the spectrum, a typical width is
            // 24 cm^-1. By integrating over the width of the peak, we
            // get the units (almost) correct. nm^2 g / cm mol. This means
            // (10^-18/10^-2) m g / mol = 10^-16 m g / mol or
            // 10^-13 g km/mol. Not clear yet how to get rid of the gram.
            In[icol] += gmx::square(dpdq[atomI][k]*myev*
                                    gmx::invsqrt(atoms[ai].mass()));
        }
    }
    intensities->clear();
    for(int i = rot_trans; i < matrixSide; i++)
    {
        intensities->push_back(In[i]);
    }
    std::vector<double> eigenvalues(matrixSide, 0.0);
    for(int col = 0; col < matrixSide; col++)
    {
        eigenvalues[col]    = iv[col].value;
    }
    computeFrequencies(eigenvalues, rot_trans, frequencies, output);
}

static void solveLapack(const MyMol                  *mol,
                        const MatrixWrapper          &hessian,
                        const std::vector<int>       &atomIndex,
                        const std::vector<ActAtom>   &atoms,
                        const std::vector<gmx::RVec> &dpdq,
                        int                           rot_trans,
                        std::vector<double>          *frequencies,
                        std::vector<double>          *intensities,
                        std::vector<std::string>     *output)
{
    if (output)
    {
        output->push_back("Using LAPACK to solve eigenvalue problem.");
    }
    int matrixSide = hessian.nColumn();
    std::vector<double> eigenvalues(matrixSide);
    std::vector<double> eigenvectors(matrixSide*matrixSide);
    auto                hessianFlat = hessian.flatten();
    eigensolver(hessianFlat.data(), matrixSide, 0, matrixSide - 1,
                eigenvalues.data(), eigenvectors.data());
    // Scale the output eigenvectors
    for (int i = 0; i < matrixSide; i++)
    {
        for (size_t j = 0; j < atomIndex.size(); j++)
        {
            size_t aj      = atomIndex[j];
            double massFac = gmx::invsqrt(atoms[aj].mass());
            for (size_t k = 0; (k < DIM); k++)
            {
                eigenvectors[i * matrixSide + j * DIM + k] *= massFac;
            }
        }
    }
    computeIntensities(eigenvalues, eigenvectors, atomIndex,
                       atoms, dpdq, intensities);
    computeFrequencies(eigenvalues, rot_trans, frequencies, output);
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
        // Get the eigenvectors into a MatrixWrapper
        MatrixWrapper eigenvecMat(eigenvectors, matrixSide);
        
        fprintf(debug, "]\nHessian eigenvectors:\n%s\n", eigenvecMat.toString().c_str());
    }
}

void MolHandler::nma(const MyMol              *mol,
                     const ForceComputer      *forceComp,
                     std::vector<gmx::RVec>   *coords,
                     std::vector<double>      *frequencies,
                     std::vector<double>      *intensities,
                     std::vector<std::string> *output,
                     bool                      useLapack,
                     bool                      debugNMA) const
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
    const int              matrixSide = DIM*atomIndex.size();
    MatrixWrapper          hessian(matrixSide, matrixSide);
    std::vector<double>    f0;
    std::vector<gmx::RVec> dpdq;
    std::map<InteractionType, double> energies;
    computeHessian(mol, forceComp, coords, atomIndex, &hessian, &f0, &energies, &dpdq);
    hessian.averageTriangle();

    if (output && debugNMA)
    {
        output->push_back("Symmetrized Hessian:");
        output->push_back(hessian.toString());
        for(size_t i = 0; i < atomIndex.size(); i++)
        {
            output->push_back(gmx::formatString("dpdq[%2zu] =  %10g  %10g  %10g",
                                                i, dpdq[i][XX], dpdq[i][YY], dpdq[i][ZZ]));
        }
    }
    // Divide the elements hessian[i][j] by sqrt(mass[i])*sqrt(mass[j])
    for (size_t i = 0; (i < atomIndex.size()); i++)
    {
        size_t ai = atomIndex[i];
        for (size_t j = 0; (j < DIM); j++)
        {
            for (size_t k = 0; (k < atomIndex.size()); k++)
            {
                size_t ak      = atomIndex[k];
                double massFac = gmx::invsqrt(atoms[ai].mass() * atoms[ak].mass());
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
    if (output)
    {
        output->push_back("Diagonalizing Hessian to find eigenvectors.");
    }
    if (useLapack)
    {
        solveLapack(mol, hessian, atomIndex, atoms, dpdq, rot_trans,
                    frequencies, intensities, output);
    }
    else
    {
        solveEigen(hessian, atomIndex, atoms, dpdq, rot_trans,
                   frequencies, intensities, output, debugNMA);
    }

    if (output)
    {
        outputFreqInten(mol, *frequencies, *intensities, output);
    }

}

static void printEnergies(FILE *logFile, int myIter, double msAtomForce,
                          const std::map<InteractionType, double> &energies)
{
    if (nullptr == logFile)
    {
        return;
    }
    fprintf(logFile, "Iter %5d rmsForce %10g", myIter, std::sqrt(msAtomForce));
    for(const auto &ee : energies)
    {
        fprintf(logFile, "  %s %8g", interactionTypeToString(ee.first).c_str(), ee.second);
    }
    fprintf(logFile, "\n");
}

static double msForce(const std::vector<int>       &theAtoms,
                      const std::vector<gmx::RVec> &forces)
{
    double msAtomForce  = 0;
    for(size_t kk = 0; kk < theAtoms.size(); kk++)
    {
        int atomI = theAtoms[kk];
        msAtomForce  += iprod(forces[atomI], forces[atomI]);
    }
    return msAtomForce /= theAtoms.size();
}

eMinimizeStatus MolHandler::minimizeCoordinates(const MyMol                       *mol,
                                                const ForceComputer               *forceComp,
                                                const SimulationConfigHandler     &simConfig,
                                                std::vector<gmx::RVec>            *coords,
                                                std::map<InteractionType, double> *energies,
                                                FILE                              *logFile) const
{
    bool      converged    = false;
    int       myIter       = 0;
    // Check for meaningful convergence tolerance
    double    shellToler2  = gmx::square(forceComp->convergenceTolerance());
    double    msForceToler = simConfig.forceTolerance();
    if (msForceToler < shellToler2)
    {
        if (logFile && msForceToler != 0)
        {
            fprintf(logFile, "Atom mean square force tolerance (%g) more strict than shell force tolerance (%g).\n",
                    msForceToler, shellToler2);
            fprintf(logFile, "Increasing atom mean square force tolerance to hundred times that for shells.\n");
        }
        msForceToler = 100*shellToler2;
    }
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
    std::vector<gmx::RVec> forces(myatoms.size());
    std::vector<double>    f0, f00;
    bool                   firstStep = true;
    // Now start the minimization loop.
    if (logFile)
    {
        fprintf(logFile, "Starting minimization of '%s' using %s algorithm.\n",
                mol->getMolname().c_str(),
                eMinimizeAlgorithmToString(simConfig.minAlg()).c_str());
    }
    std::vector<double> deltaX(DIM*theAtoms.size(), 0.0);
    // Set a maximum displacement to prevent exploding molecules
    double              deltaXTolerance = 0.002; // nm
    double              epotMin = 0;
    // Two sets of coordinates
    std::vector<gmx::RVec> newCoords[2];
    newCoords[0] = *coords;
    newCoords[1] = *coords;
    std::map<InteractionType, double> newEnergies[2];
    int current = 0;
#define next (1-current)
    do
    {
        auto eMin = eMinimizeStatus::OK;
        switch (simConfig.minAlg())
        {
        case eMinimizeAlgorithm::Newton:
            {
                // Below is a Newton-Rhapson algorithm
                // https://en.wikipedia.org/wiki/Newton%27s_method_in_optimization
                // Create Hessian, just containings atoms, not shells
                MatrixWrapper Hessian(DIM*theAtoms.size(), DIM*theAtoms.size());
                computeHessian(mol, forceComp, &newCoords[current],
                               theAtoms, &Hessian, &f0, &newEnergies[current]);
                // Solve H delta X = -grad (E) = force(E)
                int result = Hessian.solve(f0, &deltaX);
                if (0 != result)
                {
                    eMin = eMinimizeStatus::Solver;
                }
            }
            break;
        case eMinimizeAlgorithm::Steep:
            {
                (void) mol->calculateEnergy(forceComp, &newCoords[current],
                                            &forces, &newEnergies[current]);
                int    i      = 0;
                double factor = 0.001;
                for(auto atomI : theAtoms)
                {
                    for(int m = 0; m < DIM; m++)
                    {
                        deltaX[i++] = factor*forces[atomI][m];
                    }
                }
            }
            break;
        }
        if (firstStep)
        {
            // Store energy
            if (energies)
            {
                *energies = newEnergies[current];
            }
            epotMin   = newEnergies[current][InteractionType::EPOT];
            // Store forces
            f00       = f0;
            // Write stuff
            double msf = 0;
            for (const auto &f : f0)
            {
                msf += f*f;
            }
            printEnergies(logFile, myIter, (msf/theAtoms.size()), newEnergies[current]);
            firstStep = false;
        }
        if (eMinimizeStatus::OK != eMin)
        {
            return eMin;
        }
        bool acceptStep = false;
        do
        {
            // Update coordinates and check energy
            double maxDeltaX       = *std::max_element(deltaX.begin(), deltaX.end());
            double scaleDeltaX     = simConfig.overRelax();
            if (maxDeltaX > deltaXTolerance)
            {
                scaleDeltaX = deltaXTolerance/maxDeltaX;
            }
            int i = 0;
            for (auto &atomI : theAtoms)
            {
                for (int m = 0; m < DIM; m++)
                {
                    newCoords[next][atomI][m] = newCoords[current][atomI][m] + scaleDeltaX*deltaX[i++];
                }
            }
            
            // Do an energy and force calculation with the new coordinates
            (void) mol->calculateEnergy(forceComp, &newCoords[next], &forces, &newEnergies[next]);
            acceptStep = (newEnergies[next][InteractionType::EPOT] <= epotMin);
            if (!acceptStep)
            {
                deltaXTolerance *= 0.8;
            }
            else
            {
                deltaXTolerance *= 1.25;
                epotMin = newEnergies[next][InteractionType::EPOT];
                current = next;
            }
            myIter += 1;
        } while (!acceptStep && (myIter < simConfig.maxIter() || 0 == simConfig.maxIter()));
        
        // Now check force convergencs
        double msAtomForce  = msForce(theAtoms, forces);
        converged = msAtomForce <= msForceToler;
        
        printEnergies(logFile, myIter, msAtomForce, newEnergies[current]);
        if (debug)
        {
            for(size_t kk = 0; kk < mol->topology()->nAtoms(); kk++)
            {
                fprintf(debug, "f[%2zu] =  %10g  %10g  %10g x[%2zu] = %10g  %10g  %10g\n", kk,
                        forces[kk][XX], forces[kk][YY], forces[kk][ZZ], kk,
                        (*coords)[kk][XX], (*coords)[kk][YY], (*coords)[kk][ZZ]);
            }
            for(const auto &ee : newEnergies[current])
            {
                fprintf(debug, "%-20s  %10g\n",
                        interactionTypeToString(ee.first).c_str(), 
                        ee.second);
            }
        }
    }
    while (!converged && (myIter < simConfig.maxIter() || 0 == simConfig.maxIter()));
    if (converged)
    {
        *coords   = newCoords[current];
        if (energies && !newEnergies[current].empty())
        {
            *energies = newEnergies[current];
        }
        // Re-compute the energy one last time.
        // TODO: is this really needed?
        // (void) mol->calculateEnergy(forceComp, coords, &forces, energies);
        return eMinimizeStatus::OK;
    }
    else
    {
        return eMinimizeStatus::TooManySteps;
    }
}

double MolHandler::coordinateRmsd(const MyMol                  *mol,
                                  const std::vector<gmx::RVec> &xref,
                                  std::vector<gmx::RVec>       *xfit) const
{
    // Compute RMSD
    auto   myatoms = mol->atomsConst();
    double tmass   = mol->totalMass();
    std::vector<real> myMass;
    rvec ref_com = { 0, 0, 0 };
    rvec fit_com = { 0, 0, 0 };
    for (size_t i = 0; i < xref.size(); ++i)
    {
        myMass.push_back(myatoms[i].mass());
        double relativeMass = myatoms[i].mass()/tmass;
        for(int m = 0; m < DIM; m++)
        {
            ref_com[m] += relativeMass*xref[i][m];
            fit_com[m] += relativeMass*(*xfit)[i][m];
        }
    }
    std::vector<gmx::RVec> xref_com(xref.size());
    for (size_t i = 0; i < xref.size(); ++i)
    {
        rvec_sub(xref[i], ref_com, xref_com[i]);
        rvec_dec((*xfit)[i], fit_com);
    }
    do_fit(myatoms.size(), myMass.data(),
           as_rvec_array(xref_com.data()),
           as_rvec_array(xfit->data()));
    
    double msd  = 0;
    for (size_t i = 0; i < myMass.size(); i++)
    {
        if (myMass[i] > 0)
        {
            rvec dx;
            rvec_sub(xref_com[i], (*xfit)[i], dx);
            msd += iprod(dx, dx);
        }
        rvec_inc((*xfit)[i], fit_com);
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
    auto &xmol = mol->xOriginal();
    coordinates->resize(xmol.size());
    velocities->resize(xmol.size());
    forcesPrev->resize(xmol.size());
    forcesCur->resize(xmol.size());
    gmx::RVec zero = { 0, 0, 0 };
    auto &atoms = mol->atomsConst();
    for (size_t i = 0; i < xmol.size(); i++)
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
    auto xmol    = mol->xOriginal();
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
        for(size_t i = 0; i < xmol.size(); i++)
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
        for(size_t i = 0; i < xmol.size(); i++)
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

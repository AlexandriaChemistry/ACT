/*
 * This source file is part of the Alexandria Chemistry Toolkit.
 *
 * Copyright (C) 2022-2024
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

// This include has to come first, to prevent the GROMACS definition of
// "real" from messing up the library.
#pragma clang diagnostic ignored "-Wdeprecated-copy-with-dtor"
#pragma GCC diagnostic ignored "-Wdeprecated-copy-dtor"
#include "Eigen/Eigenvalues"
#pragma GCC diagnostic pop
#pragma clang diagnostic pop

#include "molhandler.h"

#include <numeric>
#include <random>

#include "act/alexandria/pdbwriter.h"
#include "act/alexandria/velocityhandler.h"
#include "act/molprop/molpropobservable.h"
#include "act/utility/units.h"
#include "external/stlbfgs/stlbfgs.h"
#include "gromacs/fileio/xvgr.h"
#include "gromacs/math/do_fit.h"
#include "gromacs/math/units.h"
#include "gromacs/math/vec.h"
#include "gromacs/math/vectypes.h"
#include "gromacs/mdtypes/enerdata.h"
#include "gromacs/mdtypes/inputrec.h"
#include "gromacs/utility/futil.h"

namespace alexandria
{

std::map<eMinimizeStatus, std::string> eMinStat2String = {
    { eMinimizeStatus::OK,           "OK"             },
    { eMinimizeStatus::TooManySteps, "Too many steps" },
    { eMinimizeStatus::Solver,       "Solver failed"  },
    { eMinimizeStatus::NoMinimum,    "There is no minimum"  }
};

const std::string &eMinimizeStatusToString(eMinimizeStatus e)
{
    return eMinStat2String[e];
}

void MolHandler::computeHessian(const ForceField                  *pd,
                                const ACTMol                      *mol,
                                const ForceComputer               *forceComp,
                                std::vector<gmx::RVec>            *coords,
                                const std::vector<int>            &atomIndex,
                                MatrixWrapper                     *hessian,
                                std::vector<double>               *forceZero,
                                std::map<InteractionType, double> *energyZero,
                                std::vector<gmx::RVec>            *dpdq) const
{
    const auto &atoms    = mol->atomsConst();
    std::vector<gmx::RVec> fzero(atoms.size());
    (void) forceComp->compute(pd, mol->topology(), coords, &fzero, energyZero);
    forceZero->clear();
    for(auto &atom : atomIndex)
    {
        for(int m = 0; m < DIM; m++)
        {
            forceZero->push_back(fzero[atom][m]);
        }
    }
    
    double      stepSize = 1e-12; // nm
#define GMX_DOUBLE_EPS 2.2204460492503131e-16
    //    stepSize = 1.0 * std::sqrt(GMX_DOUBLE_EPS);

    if (dpdq)
    {
        gmx::RVec zero = { 0, 0, 0 };
        dpdq->resize(atomIndex.size(), zero);
    }
    std::vector<gmx::RVec> forces[2];
    forces[0].resize(atoms.size());
    forces[1].resize(atoms.size());
    for(size_t ai = 0; ai < atomIndex.size(); ai++)
    {
        gmx::RVec mu[2] = { { 0, 0, 0 }, { 0, 0, 0 } };
        auto atomI = atomIndex[ai];
        for(int atomXYZ = 0; atomXYZ < DIM; atomXYZ++)
        {
            int    column = ai*DIM+atomXYZ;
            double xyzRef = (*coords)[atomI][atomXYZ];
            double xxx[2];
            bool   stepsDone = false;
            double myStep    = stepSize;
            do {
                for(int delta = 0; delta <= 1; delta++)
                {
                    xxx[delta] = xyzRef + (2*delta-1)*myStep;
                }
                // Make sure that the coordinate is in fact chenged.
                // If the delta is too small, numeric issues may kick
                // in.
                stepsDone = (xxx[0] < xyzRef && xxx[1] > xyzRef);
                if (!stepsDone)
                {
                    myStep *= 2;
                }
            }
            while (!stepsDone);
            for(int delta = 0; delta <= 1; delta++)
            {
                (*coords)[atomI][atomXYZ] = xxx[delta];
                std::map<InteractionType, double> energies;
                (void) forceComp->compute(pd, mol->topology(), coords, &forces[delta], &energies);

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
            }
            // Restore position
            (*coords)[atomI][atomXYZ] = xyzRef;
            
            for(size_t aj = 0; aj < atomIndex.size(); aj++)
            {
                int    atomJ = atomIndex[aj];
                for(int d = 0; d < DIM; d++)
                {
                    int    row   = aj*DIM+d;
                    double value = -(forces[1][atomJ][d]-forces[0][atomJ][d])/(xxx[1]-xxx[0]);
                    if (false && debug)
                    {
                        fprintf(debug, "Setting H[%2d][%2d] = %10g\n", row, column, value); 
                    }
                    hessian->set(column, row, value);
                }
            }
            if (dpdq)
            {
                for(int d = 0; d < DIM; d++)
                {
                    (*dpdq)[ai][d] += (mu[1][d]-mu[0][d])/(xxx[1]-xxx[0]);
                }
            }
        }
    }
    (void) forceComp->compute(pd, mol->topology(), coords, &fzero, energyZero);
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

static void outputFreqInten(const ACTMol               *mol,
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
    // We need this special solver for hermitian matrices.
    // It will automatically sort the eigenvalues.
    Eigen::SelfAdjointEigenSolver<Eigen::MatrixXd> solver(mat);
    auto evecs = solver.eigenvectors();
    if (output)
    {
        for(Eigen::Index i =  0; i < matrixSide; i++)
        {
            double sum = 0;
            auto evv = evecs.col(i);
            std::string evi;
            for(Eigen::Index j = 0; j < matrixSide; j++)
            {
                auto evec = evv[j];
                sum       += gmx::square(evec);
                evi += gmx::formatString("  %10g", evec);
            }
            output->push_back(gmx::formatString("|v[%2zu]|^2 = %8g Eigenval = %g",
                                                i, sum, solver.eigenvalues()[i]));
        }
        if (debugNMA)
        {
            for(Eigen::Index i =  0; i < matrixSide; i++)
            {
                auto evv = evecs.col(i);
                std::string evi;
                for(Eigen::Index j = 0; j < matrixSide; j++)
                {
                    evi += gmx::formatString("  %10g", evv[j]);
                }
                output->push_back(gmx::formatString("EV%ld %s", i, evi.c_str()));
            }
        }
    }

    auto &evals     = solver.eigenvalues();
    int   nnegative = 0;
    for(Eigen::Index i = 0; i < matrixSide; i++)
    {
        if (evals[i] < 0)
        {
            nnegative += 1;
        }
    }
    if (output && nnegative > rot_trans)
    {
        output->push_back(gmx::formatString("There are %d negative eigenvalues.",
                                            nnegative-rot_trans));
    }

    std::vector<double> In(matrixSide, 0);
    for(Eigen::Index col = 0; col < matrixSide; col++)
    {
        auto evec = evecs.col(col);
        for (Eigen::Index row = 0; row < matrixSide; row++)
        {
            // The eigenvectors should be normalized (length 1).
            // This is indeed the case.
            auto myev  = evec[row];
            auto atomI = row / DIM;
            int k      = row % DIM;
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
            In[col] += gmx::square(dpdq[atomI][k]*myev*
                                   gmx::invsqrt(atoms[ai].mass()));
        }
    }
    intensities->clear();
    for(int i = rot_trans; i < matrixSide; i++)
    {
        intensities->push_back(In[i]);
    }
    std::vector<double> eigenvalues;
    for(int i = 0; i < matrixSide; i++)
    {
        eigenvalues.push_back(evals(i));
    }
    computeFrequencies(eigenvalues, rot_trans, frequencies, output);
}

static void resortFreqIntens(std::vector<double> *frequencies, 
                             std::vector<double> *intensities,
                             double               freq_toler)
{
    std::vector<double> freq  = *frequencies;
    std::vector<double> inten = *intensities;
    std::vector<int> renumber(freq.size(), -1);
    size_t i       = 0;
    bool   swapped = false;
    while (i < freq.size())
    {
        if (i < freq.size()-1 &&
            abs(freq[i] - freq[i+1]) < freq_toler &&
            inten[i+1] > inten[i])
        {
            renumber[i]   = i+1;
            renumber[i+1] = i;
            i += 2;
            swapped = true;
        }
        else
        {
            renumber[i]   = i;
            i += 1;
        }
    }
    if (swapped)
    {
        for (size_t i = 0; i < renumber.size(); i++)
        {
            (*frequencies)[i] = freq[renumber[i]];
            (*intensities)[i] = inten[renumber[i]];
        }
    }
}

void MolHandler::nma(const ForceField         *pd,
                     const ACTMol             *mol,
                     const ForceComputer      *forceComp,
                     std::vector<gmx::RVec>   *coords,
                     std::vector<double>      *frequencies,
                     std::vector<double>      *intensities,
                     std::vector<std::string> *output,
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
    computeHessian(pd, mol, forceComp, coords, atomIndex, &hessian, &f0, &energies, &dpdq);
    //hessian.averageTriangle();

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
    solveEigen(hessian, atomIndex, atoms, dpdq, rot_trans,
               frequencies, intensities, output, debugNMA);
    // Check whether there are very similar frequencies, then change the sorting 
    // according to intensities.
    // TODO: make ftoler a parameter.
    double ftoler = 0.001; // Frequency unit is internal unit.
    resortFreqIntens(frequencies, intensities, ftoler);

    if (output)
    {
        outputFreqInten(mol, *frequencies, *intensities, output);
    }

}

static void printEnergies(FILE *logFile, int myIter,
                          double msAtomForce,
                          double msShellForce,
                          double fcurprev,
                          const std::map<InteractionType, double> &energies,
                          double gamma)
{
    if (nullptr == logFile)
    {
        return;
    }
    fprintf(logFile, "Iter %5d msAtomForce %10g msShellForce %10g fcur.fprev %10g gamma %10g\n    ",
            myIter, msAtomForce, msShellForce, fcurprev, gamma);
    double evdw = 0;
    for(const auto &ee : energies)
    {
        fprintf(logFile, "  %s %.5f", interactionTypeToString(ee.first).c_str(), ee.second);
        if (InteractionType::DISPERSION == ee.first || 
            InteractionType::EXCHANGE == ee.first)
        {
            evdw += ee.second;
        }
    }
    fprintf(logFile," (VANDERWAALS %.5f)\n", evdw);
}

static double msForce(const std::vector<int>       &theAtoms,
                      const std::vector<gmx::RVec> &forces)
{
    double msAtomForce  = 0;
    for(auto kk : theAtoms)
    {
        msAtomForce += iprod(forces[kk], forces[kk]);
    }
    if (theAtoms.size() > 0)
    {
        return msAtomForce /= theAtoms.size();
    }
    else
    {
        return 0;
    }
}

static void updateCoords(const std::vector<int>       &theAtoms,
                         const std::vector<gmx::RVec> &xold,
                         std::vector<gmx::RVec>       *xnew,
                         double                        factor,
                         const std::vector<double>    &deltaX)
{
    int i = 0;
    for (auto &atomI : theAtoms)
    {
        for (int m = 0; m < DIM; m++)
        {
            (*xnew)[atomI][m] = xold[atomI][m] + factor*deltaX[i++];
        }
    }
}

/*! \brief Local class to provide data to the L-BFGS algorithm
 * Content is just a copy of inputs, respectively working variables.
 */
class StlbfgsHandler
{
private:
    const ForceField                  *pd_;
    const ACTMol                      *mol_;
    const ForceComputer               *forceComp_;
    std::vector<gmx::RVec>             coords_;
    std::vector<gmx::RVec>             forces_;
    std::map<InteractionType, double>  energies_;
    const std::vector<int>             theAtoms_;
public:
    StlbfgsHandler(const ForceField       *pd,
                   const ACTMol           *mol,
                   const ForceComputer    *forceComp,
                   const std::vector<int> &theAtoms) : 
        pd_(pd), mol_(mol), forceComp_(forceComp), theAtoms_(theAtoms)
    {
        gmx::RVec vzero = { 0, 0, 0 };
        coords_ = mol_->xOriginal();
        forces_.resize(mol_->atomsConst().size(), vzero);
    }

    const ACTMol *mol() const { return mol_; }

    const ForceField *pd() const { return pd_; }

    const ForceComputer *forceComp() const { return forceComp_; }

    std::vector<gmx::RVec> *coordinates() { return &coords_; }

    std::vector<gmx::RVec> *forces() { return &forces_; }

    std::map<InteractionType, double> *energies() { return &energies_; }

    double energy(InteractionType itype) { return energies_[itype]; }

    const std::vector<int> &theAtoms() const { return theAtoms_; }
};

// The return of the son of global variables. Easy but ugly.
static StlbfgsHandler *lbfgs = nullptr;

static void func(const std::vector<double> &x, double &f, std::vector<double> &g)
{
    auto theAtoms = lbfgs->theAtoms();
    if (x.size() != DIM*theAtoms.size())
    {
        GMX_THROW(gmx::InternalError(gmx::formatString("Number of coordinates %lu does not match number of real atoms %lu",
                                                       x.size(), theAtoms.size()).c_str()));
    }
    auto coords = lbfgs->coordinates();
    auto forces = lbfgs->forces();
    size_t jj   = 0;
    for (auto i : theAtoms)
    {
        for(int j = 0; j < DIM; j++)
        {
            (*coords)[i][j] = x[jj++];
        }
    }
    gmx::RVec field = { 0, 0, 0 };
    bool minimizeShells = false;
    lbfgs->forceComp()->compute(lbfgs->pd(), lbfgs->mol()->topology(), coords,
                                forces, lbfgs->energies(), field,
                                minimizeShells);
    f  = lbfgs->energy(InteractionType::EPOT);
    jj = 0;
    for (auto i : theAtoms)
    {
        for(int j = 0; j < DIM; j++)
        {
            g[jj++] = -(*forces)[i][j];
        }
    }
}

eMinimizeStatus MolHandler::minimizeCoordinates(const ForceField                  *pd,
                                                const ACTMol                      *mol,
                                                const ForceComputer               *forceComp,
                                                const SimulationConfigHandler     &simConfig,
                                                std::vector<gmx::RVec>            *coords,
                                                std::map<InteractionType, double> *energies,
                                                FILE                              *logFile,
                                                const std::vector<int>            &freeze,
                                                double                            *rmsForce) const
{
    bool      converged    = false;
    int       myIter       = 0;
    // Check for meaningful convergence tolerance
    double    shellToler2  = gmx::square(forceComp->forceTolerance());
    double    msForceToler = simConfig.forceTolerance();
    if (msForceToler < shellToler2)
    {
        if (logFile && msForceToler != 0)
        {
            fprintf(logFile, "Atom mean square force tolerance (%g) more strict than shell force tolerance (%g).\n",
                    msForceToler, shellToler2);
            fprintf(logFile, "Increasing atom mean square force tolerance to hundred times that for shells.\n");
        }
        msForceToler = shellToler2;
    }
    // List of atoms and shells but not vsites and weighting factors
    auto              myatoms = mol->atomsConst();
    std::vector<int>  theAtoms;
    for(size_t atom = 0; atom < myatoms.size(); atom++)
    {
        // Store the atom numbers for the frozen atoms in the counting incl. shells
        if ((myatoms[atom].pType() == eptAtom ||
             myatoms[atom].pType() == eptShell) &&
            freeze.end() == std::find(freeze.begin(), freeze.end(), atom))
        {
            theAtoms.push_back(atom);
        }
    }
    // Two sets of forces
    std::vector<gmx::RVec> forces[2];
    forces[0].resize(myatoms.size());
    forces[1].resize(myatoms.size());
    bool                   firstStep = true;
    // Now start the minimization loop.
    if (logFile)
    {
        fprintf(logFile, "Starting minimization of '%s' using %s algorithm.\n",
                mol->getMolname().c_str(),
                eMinimizeAlgorithmToString(simConfig.minAlg()).c_str());
    }
    std::vector<double> deltaX[2];
    deltaX[0].resize(DIM*theAtoms.size(), 0.0);
    deltaX[1].resize(DIM*theAtoms.size(), 0.0);
    double              epotMin = 1e8;
    double              msfMin  = 1e16;
    // Two sets of coordinates
    std::vector<gmx::RVec> newCoords[2];
    newCoords[0] = *coords;
    newCoords[1] = *coords;
    std::map<InteractionType, double> newEnergies[2];
    double gamma_max    = 1e-5;
    double gamma        = gamma_max;
    double msShellForce = 0;
    int current = 0;
#define next (1-current)
    if (simConfig.minAlg() == eMinimizeAlgorithm::LBFGS)
    {
        lbfgs = new StlbfgsHandler(pd, mol, forceComp, theAtoms);
        STLBFGS::Optimizer                      opt{func};
        if (logFile)
        {
            opt.setVerbose();
        }
        opt.setFtol(gmx::square(msForceToler));
        opt.setGtol(msForceToler);
        // One-dimensional array
        std::vector<double>                     sx(theAtoms.size()*DIM);
        std::random_device                      rd;
        std::mt19937                            gen(rd());
        std::uniform_real_distribution<double>  dis(std::uniform_real_distribution<>(-1.0, 1.0));

        double displacement = simConfig.minimizeDisplacement();
        int maxRetry        = simConfig.minimizeRetries();
        int retry           = 0;
        // Large number!
        double eMin         = 1e8;
        while (retry < maxRetry && (!converged || simConfig.forceReminimize()))
        {
            if (logFile && retry > 0)
            {
                fprintf(logFile, "Will retry minimization with slightly modified input coordinates.\n");
            }
            for(size_t i = 0; i < theAtoms.size(); i++)
            {
                for(int j = 0; j < DIM; j++)
                {
                    sx[DIM*i+j] = (*coords)[theAtoms[i]][j] + retry * dis(gen) * displacement;
                }
            }
            converged  = opt.run(sx);
            double msf = 0;
            for(const auto &f : *lbfgs->forces())
            {
                msf += iprod(f, f);
            }
            *rmsForce = std::sqrt(msf/lbfgs->forces()->size());
            if (debug)
            {
                fprintf(debug, "RMS force after minimization %g\n", *rmsForce);
            }
            if (converged)
            {
                double enew = lbfgs->energy(InteractionType::EPOT);
                if (logFile)
                {
                    fprintf(logFile, "Minimization iteration %d/%d energy %g rms force %g\n",
                            retry+1, maxRetry, enew, *rmsForce);
                }
                if (enew < eMin)
                {
                    // Store new structure only when it has lower energy.
                    newCoords[current]   = *lbfgs->coordinates();
                    newEnergies[current] = *lbfgs->energies();
                    forces[current]      = *lbfgs->forces();
                    eMin                 = enew;
                }
            }
            retry += 1;
        }
        delete lbfgs;
    }
    else
    {
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
                std::vector<double> f0;
                MatrixWrapper Hessian(DIM*theAtoms.size(), DIM*theAtoms.size());
                computeHessian(pd, mol, forceComp, &newCoords[current],
                               theAtoms, &Hessian, &f0, &newEnergies[current]);
                               
                if (logFile && firstStep && false)
                {
                    fprintf(logFile, "Hessian:\n%s\n", Hessian.toString().c_str());
                }
                // Averaging the Hessian completely messes up the information here
                // since it is exactly the deviations from symmetry we need to
                // address. Leave this comment in here for future reference.
                // Hessian.averageTriangle();
                MatrixWrapper H2(Hessian);
                // Solve H delta X = -grad (E) = force(E)
                int result = Hessian.solve(f0, &deltaX[current]);
                if (0 != result)
                {
                    eMin = eMinimizeStatus::Solver;
                }
                else
                {
                    // Copy f0 to forces
                    size_t jj = 0;
                    for(auto ii : theAtoms)
                    {
                        for(int m = 0; m < DIM; m++)
                        {
                            forces[current][ii][m] = f0[jj++];
                        }
                    }
                    if (debug && firstStep)
                    {
                        fprintf(debug, "Sum of forces: %g sum of displacement: %g\n",
                                std::accumulate(f0.begin(), f0.end(), 0.0),
                                std::accumulate(deltaX[current].begin(), deltaX[current].end(), 0.0));
                        fprintf(debug, "H deltaX\n");
                        for(size_t ii = 0; ii < DIM*theAtoms.size(); ii++)
                        {
                            double f1 = 0;
                            for(size_t jj = 0; jj < DIM*theAtoms.size(); jj++)
                            {
                                f1 += H2.get(ii, jj) * deltaX[current][jj];
                            }
                            fprintf(debug, "f0[%2zu] = %12g  deltaX = %12g, f0-H.deltaX = %12g\n",
                                    ii, f0[ii], deltaX[current][ii], f0[ii]-f1);
                        }
                    }
                }
                gamma = 1;
                for(size_t ii = 0; ii < theAtoms.size(); ii++)
                {
                    for(int m = 0; m < DIM; m++)
                    {
                        gamma += gmx::square(deltaX[current][DIM*ii+m]);
                    }
                }
                gamma *= simConfig.overRelax();
            }
            break;
        case eMinimizeAlgorithm::Steep:
            {
                // https://en.wikipedia.org/wiki/Gradient_descent
                msShellForce = forceComp->compute(pd, mol->topology(), &newCoords[current],
                                                  &forces[current], &newEnergies[current]);
                if (!firstStep)
                {
                    double teller = 0, noemer = 0;
                    for(auto atomI : theAtoms)
                    {
                        for(int m = 0; m < DIM; m++)
                        {
                            double df  = (forces[current][atomI][m]-forces[next][atomI][m]);
                            // TODO check which one it should be
                            teller    += (newCoords[current][atomI][m]-newCoords[next][atomI][m])*df;
                            //teller    += (deltaX[current][i]-deltaX[next][i])*df;
                            noemer    += df*df;
                        }
                    }
                    if (noemer > 0)
                    {
                        gamma = std::min(gamma_max, std::abs(teller)/noemer);
                    }
                }
                int i = 0;
                for(auto atomI : theAtoms)
                {
                    for(int m = 0; m < DIM; m++)
                    {
                        deltaX[current][i++] = forces[current][atomI][m];
                    }
                }
            }
            break;
        default:
            break;
        }
        if (firstStep)
        {
            // Store energy
            if (energies)
            {
                *energies = newEnergies[current];
            }
            epotMin = newEnergies[current][InteractionType::EPOT];
            msfMin  = msForce(theAtoms, forces[current]);
            // Write stuff
            printEnergies(logFile, myIter, msfMin, 0, 0, newEnergies[current], gamma);
            firstStep = false;
        }
        if (eMinimizeStatus::OK != eMin)
        {
            return eMin;
        }
        bool acceptStep  = false;
        // Update coordinates and check energy
        do
        {
            updateCoords(theAtoms, newCoords[current], &(newCoords[next]), gamma, deltaX[current]);
            
            // Do an energy and force calculation with the new coordinates
            msShellForce   = forceComp->compute(pd, mol->topology(), &newCoords[next], &forces[next], &newEnergies[next]);
            double epotNew = newEnergies[next][InteractionType::EPOT];
            double msfNew  = msForce(theAtoms, forces[next]);
            acceptStep     = (epotNew <= epotMin) || (msfNew < msfMin);
            double gmin    = gamma;
            if (!acceptStep)
            {
                if (logFile)
                {
                    fprintf(logFile, "Trying line minimization, step 1.\n");
                    printEnergies(logFile, myIter, msfNew, 0, 0, newEnergies[next], gamma);
                }
                updateCoords(theAtoms, newCoords[current], &(newCoords[next]), 0.5*gamma, deltaX[current]);
                msShellForce    = forceComp->compute(pd, mol->topology(), &newCoords[next], &forces[next], &newEnergies[next]);
                if (logFile)
                {
                    fprintf(logFile, "Trying line minimization, step 2.\n");
                    printEnergies(logFile, myIter, msForce(theAtoms, forces[next]), 0, 0, newEnergies[next], gamma);
                }
                double epotHalf = newEnergies[next][InteractionType::EPOT];
                // Solve line minimization. We have x = 0, 0.5, 1.0 and y epotMin, epotHalf, epotNew
                // We solve M abc = yyy using abc = Minverse yyy
                // M = { { x0^2, x0, 1.0 }, { x1^2, x1, 1.0 }, { x2^2, x2, 1.0 } }
                matrix Minverse = { {2., -4., 2.}, {-3., 4., -1.}, {1., 0., 0.} };
                rvec   abc, yyy = { epotMin, epotHalf, epotNew };
                mvmul(Minverse, yyy, abc);
                if (abc[0] > 0)
                {
                    gmin = -gamma*abc[1]/(2*abc[0]);
                    updateCoords(theAtoms, newCoords[current], &(newCoords[next]), gmin, deltaX[current]);
                    (void) forceComp->compute(pd, mol->topology(), &newCoords[next], &forces[next], &newEnergies[next]);
                    double epotGmin = newEnergies[next][InteractionType::EPOT];
                    acceptStep = (epotGmin <= epotMin);
                }
                else if (logFile)
                {
                    fprintf(logFile, "Line minimization parameters a = %g b = %g c = %g\n",
                            abc[0], abc[1], abc[2]);
                    fprintf(logFile, "EpotMin = %g, EpotHalf = %g, EpotNew = %g\n",
                            epotMin, epotHalf, epotNew);
                }
            }
            if (acceptStep)
            {
                if (epotNew < epotMin)
                {
                    epotMin = epotNew;
                }
                if (msfNew < msfMin)
                {
                    msfMin = msfNew;
                }
                gamma   = gmin;
                current = next;
            }
            else
            {
                return eMinimizeStatus::NoMinimum;
            }
            myIter += 1;
        } while (!acceptStep && (myIter < simConfig.maxIter() || 0 == simConfig.maxIter()));
        
        // Now check force convergence
        double msAtomForce = msForce(theAtoms, forces[current]);
        converged = msAtomForce <= msForceToler;
        double fcurprev = 0;
        for(size_t ii = 0; ii < theAtoms.size(); ii++)
        {
            auto i = theAtoms[ii];
            for(int m = 0; m < DIM; m++)
            {
                fcurprev += forces[current][i][m]*forces[next][i][m];
            }
        }
        printEnergies(logFile, myIter, msAtomForce, msShellForce,
                      fcurprev, newEnergies[current], gamma);
        if (debug)
        {
            for(size_t kk = 0; kk < mol->topology()->nAtoms(); kk++)
            {
                fprintf(debug, "f[%2zu] =  %10g  %10g  %10g x[%2zu] = %10g  %10g  %10g\n", kk,
                        forces[current][kk][XX], forces[current][kk][YY], forces[current][kk][ZZ], kk,
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
    }
    if (converged)
    {
        *coords   = newCoords[current];
        if (energies && !newEnergies[current].empty())
        {
            *energies = newEnergies[current];
        }
        // Re-compute the energy one last time.
        // TODO: is this really needed?
        // (void) forceComp->compute(pd, mol->topology(), coords, &forces[current], energies);
        return eMinimizeStatus::OK;
    }
    else
    {
        return eMinimizeStatus::TooManySteps;
    }
}

double MolHandler::coordinateRmsd(const ACTMol                 *mol,
                                  const std::vector<gmx::RVec> &xref,
                                  std::vector<gmx::RVec>       *xfit) const
{
    // Compute RMSD
    auto   &myatoms = mol->atomsConst();
    double  tmass   = mol->totalMass();
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
        fprintf(exvg, "  %.5f", e.second);
    }
    auto epot = energies.find(InteractionType::EPOT)->second;
    auto temp = ekin*2/(3*nDOF*BOLTZ);
    fprintf(exvg, " %.5f %.5f %.5f\n", ekin, ekin+epot, temp);
}

static void initArrays(const ACTMol           *mol,
                       std::vector<gmx::RVec> *coordinates, 
                       std::vector<gmx::RVec> *velocities, 
                       std::vector<gmx::RVec> *forcesCur,
                       std::vector<gmx::RVec> *forcesPrev,
                       std::vector<double>    *inv2mass,
                       std::vector<double>    *mass_2)
{
    auto xmol = mol->xOriginal();
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
                       
void MolHandler::simulate(const ForceField              *pd,
                          ACTMol                        *mol,
                          const ForceComputer           *forceComp,
                          const SimulationConfigHandler &simConfig,
                          FILE                          *logFile,
                          const char                    *trajectoryFile,
                          const char                    *energyFile,
                          const gmx_output_env_t        *oenv) const
{
    if (simConfig.nsteps() <= 0)
    {
        // No simulation to be done!
        return;
    }
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
    std::vector<int> trajIndex;
    {
        auto atoms = mol->atomsConst();
        for(size_t ii = 0; ii < atoms.size(); ii++)
        {
            if (simConfig.writeShells() || eptAtom == atoms[ii].pType())
            {
                trajIndex.push_back(ii);
            }
        }
    }

    int nDOF = DIM*mol->nRealAtoms();
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
        forceComp->compute(pd, mol->topology(), &coordinates,
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
            // TODO: rewrite without gromacs fluff.
            pdbWriter(traj, title, mol->atomsConst(),
                      coordinates, mol->topology()->residueNames(),
                      epbcNONE,
                      box, chain, step+1, trajIndex,
                      nullptr, true);
        }
        // Swap force arrays
        cur = prev;
    }
    
    gmx_ffclose(traj);
    xvgrclose(exvg);
    fprintf(logFile, "Simulation finished correctly.\n");
}

} // namespace alexandria

/*
 * This source file is part of the Alexandria Chemistry Toolkit.
 *
 * Copyright (C) 2022-2026
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
 * \author Julian Marrades <julian@marrad.es>
 */

// This include has to come first, to prevent the GROMACS definition of
// "real" from messing up the library.
#pragma GCC diagnostic push
#pragma GCC diagnostic ignored "-Wpragmas"
#pragma GCC diagnostic ignored "-Wunknown-warning-option"
#pragma GCC diagnostic ignored "-Wmaybe-uninitialized"
#pragma GCC diagnostic ignored "-Wdeprecated-redundant-constexpr-static-def"
#pragma GCC diagnostic ignored "-Wdeprecated-copy-with-dtor"
#pragma GCC diagnostic ignored "-Wdeprecated-enum-enum-conversion"
#include "Eigen/Eigenvalues"
#pragma GCC diagnostic pop

#include "molhandler.h"

#include <algorithm>
#include <functional>
#include <array>
#include <cmath>
#include <deque>
#include <limits>
#include <numeric>
#include <random>

#include "act/alexandria/pdbwriter.h"
#include "act/basics/coordinate.h"
#include "act/constraints/shake.h"
#include "act/properties/velocityhandler.h"
#include "act/molprop/molpropobservable.h"
#include "act/topology/actmol.h"
#include "act/utility/rotation.h"
#include "act/utility/timings.h"
#include "act/utility/units.h"
#include "external/stlbfgs/stlbfgs.h"
#include "gromacs/fileio/xvgr.h"
#include "gromacs/math/units.h"
#include "gromacs/math/vec.h"
#include "gromacs/math/vectypes.h"
#include "gromacs/mdtypes/enerdata.h"
#include "gromacs/pbcutil/pbc.h"
#include "gromacs/utility/futil.h"

namespace alexandria
{

//! Map enum to string for eMinimizeStatus
std::map<eMinimizeStatus, std::string> eMinStat2String = {
    { eMinimizeStatus::OK,           "OK"             },
    { eMinimizeStatus::TooManySteps, "Too many steps" },
    { eMinimizeStatus::Solver,       "Solver failed"  },
    { eMinimizeStatus::NoMinimum,    "Could not find a lower minimum"  },
    { eMinimizeStatus::Constraints,  "Could not satisfy the constraints"  }
};

const std::string &eMinimizeStatusToString(eMinimizeStatus e)
{
    return eMinStat2String[e];
}

std::string minimizeStatusMessage(eMinimizeStatus e, double energy)
{
    if (eMinimizeStatus::NoMinimum == e)
    {
        return gmx::formatString("%s than %g kJ/mol", eMinimizeStatusToString(e).c_str(), energy);
    }
    return eMinimizeStatusToString(e);
}

void MolHandler::computeHessian(const ForceField                  *pd,
                                const Topology                    *topology,
                                const ForceComputer               *forceComp,
                                std::vector<gmx::RVec>            *coords,
                                const std::vector<int>            &atomIndex,
                                MatrixWrapper                     *hessian,
                                std::vector<double>               *forceZero,
                                std::map<InteractionType, double> *energyZero,
                                std::vector<gmx::RVec>            *dpdq) const
{
    Coordinates x(*coords);
    computeHessian(pd, topology, forceComp, &x, atomIndex, hessian, forceZero, energyZero, dpdq);
    *coords = x.toRVec();
}

void MolHandler::computeHessian(const ForceField                  *pd,
                                const Topology                    *topology,
                                const ForceComputer               *forceComp,
                                Coordinates                       *coords,
                                const std::vector<int>            &atomIndex,
                                MatrixWrapper                     *hessian,
                                std::vector<double>               *forceZero,
                                std::map<InteractionType, double> *energyZero,
                                std::vector<gmx::RVec>            *dpdq) const
{
    const auto &atoms    = topology->atoms();
    /* The derivatives of the forces are central differences. In double
     * precision of second order, with a step in nm of the square root of
     * the machine precision, which balances the truncation error against
     * the rounding error of the forces. In single precision these rounding
     * errors would change the squares of the frequencies by about 100
     * cm^-2, so they are of fourth order there, with a larger step that is
     * a power of two.
     */
    const std::vector<int> offsets  = GMX_DOUBLE ? std::vector<int>{ -1, 1 } : std::vector<int>{ -2, -1, 1, 2 };
    const size_t           noffset  = offsets.size();
    double                 stepSize = GMX_DOUBLE ? std::sqrt(GMX_REAL_EPS) : 1.0/512;

    // The input structure first, then for each coordinate of each atom
    // the structures with that coordinate shifted. They are all computed
    // together.
    std::vector<Coordinates>         xs = { *coords };
    std::vector<std::vector<double>> xxx;
    for(size_t ai = 0; ai < atomIndex.size(); ai++)
    {
        auto atomI = atomIndex[ai];
        for(int atomXYZ = 0; atomXYZ < DIM; atomXYZ++)
        {
            // The displacements, rounded to the unit of the integer
            // coordinates, such that the displaced coordinates are exact,
            // see coordinate.h
            std::vector<double> xd(noffset);
            for(size_t k = 0; k < noffset; k++)
            {
                xd[k] = fromCoordinate(toCoordinate(offsets[k]*stepSize));
                xs.push_back(*coords);
                xs.back()[atomI][atomXYZ] += xd[k];
            }
            xxx.push_back(xd);
        }
    }
    ConformationSet set;
    set.top = topology;
    set.setCoordinates(xs);
    // The derivatives are differences of nearly equal forces, which are
    // computed in double precision also in a single precision build
    set.data.doublePrecision = true;
    std::vector<std::map<InteractionType, double>> energies;
    forceComp->computeTotal(nullptr, pd, &set, ShellOptions(), &energies);
    const auto &cdata = set.data;
    // The input structure with the relaxed shells and the virtual sites
    *coords     = Coordinates(cdata.x(0));
    *energyZero = energies[0];
    forceZero->clear();
    for(auto &atom : atomIndex)
    {
        for(int m = 0; m < DIM; m++)
        {
            forceZero->push_back(cdata.f(0)[atom][m]);
        }
    }
    if (dpdq)
    {
        const gmx::RVec zero = { 0, 0, 0 };
        dpdq->resize(DIM*atomIndex.size(), zero);
    }
    // The derivative from the values at the shifted coordinates xd
    const auto derivative = [noffset](const std::vector<double> &xd, auto value)
    {
        if (noffset == 2)
        {
            return (value(1) - value(0))/(xd[1] - xd[0]);
        }
        double h = (xd[2] - xd[1])/2;
        return (8*(value(2) - value(1)) - (value(3) - value(0)))/(12*h);
    };
    for(size_t ai = 0; ai < atomIndex.size(); ai++)
    {
        for(int atomXYZ = 0; atomXYZ < DIM; atomXYZ++)
        {
            int         column = ai*DIM+atomXYZ;
            const auto &xd     = xxx[column];
            auto        conf   = [&](size_t k) { return 1+noffset*column+k; };
            for(size_t aj = 0; aj < atomIndex.size(); aj++)
            {
                int    atomJ = atomIndex[aj];
                for(int d = 0; d < DIM; d++)
                {
                    int    row   = aj*DIM+d;
                    double value = -derivative(xd, [&](size_t k) { return cdata.f(conf(k))[atomJ][d]; });
                    hessian->set(column, row, value);
                }
            }
            if (dpdq)
            {
                // To compute dipole we need to take shells into account as well!
                // The same goes for virtual sites.
                std::vector<gmx::DVec> mu(noffset, { 0, 0, 0 });
                for(size_t k = 0; k < noffset; k++)
                {
                    const auto x = cdata.x(conf(k));
                    for(size_t ii = 0; ii < atoms.size(); ii++)
                    {
                        for(int m = 0; m < DIM; m++)
                        {
                            mu[k][m] += atoms[ii].totalCharge()*x[ii][m];
                        }
                    }
                }
                for(int d = 0; d < DIM; d++)
                {
                    (*dpdq)[column][d] = derivative(xd, [&](size_t k) { return mu[k][d]; });
                }
            }
        }
    }
}

/*! \brief Orthonormal directions that a step of the coordinates should not take
 * These are the translations and rotations of all moving particles, if
 * they may move freely, and the gradients of the constraints. With
 * weights, the square roots of the masses, the translations and rotations
 * are those of mass-weighted coordinates, about the center of mass, as for
 * a normal mode analysis in the Eckart frame.
 * \param[in] theAtoms    The moving particles
 * \param[in] x           The coordinates of all particles
 * \param[in] rigidBody   Whether to include the translations and rotations
 * \param[in] constrainer The constraints, or nullptr if there are none
 * \param[in] weights     A weight for each moving particle, or empty for one
 * \return the directions, each with the coordinates of the moving particles
 */
static std::vector<std::vector<double>> fixedDirections(const std::vector<int>       &theAtoms,
                                                        const std::vector<gmx::RVec> &x,
                                                        bool                          rigidBody,
                                                        const Constrainer            *constrainer,
                                                        const std::vector<double>    &weights = {})
{
    size_t                           n = DIM*theAtoms.size();
    std::vector<std::vector<double>> raw;
    const auto weight = [&weights](size_t i) { return weights.empty() ? 1.0 : weights[i]; };
    if (rigidBody && !theAtoms.empty())
    {
        rvec   center = { 0, 0, 0 };
        double wtot   = 0;
        for (size_t i = 0; i < theAtoms.size(); i++)
        {
            double w2 = gmx::square(weight(i));
            for (int m = 0; m < DIM; m++)
            {
                center[m] += w2*x[theAtoms[i]][m];
            }
            wtot += w2;
        }
        svmul(1.0/wtot, center, center);
        for (int m = 0; m < DIM; m++)
        {
            std::vector<double> translation(n, 0.0), rotation(n, 0.0);
            for (size_t i = 0; i < theAtoms.size(); i++)
            {
                translation[DIM*i+m] = weight(i);
                // The axis m crossed with the position from the center
                rvec axis = { 0, 0, 0 }, r, v;
                axis[m] = 1;
                rvec_sub(x[theAtoms[i]], center, r);
                cprod(axis, r, v);
                for (int k = 0; k < DIM; k++)
                {
                    rotation[DIM*i+k] = weight(i)*v[k];
                }
            }
            raw.push_back(translation);
            raw.push_back(rotation);
        }
    }
    size_t nrigid = raw.size();
    if (constrainer)
    {
        std::map<int, size_t> index;
        for (size_t i = 0; i < theAtoms.size(); i++)
        {
            index.insert({ theAtoms[i], i });
        }
        for (const auto &pair : constrainer->pairs())
        {
            std::vector<double> gradient(n, 0.0);
            rvec dx;
            rvec_sub(x[pair[0]], x[pair[1]], dx);
            for (int k = 0; k < 2; k++)
            {
                const auto ai = index.find(pair[k]);
                if (ai == index.end())
                {
                    // A frozen particle
                    continue;
                }
                for (int m = 0; m < DIM; m++)
                {
                    gradient[DIM*ai->second+m] = (0 == k ? 1 : -1)*dx[m];
                }
            }
            raw.push_back(gradient);
        }
    }
    // The largest rotation, the rotation about the axis of a linear molecule
    // is small compared to it, also if the molecule is slightly bent
    double maxRotation = 0;
    for (size_t r = 1; r < nrigid; r += 2)
    {
        maxRotation = std::max(maxRotation, std::sqrt(std::inner_product(raw[r].begin(), raw[r].end(), raw[r].begin(), 0.0)));
    }
    // Gram-Schmidt, twice for precision, leaving out dependent directions
    // such as the third rotation of a linear molecule
    std::vector<std::vector<double>> directions;
    for (size_t r = 0; r < raw.size(); r++)
    {
        auto  &v         = raw[r];
        bool   rotation  = r < nrigid && r % 2 == 1;
        double threshold = (rotation ? 1e-3*maxRotation :
                            1e-8*std::sqrt(std::inner_product(v.begin(), v.end(), v.begin(), 0.0)));
        for (int pass = 0; pass < 2; pass++)
        {
            for (const auto &d : directions)
            {
                double c = std::inner_product(d.begin(), d.end(), v.begin(), 0.0);
                for (size_t k = 0; k < n; k++)
                {
                    v[k] -= c*d[k];
                }
            }
        }
        double norm = std::sqrt(std::inner_product(v.begin(), v.end(), v.begin(), 0.0));
        if (norm > threshold)
        {
            for (auto &vk : v)
            {
                vk /= norm;
            }
            directions.push_back(v);
        }
    }
    return directions;
}

/*! \brief Compute frequencies
 * \param[in]  eigenvalues The eigenvalues of the mass-weighted Hessian
 * \param[in]  vibration   Whether each eigenvector is a vibration rather
 *                         than a translation or rotation
 * \param[out] frequencies The frequencies of the vibrations
 * \param[out] output      Warnings, if not nullptr
 */
static void computeFrequencies(const std::vector<double> &eigenvalues,
                               const std::vector<bool>   &vibration,
                               std::vector<double>       *frequencies,
                               std::vector<std::string>  *output)
{
    frequencies->clear();
    auto        mpo  = MolPropObservable::FREQUENCY;
    const char *unit = mpo_unit2(mpo);
    for (size_t i = 0; i < eigenvalues.size(); i++)
    {
        if (!vibration[i])
        {
            continue;
        }
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

//! \brief Generate output of frequencies and intensities
static void outputFreqInten(const ACTMol               *mol,
                            const std::vector<double> &frequencies,
                            const std::vector<double> &intensities,
                            std::vector<std::string>  *output)
{
    auto mpo = MolPropObservable::FREQUENCY;
    const char *unit = mpo_unit2(mpo);
    double delta = 0;
    for (auto &ee : mol->experimentConst())
    {
        if (ee.hasMolPropObservable(mpo))
        {
            auto &harm = ee.propertyConst(mpo);
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
            break;
        }
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

//! \brief Solve the eigenvalue equation using the Eigen package
static void solveEigen(const MatrixWrapper          &hessian,
                       const std::vector<int>       &atomIndex,
                       const std::vector<ActAtom>   &atoms,
                       const std::vector<gmx::RVec> &dpdq,
                       const std::vector<std::vector<double>> &rigidBody,
                       std::vector<double>          *frequencies,
                       std::vector<double>          *intensities,
                       std::vector<std::string>     *output,
                       bool                          debugNMA)
{
    if (output)
    {
        output->emplace_back("Using Eigen to solve eigenvalue problem.");
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
    // The translations and rotations are projected out, such that they are
    // eigenvectors with eigenvalue zero, also if the structure is not
    // exactly at a minimum, and they are not mixed with vibrations.
    Eigen::MatrixXd projection = Eigen::MatrixXd::Identity(matrixSide, matrixSide);
    for (const auto &d : rigidBody)
    {
        Eigen::Map<const Eigen::VectorXd> dv(d.data(), matrixSide);
        projection -= dv*dv.transpose();
    }
    mat = projection*mat*projection;
    // We need this special solver for hermitian matrices.
    // It will automatically sort the eigenvalues.
    const Eigen::SelfAdjointEigenSolver<Eigen::MatrixXd> solver(mat);
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
    // The eigenvectors that lie in the space of translations and rotations
    std::vector<bool> vibration(matrixSide, true);
    int   nnegative = 0;
    for(Eigen::Index i = 0; i < matrixSide; i++)
    {
        double overlap = 0;
        for (const auto &d : rigidBody)
        {
            const Eigen::Map<const Eigen::VectorXd> dv(d.data(), matrixSide);
            overlap += std::pow(dv.dot(evecs.col(i)), 2);
        }
        vibration[i] = overlap < 0.5;
        if (vibration[i] && evals[i] < 0)
        {
            nnegative += 1;
        }
    }
    if (output && nnegative > 0)
    {
        output->push_back(gmx::formatString("There are %d negative eigenvalues.",
                                            nnegative));
    }

    // The intensity of a mode is the square of the derivative of the
    // dipole moment along the mode, which is the sum of the derivatives
    // with respect to the coordinates, weighted by the displacements.
    std::vector<double> In(matrixSide, 0);
    for(Eigen::Index col = 0; col < matrixSide; col++)
    {
        auto      evec = evecs.col(col);
        gmx::RVec dmudQ = { 0, 0, 0 };
        for (Eigen::Index row = 0; row < matrixSide; row++)
        {
            // The eigenvectors should be normalized (length 1).
            // This is indeed the case.
            auto myev  = evec[row];
            auto atomI = row / DIM;
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
            for(int d = 0; d < DIM; d++)
            {
                dmudQ[d] += dpdq[row][d]*myev*gmx::invsqrt(atoms[ai].mass());
            }
        }
        In[col] = iprod(dmudQ, dmudQ);
    }
    intensities->clear();
    for(int i = 0; i < matrixSide; i++)
    {
        if (vibration[i])
        {
            intensities->push_back(In[i]);
        }
    }
    std::vector<double> eigenvalues;
    for(int i = 0; i < matrixSide; i++)
    {
        eigenvalues.push_back(evals(i));
    }
    computeFrequencies(eigenvalues, vibration, frequencies, output);
}

//! \brief Sort Frequencies and Intensities
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
            std::abs(freq[i] - freq[i+1]) < freq_toler &&
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
        if (atoms[atom].pType() == ActParticle::Atom)
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
    computeHessian(pd, mol->topology(), forceComp, coords,
                   atomIndex, &hessian, &f0, &energies, &dpdq);

    if (output && debugNMA)
    {
        output->emplace_back("Symmetrized Hessian:");
        output->push_back(hessian.toString());
        for(size_t i = 0; i < dpdq.size(); i++)
        {
            output->push_back(gmx::formatString("dpdq[%2zu][%c] =  %10g  %10g  %10g",
                                                i/DIM, 'x' + static_cast<char>(i % DIM),
                                                dpdq[i][XX], dpdq[i][YY], dpdq[i][ZZ]));
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

    // The translations and rotations of the mass-weighted coordinates, five
    // for a linear molecule, since its rotation about its axis is left out
    std::vector<double> sqrtMass;
    for (auto ai : atomIndex)
    {
        sqrtMass.push_back(std::sqrt(atoms[ai].mass()));
    }
    const auto rigidBody = fixedDirections(atomIndex, *coords, true, nullptr, sqrtMass);
    if (output && rigidBody.size() != (mol->linearMolecule() ? 5u : 6u))
    {
        output->push_back(gmx::formatString("Warning: %zu translations and rotations for a %slinear molecule",
                                            rigidBody.size(), mol->linearMolecule() ? "" : "non-"));
    }
    // Call diagonalization routine
    if (output)
    {
        output->emplace_back("Diagonalizing Hessian to find eigenvectors.");
    }
    solveEigen(hessian, atomIndex, atoms, dpdq, rigidBody,
               frequencies, intensities, output, debugNMA);
    // Check whether there are very similar frequencies, then change the sorting 
    // according to intensities.
    //! \todo make ftoler a parameter.
    double ftoler = 0.001; // Frequency unit is internal unit.
    resortFreqIntens(frequencies, intensities, ftoler);

    if (output)
    {
        outputFreqInten(mol, *frequencies, *intensities, output);
    }

}

//! \brief Yet another routine to print energies
static void printEnergies(gmx::TextWriter *tw, int myIter,
                          double msAtomForce,
                          double msShellForce,
                          double fcurprev,
                          const std::map<InteractionType, double> &energies,
                          double gamma)
{
    if (nullptr == tw)
    {
        return;
    }
    tw->writeLineFormatted("Iter %5d msAtomForce %10g msShellForce %10g fcur.fprev %10g gamma %10g",
                           myIter, msAtomForce, msShellForce, fcurprev, gamma);
    double evdw = 0;
    for(const auto &ee : energies)
    {
        tw->writeStringFormatted("  %s %.5f", interactionTypeToString(ee.first).c_str(), ee.second);
        if (InteractionType::DISPERSION == ee.first || 
            InteractionType::EXCHANGE == ee.first)
        {
            evdw += ee.second;
        }
    }
    tw->writeLineFormatted(" (VANDERWAALS %.5f)", evdw);
}

//! \return the mean square force on atoms
static double msForce(const std::vector<int>       &theAtoms,
                      const std::vector<gmx::RVec> &forces)
{
    double msAtomForce  = 0;
    for(auto kk : theAtoms)
    {
        msAtomForce += iprod(forces[kk], forces[kk]);
    }
    if (!theAtoms.empty())
    {
        return msAtomForce /= theAtoms.size();
    }
    else
    {
        return 0;
    }
}

/*! \brief The energy change from x to xt
 * The energies are sums of many terms, with a rounding error of the
 * machine precision times the size of the terms, which in single
 * precision is larger than the changes close to a minimum. Smaller
 * changes are computed from the forces with the trapezoidal rule
 * instead, which is accurate to third order in the step. For large steps
 * it is not, which shows when it differs from the energy change by more
 * than the rounding error, then the energy change is used.
 * The displacements are computed from the integers, which is exact, see
 * coordinate.h.
 * \param[in] theAtoms The atoms that move
 * \param[in] e        The energies at x
 * \param[in] x        The old coordinates
 * \param[in] f        The forces at x
 * \param[in] enew     The potential energy at xt
 * \param[in] xt       The new coordinates
 * \param[in] ft       The forces at xt
 * \param[in] noise    The rounding error of the energy, see energyNoise
 * \return the change of the potential energy
 */
static double energyChange(const std::vector<int>                  &theAtoms,
                           const std::map<InteractionType, double> &e,
                           const Coordinates                       &x,
                           const std::vector<gmx::RVec>            &f,
                           double                                   enew,
                           const Coordinates                       &xt,
                           const std::vector<gmx::RVec>            &ft,
                           double                                   noise)
{
    double eold  = e.find(InteractionType::EPOT)->second;
    double scale = 0;
    for (const auto &term : e)
    {
        if (term.first != InteractionType::EPOT)
        {
            scale += std::abs(term.second);
        }
    }
    double de       = enew - eold;
    double rounding = std::max(10*GMX_REAL_EPS*scale, noise);
    if (std::abs(de) > rounding)
    {
        return de;
    }
    double work = 0;
    for (auto atom : theAtoms)
    {
        for (int m = 0; m < DIM; m++)
        {
            work += (f[atom][m] + ft[atom][m])*fromCoordinate(coordinateDifference(xt[atom].units(m), x[atom].units(m)));
        }
    }
    double trapezoid = -0.5*work;
    return std::abs(trapezoid - de) <= rounding ? trapezoid : de;
}

/*! \brief Estimate the rounding error of the energy
 * The energy terms are sums of terms that cancel, such that the rounding
 * error can be much larger than the machine precision times the terms,
 * in particular with shells with large charges. It is estimated from
 * small displacements of the particles, for which the trapezoidal rule
 * gives the energy change to well within the rounding error, the
 * estimate is four times the largest of a few samples.
 * \param[in] theAtoms The particles that move
 * \param[in] x        The coordinates
 * \param[in] evaluate Computes the energy and forces as the minimizer does
 * \return the largest rounding error found
 */
static double energyNoise(const std::vector<int>       &theAtoms,
                          const std::vector<gmx::RVec> &x,
                          const std::function<double(std::vector<gmx::RVec> *,
                                                     std::vector<gmx::RVec> *)> &evaluate)
{
    const double           step = 1e-5;
    std::vector<gmx::RVec> x0   = x, f0(x.size());
    double                 e0   = evaluate(&x0, &f0);
    double                 noise = 0;
    unsigned int           seed  = 1993;
    for (int probe = 0; probe < 5; probe++)
    {
        std::vector<gmx::RVec> xt = x0, ft(x.size());
        for (auto a : theAtoms)
        {
            for (int m = 0; m < DIM; m++)
            {
                seed      = seed*1103515245u + 12345u;
                xt[a][m] += ((seed >> 16) & 1) ? step : -step;
            }
        }
        double et   = evaluate(&xt, &ft);
        double work = 0;
        for (auto a : theAtoms)
        {
            for (int m = 0; m < DIM; m++)
            {
                work += (f0[a][m] + ft[a][m])*(xt[a][m] - x0[a][m]);
            }
        }
        if (std::isfinite(et) && std::isfinite(e0))
        {
            noise = std::max(noise, std::abs(et - e0 + 0.5*work));
        }
    }
    // The maximum of a few samples underestimates the largest errors
    return 4*noise;
}

/*! \brief Local class to provide data to the L-BFGS algorithm
 * Content is just a copy of inputs, respectively working variables.
 */
class StlbfgsHandler
{
private:
    //! Pointer to force field
    const ForceField                  *pd_;
    //! Pointer to message handler
    MsgHandler                        *mh_;
    //! Pointer to topology
    const Topology                    *topology_;
    //! Pointer to force computer
    const ForceComputer               *forceComp_;
    //! Atomic coordinates, integers such that steps are exact, see coordinate.h
    Coordinates                        coords_;
    //! Atomic forces
    std::vector<gmx::RVec>             forces_;
    //! Energies per interaction type
    std::map<InteractionType, double>  energies_;
    //! Vector of real atoms
    const std::vector<int>             theAtoms_;
    //! Constraints, or nullptr if there are none
    Constrainer                       *constrainer_;
    //! The lowest energy found
    double                             bestEnergy_ = std::numeric_limits<double>::max();
    //! The coordinates with the lowest energy, empty if there are none
    Coordinates                        bestCoords_;
    //! Whether there is a previous evaluation with a finite energy
    bool                               havePrevious_ = false;
    //! The coordinates of the previous evaluation
    Coordinates                        prevCoords_;
    //! The forces of the previous evaluation
    std::vector<gmx::RVec>             prevForces_;
    //! The energies of the previous evaluation
    std::map<InteractionType, double>  prevEnergies_;
    //! The energy along the path of the evaluations
    double                             pathEnergy_ = 0;
    //! The rounding error of the energy, see energyNoise
    double                             energyNoise_ = 0;
public:
    //! Constructor
    StlbfgsHandler(const ForceField       *pd,
                   MsgHandler             *mh,
                   const Topology         *topology,
                   const ForceComputer    *forceComp,
                   const std::vector<int> &theAtoms,
                   const std::vector<gmx::RVec> &coords,
                   Constrainer            *constrainer) :
        pd_(pd), mh_(mh), topology_(topology), forceComp_(forceComp), coords_(coords), theAtoms_(theAtoms),
        constrainer_(constrainer)
    {
        const gmx::RVec vzero = { 0, 0, 0 };
        forces_.resize(topology->atoms().size(), vzero);
    }

    //! \return the molecule
    const Topology *topology() const { return topology_; }

    //! \return the force field
    const ForceField *pd() const { return pd_; }

    //! \return the message handler
    MsgHandler *mh() const { return mh_; }

    //! \return the force computer
    const ForceComputer *forceComp() const { return forceComp_; }

    //! \return the coordinates
    Coordinates *coordinates() { return &coords_; }

    //! \return the forces
    std::vector<gmx::RVec> *forces() { return &forces_; }

    //! \return the energy mao
    std::map<InteractionType, double> *energies() { return &energies_; }

    //! \return a particular energy
    double energy(InteractionType itype) { return energies_[itype]; }

    //! \return the atom indices
    const std::vector<int> &theAtoms() const { return theAtoms_; }

    //! \return the constraints, or nullptr if there are none
    Constrainer *constrainer() const { return constrainer_; }

    /*! \brief Store the current coordinates if their energy is the lowest
     * \param[in] epot The potential energy of the current coordinates
     */
    void storeIfBest(double epot)
    {
        if (std::isfinite(epot) && epot < bestEnergy_)
        {
            bestEnergy_ = epot;
            bestCoords_ = coords_;
        }
    }

    //! \return the coordinates with the lowest energy, empty if there are none
    const Coordinates &bestCoordinates() const { return bestCoords_; }

    //! \return the lowest energy found
    double bestEnergy() const { return bestEnergy_; }

    /*! \brief Continue from the coordinates with the lowest energy
     * The energy along the path starts anew, see pathEnergy.
     */
    void restartFromBest()
    {
        coords_       = bestCoords_;
        havePrevious_ = false;
    }

    /*! \brief The energy along the path of the evaluations
     * The energy of the current coordinates is that of the previous
     * evaluation plus the change, which close to a minimum is computed
     * from the forces, see energyChange. The line search then also finds
     * decreases of the energy that are smaller than its rounding error.
     * \return the energy at the current coordinates
     */
    double pathEnergy()
    {
        double enew = energies_[InteractionType::EPOT];
        if (!std::isfinite(enew))
        {
            return enew;
        }
        if (havePrevious_)
        {
            pathEnergy_ += energyChange(theAtoms_, prevEnergies_, prevCoords_, prevForces_,
                                        enew, coords_, forces_, energyNoise_);
        }
        else
        {
            pathEnergy_ = enew;
        }
        havePrevious_ = true;
        prevCoords_   = coords_;
        prevForces_   = forces_;
        prevEnergies_ = energies_;
        return pathEnergy_;
    }

    //! \return the energy along the path at the previous evaluation
    double previousPathEnergy() const { return pathEnergy_; }

    //! \brief Set the rounding error of the energy
    void setEnergyNoise(double noise) { energyNoise_ = noise; }

    //! \return the rounding error of the energy, see energyNoise
    double energyNoise() const { return energyNoise_; }

    //! \brief print force on frozen atom
    void printFrozenForce(gmx::TextWriter *tw, const std::vector<int> &myFreeze)
    {
        if (tw && myFreeze.size() == 2)
        {
            // If minimization with constraints worked well, the resulting
            // force should be along the vector between the two frozen atoms.
            // Compute this vector first:
            gmx::RVec dx = coords_[myFreeze[0]] - coords_[myFreeze[1]];
            copy_rvec(forces_[myFreeze[0]], dx);
            // Now determine the rotation of the force vector to the Z axis
            rvec frot;
            // If there is a component away from the Z-axis, rotate around
            // Z first, and then Y
            real xy     = std::sqrt(dx[XX]*dx[XX] + dx[YY]*dx[YY]);
            if (xy > 0)
            {
                real sinphi = dx[YY]/xy;
                real cosphi = dx[XX]/xy;
                // Rotation matrix around Z in negative phi direction
                matrix m1 = { { cosphi, sinphi, 0 }, { -sinphi, cosphi, 0 }, { 0, 0, 1 } };
                real xnorm  = norm(dx);
                real costh  = dx[ZZ]/xnorm;
                real sinth  = std::sqrt(1 - costh*costh);
                // Rotation matrix around Y in negative theta direction
                matrix m2 = { { costh, 0, sinth }, { 0, 1, 0 }, { -sinth, 0, costh } };
                matrix m12;
                mmul(m1, m2, m12);
                mvmul(m12, forces_[myFreeze[0]], frot);
            }
            else
            {
                // Otherwise just copy the force on the first atom
                copy_rvec(forces_[myFreeze[0]], frot);
            }
            tw->writeLineFormatted("Force on frozen atom %d = %12f %12f %12f",
                                   myFreeze[0], frot[XX], frot[YY], frot[ZZ]);
        }
    }
};

//! The return of the son of global variables. Easy but ugly.
static StlbfgsHandler *lbfgs = nullptr;

//! \brief callback function for lbfgs minimizer
static void func(const std::vector<double> &x, double &f, std::vector<double> &g)
{
    const auto theAtoms = lbfgs->theAtoms();
    if (x.size() != DIM*theAtoms.size())
    {
        GMX_THROW(gmx::InternalError(gmx::formatString("Number of coordinates %lu does not match number of real atoms %lu",
                                                       x.size(), theAtoms.size()).c_str()));
    }
    auto coords = lbfgs->coordinates();
    auto forces = lbfgs->forces();
    // The previous coordinates, that satisfy the constraints and are
    // within the range of the integer coordinates, see coordinate.h
    const Coordinates      reference = *coords;
    bool                   inRange   = std::all_of(x.begin(), x.end(),
                                                   [](double xi) { return std::abs(xi) < 0.5*coordinateMax; });
    size_t jj   = 0;
    for (auto i : theAtoms)
    {
        for(int j = 0; j < DIM && inRange; j++)
        {
            // Rounded to real, as before the coordinates were integers.
            // With the energies rounded to single precision, the line
            // search only converges if small steps can leave the
            // coordinates and the energy unchanged.
            (*coords)[i][j] = static_cast<real>(x[jj++]);
        }
    }
    // Energy penalty to make the line search back off
    double penalty = 0;
    if (!inRange)
    {
        // A step that takes atoms too far from the origin
        *coords = reference;
        penalty = 1e6;
    }
    else if (lbfgs->constrainer())
    {
        // The energy is that of the closest coordinates that satisfy
        // the constraints, and the gradient is along the constraints.
        // The constraint directions are taken from x itself, such that
        // the energy only depends on x, as the line search requires.
        // The constraints work with coordinates in nm.
        auto xc = coords->toRVec();
        if (lbfgs->constrainer()->shake(lbfgs->mh(), &xc) != 0)
        {
            // A step that is too large for the constraints
            *coords = reference;
            penalty = 1e6;
        }
        else
        {
            coords->assign(xc);
        }
    }
    const gmx::RVec field = { 0, 0, 0 };
    lbfgs->forceComp()->constructVsiteCoordinates(lbfgs->topology(), coords);
    lbfgs->forceComp()->computeOnce(lbfgs->mh(),
                                    lbfgs->pd(), lbfgs->topology(), coords,
                                    forces, lbfgs->energies(), field);
    lbfgs->forceComp()->spreadVsiteForces(lbfgs->topology(), coords,
                                          forces);
    if (lbfgs->constrainer())
    {
        lbfgs->constrainer()->projectForces(lbfgs->mh(), coords->toRVec(), forces);
    }
    // With a penalty the coordinates are those of the previous evaluation
    f  = 0 == penalty ? lbfgs->pathEnergy() : lbfgs->previousPathEnergy() + penalty;
    if (0 == penalty)
    {
        lbfgs->storeIfBest(f);
    }
    jj = 0;
    rvec   fsum = { 0, 0, 0 };
    double fmax = 0;
    int    imax = -1;
    for (auto i : theAtoms)
    {
        for(int j = 0; j < DIM; j++)
        {
            double ff = (*forces)[i][j];
            g[jj++]   = -ff;
            if (std::abs(ff) > fmax)
            {
                fmax = std::abs(ff);
                imax = i;
            }
        }
        if (lbfgs->mh()->debug())
        {
            lbfgs->mh()->writeDebug(gmx::formatString("Force %03d  %12.3f  %12.3f  %12.3f\n", i,
                                                      (*forces)[i][XX], (*forces)[i][YY], (*forces)[i][ZZ]));
        }
        rvec_inc(fsum, (*forces)[i]);
    }
    if (lbfgs->mh()->debug())
    {
        lbfgs->mh()->writeDebug(gmx::formatString("Energy %.3f sum of forces is %g %g %g fmax %g on atom %d out of %zu\n",
                                                  f, fsum[XX], fsum[YY], fsum[ZZ], fmax, imax, theAtoms.size()));
    }
}

//! \return the largest displacement of a particle in a step of all coordinates
static double maxDisplacement(const std::vector<double> &dx)
{
    double dmax2 = 0;
    for (size_t i = 0; i + DIM <= dx.size(); i += DIM)
    {
        dmax2 = std::max(dmax2, dx[i]*dx[i] + dx[i+1]*dx[i+1] + dx[i+2]*dx[i+2]);
    }
    return std::sqrt(dmax2);
}

/*! \brief Newton steps from the eigenvectors of the Hessian
 * The Hessian is symmetrized and the fixed directions are projected out.
 * Along each eigenvector u with eigenvalue lambda the step is the force
 * along u divided by |lambda| + mu: the absolute value goes downhill also
 * where the curvature is negative and mu damps the step, as in the method
 * of Levenberg and Marquardt. Without damping near a minimum this is the
 * Newton step. Directions with an eigenvalue that is zero compared to the
 * largest, such as those of shells, which are relaxed in each energy
 * calculation, are left out.
 */
class NewtonStep
{
private:
    //! The eigenvalues
    Eigen::VectorXd values_;
    //! The eigenvectors
    Eigen::MatrixXd vectors_;
    //! The force along each eigenvector
    Eigen::VectorXd force_;
    //! Whether to use each eigenvector
    std::vector<bool> use_;
    //! The largest absolute eigenvalue
    double            largest_ = 0;
public:
    /*! \brief Decompose a Hessian
     * \param[in] hessian    The Hessian of the moving coordinates
     * \param[in] force      The forces on them
     * \param[in] directions The fixed directions, see fixedDirections
     */
    void decompose(const MatrixWrapper                    &hessian,
                   const std::vector<double>              &force,
                   const std::vector<std::vector<double>> &directions)
    {
        Eigen::Index    n = force.size();
        Eigen::MatrixXd h(n, n);
        for (Eigen::Index i = 0; i < n; i++)
        {
            for (Eigen::Index j = 0; j < n; j++)
            {
                h(i, j) = 0.5*(hessian.get(i, j) + hessian.get(j, i));
            }
        }
        Eigen::MatrixXd projection = Eigen::MatrixXd::Identity(n, n);
        for (const auto &d : directions)
        {
            Eigen::Map<const Eigen::VectorXd> dv(d.data(), n);
            projection -= dv*dv.transpose();
        }
        h = projection*h*projection;
        const Eigen::SelfAdjointEigenSolver<Eigen::MatrixXd> solver(h);
        values_  = solver.eigenvalues();
        vectors_ = solver.eigenvectors();
        const Eigen::Map<const Eigen::VectorXd> fv(force.data(), n);
        force_   = vectors_.transpose()*(projection*fv);
        largest_ = values_.cwiseAbs().maxCoeff();
        use_.assign(n, true);
        for (Eigen::Index k = 0; k < n; k++)
        {
            // The fixed directions are eigenvectors with eigenvalue zero
            double fixed = 0;
            for (const auto &d : directions)
            {
                const Eigen::Map<const Eigen::VectorXd> dv(d.data(), n);
                fixed += std::pow(dv.dot(vectors_.col(k)), 2);
            }
            use_[k] = fixed < 0.5 && std::abs(values_(k)) > 1e-10*largest_;
        }
    }

    //! \return the largest absolute eigenvalue
    double largestEigenvalue() const { return largest_; }

    /*! \brief The step
     * \param[in] mu    The damping
     * \param[in] scale Factor for the step, such as an overrelaxation
     * \return the step of the moving coordinates
     */
    std::vector<double> step(double mu, double scale) const
    {
        Eigen::VectorXd dx = Eigen::VectorXd::Zero(force_.size());
        for (Eigen::Index k = 0; k < force_.size(); k++)
        {
            if (use_[k])
            {
                dx += scale*force_(k)/(std::abs(values_(k)) + mu)*vectors_.col(k);
            }
        }
        return std::vector<double>(dx.data(), dx.data() + dx.size());
    }
};

eMinimizeStatus MolHandler::minimizeCoordinates(MsgHandler                        *msghandler,
                                                const ForceField                  *pd,
                                                const Topology                    *topology,
                                                const ForceComputer               *forceComp,
                                                const SimulationConfigHandler     &simConfig,
                                                std::vector<gmx::RVec>            *coords,
                                                std::map<InteractionType, double> *energies,
                                                const std::vector<int>            &freeze) const
{
    bool      converged    = false;
    int       myIter       = 0;
    // Check for meaningful convergence tolerance, both are on the mean
    // square force. The warning is given once per configuration.
    double    shellToler   = forceComp->forceTolerance();
    double    msForceToler = simConfig.forceTolerance();
    if (msForceToler < shellToler)
    {
        if (msForceToler != 0 && msghandler && simConfig.toleranceWarning())
        {
            msghandler->msg(ACTStatus::Warning,
                            gmx::formatString("Atom mean square force tolerance (%g) more strict than shell mean square force tolerance (%g). Increasing atom mean square force tolerance to hundred times that for shells.",
                                              msForceToler, shellToler));
        }
        msForceToler = 100*shellToler;
    }
    // List of atoms (not shells) and weighting factors
    const auto             &myatoms = topology->atoms();
    // Copy and sort the freeze array
    std::vector<int>  myFreeze = freeze;
    std::sort(myFreeze.begin(), myFreeze.end());
    // Renumber the myFreeze atoms
    for(size_t i = 0; i < myFreeze.size(); i++)
    {
        int count  = 0;
        bool found = false;
        for(size_t j = 0; !found && j < myatoms.size(); j++)
        {
            if (myatoms[j].pType() == ActParticle::Atom)
            {
                if (myFreeze[i] == count)
                {
                    myFreeze[i] = j;
                    found = true;
                }
                count += 1;
            }
        }
    }
    std::vector<int>  theAtoms;
    for(size_t atom = 0; atom < myatoms.size(); atom++)
    {
        // Store the atom numbers for the frozen atoms in the counting incl. shells
        if ((myatoms[atom].pType() == ActParticle::Atom || myatoms[atom].pType() == ActParticle::Shell) &&
            myFreeze.end() == std::find(myFreeze.begin(), myFreeze.end(), atom))
        {
            theAtoms.push_back(atom);
        }
    }
    // Constraints, the starting coordinates have to satisfy them
    std::vector<gmx::RVec> xstart = *coords;
    // Near the minimum the steps are much smaller than the errors that
    // the tolerance for simulations allows, which would make the energy
    // change irregularly. The Newton iterations in SHAKE converge
    // quadratically, such that a much smaller tolerance is cheap.
    Constrainer            constrainer(simConfig.constraintMaxIter(),
                                       std::min<double>(simConfig.constraintTolerance(),
                                                        1000*GMX_REAL_EPS), 1.0);
    bool                   constrain = simConfig.constrain();
    if (constrain)
    {
        constrainer.makeConstraints(msghandler, simConfig.constraintType(),
                                    pd, topology, xstart);
        constrainer.setFrozen(myFreeze);
        if (msghandler->verbose())
        {
            msghandler->write(gmx::formatString("Constraining %zu bonds and %zu angles in minimization.",
                                                constrainer.numBondConstraints(),
                                                constrainer.numAngleConstraints()));
        }
        if (constrainer.shake(msghandler, &xstart) != 0)
        {
            return eMinimizeStatus::Constraints;
        }
    }
    // Remove the force components along the constraints
    auto constrainForces = [&](const std::vector<gmx::RVec> &x,
                               std::vector<gmx::RVec>       *f)
    {
        if (constrain)
        {
            constrainer.projectForces(msghandler, x, f);
        }
    };
    /* The mean square change of the forces when each coordinate of the
     * atoms changes by one unit in the last place, in directions that
     * differ between the coordinates. Since the coordinates are rounded,
     * the forces at a minimum are not smaller than this, which in single
     * precision is much larger than in double precision.
     */
    auto roundingMsForce = [&](const std::vector<gmx::RVec> &x,
                               const std::vector<gmx::RVec> &f)
    {
        // The differences of the coordinates are computed from the integer
        // coordinates and rounded to real, see coordinate.h, by up to the
        // rounding of the largest difference. Each coordinate is moved by
        // that much and one unit in the last place.
        double extent = 0;
        for (auto a : theAtoms)
        {
            for (int m = 0; m < DIM; m++)
            {
                extent = std::max(extent, 2*std::abs(static_cast<double>(x[a][m])));
            }
        }
        double       step = 0.5*std::numeric_limits<real>::epsilon()*extent;
        auto         xr   = x;
        unsigned int seed = 1993;
        for (auto a : theAtoms)
        {
            for (int m = 0; m < DIM; m++)
            {
                seed     = seed*1103515245u + 12345u;
                bool up  = (seed >> 16) & 1;
                real to  = up ? std::numeric_limits<real>::max() : std::numeric_limits<real>::lowest();
                xr[a][m] = std::nextafter(static_cast<real>(xr[a][m] + (up ? step : -step)), to);
            }
        }
        std::vector<gmx::RVec>            fr(f.size());
        std::map<InteractionType, double> er;
        forceComp->compute(msghandler, pd, topology, &xr, &fr, &er);
        constrainForces(xr, &fr);
        std::vector<gmx::RVec> df(f.size(), { 0, 0, 0 });
        for (auto a : theAtoms)
        {
            rvec_sub(fr[a], f[a], df[a]);
        }
        return msForce(theAtoms, df);
    };
    // Constrain coordinates, using the reference coordinates for the directions
    auto constrainCoords = [&](const std::vector<gmx::RVec> &xref,
                               std::vector<gmx::RVec>       *x)
    {
        if (constrain)
        {
            constrainer.shake(msghandler, xref, x, nullptr);
        }
    };
    // Two sets of coordinates
    std::vector<gmx::RVec> newCoords[2];
    newCoords[0] = xstart;
    newCoords[1] = xstart;
    std::map<InteractionType, double> newEnergies[2];
    int current = 0;
    // Why the minimization did not converge
    auto failure = eMinimizeStatus::TooManySteps;
    if (simConfig.minAlg() == eMinimizeAlgorithm::LBFGS)
    {
        lbfgs = new StlbfgsHandler(pd, msghandler, topology, forceComp, theAtoms, xstart,
                                   constrain ? &constrainer : nullptr);
        // The energies and forces as in func, the shells are minimized
        // with the atoms
        lbfgs->setEnergyNoise(energyNoise(theAtoms, xstart,
                                          [&](std::vector<gmx::RVec> *xx, std::vector<gmx::RVec> *ff)
                                          {
                                              std::map<InteractionType, double> ee;
                                              const gmx::RVec field = { 0, 0, 0 };
                                              forceComp->constructVsiteCoordinates(topology, xx);
                                              forceComp->computeOnce(msghandler, pd, topology, xx, ff, &ee, field);
                                              forceComp->spreadVsiteForces(topology, xx, ff);
                                              constrainForces(*xx, ff);
                                              return ee[InteractionType::EPOT];
                                          }));

        STLBFGS::Optimizer                      opt{func, 1000};
        // The optimizer stops when the norm of the gradient divided by that
        // of the coordinates is below gtol, the norm of the gradient is the
        // square root of the number of particles times the mean square
        // force. The tolerance is the one requested, also if it is below
        // that of the shells, such that a normal mode analysis can get as
        // close to the minimum as possible. The relative decrease of the
        // energy only stops it when the energy cannot decrease within the
        // precision of the energy, a sum of many terms, since the energy
        // can decrease little while the forces are still large.
        double xnorm = 0;
        for (auto a : theAtoms)
        {
            xnorm += iprod(xstart[a], xstart[a]);
        }
        opt.gtol    = 0.9*std::sqrt(simConfig.forceTolerance()*theAtoms.size())/std::max(1.0, std::sqrt(xnorm));
        opt.ftol    = 1e-12;
        opt.verbose = msghandler->verbose();
        // Zero is until convergence
        opt.maxiter = simConfig.maxIter() > 0 ? simConfig.maxIter() : std::numeric_limits<int>::max();
        bool reachedMaxIter = false;

        // One-dimensional array
        std::vector<double>                     sx(theAtoms.size()*DIM);
        std::random_device                      rd;
        std::mt19937                            gen(rd());
        std::uniform_real_distribution<double>  dis(std::uniform_real_distribution<>(-1.0, 1.0));

        double displacement = simConfig.minimizeDisplacement();
        int maxRetry        = simConfig.minimizeRetries();
        int retry           = 0;
        // The lowest energy found, a structure with an energy that is not
        // finite is never stored
        double eMin         = std::numeric_limits<double>::max();
        bool   haveMinimum  = false;
        while (retry < maxRetry && (!converged || simConfig.forceReminimize()))
        {
            if (retry > 0)
            {
                msghandler->msg(ACTStatus::Info,
                                "Will retry minimization with slightly modified input coordinates.\n");
            }
            for(size_t i = 0; i < theAtoms.size(); i++)
            {
                for(int j = 0; j < DIM; j++)
                {
                    sx[DIM*i+j] = xstart[theAtoms[i]][j] + retry * dis(gen) * displacement;
                }
            }
            // The line search fails when the energy cannot decrease within
            // the precision, then the structure is a minimum if the forces
            // are within the tolerance. It also fails along very soft
            // directions, for instance near a saddle point, where the
            // energy decreases by about its rounding error while the
            // forces are large. Then the minimization continues from the
            // structure with the lowest energy, with steepest descent, as
            // long as that lowers the energy by more than its rounding
            // error, a limited number of times.
            const int maxRestart = 10;
            int    maxIter  = opt.maxiter;
            bool   done     = opt.run(sx);
            int    niter    = opt.niter;
            double eRestart = std::numeric_limits<double>::max();
            for (int restart = 0;
                 !done && restart < maxRestart && niter < maxIter &&
                 lbfgs->bestEnergy() < eRestart - lbfgs->energyNoise();
                 restart++)
            {
                eRestart = lbfgs->bestEnergy();
                lbfgs->restartFromBest();
                lbfgs->forceComp()->compute(msghandler, lbfgs->pd(), topology,
                                            lbfgs->coordinates(),
                                            lbfgs->forces(), lbfgs->energies());
                constrainForces(lbfgs->coordinates()->toRVec(), lbfgs->forces());
                double msf = msForce(theAtoms, *lbfgs->forces());
                if (msf <= msForceToler ||
                    msf <= 4*roundingMsForce(lbfgs->coordinates()->toRVec(), *lbfgs->forces()))
                {
                    break;
                }
                msghandler->msg(ACTStatus::Verbose,
                                gmx::formatString("Line search failed after %d iterations with mean square force %g, continuing from the lowest energy %g\n",
                                                  niter, msf, eRestart));
                for(size_t i = 0; i < theAtoms.size(); i++)
                {
                    for(int j = 0; j < DIM; j++)
                    {
                        sx[DIM*i+j] = (*lbfgs->coordinates())[theAtoms[i]][j];
                    }
                }
                opt.invH.reset();
                opt.maxiter = maxIter - niter;
                done        = opt.run(sx);
                niter      += opt.niter;
            }
            opt.maxiter    = maxIter;
            reachedMaxIter = reachedMaxIter || niter >= maxIter;
            if (!done && !lbfgs->bestCoordinates().empty())
            {
                // The last evaluation of a failed line search is not the
                // structure with the lowest energy
                lbfgs->restartFromBest();
            }
            // Final minimization for shells.
            lbfgs->forceComp()->compute(msghandler, lbfgs->pd(), topology,
                                        lbfgs->coordinates(),
                                        lbfgs->forces(), lbfgs->energies());
            constrainForces(lbfgs->coordinates()->toRVec(), lbfgs->forces());
            converged = done || msForce(theAtoms, *lbfgs->forces()) <= msForceToler;
            // The line search also fails when the energy cannot decrease
            // within the precision of the coordinates, then the forces are
            // as small as their change for one unit in the last place
            if (!converged && niter < maxIter)
            {
                converged = (msForce(theAtoms, *lbfgs->forces()) <=
                             4*roundingMsForce(lbfgs->coordinates()->toRVec(), *lbfgs->forces()));
            }
            if (converged)
            {
                double enew = lbfgs->energy(InteractionType::EPOT);
                if (msghandler->verbose())
                {
                    double msf = 0;
                    for(const auto &f : *lbfgs->forces())
                    {
                        msf += iprod(f,f);
                    }
                    msf /= theAtoms.size();
                    msghandler->write(gmx::formatString("Minimization of %zu atoms, iteration %d/%d energy %g RMS force %g\n",
                                                        theAtoms.size(), retry+1, maxRetry, enew, std::sqrt(msf)));
                }
                if (enew < eMin)
                {
                    // Store new structure only when it has lower energy.
                    newCoords[current]   = lbfgs->coordinates()->toRVec();
                    newEnergies[current] = *lbfgs->energies();
                    eMin                 = enew;
                    haveMinimum          = true;
                }
            }
            retry += 1;
        }
        // Only a stored structure is a minimum
        converged = haveMinimum;
        if (!converged)
        {
            failure = reachedMaxIter ? eMinimizeStatus::TooManySteps : eMinimizeStatus::NoMinimum;
            if (!lbfgs->bestCoordinates().empty())
            {
                // The structure with the lowest energy, although it is not
                // a minimum, with relaxed shells
                newCoords[current] = lbfgs->bestCoordinates().toRVec();
                std::vector<gmx::RVec> forces(newCoords[current].size());
                forceComp->compute(msghandler, pd, topology, &newCoords[current],
                                   &forces, &newEnergies[current]);
            }
        }
        // If applicable, print the force between the frozen two atoms
        lbfgs->printFrozenForce(msghandler->tw(), myFreeze);

        delete lbfgs;
    }
    else
    {
        // Newton or steepest descent, from the coordinates that satisfy
        // the constraints
        // The coordinates are integers, such that the steps are exact, see
        // coordinate.h, the constraints work with coordinates in nm
        Coordinates                       x(xstart), xt(xstart);
        std::vector<gmx::RVec>            f(myatoms.size()), ft(myatoms.size());
        std::map<InteractionType, double> e, et;
        auto evaluate = [&](Coordinates                       *xx,
                            std::vector<gmx::RVec>            *ff,
                            std::map<InteractionType, double> *ee)
        {
            forceComp->compute(msghandler, pd, topology, xx, ff, ee);
            constrainForces(xx->toRVec(), ff);
            return (*ee)[InteractionType::EPOT];
        };
        // The energy at x plus scale times the step, in xt, ft and et
        const auto tryStep = [&](const std::vector<double> &dx, double scale)
        {
            xt = x;
            size_t k = 0;
            for (auto atom : theAtoms)
            {
                for (int m = 0; m < DIM; m++)
                {
                    xt[atom][m] += scale*dx[k++];
                }
            }
            if (constrain)
            {
                auto xr = xt.toRVec();
                constrainCoords(x.toRVec(), &xr);
                xt.assign(xr);
            }
            return evaluate(&xt, &ft, &et);
        };
        // The rounding error of the energy, the shells are relaxed
        double noise = energyNoise(theAtoms, x.toRVec(),
                                   [&](std::vector<gmx::RVec> *xx, std::vector<gmx::RVec> *ff)
                                   {
                                       std::map<InteractionType, double> ee;
                                       forceComp->compute(msghandler, pd, topology, xx, ff, &ee);
                                       constrainForces(*xx, ff);
                                       return ee[InteractionType::EPOT];
                                   });
        const auto accept = [&]()
        {
            std::swap(x, xt);
            std::swap(f, ft);
            std::swap(e, et);
        };
        // The energy change from x to xt
        const auto energyChange = [&](double enew)
        {
            return alexandria::energyChange(theAtoms, e, x, f, enew, xt, ft, noise);
        };
        // The energy along the path of the minimization, the sum of the
        // energy changes
        double epot = evaluate(&x, &f, &e);
        double msf  = msForce(theAtoms, f);
        // The structure with the lowest energy, steepest descent may
        // accept steps that increase the energy
        double                            bestEpot     = epot;
        Coordinates                       bestX        = x;
        std::map<InteractionType, double> bestEnergies = e;
        // Steps since the lowest energy decreased more than the precision
        // of the energy, a sum of many terms, as with L-BFGS
        int  noProgress = 0;
        const auto keepBest = [&]()
        {
            if (bestEpot - epot > 1e-12*std::max({ std::abs(bestEpot), std::abs(epot), 1.0 }))
            {
                noProgress = 0;
            }
            else
            {
                noProgress += 1;
            }
            if (epot < bestEpot)
            {
                bestEpot     = epot;
                bestX        = x;
                bestEnergies = e;
            }
        };
        // The number of steps without progress after which there is no
        // minimum. Steepest descent may accept steps that increase the
        // energy, and needs more of them.
        const int maxNoProgress = simConfig.minAlg() == eMinimizeAlgorithm::Newton ? 10 : 100;
        if (energies)
        {
            *energies = e;
        }
        printEnergies(msghandler->tw(), myIter, msf, 0, 0, e, 0);
        // The forces do not become smaller than their change for one unit
        // in the last place of the coordinates, which in single precision
        // is larger than the tolerance. It is estimated at the start and
        // computed again when the forces are that small.
        double floorMsf    = roundingMsForce(x.toRVec(), f);
        const auto isConverged = [&]()
        {
            if (msf <= msForceToler)
            {
                return true;
            }
            if (msf <= 4*floorMsf)
            {
                floorMsf = roundingMsForce(x.toRVec(), f);
                return msf <= 4*floorMsf;
            }
            return false;
        };
        converged = isConverged();
        // The largest step of a particle and the starting one, in nm. The
        // step is limited to the trust radius, which grows after a step that
        // lowered the energy and shrinks after one that did not.
        const double maxRadius = 0.1;
        double       radius    = 0.01;
        size_t       n         = DIM*theAtoms.size();
        const auto moreIterations = [&]()
        {
            return myIter < simConfig.maxIter() || 0 == simConfig.maxIter();
        };
        if (simConfig.minAlg() == eMinimizeAlgorithm::Newton)
        {
            // Without frozen particles the molecule can move and rotate freely
            bool        rigidBody = myFreeze.empty();
            NewtonStep  newton;
            bool        haveHessian = false;
            double      mu          = 0;
            while (!converged && moreIterations())
            {
                if (!haveHessian)
                {
                    MatrixWrapper                     hessian(n, n);
                    std::vector<double>               f0;
                    std::map<InteractionType, double> e0;
                    computeHessian(pd, topology, forceComp, &x, theAtoms, &hessian, &f0, &e0);
                    newton.decompose(hessian, f0,
                                     fixedDirections(theAtoms, x.toRVec(), rigidBody,
                                                     constrain ? &constrainer : nullptr));
                    haveHessian = true;
                }
                const auto dx    = newton.step(mu, simConfig.overRelax());
                double dmax  = maxDisplacement(dx);
                double scale = (dmax > radius) ? radius/dmax : 1.0;
                double enew  = tryStep(dx, scale);
                double de    = std::isfinite(enew) ? energyChange(enew) : enew;
                if (std::isfinite(de) && de <= 0)
                {
                    epot       += de;
                    accept();
                    keepBest();
                    if (noProgress >= maxNoProgress)
                    {
                        failure = eMinimizeStatus::NoMinimum;
                        break;
                    }
                    haveHessian = false;
                    radius      = std::min(maxRadius, 2*radius);
                    mu         /= 10;
                    if (mu < 1e-12*newton.largestEigenvalue())
                    {
                        mu = 0;
                    }
                }
                else
                {
                    // The same Hessian with a shorter, more damped step
                    radius /= 2;
                    mu      = std::max(10*mu, 1e-6*newton.largestEigenvalue());
                    if (radius < 1e-9)
                    {
                        failure = eMinimizeStatus::NoMinimum;
                        break;
                    }
                }
                myIter += 1;
                msf       = msForce(theAtoms, f);
                converged = isConverged();
                printEnergies(msghandler->tw(), myIter, msf, 0, 0, e, scale);
            }
        }
        else
        {
            // Steepest descent with the step sizes of Barzilai and Borwein,
            // which only work well if the energy may increase in some steps.
            // A step is accepted if the energy is lower than the highest of
            // the last ones, with the rule of Grippo, Lampariello and Lucidi,
            // otherwise the step is halved.
            const size_t           nlast = 10;
            std::deque<double>     last  = { epot };
            std::vector<double>    d(n), dprev(n), xprev(n), xcur(n);
            double                 gamma = 0;
            const auto flatten = [&](const auto &v, std::vector<double> *flat)
            {
                size_t k = 0;
                for (auto atom : theAtoms)
                {
                    for (int m = 0; m < DIM; m++)
                    {
                        (*flat)[k++] = v[atom][m];
                    }
                }
            };
            while (!converged && moreIterations())
            {
                flatten(f, &d);
                flatten(x, &xcur);
                double dmax = maxDisplacement(d);
                if (0 == myIter)
                {
                    gamma = radius/dmax;
                }
                else
                {
                    double sy = 0, yy = 0;
                    for (size_t k = 0; k < n; k++)
                    {
                        double y = d[k] - dprev[k];
                        sy += (xcur[k] - xprev[k])*y;
                        yy += y*y;
                    }
                    if (yy > 0 && sy != 0)
                    {
                        gamma = std::abs(sy)/yy;
                    }
                }
                gamma = std::min(gamma, maxRadius/dmax);
                double eref = *std::max_element(last.begin(), last.end());
                double dd   = std::inner_product(d.begin(), d.end(), d.begin(), 0.0);
                bool   ok   = false;
                for (int halving = 0; !ok && halving < 40; halving++)
                {
                    double enew = tryStep(d, gamma);
                    double de   = std::isfinite(enew) ? energyChange(enew) : enew;
                    ok = std::isfinite(de) && epot + de <= eref - 1e-4*gamma*dd;
                    if (ok)
                    {
                        epot += de;
                    }
                    else
                    {
                        gamma /= 2;
                    }
                }
                if (!ok)
                {
                    failure = eMinimizeStatus::NoMinimum;
                    break;
                }
                xprev = xcur;
                dprev = d;
                accept();
                keepBest();
                if (noProgress >= maxNoProgress)
                {
                    failure = eMinimizeStatus::NoMinimum;
                    break;
                }
                last.push_back(epot);
                if (last.size() > nlast)
                {
                    last.pop_front();
                }
                myIter += 1;
                msf       = msForce(theAtoms, f);
                converged = isConverged();
                printEnergies(msghandler->tw(), myIter, msf, 0, 0, e, gamma);
            }
        }
        // Without progress the energy may not decrease within the precision
        // of the coordinates, then the forces are as small as their change
        // for one unit in the last place
        if (!converged && eMinimizeStatus::NoMinimum == failure)
        {
            msf       = msForce(theAtoms, f);
            floorMsf  = roundingMsForce(x.toRVec(), f);
            converged = msf <= 4*floorMsf;
        }
        // A converged structure is the last one, which may not have the
        // lowest energy with steepest descent
        newCoords[current]   = (converged ? x : bestX).toRVec();
        newEnergies[current] = converged ? e : bestEnergies;
    }
    // Also if it is not a minimum, the structure with the lowest energy
    *coords = newCoords[current];
    if (energies && !newEnergies[current].empty())
    {
        *energies = newEnergies[current];
    }
    return converged ? eMinimizeStatus::OK : failure;
}

double MolHandler::coordinateRmsd(const std::vector<ActAtom>   &myatoms,
                                  const std::vector<gmx::RVec> &xref,
                                  std::vector<gmx::RVec>       *xfit) const
{
    GMX_RELEASE_ASSERT(myatoms.size() == xref.size(), "Size mismatch in arrays");
    GMX_RELEASE_ASSERT(myatoms.size() == xfit->size(), "Size mismatch in arrays");
    // Compute RMSD
    std::vector<real> myMass;
    double            tmass      = 0;
    int               nRealAtoms = 0;
    for(const auto &ma: myatoms)
    {
        auto mass = ma.mass();
        myMass.push_back(mass);
        if (mass > 0)
        {
            tmass      += mass;
            nRealAtoms += 1;
        }
    }
    rvec ref_com = { 0, 0, 0 };
    rvec fit_com = { 0, 0, 0 };
    for (size_t i = 0; i < xref.size(); ++i)
    {
        double relativeMass = myMass[i]/tmass;
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
    do_fit(myMass, xref_com, xfit);
    
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
    return std::sqrt(msd/nRealAtoms);
}

//! \brief Write energy legend to xvg file
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

//! \brief Compute square norm of vector
static double rvecNormSquared(const gmx::RVec &v)
{
    return v[XX]*v[XX]+v[YY]*v[YY]+v[ZZ]*v[ZZ];
}

//! \brief Print energies to a file
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
    double temp = 0;
    if (nDOF > 0)
    {
        temp = 2*ekin/(nDOF*BOLTZ);
    }
    fprintf(exvg, " %.5f %.5f %.5f\n", ekin, ekin+epot, temp);
}

/*! \brief Positions relative to an atom of each group of constrained atoms
 * The constraints are computed with coordinates in double precision in
 * nm, the differences of the integer coordinates, see coordinate.h,
 * relative to one atom of the atoms that constraints connect. They are
 * as precise as double allows for the size of a molecule, in either
 * precision. The constrained particles are put at their anchors plus
 * their relative positions.
 */
class ConstraintFrame
{
public:
    /*! \brief Constructor
     * \param[in] pairs The constrained atom pairs
     * \param[in] natom The number of particles
     */
    ConstraintFrame(const std::vector<std::array<int, 2>> &pairs, size_t natom) : anchor_(natom)
    {
        // The atom with the lowest index of each group of connected atoms
        for (size_t i = 0; i < natom; i++)
        {
            anchor_[i] = i;
        }
        const auto root = [this](size_t i)
        {
            while (anchor_[i] != i)
            {
                i = anchor_[i];
            }
            return i;
        };
        for (const auto &p : pairs)
        {
            size_t a = root(p[0]);
            size_t b = root(p[1]);
            anchor_[std::max(a, b)] = std::min(a, b);
        }
        for (size_t i = 0; i < natom; i++)
        {
            anchor_[i] = root(i);
        }
    }

    /*! \brief The positions relative to their anchors
     * \param[in] x The coordinates
     * \return the positions in nm
     */
    std::vector<gmx::DVec> relative(const Coordinates &x) const
    {
        std::vector<gmx::DVec> r(x.size());
        for (size_t i = 0; i < x.size(); i++)
        {
            for (int m = 0; m < DIM; m++)
            {
                r[i][m] = fromCoordinate(coordinateDifference(x[i].units(m), x[anchor_[i]].units(m)));
            }
        }
        return r;
    }

    /*! \brief Set the coordinates from relative positions
     * The particles that are their own anchors keep their coordinates,
     * the others are put at the anchors plus their relative positions.
     * The anchors move by their relative positions as well.
     * \param[in]    r          The relative positions
     * \param[inout] x          The coordinates
     * \param[inout] velocities If not nullptr, these change by the displacements over dt
     * \param[in]    invdt      The inverse time step
     */
    void apply(const std::vector<gmx::DVec> &r,
               Coordinates                  *x,
               std::vector<gmx::RVec>       *velocities = nullptr,
               double                        invdt      = 0) const
    {
        Coordinates x0 = *x;
        for (size_t i = 0; i < x->size(); i++)
        {
            for (int m = 0; m < DIM; m++)
            {
                (*x)[i][m]  = x0[anchor_[i]][m];
                (*x)[i][m] += r[i][m];
            }
        }
        if (velocities)
        {
            for (size_t i = 0; i < x->size(); i++)
            {
                gmx::RVec dx = (*x)[i] - x0[i];
                for (int m = 0; m < DIM; m++)
                {
                    (*velocities)[i][m] += dx[m]*invdt;
                }
            }
        }
    }

private:
    //! The anchor of each particle
    std::vector<size_t> anchor_;
};

//! \brief Make arrays the length of the number of atoms
static void initArrays(const ACTMol                 *mol,
                       const std::vector<gmx::RVec> &xmol,
                       Coordinates                  *coordinates,
                       std::vector<gmx::RVec>       *velocities,
                       std::vector<gmx::RVec>       *forcesCur,
                       std::vector<gmx::RVec>       *forcesPrev,
                       std::vector<double>          *inv2mass,
                       std::vector<double>          *mass_2)
{
    *coordinates = Coordinates(xmol);
    velocities->resize(xmol.size());
    forcesPrev->resize(xmol.size());
    forcesCur->resize(xmol.size());
    gmx::RVec zero = { 0, 0, 0 };
    auto &atoms = mol->atomsConst();
    for (size_t i = 0; i < xmol.size(); i++)
    {
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
                       
void MolHandler::simulate(MsgHandler                    *msghandler,
                          const ForceField              *pd,
                          ACTMol                        *mol,
                          const std::vector<gmx::RVec>  &xmol,
                          const ForceComputer           *forceComp,
                          const SimulationConfigHandler &simConfig,
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
    // The coordinates are integers, such that the steps are exact, see
    // coordinate.h
    Coordinates                        coordinates, xprev;
    std::vector<gmx::RVec>             velocities, forces[2];
    std::map<InteractionType, double>  energies;
    std::vector<double>                inv2mass;
    std::vector<double>                mass_2;

    // To swap between forces efficiently
    int cur = 0;
#define prev (1-cur)

    // Initiate coordinates, etc.
    initArrays(mol, xmol, &coordinates, &velocities,
               &forces[cur], &forces[prev], &inv2mass, &mass_2);

    auto deltat  = simConfig.deltat();
    Constrainer constrainer(simConfig.constraintMaxIter(),
                            simConfig.constraintTolerance(),
                            1.0/deltat);
    // The motion of the center of mass that is removed from the initial
    // velocities, and every nstcomm steps, removes three degrees of
    // freedom for the translation and two or three for the rotation.
    // Without external forces the center of mass does not move in
    // between. Each constraint removes one more degree of freedom.
    int nDOF           = DIM*mol->nRealAtoms();
    // The particles that move, for counting the floating point operations
    double nmoving     = std::count_if(inv2mass.begin(), inv2mass.end(),
                                       [](double m) { return m > 0; });
    if (simConfig.constrain())
    {
        constrainer.makeConstraints(msghandler, simConfig.constraintType(),
                                    pd, mol->topology(), coordinates.toRVec());
        msghandler->write(gmx::formatString("Constraining %zu bonds and %zu angles.",
                                            constrainer.numBondConstraints(),
                                            constrainer.numAngleConstraints()));
        if (constrainer.numRedundantAngles() > 0)
        {
            msghandler->write(gmx::formatString("Left out %zu angle constraints that are redundant.",
                                                constrainer.numRedundantAngles()));
        }
        nDOF -= constrainer.numConstraints();
    }
    const ConstraintFrame frame(constrainer.pairs(), coordinates.size());
    if (simConfig.constrain())
    {
        // Start from coordinates that satisfy the constraints
        auto r     = frame.relative(coordinates);
        int  error = constrainer.shakeDouble(msghandler, r, &r);
        frame.apply(r, &coordinates);
        if (error != 0)
        {
            // With angles the reference geometry may not be possible,
            // otherwise the tolerance may be smaller than the precision
            // allows
            std::string reason;
            if (constrainer.numAngleConstraints() > 0)
            {
                reason = " The reference bond lengths and angles may be inconsistent, for instance in rings, then try -constrain bonds.";
            }
            msghandler->fatal(gmx::formatString("Could not constrain the starting coordinates within the relative tolerance %g.%s Try a larger tolerance with -constr_toler.",
                                                simConfig.constraintTolerance(), reason.c_str()));
        }
    }
    const COMMRemover commRemover(simConfig.commRemoval(), mol->atomsConst(), coordinates.toRVec());
    nDOF -= commRemover.nDOF();
    if (COMMRemoval::None != commRemover.mode())
    {
        msghandler->write(gmx::formatString("Removing the %s the center of mass %s, which removes %d degrees of freedom.",
                                            COMMRemoval::Angular == commRemover.mode() ? "translation of and rotation around" : "translation of",
                                            simConfig.nstcomm() > 0 ? gmx::formatString("every %d steps", simConfig.nstcomm()).c_str() : "from the initial velocities",
                                            commRemover.nDOF()));
    }
    // Generate velocities
    if (simConfig.temperature() > 0)
    {
        maxwell_speed(simConfig.temperature(), simConfig.seed(),
                      mol->atomsConst(), &velocities, msghandler->tw());
        commRemover.remove(coordinates.toRVec(), &velocities);
        if (simConfig.constrain())
        {
            constrainer.rattleDouble(msghandler, frame.relative(coordinates), &velocities);
        }
        // Scale the velocities to the requested temperature for the
        // remaining degrees of freedom.
        double ekin = 0;
        for(size_t i = 0; i < xmol.size(); i++)
        {
            ekin += mass_2[i]*rvecNormSquared(velocities[i]);
        }
        if (ekin > 0 && nDOF > 0)
        {
            double scale = std::sqrt(0.5*nDOF*BOLTZ*simConfig.temperature()/ekin);
            for(auto &v : velocities)
            {
                svmul(scale, v, v);
            }
        }
    }
    auto deltat2 = deltat*deltat;
    // For trajectory writing
    matrix      box   = { { 4, 0, 0 }, { 0, 4, 0 }, { 0, 0, 4 } };
    char        chain = ' ';
    const char *title = "ACT trajectory";
    std::vector<int> trajIndex;
    {
        const auto &atoms = mol->atomsConst();
        for(size_t ii = 0; ii < atoms.size(); ii++)
        {
            if (simConfig.writeShells() || ActParticle::Atom == atoms[ii].pType())
            {
                trajIndex.push_back(ii);
            }
        }
    }
    msghandler->write("\nWill start simulation.\n");
    // The forces are timed where they are computed
    auto               &timings          = Timings::instance();
    static const size_t integrationTimer = timings.index("MD integration");
    static const size_t constraintTimer  = timings.index("MD constraints");
    static const size_t energyTimer      = timings.index("MD energy output");
    static const size_t trajectoryTimer  = timings.index("MD trajectory output");
    // Compute forces and energies for the starting coordinates
    forceComp->compute(msghandler, pd, mol->topology(), &coordinates,
                       &forces[prev], &energies);
    // After that the shells start from a prediction of their relaxed
    // positions, extrapolated from their positions relative to their
    // parent particles in the two previous steps, rather than from their
    // parent particles. Since the induction is the change of the energy
    // from the shells on their parents, these energies are computed as
    // well in the steps of which the energies are written.
    std::vector<std::pair<size_t, size_t>> shellParent;
    {
        auto &atoms = mol->atomsConst();
        for (size_t i = 0; i < atoms.size(); i++)
        {
            if (ActParticle::Shell == atoms[i].pType() && atoms[i].parents().size() == 1)
            {
                shellParent.emplace_back(i, static_cast<size_t>(atoms[i].parents()[0]));
            }
        }
    }
    // The shell positions relative to their parents in this and the previous step
    std::vector<gmx::RVec> shellOffset(shellParent.size()), shellOffsetPrev(shellParent.size());
    bool                   haveOffsetPrev = false;
    for (int step = 0; step < simConfig.nsteps(); step++)
    {
        for (size_t k = 0; k < shellParent.size(); k++)
        {
            shellOffsetPrev[k] = shellOffset[k];
            shellOffset[k] = coordinates[shellParent[k].first] - coordinates[shellParent[k].second];
        }
        // Print the energies at the current time, i.e. before integrating
        if (simConfig.nstener() > 0 && step % simConfig.nstener() == 0)
        {
            const ScopedTimer timer(energyTimer);
            double ekin = 0;
            for(size_t i = 0; i < xmol.size(); i++)
            {
                ekin += mass_2[i]*rvecNormSquared(velocities[i]);
            }
            printEnergies(exvg, step, deltat, energies, ekin, nDOF, oenv);
        }
        // Write the coordinates of the same time, numbering the frames from 1
        if (simConfig.nstxout() > 0 && step % simConfig.nstxout() == 0)
        {
            const ScopedTimer timer(trajectoryTimer);
            //! \todo rewrite without gromacs fluff.
            pdbWriter(traj, title, mol->atomsConst(),
                      coordinates.toRVec(), mol->topology()->residueNames(),
                      epbcNONE,
                      box, chain, step/simConfig.nstxout()+1, trajIndex,
                      nullptr, true);
        }
        // Integrate positions using velocity Verlet
        timings.start(integrationTimer);
        if (simConfig.constrain())
        {
            xprev = coordinates;
        }
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
            }
        }
        timings.addFlops(13*nmoving);
        timings.stop(integrationTimer);
        if (simConfig.constrain())
        {
            const ScopedTimer timer(constraintTimer);
            // First part of rattle, constrain the positions along the
            // previous bonds and correct the velocities accordingly.
            // The velocities change by the displacements of the integer
            // coordinates, such that they stay consistent
            auto r     = frame.relative(coordinates);
            int  error = constrainer.shakeDouble(msghandler, frame.relative(xprev), &r);
            frame.apply(r, &coordinates, &velocities, 1.0/deltat);
            if (error != 0)
            {
                msghandler->msg(ACTStatus::Warning,
                                gmx::formatString("SHAKE did not converge at step %d.", step));
            }
        }
        // Predict the shell positions
        for (size_t k = 0; k < shellParent.size(); k++)
        {
            gmx::RVec offset = shellOffset[k];
            if (haveOffsetPrev)
            {
                for (int m = 0; m < DIM; m++)
                {
                    offset[m] = 2*shellOffset[k][m] - shellOffsetPrev[k][m];
                }
            }
            for (int m = 0; m < DIM; m++)
            {
                coordinates[shellParent[k].first][m]  = coordinates[shellParent[k].second][m];
                coordinates[shellParent[k].first][m] += offset[m];
            }
        }
        haveOffsetPrev = true;
        // Compute forces and energies, for the induction energy also with
        // the shells on their parents if the energies are written
        bool writeEnergies = simConfig.nstener() > 0 && (step + 1) % simConfig.nstener() == 0;
        forceComp->compute(msghandler, pd, mol->topology(), &coordinates,
                           &forces[cur], &energies, { 0, 0, 0 }, writeEnergies, {}, writeEnergies);
        // Integrate velocities using velocity Verlet
        timings.start(integrationTimer);
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
        timings.addFlops(12*nmoving);
        timings.stop(integrationTimer);
        if (simConfig.constrain())
        {
            const ScopedTimer timer(constraintTimer);
            if (constrainer.rattleDouble(msghandler, frame.relative(coordinates), &velocities) != 0)
            {
                msghandler->msg(ACTStatus::Warning,
                                gmx::formatString("RATTLE did not converge at step %d.", step));
            }
        }
        // Remove the motion of the center of mass, which leaves the
        // constraints satisfied, since it is a rigid body motion.
        if (simConfig.nstcomm() > 0 && (step + 1) % simConfig.nstcomm() == 0)
        {
            const ScopedTimer timer(integrationTimer);
            commRemover.remove(coordinates.toRVec(), &velocities);
            timings.addFlops(60*nmoving);
        }
        // Swap force arrays
        cur = prev;
    }
    
    gmx_ffclose(traj);
    xvgrclose(exvg);
    msghandler->write("Simulation finished correctly.");
}

} // namespace alexandria

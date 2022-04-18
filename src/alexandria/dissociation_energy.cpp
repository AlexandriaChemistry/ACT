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
 */

#include "dissociation_energy.h"

#include <cstdio>

#include <random>

#include "gromacs/statistics/statistics.h"
#include "gromacs/utility/futil.h"

#include "act/utility/regression.h"

namespace alexandria
{
    
/*! \brief Write a csv file containing molecule names and bond energy
 *
 * Writes the whole bond energy matrix.
 * \param[in] csvFile             Filename for csv file for debugging, may be a nullptr
 * \param[in] used                Gives a key, the index in molecule vector of those compounds that
 *                                are used, and the value is the number of times it is used
 *                                in the bootstrap.
 * \param[in] bondIdToIndex       Maps the bond identifiers to the column in the matrix
 * \param[in] mm                  The molecules
 * \param[in] a                   The matrix
 * \param[in] edissoc             Vector of bond energies
 * \param[in] ntrain              Data structure for counting the amount of test data per bond
 */
static void dump_csv(const char                      *csvFile,
                     const std::map<int, int>        &used,
                     const std::map<Identifier, int> &bondIdToIndex,
                     const std::vector<MyMol>        &mm,
                     const MatrixWrapper             &a,
                     const std::vector<double>       &edissoc,
                     const std::map<Identifier, int> *ntrain)
{
    std::map<int, Identifier> indexToBondId;
    for (auto &j : bondIdToIndex)
    {
        indexToBondId.insert(std::pair<int, Identifier>(j.second, j.first));
    }
    FILE *csv = gmx_ffopen(csvFile, "w");
    fprintf(csv, ",");
    for (auto &j : indexToBondId)
    {
        fprintf(csv, "%s,", j.second.id().c_str());
    }
    fprintf(csv, "DeltaE0\n");
    int row = 0;
    for (const auto &j : used)
    {
        auto mymol = &(mm[j.first]);
        fprintf(csv, "%s,", mymol->getMolname().c_str());
        for (size_t i = 0; i < edissoc.size(); i++)
        {
            fprintf(csv, "%g,", a.get(i, row));
        }
        double deltaE0;
        GMX_RELEASE_ASSERT(mymol->energy(MolPropObservable::DELTAE0, &deltaE0),
                           gmx::formatString("No DeltaE0 for %s",
                                             mymol->getMolname().c_str()).c_str());
        fprintf(csv, "%.3f\n", -deltaE0*j.second);
        row++;
    }
    if (ntrain)
    {
        fprintf(csv, "Total,");
        for (auto &j : indexToBondId)
        {
            auto nt = ntrain->find(j.second);
            GMX_RELEASE_ASSERT(ntrain->end() != nt, gmx::formatString("Cannot find %s in training counter", j.second.id().c_str()).c_str());
            fprintf(csv, "%d,", nt->second);
        }
        fprintf(csv, "\n");
    }
    fprintf(csv, "Edissoc,");
    for (auto &j : edissoc)
    {
        fprintf(csv, "%.3f,", j);
    }
    fprintf(csv, "\n");
    fclose(csv);
}

/*! \brief Calculate the dissociation energies once
 * \param[in] fplog               File pointer for logging information
 * \param[in] pd                  The input force field
 * \param[in] molset              The molecules
 * \param[in] pickRandomMolecules Whether or not to pick random molecules or just everything
 * \param[in] hasExpData          Vector of indices pointing to those molecules for which 
 *                                there is data 
 * \param[inout] edissoc          Map from the bond identifier to a statistics container
 * \param[in] gen                 Random number generator. Must be a pointer, otherwise the 
 *                                internal data structure is not updated and the code will
 *                                repeat the same random sequence for each invocation
 * \param[in] uniform             Uniform distribution between 0 and 1
 * \param[in] csvFile             Filename for csv file for debugging, may be a nullptr
 * \param[in] ntrain              Data structure for counting the amount of test data per bond
 * \param[out] rmsd               Root mean square deviation, if not nullptr
 * \return true if successful.
 */
static bool calcDissoc(FILE                              *fplog,
                       const Poldata                     *pd,
                       const std::vector<MyMol>          &molset,
                       bool                               pickRandomMolecules,
                       const std::vector<int>            &hasExpData,
                       std::map<Identifier, gmx_stats>   *edissoc,
                       std::mt19937                      *gen,  
                       std::uniform_real_distribution<>   uniform,
                       const char                        *csvFile,
                       std::map<Identifier, int>         *ntrain,
                       double                            *rmsd)
{
    // Determine which molecules to use
    std::map<int, int>  used;
    size_t              nExpData = hasExpData.size();
    for (size_t i = 0; i < nExpData; i++)
    {
        // By default, copy the input array
        int mytry = hasExpData[i];
        if (pickRandomMolecules)
        {
            // Pick nMol random compounds from the set for which we have exp data
            mytry = hasExpData[int(nExpData*uniform(*gen)) % nExpData];
        }
        if (used.find(mytry) == used.end())
        {
            used.insert(std::pair<int, int>(mytry, 1));
        }
        else
        {
            used.find(mytry)->second += 1;
        }
    }
    // Total number of compounds in the matrix is equalt to the number of rows
    int nRow = used.size();
    
    // Now time to find out which bonds are present in this subset of compounds
    std::map<Identifier, int> bondIdToIndex;
    int                       nColumn = 0;
    for (const auto &uu : used)
    {
        auto mymol   = &(molset[uu.first]);
        auto myatoms = mymol->atomsConst();
        for (auto &b : mymol->bondsConst())
        {
            auto atypeI = *myatoms.atomtype[b.aI()];
            auto atypeJ = *myatoms.atomtype[b.aJ()];
            std::string btypeI, btypeJ;
            if (pd->atypeToBtype(atypeI, &btypeI) &&
                pd->atypeToBtype(atypeJ, &btypeJ))
            {
                Identifier bondId({btypeI, btypeJ}, { b.bondOrder() }, CanSwap::Yes);
                if (bondIdToIndex.find(bondId) == bondIdToIndex.end())
                {
                    bondIdToIndex.insert(std::pair<Identifier, int>(bondId, nColumn++));
                }
            }
            else
            {
                gmx_fatal(FARGS, "No parameters for bond in the force field, atoms %s-%s mol %s",
                          atypeI, atypeJ,
                          mymol->getIupac().c_str());
            }
        }
    }
    if (fplog)
    {
        fprintf(fplog, "There are %d different bondtypes and %d reference datapoints to optimize the heats of formation\n", nColumn,  nRow);
    }
    if (nColumn > nRow)
    {
        if (fplog)
        {
            fprintf(fplog, "Not enough data. Try again.\n");
        }
        return false;
    }
    // Now we can allocate our matrices and the right-hand-side
    MatrixWrapper       a(nColumn, nRow);
    MatrixWrapper       a_copy(nColumn, nRow);
    std::vector<double> rhs;
    int                 row = 0;

    // Now it is time to fill the matrices
    for (const auto &uu : used)
    {
        auto mymol = &(molset[uu.first]);
        auto myatoms = mymol->atomsConst();
        for (auto &b : mymol->bondsConst())
        {
            const char *atypeI = *myatoms.atomtype[b.aI()];
            const char *atypeJ = *myatoms.atomtype[b.aJ()];
            std::string btypeI, btypeJ;
            if (pd->atypeToBtype(atypeI, &btypeI) &&
                pd->atypeToBtype(atypeJ, &btypeJ))
            {
                Identifier bondId({btypeI, btypeJ}, { b.bondOrder() }, CanSwap::Yes);
                int column = bondIdToIndex[bondId];
                GMX_RELEASE_ASSERT(column < nColumn && column >= 0, gmx::formatString("Column %d should be within 0..%d", column, nColumn).c_str());
                GMX_RELEASE_ASSERT(row < nRow && row >= 0, gmx::formatString("Row %d should be within 0..%d", row, nRow).c_str());
                a.set(column, row, a.get(column, row) + uu.second);
                a_copy.set(column, row, a.get(column, row));
                if (ntrain)
                {
                    auto nnn = ntrain->find(bondId);
                    if (ntrain->end() == nnn)
                    {
                        ntrain->insert(std::pair<Identifier, int>(bondId, 0));
                        nnn = ntrain->find(bondId);
                    }
                    nnn->second += uu.second;
                }
            }
        }
        double deltaE0;
        GMX_RELEASE_ASSERT(mymol->energy(MolPropObservable::DELTAE0, &deltaE0),
                           gmx::formatString("No molecular energy for %s",
                                             mymol->getMolname().c_str()).c_str());
        rhs.push_back(-deltaE0 * uu.second);
        row += 1;
    }

    // Let's try and solve the equation.
    std::vector<double> Edissoc(nColumn, 0.0);
    if (0 == a.solve(rhs, &Edissoc))
    {
        // Check for large numbers:
        bool bLargeNumbers = false;
        for(auto &E : Edissoc)
        {
            if (fabs(E) > 1000)
            {
                bLargeNumbers = true;
            }
        }
        if (bLargeNumbers)
        {
            if (csvFile)
            {
                dump_csv(csvFile, used, bondIdToIndex, molset,
                         a_copy, Edissoc, ntrain);
            }
            return false;
        }
        // Compute rmsd
        if (rmsd)
        {
            double msd = 0;
            for(int i = 0; i < nRow; i++)
            {
                double result = 0;
                for(int j = 0; j < nColumn; j++)
                {
                    result += a_copy.get(j, i)*Edissoc[j];
                }
                msd += gmx::square(result-rhs[i]);
            }
            *rmsd = std::sqrt(msd/nRow);
        }
        // Copy to the output map.
        for (const auto &b : bondIdToIndex)
        {
            auto ed = edissoc->find(b.first);
            if (edissoc->end() == ed)
            {
                // New bond type!
                gmx_stats gs;
                edissoc->insert(std::pair<Identifier, gmx_stats>(b.first, std::move(gs)));
            }
            auto  &gs = edissoc->find(b.first)->second;
            size_t N  = gs.get_npoints();
            gs.add_point(N, Edissoc[b.second], 0, 0);
            if (fplog && fabs(Edissoc[b.second]) > 1000)
            {
                fprintf(fplog, "Adding energy %g for %s\n", Edissoc[b.second],
                        b.first.id().c_str());
            }
        }
        return true;
    }
    else
    {
        return false;
    }
}
                           
double getDissociationEnergy(FILE               *fplog,
                             Poldata            *pd,
                             std::vector<MyMol> *molset,
                             iqmType             iqm,
                             const char         *csvFile,
                             const std::string  &method,
                             const std::string  &basis,
                             int                 nBootStrap)
{
    std::random_device               rd;
    std::mt19937                     gen(rd());  
    std::uniform_real_distribution<> uniform(0.0, 1.0);
    std::map<Identifier, int>        bondIdToIndex;
    std::vector<int>                 hasExpData;
    std::map<MolPropObservable, iqmType> myprops = {
        { MolPropObservable::DELTAE0, iqm }
    };
    std::map<iqmType, double> tmap = {
        { iqmType::QM, 0 }, { iqmType::Both, -1 }
    };
    // Loop over molecules to find the ones with experimental DeltaHform
    for (size_t i = 0; i < molset->size(); i++)
    {
        auto mymol = &((*molset)[i]);
        if (immStatus::OK == mymol->getExpProps(myprops, method, basis, pd, 
                                                tmap[iqm]))
        {
            double deltaE0;
            if (mymol->energy(MolPropObservable::DELTAE0, &deltaE0))
            {
                hasExpData.push_back(i);
            }
        }
    }
    if (hasExpData.size() < 2)
    {
        fprintf(fplog, "Not enough molecules with experimental data to determine dissocation energy.\n");
        return -1;
    }
    // Call the low level routine once to get optimal values and to
    // establish all the entries in the edissoc map.
    std::map<Identifier, gmx_stats> edissoc;
    std::map<Identifier, int>       ntrain;
    double                          rmsd = 0;
    if (!calcDissoc(fplog, pd, *molset, false, hasExpData, &edissoc, &gen, uniform, csvFile, &ntrain, &rmsd))
    {
        gmx_fatal(FARGS, "Cannot solve the matrix equations for determining the dissociation energies");
    }
    // Now run the bootstrapping
    std::map<Identifier, gmx_stats> edissoc_bootstrap;
    int maxBootStrap = 2*nBootStrap;
    int nBStries = 0;
    int iter;
    for(iter = 0; iter < nBootStrap && nBStries < maxBootStrap; )
    {
        std::map<Identifier, int> nnn;
        if (calcDissoc(nullptr, pd, *molset, true, hasExpData,
                       &edissoc_bootstrap, &gen,
                       uniform, csvFile, &nnn, nullptr))
        {
            iter++;
        }
        nBStries++;
    }
    if (nBStries == maxBootStrap)
    {
        fprintf(fplog, "Maximum number of tries %d for running bootstraps reached.\n", maxBootStrap);
    }
    
    if (fplog)
    {
        fprintf(fplog, "Optimized dissociation energy based on %4d bootstraps.\n", iter);
        fprintf(fplog, "%-14s  %6s  %10s  %10s\n", "Bond", "N", "Edissoc", "Std.Dev.");
    }
    // Finally copy the new dissociation energies to the force field.
    auto iBonds = InteractionType::BONDS;
    GMX_RELEASE_ASSERT(pd->interactionPresent(iBonds), "No bonds in force field file");
    auto fs  = pd->findForces(InteractionType::BONDS);
    for (auto &bi : edissoc)
    {
        double average, error = 0;
        auto estats = bi.second.get_average(&average);
        if (eStats::OK != estats)
        {
            if (fplog)
            {
                fprintf(fplog, "%s: %s\n", bi.first.id().c_str(), 
                        gmx_stats_message(estats));
            }
            continue;
        }
        if (nBootStrap > 0)
        {
            auto ed = edissoc_bootstrap.find(bi.first);
            // We have to check whether this particular bond exists.
            // If there are few bootstraps, a rare bond may not be there.
            if (edissoc_bootstrap.end() != ed)
            {
                estats = edissoc_bootstrap[bi.first].get_sigma(&error);
                GMX_RELEASE_ASSERT(eStats::OK == estats, gmx_stats_message(estats));
            }
        }
        double delta = error;
        if (error == 0)
        {
            // Can happen when there is one data point only
            delta = std::abs(average*0.1);
        }
        // Fetch the parameter from the force field
        auto fp  = fs->findParameterType(bi.first, "De");
        int  ntr = ntrain.find(bi.first)->second;
        // Print to the log file    
        if (fplog)
        {
            fprintf(fplog, "%-14s  %6d  %10.1f  %10.1f\n", 
                    bi.first.id().c_str(), ntr, average, error);
        }
        // Now add the new parameter to the force field
        if (fp->mutability() == Mutability::Free ||
            fp->mutability() == Mutability::Bounded)
        {
            fp->setMinimum(average-2*delta);
            fp->setMaximum(average+2*delta);
            fp->setValue(average);
            fp->setUncertainty(error);
            fp->setNtrain(ntr);
        }
        else
        {
            if (fplog)
            {
                fprintf(fplog, "Dissociation energy for %s estimated to be %g, but the parameter is not mutable.\n",
                        bi.first.id().c_str(), average);
            }
        }
    }
    // TODO free the gmx_stats_t
    return rmsd;
}

} // namespace alexandria

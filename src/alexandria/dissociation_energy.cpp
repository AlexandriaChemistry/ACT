/*
 * This source file is part of the Alexandria Chemistry Toolkit.
 *
 * Copyright (C) 2014-2021
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

#include "gromacs/utility/futil.h"

#include "regression.h"
#include "tune_fc_utils.h"
    
/*! \brief Write a csv file containing molecule names and bond energy
 *
 * Writes the whole bond energy matrix.
 */
static void dump_csv(const char                            *csvFile,
                     const std::map<Identifier, int>       &bondIdToIndex,
                     const std::vector<alexandria::MyMol>  &mm,
                     const std::vector<int>                &ntest,
                     const std::vector<double>             &Edissoc,
                     const MatrixWrapper                   &a,
                     const double                           x[])
{
    std::vector<const char*> b2i;
    b2i.resize(bondIdToIndex.size(), nullptr);
    for (auto &j : bondIdToIndex)
    {
        b2i[j.second] = j.first.id().c_str();
    }
    FILE *csv = gmx_ffopen(csvFile, "w");
    fprintf(csv, ",");
    for (auto &j : b2i)
    {
        fprintf(csv, "%s,", j);
    }
    fprintf(csv, "Emol,DeltaHf\n");
    int i = 0;
    for (auto &mymol : mm)
    {
        if (mymol.Hform_ != 0)
        {
            fprintf(csv, "%s,", mymol.getMolname().c_str());
            for (size_t j = 0; (j < bondIdToIndex.size()); j++)
            {
                fprintf(csv, "%g,", a.get(j, i));
            }
            fprintf(csv, "%.3f,%.3f\n", x[i], mymol.Hform_);
            i++;
        }
    }
    fprintf(csv, "Total,");
    for (auto j : ntest)
    {
        fprintf(csv, "%d,", j);
    }
    fprintf(csv, "\n");
    fprintf(csv, "Edissoc,");
    for (auto j : Edissoc)
    {
        fprintf(csv, "%.3f,", j);
    }
    fprintf(csv, "\n");
    fclose(csv);
}

void getDissociationEnergy(FILE                     *fplog,
                           Poldata                  *pd,
                           std::vector<MyMol>       *molset,
                           const char               *csvFile,
                           const std::string        &method,
                           const std::string        &basis)
{
    iqmType iqm = iqmType::Exp;
    std::map<Identifier, int> bondIdToIndex;
    int nColumn = 0;
    int nRow    = 0;
    // Loop over molecules to count number of dissociation energies first
    for (auto mymol = molset->begin(); mymol < molset->end(); ++mymol)
    {
        if (immStatus::OK == mymol->getExpProps(iqm, true, false, true,
                                                method, basis, pd))
        {
            auto myatoms = mymol->atomsConst();
            for (auto &b : mymol->bondsConst())
            {
                auto atypeI = *myatoms.atomtype[b.getAi()-1];
                auto atypeJ = *myatoms.atomtype[b.getAj()-1];
                std::string btypeI, btypeJ;
                if (pd->atypeToBtype(atypeI, &btypeI) &&
                    pd->atypeToBtype(atypeJ, &btypeJ))
                {
                    Identifier bondId({btypeI, btypeJ}, b.getBondOrder(), CanSwap::Yes);
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
            nRow += 1;
        }
    }
    // Now fill the matrix
    std::vector<double> rhs;
    std::vector<int>    ntest;
    int                 row  = 0;

    if ((0 == nColumn) || (0 == nRow))
    {
        gmx_fatal(FARGS, "Number of dissociation energies is %d and number of molecules is %d",
                  nColumn, nRow);
    }

    MatrixWrapper a(nColumn, nRow);
    MatrixWrapper a_copy(nColumn, nRow);
    ntest.resize(nColumn, 0);

    fprintf(fplog, "There are %d different bondtypes to optimize the heat of formation\n", nColumn);
    fprintf(fplog, "There are %d (experimental) reference heat of formation.\n", nRow);

    for (auto mymol = molset->begin(); mymol < molset->end(); ++mymol)
    {
        if (immStatus::OK == mymol->getExpProps(iqm, true, false, true,
                                                method, basis, pd))
        {
            auto myatoms = mymol->atomsConst();
            for (auto &b : mymol->bondsConst())
            {
                const char *atypeI = *myatoms.atomtype[b.getAi()-1];
                const char *atypeJ = *myatoms.atomtype[b.getAj()-1];
                std::string btypeI, btypeJ;
                if (pd->atypeToBtype(atypeI, &btypeI) &&
                    pd->atypeToBtype(atypeJ, &btypeJ))
                {
                    Identifier bondId({btypeI, btypeJ}, b.getBondOrder(), CanSwap::Yes);
                    int column = bondIdToIndex[bondId];
                    
                    a.set(column, row, a.get(column, row) + 1);
                    a_copy.set(column, row, a.get(column, row));
                    ntest[column]++;
                }
            }
            rhs.push_back(-mymol->Emol_);
            row += 1;
        }
    }

    std::string buf = gmx::formatString("Inconsistency in number of energies nRow %d != #rhs %zu", nRow, rhs.size());
    GMX_RELEASE_ASSERT(static_cast<int>(rhs.size()) == nRow, buf.c_str());

    auto nzero = std::count_if(ntest.begin(), ntest.end(), [](const int n)
                               {
                                   return n == 0;
                               });

    GMX_RELEASE_ASSERT(nzero == 0, "Inconsistency in the number of bonds in poldata and ForceConstants_");

    std::vector<double> Edissoc(nColumn);
    a.solve(rhs, &Edissoc);
    if (csvFile)
    {
        dump_csv(csvFile, bondIdToIndex,  *molset, ntest, Edissoc, a_copy, rhs.data());
    }
    for (auto &bi : bondIdToIndex)
    {
        if (fplog)
        {
            fprintf(fplog, "Optimized dissociation energy for %8s with %4d copies to %g\n",
                    bi.first.id().c_str(), ntest[bi.second], Edissoc[bi.second]);
        }
    }

    auto iBonds = InteractionType::BONDS;
    GMX_RELEASE_ASSERT(pd->interactionPresent(iBonds), "No bonds in force field file");
    auto fs  = pd->findForces(InteractionType::BONDS);
    for (auto &b : bondIdToIndex)
    {
        auto fp = fs->findParameterType(b.first, "Dm");
        if (fp->mutability() == Mutability::Free ||
            fp->mutability() == Mutability::Bounded)
        {
            fp->setValue(std::max(100.0, Edissoc[b.second]));
        }
        else
        {
            if (fplog)
            {
                fprintf(fplog, "Dissociation energy for %s estimated to be %g, but the parameter is not mutable.\n",
                        b.first.id().c_str(), Edissoc[b.second]);
            }
        }
    }
}


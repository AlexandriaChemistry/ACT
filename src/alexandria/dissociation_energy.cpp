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
#include "tune_fc_utils.h"
    
/*! \brief Write a csv file containing molecule names and bond energy
 *
 * Writes the whole bond energy matrix.
 */
static void dump_csv(const std::vector<std::string>        &ctest,
                     const std::vector<alexandria::MyMol>  &mm,
                     const std::vector<int>                &ntest,
                     const std::vector<double>             &Edissoc,
                     const MatrixWrapper                   &a,
                     const double                           x[])
{
    FILE *csv = gmx_ffopen("tune_fc.csv", "w");
    fprintf(csv, ",");
    for (auto j : ctest)
    {
        fprintf(csv, "%s,", j.c_str());
    }
    fprintf(csv, "\n");
    int i = 0;
    for (auto &mymol : mm)
    {
        fprintf(csv, "%s,", mymol.getMolname().c_str());
        for (size_t j = 0; (j < ctest.size()); j++)
        {
            fprintf(csv, "%g,", a.get(j, i));
        }
        fprintf(csv, "%.3f\n", x[i]);
        i++;
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
                           const std::vector<MyMol> &molset,
                           int                       nDissociation)
{
    std::vector<double>         rhs;
    std::vector<int>            ntest;
    std::vector<std::string>    ctest;

    int nMol = molset().size();

    if ((0 == nDissociation) || (0 == nMol))
    {
        gmx_fatal(FARGS, "Number of variables is %d and number of molecules is %d",
                  nDissociation, nMol);
    }

    MatrixWrapper a(nDissociation, nMol);
    MatrixWrapper a_copy(nDissociation, nMol);
    ntest.resize(nDissociation, 0);
    ctest.resize(nDissociation);

    fprintf(fplog, "There are %d different bondtypes to optimize the heat of formation\n", nDissociation);
    fprintf(fplog, "There are %d (experimental) reference heat of formation.\n", nMol);

    auto fs  = poldata()->findForces(InteractionType::BONDS);
    auto j   = 0;

    for (auto &mymol :  molset)
    {
        auto myatoms = mymol.atoms();
        for (auto &b : mymol.bondsConst())
        {
            const char *atypeI = *myatoms->atomtype[b.getAi()];
            const char *atypeJ = *myatoms->atomtype[b.getAj()];
            std::string btypeI, btypeJ;
            if (poldata()->atypeToBtype(atypeI, &btypeI) &&
                poldata()->atypeToBtype(atypeJ, &btypeJ))
            {
                Identifier bondId({btypeI, btypeJ}, b.getBondOrder(), CanSwap::Yes);
                auto f   = fs->findParameterTypeConst(bondId, "Dm");
                auto gt  = f.index();
                auto gti = ForceConstants_[InteractionType::BONDS].reverseIndex(gt);
                a.set(gti, j, a.get(gti, j) + 1);
                a_copy.set(gti, j, a.get(gti, j));
                ntest[gti]++;
                if (ctest[gti].empty())
                {
                    ctest[gti].assign(bondId.id());
                }
            }
            else
            {
                gmx_fatal(FARGS, "No parameters for bond in the force field, atoms %s-%s mol %s",
                          atypeI, atypeJ,
                          mymol.getIupac().c_str());
            }
        }
        rhs.push_back(-mymol.Emol_);
    }

    char buf[STRLEN];
    snprintf(buf, sizeof(buf), "Inconsistency in number of energies nMol %d != #rhs %zu", nMol, rhs.size());
    GMX_RELEASE_ASSERT(static_cast<int>(rhs.size()) == nMol, buf);

    auto nzero = std::count_if(ntest.begin(), ntest.end(), [](const int n)
                               {
                                   return n == 0;
                               });

    GMX_RELEASE_ASSERT(nzero == 0, "Inconsistency in the number of bonds in poldata and ForceConstants_");

    std::vector<double> Edissoc(nDissociation);
    a.solve(rhs, &Edissoc);
    if (debug)
    {
        dump_csv(ctest,  molset(), ntest, Edissoc, a_copy, rhs.data());
    }
    for (size_t i = 0; i < ctest.size(); i++)
    {
        if (fplog)
        {
            fprintf(fplog, "Optimized dissociation energy for %8s with %4d copies to %g\n",
                    ctest[i].c_str(), ntest[i], Edissoc[i]);
        }
    }

    int i = 0;
    for (auto &b : ForceConstants_[InteractionType::BONDS].bondNames())
    {
        auto fs = poldata()->findForces(InteractionType::BONDS);
        for(auto &fp : *(fs->findParameters(b.first)))
        {
            if (fp.second.mutability() == Mutability::Free ||
                fp.second.mutability() == Mutability::Bounded)
            {
                if (fp.first == "De")
                {
                    fp.second.setValue(std::max(100.0, Edissoc[i++]));
                }
            }
        }
    }
}


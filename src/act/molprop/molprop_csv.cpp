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


#include "actpre.h"

#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "gromacs/commandline/pargs.h"
#include "gromacs/topology/atomprop.h"
#include "gromacs/topology/topology.h"
#include "gromacs/utility/futil.h"

#include "alexandria/alex_modules.h"
#include "act/molprop/molprop.h"
#include "act/molprop/molprop_util.h"
#include "act/molprop/molprop_xml.h"
#include "act/poldata/poldata_xml.h"

namespace alexandria
{

static void gmx_molprop_csv(const char                 *fn,
                            const std::vector<MolProp> &mp,
                            const char                 *dip_str,
                            const char                 *pol_str,
                            const char                 *ener_str)
{
    MolPropConstIterator  mpi;
    FILE                 *fp;
    double                T = 0.0;
#define NEMP 3
    MolPropObservable mpo[NEMP]   = { MolPropObservable::DIPOLE, MolPropObservable::POLARIZABILITY, MolPropObservable::DHFORM  };
    QmCount           qmc[NEMP];

    find_calculations(mp, mpo[0], dip_str, &qmc[0]);
    find_calculations(mp, mpo[1], pol_str, &qmc[1]);
    find_calculations(mp, mpo[2], ener_str, &qmc[2]);
    for (int k = 0; (k < NEMP); k++)
    {
        printf("--------------------------------------------------\n");
        printf("      Some statistics for %s\n", mpo_name(mpo[k]));
        for (auto q = qmc[k].beginCalc(); q < qmc[k].endCalc(); ++q)
        {
            printf("There are %d calculation results using %s type %s\n",
                   q->count(), q->lot().c_str(), q->type().c_str());
        }
    }
    printf("--------------------------------------------------\n");

    fp = gmx_ffopen(fn, "w");
    fprintf(fp, "\"%s\",\"%s\",\"%s\",\"%s\",\"%s\",\"%s\"",
            "Molecule", "Formula", "InChi", "Charge", "Multiplicity", "Mass");
    for (int k = 0; (k < NEMP); k++)
    {
        for (size_t j = 0; (j < qmc[k].nCalc()+2); j++)
        {
            fprintf(fp, ",\"%s\"", mpo_name(mpo[k]));
        }
    }
    fprintf(fp, "\n");
    for (int ll = 0; (ll < NEMP); ll++)
    {
        fprintf(fp, "\"%s\",\"%s\",\"%s\",\"%s\",\"%s\",\"%s\"", "", "", "", "", "", "");
        for (int k = 0; (k < 3); k++)
        {
            if (ll == 0)
            {
                fprintf(fp, ",\"Experiment\",\"Reference\"");
            }
            else
            {
                fprintf(fp, ",\"\",\"\"");
            }
            for (auto j = qmc[k].beginCalc(); j < qmc[k].endCalc(); ++j)
            {
                switch (ll)
                {
                    case 0:
                        fprintf(fp, ",\"%s\"", j->method().c_str());
                        break;
                    case 1:
                        fprintf(fp, ",\"%s\"", j->basis().c_str());
                        break;
                    case 2:
                        fprintf(fp, ",\"%s\"", j->type().c_str());
                        break;
                    default:
                        fprintf(stderr, "BOE\n");
                        exit(1);
                }
            }
        }
        fprintf(fp, "\n");
    }
    for (mpi = mp.begin(); (mpi < mp.end()); mpi++)
    {
        fprintf(fp, "\"%s\",\"%s\"",
                mpi->getIupac().c_str(),
                mpi->getInchi().c_str());
        for (int k = 0; (k < NEMP); k++)
        {
            std::vector<iqmType> iqms = { iqmType::Exp, iqmType::QM };
            for (auto &iqm : iqms)
            {
                auto gp = mpi->findProperty(mpo[k], iqm, T, "", "", "");
                if (gp)
                {
                    fprintf(fp, ",\"%.4f\",\"FIXME\"", gp->getValue());
                }
                else
                {
                    fprintf(fp, ",\"\",\"\"");
                }
            }
        }
        fprintf(fp, "\n");
    }
    fclose(fp);
}

int mp2csv(int argc, char*argv[])
{
    static const char               *desc[] = {
        "mp2csv converts a molprop database into a spreadsheet"
    };
    t_filenm                         fnm[] =
    {
        { efXML, "-f",  "allmols",  ffREAD },
        { efCSV, "-o",  "csvout",   ffWRITE }
    };
    int                              NFILE   = (sizeof(fnm)/sizeof(fnm[0]));
    static const char               *sort[]  = { nullptr, "molname", "formula", "composition", nullptr };
    static const char               *dip_str = "", *pol_str = "", *ener_str = "";
    static gmx_bool                  bMerge  = FALSE;
    t_pargs                          pa[]    =
    {
        { "-sort",   FALSE, etENUM, {sort},
          "Key to sort the final data file on." },
        { "-merge",  FALSE, etBOOL, {&bMerge},
          "Merge molecule records in the input file" },
        { "-dip_str", FALSE, etSTR, {&dip_str},
          "Selection of the dipole stuff you want in the tables, given as a single string with spaces like: method1/basis1/type1:method2/basis2/type2 (you may have to put quotes around the whole thing in order to prevent the shell from interpreting it)." },
        { "-pol_str", FALSE, etSTR, {&pol_str},
          "Same but for polarizabilities" },
        { "-ener_str", FALSE, etSTR, {&ener_str},
          "Same but for energies" }
    };
    std::vector<alexandria::MolProp> mp;
    gmx_atomprop_t                   ap;
    gmx_output_env_t                *oenv;

    if (!parse_common_args(&argc, argv, PCA_NOEXIT_ON_ARGS, NFILE, fnm,
                           sizeof(pa)/sizeof(pa[0]), pa,
                           sizeof(desc)/sizeof(desc[0]), desc,
                           0, nullptr, &oenv))
    {
        return 0;
    }
    MolPropRead(opt2fn("-f", NFILE, fnm), &mp);
    ap = gmx_atomprop_init();

    alexandria::MolSelect gms;
    MolPropSort(&mp, MPSA_COMPOSITION, ap, gms);

    gmx_molprop_csv(opt2fn("-o", NFILE, fnm), mp,
                    strlen(dip_str) > 0  ? dip_str : nullptr,
                    strlen(pol_str) > 0  ? pol_str : nullptr,
                    strlen(ener_str) > 0 ? ener_str : nullptr);

    return 0;
}

} // namespace alexandria

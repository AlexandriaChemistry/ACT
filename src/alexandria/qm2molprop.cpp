/*
 * This source file is part of the Alexandria Chemistry Toolkit.
 *
 * Copyright (C) 2014-2020 
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

#include <map>

#include "gromacs/commandline/filenm.h"
#include "gromacs/commandline/pargs.h"
#include "gromacs/utility/arraysize.h"
#include "gromacs/utility/exceptions.h"

#include "alex_modules.h"
#include "atype_mapping.h"
#include "babel_io.h"
#include "molprop/molprop.h"
#include "molprop/molprop_util.h"
#include "molprop/molprop_xml.h"
#include "poldata/poldata.h"
#include "poldata/poldata_xml.h"
#include "readpsi4.h"

namespace alexandria
{

int qm2molprop(int argc, char *argv[])
{
    static const char               *desc[] = 
        {
         "qm2molprop reads a series of output files from either",
         "Gaussian ([TT]-g03[tt] option) or Psi4 ([TT]-psi4[tt] option),",
         "collects useful information and saves it to molprop file.[PAR]",
         "The program can optionally map atom type names from an external",
         "source to alexandria types. In the case supply a mapping file",
         "with the [TT]-map[tt] option. The format of that file is:[PAR]",
         "ha h[PAR]",
         "cp c3[PAR]",
         "os o3[PAR]",
         "oh o3[PAR]",
         "etc. where the first atom type maps onto the second atom type.",
         "The first atom type must be unique and basic error checking is",
         "performed."
    };

    t_filenm fnm[] = {
        { efLOG, "-g03",  "gauss",   ffOPTRDMULT },
        { efOUT, "-psi4", "psi4",    ffOPTRDMULT },
        { efSDF, "-sdf",  "mol",     ffOPTRDMULT },
        { efXML, "-d",    "gentop",  ffREAD },
        { efDAT, "-map",  "mapping", ffOPTRD },
        { efXML, "-o",    "molprop", ffWRITE }
    };
#define NFILE sizeof(fnm)/sizeof(fnm[0])

    static int                       maxpot     = 100;
    static int                       nsymm      = 0;
    static char                     *molnm      = nullptr;
    static char                     *iupac      = nullptr;
    static char                     *basis      = nullptr;
    static char                     *jobtype    = (char *)"Opt";
    static char                     *conf       = (char *)"minimum";
    static gmx_bool                  bVerbose   = false;
    static gmx_bool                  compress   = false;

    t_pargs                          pa[]       = {
        { "-v",      FALSE, etBOOL, {&bVerbose},
          "Generate verbose terminal output." },
        { "-compress", FALSE, etBOOL, {&compress},
          "Compress output XML files" },
        { "-molnm", FALSE, etSTR, {&molnm},
          "Name of the molecule in *all* input files. Do not use if you have different molecules in the input files." },
        { "-nsymm", FALSE, etINT, {&nsymm},
          "Symmetry number of the molecule can be supplied here if you know there is an error in the input file" },
        { "-iupac", FALSE, etSTR, {&iupac},
          "IUPAC name of the molecule in *all* input files. Do not use if you have different molecules in the input files." },
        { "-conf",  FALSE, etSTR, {&conf},
          "Conformation of the molecule" },
        { "-basis",  FALSE, etSTR, {&basis},
          "Basis-set used in this calculation for those case where it is difficult to extract from a QM file" },
        { "-jobtype",  FALSE, etSTR, {&jobtype},
          "The job type used in the QM calculation: Opt, Polar, SP, and etc." },
        { "-maxpot", FALSE, etINT, {&maxpot},
          "Maximum percent of the electrostatic potential points that will be added to the molprop file." }
    };
    
    gmx_output_env_t                *oenv;
    alexandria::Poldata              pd;
    std::vector<alexandria::MolProp> mp;

    if (!parse_common_args(&argc, argv, 0, NFILE, fnm, asize(pa), pa, 
                           asize(desc), desc, 0, nullptr, &oenv))
    {
        return 0;
    }

    /* Read force field stuff */
    try
    {
        readPoldata(opt2fn_null("-d", NFILE, fnm), &pd);
    }
    GMX_CATCH_ALL_AND_EXIT_WITH_FATAL_ERROR;

    std::map<std::string, std::string> g2a;
    gaffToAlexandria("", &g2a);
    if (g2a.empty())
    {
        fprintf(stderr, "Don't know how to map GAFF atomtypes to Alexandria\n");
        return 0;
    }
    // Read Gaussian files
    if (opt2bSet("-g03", NFILE, fnm) || opt2bSet("-sdf", NFILE, fnm))
    {
        gmx::ArrayRef<const std::string> gfns;
        if (opt2bSet("-g03", NFILE, fnm))
        {
            gfns = ftp2fns(efLOG, NFILE, fnm);
        }
        else
        {
            gfns = ftp2fns(efSDF, NFILE, fnm);
        }
        int nread = 0;
        for (auto &i : gfns)
        {
            alexandria::MolProp mmm;
            double qtot = 0;
            if (readBabel(i.c_str(), 
                          &mmm, 
                          molnm, 
                          iupac, 
                          conf, 
                          basis,
                          maxpot, 
                          nsymm, 
                          jobtype,
                          &qtot,
                          false))
            {
                nread += 1;
                if (renameAtomTypes(&mmm, g2a))
                {
                    mmm.SetTotalCharge(qtot);
                    mp.push_back(std::move(mmm));
                }
            }
        }
        printf("Read %d molprops from %d Gaussian files.\n", 
               static_cast<int>(mp.size()), static_cast<int>(gfns.size()));
    }
    auto mpsize = mp.size();
    
    // Read Psi4 files
    if (opt2bSet("-psi4", NFILE, fnm))
    {
        gmx::ArrayRef<const std::string> pfns = ftp2fns(efOUT, NFILE, fnm);    
        for (auto &i : pfns)
        {
            alexandria::MolProp mmm;
            if (alexandria::readPsi4(i, &mmm))
            {
                if (SetMolpropAtomTypesAndBonds(&mmm))
                {
                    mp.push_back(std::move(mmm));
                }
            }
        }
        printf("Read %d molprops from %d Psi4 files.\n", 
               static_cast<int>(mp.size()-mpsize), static_cast<int>(pfns.size()));
    }
    alexandria::MolSelect gms;
    MolPropSort(&mp, MPSA_MOLNAME, nullptr, gms);
    MergeDoubleMolprops(&mp, nullptr, TRUE);
    if (mp.size() > 0)
    {
        MolPropWrite(opt2fn("-o", NFILE, fnm), mp, (int)compress);
    }

    return 0;
}

} // namespace alexandria

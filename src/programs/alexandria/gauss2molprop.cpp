/*
 * This source file is part of the Alexandria Chemistry Toolkit.
 *
 * Copyright (C) 2014-2019 
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
 
#include "gmxpre.h"

#include <ctype.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "gromacs/commandline/filenm.h"
#include "gromacs/commandline/pargs.h"
#include "gromacs/fileio/confio.h"
#include "gromacs/math/units.h"
#include "gromacs/math/vec.h"
#include "gromacs/mdtypes/md_enums.h"
#include "gromacs/pbcutil/pbc.h"
#include "gromacs/topology/atomprop.h"
#include "gromacs/utility/arraysize.h"
#include "gromacs/utility/cstringutil.h"
#include "gromacs/utility/exceptions.h"
#include "gromacs/utility/real.h"

#include "alex_modules.h"
#include "babel_io.h"
#include "molprop.h"
#include "molprop_util.h"
#include "molprop_xml.h"
#include "poldata.h"
#include "poldata_xml.h"

int alex_gauss2molprop(int argc, char *argv[])
{
    static const char               *desc[] = {
        "gauss2molprop reads a series of Gaussian output files, and collects",
        "useful information, and saves it to molprop file."
    };

    t_filenm                         fnm[] = {
        { efLOG, "-g03",  "gauss",   ffRDMULT },
        { efDAT, "-d",    "gentop",  ffREAD },
        { efDAT, "-o",    "molprop", ffWRITE }
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
          "Basis-set used in this calculation for those case where it is difficult to extract from a Gaussian file" },
        { "-jobtype",  FALSE, etSTR, {&jobtype},
          "The job type used in the Gaussian calculation: Opt, Polar, SP, and etc." },
        { "-maxpot", FALSE, etINT, {&maxpot},
          "Maximum percent of the electrostatic potential points that will be added to the molprop file." }
    };
    
    gmx_output_env_t                *oenv;
    gmx_atomprop_t                   aps;
    alexandria::Poldata              pd;
    std::vector<alexandria::MolProp> mp;

    if (!parse_common_args(&argc, argv, 0, NFILE, fnm, asize(pa), pa, 
                           asize(desc), desc, 0, nullptr, &oenv))
    {
        return 0;
    }

    /* Read standard atom properties */
    aps = gmx_atomprop_init();

    /* Read force field stuff */
    try
    {
        readPoldata(opt2fn_null("-d", NFILE, fnm), pd, aps);
    }
    GMX_CATCH_ALL_AND_EXIT_WITH_FATAL_ERROR;

    gmx::ArrayRef<const std::string> fns = ftp2fns(efLOG, NFILE, fnm);    
    for (auto &i : fns)
    {
        alexandria::MolProp mmm;
        readBabel(i.c_str(), 
                  mmm, 
                  molnm, 
                  iupac, 
                  conf, 
                  basis,
                  maxpot, 
                  nsymm, 
                  jobtype,
                  0.0);
        mp.push_back(std::move(mmm));
    }

    printf("Succesfully read %d molprops from %d Gaussian files.\n", 
           static_cast<int>(mp.size()), static_cast<int>(fns.size()));
    alexandria::MolSelect gms;
    MolPropSort(mp, MPSA_MOLNAME, nullptr, gms);
    merge_doubles(mp, nullptr, TRUE);
    if (mp.size() > 0)
    {
        MolPropWrite(opt2fn("-o", NFILE, fnm), mp, (int)compress);
    }

    return 0;
}

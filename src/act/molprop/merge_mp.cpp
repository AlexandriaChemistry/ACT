/*
 * This source file is part of the Alexandria Chemistry Toolkit.
 *
 * Copyright (C) 2014-2020 
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
 */
 
#include "actpre.h"

#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "gromacs/commandline/pargs.h"
#include "gromacs/math/units.h"
#include "gromacs/math/vec.h"
#include "gromacs/topology/topology.h"
#include "gromacs/utility/arraysize.h"
#include "gromacs/utility/cstringutil.h"
#include "gromacs/utility/exceptions.h"
#include "gromacs/utility/fatalerror.h"
#include "gromacs/utility/futil.h"
#include "gromacs/utility/smalloc.h"

#include "alexandria/alex_modules.h"
#include "act/molprop/molprop.h"
#include "act/molprop/molprop_sqlite3.h"
#include "act/molprop/molprop_util.h"
#include "act/molprop/molprop_xml.h"
#include "act/poldata/poldata.h"
#include "act/poldata/poldata_xml.h"
#include "act/utility/stringutil.h"

namespace alexandria
{

int merge_mp(int argc, char *argv[])
{
    static const char               *desc[] =
    {
        "merge_mp reads multiple molprop files and merges the molecule descriptions",
        "into a single new file. By specifying the [TT]-db[TT] option additional experimental",
        "information will be read from a SQLite3 database.[PAR]",
    };
    t_filenm                         fnm[] =
    {
        { efXML, "-mp",  "data",      ffOPTRDMULT },
        { efXML, "-o",   "allmols",   ffWRITE },
        { efXML, "-ff",  "gentop",    ffOPTRD },
        { efDAT, "-db",  "sqlite",    ffOPTRD }
    };
    int     NFILE       = asize(fnm);
    int     compress    = 1;
    real    temperature = 298.15;
    bool    forceMerge  = false;
    int     maxwarn     = 0;
    t_pargs pa[]        =
    {
        { "-compress", FALSE, etBOOL, {&compress},
          "Compress output XML files" },
        { "-maxwarn", FALSE, etINT, {&maxwarn},
          "Will only write output if number of warnings is at most this." },
        { "-temperature", FALSE, etREAL, {&temperature},
          "Temperature for properties to extract from the SQLite database" },
        { "-force", FALSE, etBOOL, {&forceMerge},
          "Force merging compounds with the same name even if not the formula matches" }
    };
    std::vector<MolProp>  mp;
    Poldata               pd;
    gmx_output_env_t     *oenv;

    if (!parse_common_args(&argc, argv, PCA_NOEXIT_ON_ARGS, NFILE, fnm,
                           asize(pa), pa, asize(desc), desc, 0, nullptr, &oenv))
    {
        return 0;
    }

    try 
    {
        alexandria::readPoldata(opt2fn("-ff", NFILE, fnm), &pd);
    }
    GMX_CATCH_ALL_AND_EXIT_WITH_FATAL_ERROR;

    auto fns = opt2fns("-mp", NFILE, fnm);
    int nwarn = merge_xml(fns, &mp, nullptr, nullptr, nullptr, forceMerge);

    if (nwarn <= maxwarn)
    {
        ReadSqlite3(opt2fn_null("-db", NFILE, fnm), &mp, temperature);

        MolPropWrite(opt2fn("-o", NFILE, fnm), mp, compress);
    }
    else
    {
        printf("Too many warnings (%d), not generating output\n", nwarn);
    }
    return 0;
}

} // namespace alexandria

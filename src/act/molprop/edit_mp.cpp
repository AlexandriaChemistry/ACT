/*
 * This source file is part of the Alexandria Chemistry Toolkit.
 *
 * Copyright (C) 2014-2023
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

#include "act/alexandria/alex_modules.h"
#include "act/molprop/molprop.h"
#include "act/molprop/molprop_sqlite3.h"
#include "act/molprop/molprop_util.h"
#include "act/molprop/molprop_xml.h"
#include "act/poldata/poldata.h"
#include "act/poldata/poldata_xml.h"
#include "act/utility/stringutil.h"
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

namespace alexandria
{

int edit_mp(int argc, char *argv[])
{
    static const char               *desc[] =
    {
        "edit_mp manipulates molprop files. It can read multiple molprop files and merges",
        "the molecule descriptions into a single new file. By specifying the [TT]-db[TT]", 
        "option additional experimental information will be read from a SQLite3 database.[PAR]",
        "The program can run in paralllel to read a file on one processer, then send it over",
        "an MPI connection to one or more other processors to write. In this manner the MPI transfer",
        "software in ACT can be tested."
    };
    t_filenm                         fnm[] =
    {
        { efXML, "-mp",  "data",    ffOPTRDMULT },
        { efXML, "-o",   "allmols", ffWRITE     },
        { efDAT, "-db",  "sqlite",  ffOPTRD     }
    };
    int      NFILE       = asize(fnm);
    int      compress    = 1;
    real     temperature = 298.15;
    bool     forceMerge  = false;
    gmx_bool bcast       = false;
    int      maxwarn     = 0;
    int      writeNode   = 0;
    t_pargs pa[] =
    {
        { "-compress", FALSE, etBOOL, {&compress},
          "Compress output XML files" },
        { "-maxwarn", FALSE, etINT, {&maxwarn},
          "Will only write output if number of warnings is at most this." },
        { "-temperature", FALSE, etREAL, {&temperature},
          "Temperature for properties to extract from the SQLite database" },
        { "-force", FALSE, etBOOL, {&forceMerge},
          "Force merging compounds with the same name even if not the formula matches" },
        { "-bcast", FALSE, etBOOL, {&bcast},
          "Use broadcast instead of send/receive when running in parallel" },
        { "-wn", FALSE, etINT, {&writeNode},
          "Processor ID to write from if in parallel." }
    };
    std::vector<MolProp>  mpt;
    Poldata               pd;
    gmx_output_env_t     *oenv;

    if (!parse_common_args(&argc, argv, PCA_NOEXIT_ON_ARGS, NFILE, fnm,
                           asize(pa), pa, asize(desc), desc, 0, nullptr, &oenv))
    {
        return 0;
    }

    auto fns = opt2fns("-mp", NFILE, fnm);

    CommunicationRecord cr;
    cr.init(cr.size());
    auto comm = MPI_COMM_WORLD;
    int root  = 0;
    if (cr.isMaster())
    {
        int nwarn   = merge_xml(fns, &mpt, nullptr, nullptr, nullptr, forceMerge);
        int mptsize = mpt.size();
        if (nwarn <= maxwarn)
        {
            ReadSqlite3(opt2fn_null("-db", NFILE, fnm), &mpt, temperature);

            printf("Read %d molecules\n", mptsize);
        }
        else
        {
            printf("Too many warnings (%d), not generating output\n", nwarn);
            mpt.clear();
        }
            
        if (bcast)
        {
            cr.bcast(&mptsize, comm);
            for(auto &mm : mpt)
            {
                mm.BroadCast(&cr, root, comm);
            }
        }
        else
        {
            for(int dest = 1; dest < cr.size(); dest++)
            {
                cr.send_int(dest, mptsize);
                for(const auto &mm : mpt)
                {
                    mm.Send(&cr, dest);
                }
            }
        }
    }
    else
    {
        if (bcast)
        {
            int nmpt;
            cr.bcast(&nmpt, comm);
            for(int i = 0; i < nmpt; i++)
            {
                MolProp mp;
                mp.BroadCast(&cr, root, comm);
                mpt.push_back(mp);
            }
        }
        else
        {
            int nmpt = cr.recv_int(cr.superior());
            for(int i = 0; i < nmpt; i++)
            {
                MolProp mp;
                mp.Receive(&cr, cr.superior());
                mpt.push_back(mp);
            }
        }
    }
    if (writeNode == cr.rank())
    {
        MolPropWrite(opt2fn("-o", NFILE, fnm), mpt, compress);
    }
    return 0;
}

} // namespace alexandria

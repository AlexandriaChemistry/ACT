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

#include <stdlib.h>

#include "gromacs/commandline/pargs.h"
#include "gromacs/fileio/oenv.h"
#include "gromacs/utility/arraysize.h"

#include "alex_modules.h"
#include "poldata.h"
#include "poldata_xml.h"

int alex_poldata_test(int argc, char*argv[])
{
    static const char               *desc[] = {
        "poldata_test reads a poldata (force field) file and writes a new one.",
    };
    gmx_output_env_t                *oenv;
    t_filenm                         fnm[] = {
        { efDAT, "-f", "pdin", ffREAD },
        { efDAT, "-o", "pdout", ffWRITE }
    };
    
    int                 NFILE = asize(fnm);
    alexandria::Poldata pd;
    
    auto cr = init_commrec();
    if (MASTER(cr))
    {
        if (!parse_common_args(&argc, argv, 0, NFILE, fnm, 0, nullptr,
                               1, desc, 0, nullptr, &oenv))
        {
            return 0;
        }
        
        try 
        {
            alexandria::readPoldata(opt2fn("-f", NFILE, fnm), &pd);
        }
        GMX_CATCH_ALL_AND_EXIT_WITH_FATAL_ERROR;
    }
    if (PAR(cr))
    {
        if (MASTER(cr))
        {
            pd.Send(cr, 1);
            std::string outfile(opt2fn("-o", NFILE, fnm));
            gmx_send_str(cr, 1, &outfile);
        }
        else if (cr->nodeid == 1)
        {
            pd.Receive(cr, 0);
            std::string outfile;
            gmx_recv_str(cr, 0, &outfile);
            alexandria::writePoldata(outfile, &pd, 0);
        }
    }
    else
    {
        alexandria::writePoldata(opt2fn("-o", NFILE, fnm), &pd, 0);
    }
    return 0;
}

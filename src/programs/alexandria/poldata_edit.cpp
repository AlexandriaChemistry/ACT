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
#include "interactiontype.h"
#include "poldata.h"
#include "poldata_xml.h"

static void modifyPoldata(alexandria::Poldata *pd,
                          const char *ptype,
                          const char *particle,
                          bool bSetMin, real pmin,
                          bool bSetMax, real pmax,
                          bool force)
{
    if (!(bSetMin || bSetMax))
    {
        return;
    }
    auto itype = pd->typeToInteractionType(ptype);
    printf("Will change parameter type %s in %s\n", ptype,
           interactionTypeToString(itype).c_str());
    printf("for particle type %s.\n", strlen(particle) > 0 ? particle : "all");
    if (bSetMin)
    {
        printf("Minimum will be set to %g\n", pmin);
    }
    if (bSetMax)
    {
        printf("Maximum will be set to %g\n", pmax);
    }
    auto fs = pd->findForces(itype)->parameters();
    for(auto &ffs : *fs)
    {
        if (strlen(particle) == 0 || particle == ffs.first.id())
        {
            for(auto &pp : ffs.second)
            {
                if (pp.first == ptype)
                {
                    if (force || pp.second.isMutable())
                    {
                        if (bSetMin)
                        {
                            pp.second.setMinimum(pmin);
                            pp.second.setValue(std::max(pmin, pp.second.value()));
                            pp.second.setMaximum(std::max(pmin, pp.second.maximum()));
                        }
                        if (bSetMax)
                        {
                            pp.second.setMaximum(pmax);
                            pp.second.setValue(std::min(pmax, pp.second.value()));
                            pp.second.setMinimum(std::min(pmax, pp.second.minimum()));
                        }
                    }
                }
            }
        }
    }
}

int alex_poldata_edit(int argc, char*argv[])
{
    static const char               *desc[] = {
        "poldata_edit reads a poldata (force field) file and writes a new one.[PAR]",
        "If run in parallel, this utility will read on the first processor,",
        "then process the data and send it to the second processor. That one",
        "will then write the new file. This is for debugging poldata",
        "communication.[PAR]",
        "Before saving the file, certain changes can be made. If the minimum",
        "or the maximum is to be changed, the actual value may be set to the",
        "new minmum or maximum if it falls outside the new bounds."
    };
    gmx_output_env_t                *oenv;
    t_filenm                         fnm[] = {
        { efDAT, "-f", "pdin", ffREAD },
        { efDAT, "-o", "pdout", ffWRITE }
    };
    static char *parameter = (char *)"";
    static char *particle  = (char *)"";
    real         pmin      = 0;
    real         pmax      = 0;
    gmx_bool     force     = false;
    t_pargs                          pa[]     = 
    {
        { "-p",      FALSE, etSTR,  {&parameter},
          "Type of parameter to change, e.g. zeta." },
        { "-a",      FALSE, etSTR,  {&particle},
          "Particle type to change, if not set all particles of the correct type will be changed." },
        { "-min",    FALSE, etREAL, {&pmin},
          "Minimum value of parameter." },
        { "-max",    FALSE, etREAL,  {&pmax},
          "Maximum value of parameter." },
        { "-force",  FALSE, etBOOL, {&force},
          "Will change also non-mutable parameters. Use with care!" }
    };
    int                 npargs = asize(pa);
    int                 NFILE  = asize(fnm);
    alexandria::Poldata pd;
    
    auto cr = init_commrec();
    if (MASTER(cr))
    {
        if (!parse_common_args(&argc, argv, 0, NFILE, fnm, npargs, pa,
                               1, desc, 0, nullptr, &oenv))
        {
            return 0;
        }
        
        try 
        {
            alexandria::readPoldata(opt2fn("-f", NFILE, fnm), &pd);
        }
        GMX_CATCH_ALL_AND_EXIT_WITH_FATAL_ERROR;
        
        modifyPoldata(&pd, parameter, particle,
                      opt2parg_bSet("-min", npargs, pa), pmin,
                      opt2parg_bSet("-max", npargs, pa), pmax,
                      force);
                      
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

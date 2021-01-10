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
#include "gmxpre.h"

#include <map>

#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "gromacs/commandline/pargs.h"
#include "gromacs/math/units.h"
#include "gromacs/math/vec.h"
#include "gromacs/statistics/statistics.h"
#include "gromacs/topology/topology.h"
#include "gromacs/utility/arraysize.h"
#include "gromacs/utility/cstringutil.h"
#include "gromacs/utility/exceptions.h"
#include "gromacs/utility/fatalerror.h"
#include "gromacs/utility/futil.h"
#include "gromacs/utility/smalloc.h"

#include "alex_modules.h"
#include "identifier.h"
#include "poldata.h"
#include "poldata_tables.h"
#include "poldata_xml.h"

static void merge_parameter(const std::vector<alexandria::Poldata> &pds,
                            alexandria::InteractionType             iType,
                            const std::string                      &parameter,
                            alexandria::Poldata                    *pdout)
{  
    std::vector<gmx_stats_t> lsq;
    std::vector<int>         ntrain;
    bool                     first = true;
    for (const auto &pd : pds)
    {
        int  j  = 0;
        auto fs = pd.findForcesConst(iType);
        for (const auto &fp : fs.parametersConst())
        {
            for (const auto &pp : fp.second)
            {
                if (parameter == pp.first)
                {
                    if (first)
                    {
                        auto newlsq = gmx_stats_init();
                        lsq.push_back(newlsq);
                        ntrain.push_back(0);
                    }
                    gmx_stats_add_point(lsq[j], 0, pp.second.value(), 0, 0);
                    ntrain[j] += pp.second.ntrain();
                    j++;
                }
            }
        }
        first = false;
    }
    auto forces = pdout->findForces(iType);
    auto param  = forces->parameters();
    int j = 0;
    for (auto &fp : *param)
    {
        for (auto &pp : fp.second)
        {
            if (parameter == pp.first)
            {
                real average = 0;
                real sigma   = 0;
                int  N       = 0;
                if ((estatsOK == gmx_stats_get_average(lsq[j], &average)) &&
                    (estatsOK == gmx_stats_get_sigma(lsq[j], &sigma)) &&
                    (estatsOK == gmx_stats_get_npoints(lsq[j], &N)))
                {
                    pp.second.setValue(average);
                    pp.second.setUncertainty(sigma);
                    pp.second.setNtrain(ntrain[j]/pds.size());
                }
                j++;
            }
        }
    }
    
    for (auto lll = lsq.begin(); lll < lsq.end(); lll++)
    {
         gmx_stats_free(*lll);
    }
}

int alex_merge_pd(int argc, char *argv[])
{
    const char *desc[] =
    {
        "merge_pd reads multiple Alexandria force field files and merges them",
        "into a single new force field file containing average values",
        "and standard deviation.",
    };    
    t_filenm    fnm[] =
    {
        { efDAT, "-di",    "pdin",  ffRDMULT},
        { efDAT, "-do",    "pdout", ffWRITE },
        { efTEX, "-latex", "pdout", ffWRITE }
    };
    int         NFILE       = asize(fnm);;
    
    gmx_bool   bcompress   = false;    
    //! String for command line to harvest the options to fit
    char       *mergeString = nullptr;

    t_pargs     pa[]        =
    {
        { "-compress", FALSE, etBOOL, {&bcompress},
          "Compress output XML files" },
        { "-merge", FALSE, etSTR, {&mergeString},
          "Quoted list of parameters to merge,  e.g. 'alpha zeta'. An empty string means all parameters will be merged." }
    };
    std::vector<alexandria::Poldata> pds;
    alexandria::Poldata              pdout;
    gmx_output_env_t                *oenv;

    if (!parse_common_args(&argc, 
                           argv, 
                           PCA_NOEXIT_ON_ARGS, 
                           NFILE, 
                           fnm,
                           asize(pa), 
                           pa, 
                           asize(desc), 
                           desc, 
                           0, 
                           nullptr, 
                           &oenv))
    {
        return 0;
    }
    std::string allParams("alpha chi jaa zeta hardness electronegativity");
    if (nullptr == mergeString)
    {
        mergeString = strdup(allParams.c_str());
    }
    /*
      Read all the gentop files.
     */
    auto filenames = opt2fns("-di", NFILE, fnm);
    if (filenames.size() < 2)
    {
        gmx_fatal(FARGS, "At least two gentop files are needed for merging!\n");
    }
    
    for (auto &i : filenames)
    {
        alexandria::Poldata pd;
        readPoldata(i.c_str(), &pd);
        pds.push_back(std::move(pd));
    }    

    // Copy the first gentop file into pdout->
    readPoldata(filenames[0].c_str(), &pdout);
    
    // We now update different parts of pdout-> 
    for(const auto &type : gmx::splitString(mergeString))
    {
        alexandria::InteractionType itype;
        if (pdout.typeToInteractionType(type, &itype))
        {
            merge_parameter(pds, itype, type, &pdout);
        }
        else
        {
            fprintf(stderr, "No such parameter type '%s' in %s\n",
                    type.c_str(), filenames[0].c_str());
        }
    }

    writePoldata(opt2fn("-do", NFILE, fnm), &pdout, bcompress);
    if (opt2bSet("-latex", NFILE, fnm))
    {
        FILE *tp = gmx_ffopen(opt2fn("-latex", NFILE, fnm), "w");
        alexandria_subtype_table(tp, &pdout);
        alexandria_eemprops_table(tp, &pdout);
        gmx_ffclose(tp);
    }           
    
    return 0;
}


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

enum EemAtomProps {
    eEMEta   = 0, 
    eEMChi   = 1, 
    eEMZeta  = 2,  
    eEMAlpha = 3,
    eEMAll   = 4,
    eEMNR    = 5
};

typedef struct {
    EemAtomProps  eEM;
    const char   *name;
} t_eemAtom_props;

t_eemAtom_props eemAtom_props[eEMNR] = {
    {eEMEta,   "eta"},
    {eEMChi,   "chi"},
    {eEMZeta,  "zeta"},
    {eEMAlpha, "alpha"},
    {eEMAll,   "all"}
};

static EemAtomProps name2eemprop(const std::string name)
{
    for (auto i = 0; i < eEMNR; i++)
    {
        if (strcasecmp(name.c_str(), eemAtom_props[i].name) == 0)
        {
            return eemAtom_props[i].eEM;
        }
    }
    return eEMNR;
}

static void merge_parameter(const std::vector<alexandria::Poldata> &pds,
                            alexandria::InteractionType             iType,
                            const std::string                      &parameter,
                            alexandria::Poldata                    *pdout)
{  
    auto nAtypes = pdout->getNatypes();
    gmx_stats_t lsq[nAtypes];
    
    for (size_t i = 0; i < nAtypes; i++)
    {
        lsq[i] =  gmx_stats_init();
    }    
    int j = 0;
    for (const auto &atp : pdout->particleTypesConst())
    {
        for (const auto& pd : pds)
        {
            auto fs    = pd.findForcesConst(iType);
            auto ztype = atp.interactionTypeToIdentifier(iType);
            auto ei    = fs.findParametersConst(ztype);
            gmx_stats_add_point(lsq[j], 0, ei[parameter].value(), 0, 0);
        }
        j++;
    }
    
    j = 0;
    for (const auto &atp : pdout->particleTypesConst())
    {
        auto fs      = pdout->findForces(iType);
        auto ztype   = atp.interactionTypeToIdentifier(iType);
        auto ei      = fs->findParameters(ztype);
        real average = 0;
        real sigma   = 0;
        int  N       = 0;
        if ((estatsOK == gmx_stats_get_average(lsq[j], &average)) &&
            (estatsOK == gmx_stats_get_sigma(lsq[j], &sigma)) &&
            (estatsOK == gmx_stats_get_npoints(lsq[j], &N)))
        {
            (*ei)[parameter].setValue(average);
            (*ei)[parameter].setUncertainty(sigma);
            (*ei)[parameter].setNtrain(N);
        }
        j++;
    }
    for (size_t i = 0; i < nAtypes; i++)
    {
         gmx_stats_free(lsq[i]);
    }
}

int alex_merge_pd(int argc, char *argv[])
{
    static const char               *desc[] =
    {
        "merge_pd reads multiple gentop files and merges them",
        "into a single new gentop file.",
    };    
    t_filenm                         fnm[] =
    {
        { efDAT, "-di",    "pdin",  ffRDMULT},
        { efDAT, "-do",    "pdout", ffWRITE },
        { efTEX, "-latex", "pdout", ffWRITE }
    };
    int                              NFILE       = asize(fnm);;
    
    static gmx_bool                  bcompress   = false;    
    static const char               *eemprop[]   = {nullptr, "eta", "chi", "zeta", "alpha", "all", nullptr};
    
    t_pargs                          pa[]        =
    {
        { "-compress", FALSE, etBOOL, {&bcompress},
          "Compress output XML files" },
        { "-eemprop", FALSE, etENUM, {eemprop},
          "Atomic property used to describe molecular eletric properties." }
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
    
    //We now update different parts of pdout->    
    EemAtomProps eem = name2eemprop(eemprop[0]);        
    if (eem == eEMEta || eem == eEMAll)
    {
        merge_parameter(pds, alexandria::InteractionType::ELECTRONEGATIVITYEQUALIZATION,
                        "jaa", &pdout);
    }
    if (eem == eEMChi || eem == eEMAll)
    {
        merge_parameter(pds, alexandria::InteractionType::ELECTRONEGATIVITYEQUALIZATION,
                        "chi", &pdout);
    }
    if (eem == eEMAlpha || eem == eEMAll)
    {
        merge_parameter(pds, alexandria::InteractionType::POLARIZATION,
                        "alpha", &pdout);
    }
    if (eem == eEMZeta || eem == eEMAll)
    {
        merge_parameter(pds, alexandria::InteractionType::CHARGEDISTRIBUTION,
                        "zeta", &pdout);
    }
    if (eem == eEMNR)
    {
        gmx_fatal(FARGS, "There is no atomic electric property called %s in alexandria.\n", eemprop[0]);
    }    
    writePoldata(opt2fn("-do", NFILE, fnm), &pdout, bcompress);
    if (opt2bSet("-latex", NFILE, fnm))
    {
        FILE        *tp;
        tp = gmx_ffopen(opt2fn("-latex", NFILE, fnm), "w");
        alexandria_eemprops_table(tp, &pdout);
        gmx_ffclose(tp);
    }           
    
    return 0;
}


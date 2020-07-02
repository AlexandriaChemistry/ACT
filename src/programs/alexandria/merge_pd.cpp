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

static void merge_Chi(std::vector<alexandria::Poldata>     pds,
                      alexandria::Poldata                 &pdout)
{  
    auto nAtypes = pdout.getNatypes();
    gmx_stats_t lsq[nAtypes];
    
    for (size_t i = 0; i < nAtypes; i++)
    {
        lsq[i] =  gmx_stats_init();
    }    
    int j = 0;
    for (auto atp = pdout.getAtypeBegin(); 
         atp < pdout.getAtypeEnd(); atp++, j++)
    {
        for (const auto& pd : pds)
        {
            auto ei = pd.atype2Eem(atp->getType());
            if (ei != pd.EndEemprops())
            {
                gmx_stats_add_point(lsq[j], 0, ei->getChi0(), 0, 0);

            }
            else
            {
                gmx_fatal(FARGS, "%s atomtype does not exists in Eeemprops.\n", 
                          atp->getType().c_str());
            }
        }
    }
    
    j = 0;
    for (auto atp = pdout.getAtypeBegin(); 
         atp < pdout.getAtypeEnd(); atp++, j++)
    {
        auto ei = pdout.atype2Eem(atp->getType());
        if (ei != pdout.EndEemprops())
        {
            real average = 0;
            real sigma   = 0;
            if ((estatsOK == gmx_stats_get_average(lsq[j], &average))&&
                (estatsOK == gmx_stats_get_sigma(lsq[j], &sigma)))
            {
                ei->setChi0(average);
                ei->setChi0_sigma(sigma);
                
            }
            else
            {
                gmx_fatal(FARGS, "estats is not OK for %s for Chi.\n", 
                          atp->getType().c_str());
            }
        }
        else
        {
            gmx_fatal(FARGS, "%s atomtype does not exists in Eeemprops.\n", 
                      atp->getType().c_str());
        }
    }    
    for (size_t i = 0; i < nAtypes; i++)
    {
         gmx_stats_free(lsq[i]);
    }
}

static void merge_Eta(std::vector<alexandria::Poldata>     pds,
                      alexandria::Poldata                 &pdout)
{
    auto nAtypes = pdout.getNatypes();
    gmx_stats_t lsq[nAtypes];
    
    for (size_t i = 0; i < nAtypes; i++)
    {
        lsq[i] =  gmx_stats_init();
    }    
    int j = 0;
    for (auto atp = pdout.getAtypeBegin(); 
         atp < pdout.getAtypeEnd(); atp++, j++)
    {
        for (const auto& pd : pds)
        {
            auto ei = pd.atype2Eem(atp->getType());
            if (ei != pd.EndEemprops())
            {
                gmx_stats_add_point(lsq[j], 0, ei->getJ0(), 0, 0);

            }
            else
            {
                gmx_fatal(FARGS, "%s atomtype does not exists in Eeemprops.\n", 
                          atp->getType().c_str());
            }
        }
    }
    
    j = 0;
    for (auto atp = pdout.getAtypeBegin(); 
         atp < pdout.getAtypeEnd(); atp++, j++)
    {
        auto ei = pdout.atype2Eem(atp->getType());
        if (ei != pdout.EndEemprops())
        {
            real average = 0;
            real sigma   = 0;
            if ((estatsOK == gmx_stats_get_average(lsq[j], &average))&&
                (estatsOK == gmx_stats_get_sigma(lsq[j], &sigma)))
            {
                ei->setJ0(average);
                ei->setJ0_sigma(sigma);
                
            }
            else
            {
                gmx_fatal(FARGS, "estats is not OK for %s for Eta.\n", 
                          atp->getType().c_str());
            }
        }
        else
        {
            gmx_fatal(FARGS, "%s atomtype does not exists in Eeemprops.\n", 
                      atp->getType().c_str());
        }
    }    
    for (size_t i = 0; i < nAtypes; i++)
    {
         gmx_stats_free(lsq[i]);
    }
}

static void merge_Alpha(std::vector<alexandria::Poldata>     pds,
                        alexandria::Poldata                 &pdout)
{
    
    auto nAtypes = pdout.getNatypes();
    gmx_stats_t lsq[nAtypes];
    
    for (size_t i = 0; i < nAtypes; i++)
    {
        lsq[i] =  gmx_stats_init();
    }    
    int j = 0;
    for (auto pt = pdout.getPtypeBegin(); 
         pt < pdout.getPtypeEnd(); pt++, j++)
    {
        for (const auto& pd : pds)
        {
            double alpha = 0;
            double sigma = 0;
            if (pd.getPtypePol(pt->getType(), &alpha, &sigma))
            {                
                gmx_stats_add_point(lsq[j], 0, alpha, 0, 0);

            }
            else
            {
                gmx_fatal(FARGS, "%s atomtype does not exists in Eeemprops.\n", 
                          pt->getType().c_str());
            }
        }
    }
    
    j = 0;
    for (auto pt = pdout.getPtypeBegin(); 
         pt < pdout.getPtypeEnd(); pt++, j++)
    {
        std::string  ptype;           
        if (pt != pdout.getPtypeEnd())
        {
            real average = 0;
            real sigma   = 0;
            if ((estatsOK == gmx_stats_get_average(lsq[j], &average))&&
                (estatsOK == gmx_stats_get_sigma(lsq[j], &sigma)))
            {
                pdout.setPtypePolarizability(pt->getType(), average, sigma);                
            }
            else
            {
                gmx_fatal(FARGS, "estats is not OK for %s for Eta.\n", 
                          pt->getType().c_str());
            }
        }
        else
        {
            gmx_fatal(FARGS, "%s atomtype does not exists in Eeemprops.\n", 
                      pt->getType().c_str());
        }
    }    
    for (size_t i = 0; i < nAtypes; i++)
    {
         gmx_stats_free(lsq[i]);
    }
}

static void merge_Zeta(std::vector<alexandria::Poldata>     pds,
                       alexandria::Poldata                 &pdout)
{   
    auto nAtypes = pdout.getNatypes();   
    
    gmx_stats_t core[nAtypes];
    gmx_stats_t shell[nAtypes];
       
    for (size_t i = 0; i < nAtypes; i++)
    {
        core[i]  =  gmx_stats_init();
        shell[i] =  gmx_stats_init();
    }     
    
    // Reading Zeta from all pds   
    int j = 0;
    for (auto atp = pdout.getAtypeBegin(); 
         atp < pdout.getAtypeEnd(); atp++, j++)
    {
        for (const auto& pd : pds)
        {
            auto ei = pd.ztype2Eem(atp->getZtype());
            if (ei != pd.EndEemprops())
            {
                auto nzeta = ei->getNzeta();
                if (nzeta == 1)
                {
                    gmx_stats_add_point(core[j], 0, ei->getZeta(0), 0, 0);
                }
                else if (nzeta == 2)
                {
                    gmx_stats_add_point(core[j],  0, ei->getZeta(0), 0, 0);
                    gmx_stats_add_point(shell[j], 0, ei->getZeta(1), 0, 0);
                }
                else
                {
                    gmx_fatal(FARGS, "The number of Zeta is wrong for %s atomtype.\n", 
                              ei->getName());
                }
                
            }   
        }
    }
    
    // Merging Zeta
    j = 0;
    for (auto atp = pdout.getAtypeBegin(); 
         atp < pdout.getAtypeEnd(); atp++, j++)
    {
        auto ei = pdout.ztype2Eem(atp->getZtype());
        if (ei != pdout.EndEemprops())
        {
            std::string zstr, z_sig;
            auto nzeta  = ei->getNzeta();            
            if (nzeta == 1)
            {
                double core_ave  = 0;
                double core_sig  = 0;
                if ((estatsOK == gmx_stats_get_average(core[j], &core_ave))&&
                    (estatsOK == gmx_stats_get_sigma(core[j],   &core_sig)))
                {
                    zstr.append(gmx::formatString("%f ", core_ave));
                    z_sig.append(gmx::formatString("%f ", core_sig));
                    ei->setZetastr(zstr);
                    ei->setZeta_sigma(z_sig);
                }
                else
                {
                    gmx_fatal(FARGS, "estats is not OK for %s.\n", 
                              ei->getName());
                }
            }
            else if (nzeta == 2)
            {
                double core_ave  = 0;
                double core_sig  = 0;
                double shell_ave = 0;
                double shell_sig = 0;
                if ((estatsOK == gmx_stats_get_average(core[j],  &core_ave))  &&
                    (estatsOK == gmx_stats_get_sigma(core[j],    &core_sig))  &&
                    (estatsOK == gmx_stats_get_average(shell[j], &shell_ave)) &&
                    (estatsOK == gmx_stats_get_sigma(shell[j],   &shell_sig)))
                {
                    zstr.append(gmx::formatString("%f ", core_ave));
                    zstr.append(gmx::formatString("%f ", shell_ave));
                    z_sig.append(gmx::formatString("%f ", core_sig));
                    z_sig.append(gmx::formatString("%f ", shell_sig));
                        
                    ei->setZetastr(zstr);
                    ei->setZeta_sigma(z_sig);
                }
                else
                {
                    gmx_fatal(FARGS, "estats is not OK for %s.\n", 
                              ei->getName());
                }
            }
            else
            {
                gmx_fatal(FARGS, "The number of Zeta is wrong for %s atomtype.\n", 
                          ei->getName());
            }
        }        
    }    
    for (size_t i = 0; i < nAtypes; i++)
    {
         gmx_stats_free(core[i]);
         gmx_stats_free(shell[i]);
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
    gmx_atomprop_t                   aps;
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
    aps = gmx_atomprop_init();    
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
        readPoldata(i.c_str(), pd, aps);
        pds.push_back(std::move(pd));
    }    

    // Copy the first gentop file into pdout.
    readPoldata(filenames[0].c_str(), pdout, aps); 
    
    //We now update different parts of pdout.    
    EemAtomProps eem = name2eemprop(eemprop[0]);        
    if (eem == eEMEta || eem == eEMAll)
    {
        merge_Eta(pds, pdout);
    }
    if (eem == eEMChi || eem == eEMAll)
    {
        merge_Chi(pds, pdout);
    }
    if (eem == eEMAlpha || eem == eEMAll)
    {
        merge_Alpha(pds, pdout);
    }
    if (eem == eEMZeta || eem == eEMAll)
    {
        merge_Zeta(pds, pdout);
    }
    if (eem == eEMNR)
    {
        gmx_fatal(FARGS, "There is no atomic electric property called %s in alexandria.\n", eemprop[0]);
    }    
    alexandria::writePoldata(opt2fn("-do", NFILE, fnm), &pdout, bcompress);
    if (opt2bSet("-latex", NFILE, fnm))
    {
        FILE        *tp;
        tp = gmx_ffopen(opt2fn("-latex", NFILE, fnm), "w");
        alexandria_eemprops_table(tp, &pdout);
        gmx_ffclose(tp);
    }           
    
    return 0;
}


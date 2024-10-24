/*
 * This source file is part of the Alexandria Chemistry Toolkit.
 *
 * Copyright (C) 2014-2024
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

#include <map>
#include <set>

#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "act/alexandria/alex_modules.h"
#include "act/basics/identifier.h"
#include "act/forcefield/forcefield.h"
#include "act/forcefield/forcefield_tables.h"
#include "act/forcefield/forcefield_xml.h"
#include "gromacs/commandline/pargs.h"
#include "gromacs/math/units.h"
#include "gromacs/math/vec.h"
#include "gromacs/statistics/statistics.h"
#include "gromacs/utility/arraysize.h"
#include "gromacs/utility/cstringutil.h"
#include "gromacs/utility/exceptions.h"
#include "gromacs/utility/fatalerror.h"
#include "gromacs/utility/futil.h"
#include "gromacs/utility/smalloc.h"

namespace alexandria
{

static void merge_parameter(const std::vector<alexandria::ForceField> &pds,
                            alexandria::InteractionType                iType,
                            const std::string                         &parameter,
                            alexandria::ForceField                    *pdout,
                            double                                     limits)
{
    std::vector<gmx_stats> lsq;
    std::vector<int>       ntrain;
    bool                   first = true;
    for (const auto &pd : pds)
    {
        // Index for lsq vector
        int  j  = 0;
        // Loop over forces
        if (pd.interactionPresent(iType))
        {
            auto fs = pd.findForcesConst(iType);
            for (const auto &fp : fs.parametersConst())
            {
                for (const auto &pp : fp.second)
                {
                    if (parameter == pp.first)
                    {
                        if (first)
                        {
                            gmx_stats newlsq;
                            lsq.push_back(newlsq);
                            ntrain.push_back(0);
                        }
                        lsq[j].add_point(0, pp.second.value(), 0, 0);
                        ntrain[j] += pp.second.ntrain();
                        j++;
                    }
                }
            }
        }
        else
        {
            // Loop over particles
            for (const auto &pp : pd.particleTypesConst())
            {
                for(const auto &ppar : pp.second.parametersConst())
                {
                    if (parameter == ppar.first)
                    {
                        if (first)
                        {
                            gmx_stats newlsq;
                            lsq.push_back(newlsq);
                            ntrain.push_back(0);
                        }
                        lsq[j].add_point(0, ppar.second.value(), 0, 0);
                        ntrain[j] += ppar.second.ntrain();
                        j++;
                    }
                }
            }
        }
        first = false;
    }
    int j = 0;
    if (pdout->interactionPresent(iType))
    {
        auto forces = pdout->findForces(iType);
        auto param  = forces->parameters();
        for (auto &fp : *param)
        {
            for (auto &pp : fp.second)
            {
                if (parameter == pp.first)
                {
                    if (pp.second.mutability() == alexandria::Mutability::Free ||
                        pp.second.mutability() == alexandria::Mutability::Bounded)
                    {
                        real   average = 0;
                        real   sigma   = 0;
                        if ((eStats::OK == lsq[j].get_average(&average)) &&
                            (eStats::OK == lsq[j].get_sigma(&sigma)))
                        {
                            pp.second.setValue(average);
                            pp.second.setUncertainty(sigma);
                            pp.second.setNtrain(ntrain[j]/pds.size());
                            if (limits > 0)
                            {
                                pp.second.setMinimum(average-limits*sigma);
                                pp.second.setMaximum(average+limits*sigma);
                            }
                        }
                    }
                    j++;
                }
            }
        }
    }
    else
    {
        // Loop over particles
        for (auto &pp : *pdout->particleTypes())
        {
            auto mypar = pp.second.parameters();
            for(auto &ppar : *mypar)
            {
                if (parameter == ppar.first)
                {
                    if (ppar.second.mutability() == alexandria::Mutability::Free ||
                        ppar.second.mutability() == alexandria::Mutability::Bounded)
                    {
                        real average = 0;
                        real sigma   = 0;
                        if ((eStats::OK == lsq[j].get_average(&average)) &&
                            (eStats::OK == lsq[j].get_sigma(&sigma)))
                        {
                            ppar.second.setValue(average);
                            ppar.second.setUncertainty(sigma);
                            ppar.second.setNtrain(ntrain[j]/pds.size());
                            if (limits > 0)
                            {
                                ppar.second.setMinimum(average-limits*sigma);
                                ppar.second.setMaximum(average+limits*sigma);
                            }
                        }
                    }
                    j++;
                }
            }
        }
    }
}

int merge_ff(int argc, char *argv[])
{
    const char *desc[] =
    {
        "merge_ff reads multiple Alexandria force field files and merges them",
        "into a single new force field file containing average values",
        "and standard deviation.",
    };    
    t_filenm    fnm[] =
    {
        { efXML, "-ff",    "aff_in",  ffRDMULT},
        { efXML, "-o",     "aff_out", ffWRITE },
        { efTEX, "-latex", "aff_out", ffWRITE }
    };
    int         NFILE       = asize(fnm);;
    
    gmx_bool    bcompress   = false;
    gmx_bool    bAtypeMap   = false;
    //! String for command line to harvest the options to fit
    char       *mergeString = nullptr;
    char       *info        = nullptr;
    int         ntrain      = 1;
    real        limits      = 0;
    t_pargs     pa[]        =
    {
        { "-compress", FALSE, etBOOL, {&bcompress},
          "Compress output XML files" },
        { "-merge", FALSE, etSTR, {&mergeString},
          "Quoted list of parameters to merge,  e.g. 'alpha zeta'. An empty string means all parameters will be merged." },
        { "-amap", FALSE, etBOOL, {&bAtypeMap},
          "Add latex table for mapping from particletype to subtypes" },
        { "-info", FALSE, etSTR, {&info},
          "Extra information to print in table captions" },
        { "-ntrain", FALSE, etINT, {&ntrain},
          "Include only variables that have their ntrain values larger or equal to this number." },
        { "-limits", FALSE, etREAL, {&limits},
          "Change the boundaries for parameters that have 1) mutability Bounded, 2) ntrain large enough (see previous option) and 3) have a standard deviation larger than zero. If limits equals zero nothing will be done, else the boundaries will be changed to +/- limits times the standard deviation, taking into account whether parameters should be non-negative." }
    };
    std::vector<alexandria::ForceField> pds;
    alexandria::ForceField              pdout;
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
    std::string allParams("alpha chi eta zeta delta_eta delta_chi charge vs2a");
    if (nullptr == mergeString)
    {
        mergeString = strdup(allParams.c_str());
    }
    /* Read all the gentop files. */
    auto filenames = opt2fns("-ff", NFILE, fnm);
    if (filenames.size() < 2)
    {
        gmx_fatal(FARGS, "At least two gentop files are needed for merging!\n");
    }
    
    for (auto &i : filenames)
    {
        alexandria::ForceField pd;
        readForceField(i.c_str(), &pd);
        pds.push_back(std::move(pd));
    }    

    // Copy the first gentop file into pdout
    readForceField(filenames[0].c_str(), &pdout);
    
    // We now update different parts of pdout
    std::set <alexandria::InteractionType> itypes;
    auto allparams = gmx::splitString(mergeString);
    if (allparams.empty())
    {
        for(const auto &fs : pds[0].forcesConst())
        {
            if (!fs.second.empty())
            {
                auto pp = fs.second.parametersConst().begin();
                for(const auto &type : pp->second)
                {
                    merge_parameter(pds, fs.first, type.first, &pdout, limits);
                    itypes.insert(fs.first);
                }
            }
        }
    }
    else
    {
        for(const auto &type : allparams)
        {
            alexandria::InteractionType itype;
            if (pdout.typeToInteractionType(type, &itype))
            {
                merge_parameter(pds, itype, type, &pdout, limits);
                itypes.insert(itype);
            }
            else
            {
                fprintf(stderr, "No such parameter type '%s' in %s\n",
                        type.c_str(), filenames[0].c_str());
            }
        }
    }

    writeForceField(opt2fn("-o", NFILE, fnm), &pdout, bcompress);
    if (opt2bSet("-latex", NFILE, fnm))
    {
        FILE *tp = gmx_ffopen(opt2fn("-latex", NFILE, fnm), "w");
        ForceFieldTable fft(tp, &pdout, ntrain);
        std::string myinfo;
        if (info != nullptr)
        {
            myinfo.assign(info);
        }
        if (bAtypeMap)
        {
            fft.subtype_table(myinfo);
        }
        //fft.zeta_table(myinfo);
        //fft.eemprops_table(myinfo);
        for(auto itype : itypes)
        {
            fft.itype_table(itype, myinfo);
        }
        gmx_ffclose(tp);
    }           
    
    return 0;
}

} // namespace alexandria

/*
 * This source file is part of the Alexandria Chemistry Toolkit.
 *
 * Copyright (C) 2014-2026
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

#include <algorithm>
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
#include "act/statistics/statistics.h"
#include "gromacs/commandline/pargs.h"
#include "gromacs/math/units.h"
#include "gromacs/math/vec.h"
#include "gromacs/utility/exceptions.h"
#include "gromacs/utility/fatalerror.h"
#include "gromacs/utility/futil.h"
#include "gromacs/utility/smalloc.h"

namespace alexandria
{

/*! \brief Merge parameters from multiple force fields into one
 * The parameters are matched by the identifiers of the interactions, or of
 * the particle types if the parameter is not of an interaction. Each
 * parameter with mutability Bounded gets the average of its values in the
 * force fields that have it, and their standard deviation as uncertainty.
 * \param[in]    pds       The force fields
 * \param[in]    iType     The interaction type of the parameter
 * \param[in]    parameter The name of the parameter
 * \param[inout] pdout     The merged force field, a copy of the first one
 * \param[in]    limits    If positive, the bounds become the average plus
 *                         and minus this number of standard deviations
 */
static void merge_parameter(const std::vector<alexandria::ForceField> &pds,
                            alexandria::InteractionType                iType,
                            const std::string                         &parameter,
                            alexandria::ForceField                    *pdout,
                            double                                     limits)
{
    // The values and the sum of the training points of each identifier
    std::map<Identifier, gmx_stats> lsq;
    std::map<Identifier, int>       ntrain;
    auto add = [&](const Identifier &id, const ForceFieldParameterMap &params)
    {
        auto pp = params.find(parameter);
        if (pp != params.end())
        {
            lsq[id].add_point(0, pp->second.value(), 0, 0);
            ntrain[id] += pp->second.ntrain();
        }
    };
    for (const auto &pd : pds)
    {
        if (pd.interactionPresent(iType))
        {
            for (const auto &fp : pd.findForcesConst(iType).parametersConst())
            {
                add(fp.first, fp.second);
            }
        }
        else
        {
            for (const auto &pt : pd.particleTypesConst())
            {
                add(pt.first, pt.second.parametersConst());
            }
        }
    }
    auto update = [&](const Identifier &id, ForceFieldParameterMap *params)
    {
        auto pp = params->find(parameter);
        auto st = lsq.find(id);
        if (pp == params->end() || st == lsq.end() ||
            pp->second.mutability() != alexandria::Mutability::Bounded)
        {
            return;
        }
        real average = 0;
        real sigma   = 0;
        if ((eStats::OK == st->second.get_average(&average)) &&
            (eStats::OK == st->second.get_sigma(&sigma)))
        {
            auto &p = pp->second;
            // The average of values within their bounds may be outside the
            // bounds of the first force field
            p.setMinimum(std::min(p.minimum(), static_cast<double>(average)));
            p.setMaximum(std::max(p.maximum(), static_cast<double>(average)));
            p.setValue(average);
            p.setUncertainty(sigma);
            p.setNtrain(ntrain[id]/st->second.get_npoints());
            if (limits > 0)
            {
                p.setMinimum(average-limits*sigma);
                p.setMaximum(average+limits*sigma);
            }
        }
    };
    if (pdout->interactionPresent(iType))
    {
        for (auto &fp : *pdout->findForces(iType)->parameters())
        {
            update(fp.first, &fp.second);
        }
    }
    else
    {
        for (auto &pt : *pdout->particleTypes())
        {
            update(pt.first, pt.second.parameters());
        }
    }
}

int merge_ff(int argc, char *argv[])
{
    std::vector<const char *> desc =
    {
        "merge_ff reads multiple Alexandria force field files and merges them",
        "into a single new force field file containing average values",
        "and standard deviation.",
    };    
    std::vector<t_filenm> fnm =
    {
        { efXML, "-ff",    "aff_in",  ffRDMULT},
        { efXML, "-o",     "aff_out", ffWRITE },
        { efTEX, "-latex", "aff_out", ffWRITE }
    };
    gmx_bool    bcompress   = false;
    gmx_bool    bAtypeMap   = false;
    //! String for command line to harvest the options to fit
    char       *mergeString = nullptr;
    char       *info        = nullptr;
    int         ntrain      = 1;
    real        limits      = 0;
    std::vector<t_pargs> pa =
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
                           fnm.size(), fnm.data(),
                           pa.size(), pa.data(), 
                           desc.size(), desc.data(), 
                           0, 
                           nullptr, 
                           &oenv))
    {
        return 0;
    }
    /* Read all the gentop files. */
    auto filenames = opt2fns("-ff", fnm.size(), fnm.data());
    
    for (auto &i : filenames)
    {
        alexandria::ForceField pd;
        readForceField(i.c_str(), &pd);
        pds.push_back(std::move(pd));
    }    

    // Copy the first gentop file into pdout
    readForceField(filenames[0].c_str(), &pdout);

    // Let's see what parameters to handle
    std::vector<std::string> allparams;
    if (nullptr != mergeString)
    {
        allparams = gmx::splitString(mergeString);
    }
    // We now update different parts of pdout
    std::set <alexandria::InteractionType> itypes;
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

    bool printSigma = filenames.size() > 1;
    writeForceField(opt2fn("-o", fnm.size(), fnm.data()), &pdout, bcompress);
    if (opt2bSet("-latex", fnm.size(), fnm.data()))
    {
        FILE *tp = gmx_ffopen(opt2fn("-latex", fnm.size(), fnm.data()), "w");
        ForceFieldTable fft(tp, &pdout, ntrain, printSigma);
        std::string myinfo;
        if (info != nullptr)
        {
            myinfo.assign(info);
        }
        if (bAtypeMap)
        {
            fft.subtype_table(myinfo);
        }
        for(auto itype : itypes)
        {
            fft.itype_table(itype, myinfo);
        }
        gmx_ffclose(tp);
    }           
    
    return 0;
}

} // namespace alexandria

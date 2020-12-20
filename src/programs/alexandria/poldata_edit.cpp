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
#include "forcefieldparameter.h"
#include "interactiontype.h"
#include "mutability.h"
#include "poldata.h"
#include "poldata_xml.h"

static void setMinMaxMut(alexandria::ForceFieldParameter *pp,
                         bool bSetMin, double pmin,
                         bool bSetVal, double pval,
                         bool bSetMax, double pmax,
                         bool bSetMut, const char *mutability)
{
    if (bSetVal && bSetMin)
    {
        GMX_RELEASE_ASSERT(pmin <= pval, "Please choose value to be at least the minimum");
    }
    if (bSetVal && bSetMax)
    {
        GMX_RELEASE_ASSERT(pmax >= pval, "Please choose value to be at most the maximum");
    }
    if (bSetMin && bSetMax)
    {
        GMX_RELEASE_ASSERT(pmin < pmax, "Please choose minimum to be at most the maximum");
    }
    if (bSetVal)
    {
        pp->setValue(pval);
    }
    if (bSetMin)
    {
        pp->setMinimum(pmin);
        pp->setValue(std::max(pmin, pp->value()));
        pp->setMaximum(std::max(pmin, pp->maximum()));
    }
    if (bSetMax)
    {
        pp->setMaximum(pmax);
        pp->setValue(std::min(pmax, pp->value()));
        pp->setMinimum(std::min(pmax, pp->minimum()));
    }
    if (bSetMut)
    {
        alexandria::Mutability mut;
        if (nameToMutability(mutability, &mut))
        {
            pp->setMutability(mut);
        }
    }
}

static void modifyPoldata(alexandria::Poldata *pd,
                          const char *ptype,
                          const char *particle,
                          bool bSetMin, double pmin,
                          bool bSetVal, double pval,
                          bool bSetMax, double pmax,
                          bool bSetMut, const char *mutability,
                          bool force)
{
    if (!(bSetMin || bSetMax || bSetMut))
    {
        return;
    }
    if (strlen(ptype) == 0)
    {
        printf("Empty parameter type. Nothing to modify.\n");
        return;
    }
    bool done = false;
    // Check whether we have a parameter type belonging to atomtypes
    if ((strlen(particle) > 0) && pd->hasParticleType(particle))
    {
        auto pt = pd->findParticleType(particle);
        if (pt->hasParameter(ptype))
        {
            setMinMaxMut(pt->parameter(ptype),
                         bSetMin, pmin, bSetVal, pval, bSetMax, pmax,
                         bSetMut, mutability);
            done = true;
        }
    }
    else
    {
        auto particleTypes = pd->particleTypes();
        for (auto pt=particleTypes->begin(); pt < particleTypes->end(); ++pt)
        {
            if (pt->hasParameter(ptype))
            {
                setMinMaxMut(pt->parameter(ptype),
                             bSetMin, pmin, bSetVal, pval, bSetMax, pmax,
                             bSetMut, mutability);
                done = true;
            }
        }
    }
    if (done)
    {
        return;
    }    
    alexandria::InteractionType itype; 
    if (!pd->typeToInteractionType(ptype, &itype))
    {
        return;
    }
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
                        setMinMaxMut(&pp.second,
                                     bSetMin, pmin, 
                                     bSetVal, pval, bSetMax, pmax,
                                     bSetMut, mutability);
                    }
                }
            }
        }
    }
}

static void analyzePoldata(alexandria::Poldata *pd)
{
    size_t mindata   = 1;
    double tolerance = 0.001;
    for(auto &fc : pd->forcesConst())
    {
        auto itype = fc.first;
        for(auto &fm : fc.second.parametersConst())
        {
            auto myid = fm.first;
            for(auto &param : fm.second)
            {
                auto ptype = param.first;
                auto ppp   = param.second;
                if (ppp.ntrain() >= mindata)
                {
                    if (fabs(ppp.value() - ppp.minimum()) < tolerance)
                    {
                        printf("%s %s %s at minimum %g\n",
                               alexandria::interactionTypeToString(itype).c_str(),
                               myid.id().c_str(),
                               ptype.c_str(), ppp.minimum());
                    }
                    else if (fabs(ppp.value() - ppp.maximum()) < tolerance)
                    {
                        printf("%s %s %s at maximum %g\n",
                               alexandria::interactionTypeToString(itype).c_str(),
                               myid.id().c_str(),
                               ptype.c_str(), ppp.maximum());
                    }
                }
            }
        }
    }
}

static void compare_pd(alexandria::Poldata *pd1,
                       alexandria::Poldata *pd2,
                       const std::string   &ptype)
{
    printf("Comparing %s and %s\n",
           pd1->filename().c_str(), pd2->filename().c_str());
    alexandria::InteractionType itype1, itype2;
    if (!pd1->typeToInteractionType(ptype, &itype1) ||
        !pd2->typeToInteractionType(ptype, &itype2))
    {
        return;
    }
    auto f1     = pd1->findForcesConst(itype1);
    auto f2     = pd2->findForcesConst(itype2);
    printf("%s %6s %9s %9s %9s   %9s %9s %9s\n",
           "Param", "Ptype", 
           "Minimum_1", "Value_1", "Maximum_1",
           "Minimum_2", "Value_2", "Maximum_2");
    double psum2 = 0;
    int    np2   = 0;
    for(auto &id1 : f1.parametersConst())
    {
        if (f2.parameterExists(id1.first))
        {
            auto fp1 = f1.findParametersConst(id1.first);
            auto fp2 = f2.findParametersConst(id1.first);
            auto p1  = fp1[ptype];
            auto p2  = fp2[ptype];
            if (p1.value() != p2.value())
            {
                printf("%s %6s %9g %9g %9g   %9g %9g %9g\n",
                       ptype.c_str(), id1.first.id().c_str(),
                       p1.minimum(), p1.value(), p1.maximum(),
                       p2.minimum(), p2.value(), p2.maximum());
                double d = p1.value()-p2.value();
                psum2   += d*d;
                np2     += 1;
            }
        }
    }
    if (np2 > 1)
    {                                           
        printf("RMSD: %9g\n", std::sqrt(psum2/np2));
    }
}

int alex_poldata_edit(int argc, char*argv[])
{
    static const char               *desc[] = {
        "poldata_edit reads a poldata (force field) file and can do ",
        "a number of things.[PAR]",
        "If the analyze flag is set it will analyze the file and print",
        "if parameters are close to their bounds.[PAR]",
        "It may write a new file that can be compared to the input.",
        "If run in parallel, this utility will read on the first processor,",
        "then process the data and send it to the second processor. That one",
        "will then write the new file. This is for debugging poldata",
        "communication.[PAR]",
        "Modifications of parameters can be made by specifying both",
        "input and output files and what parameters to change.",
        "If the value, the minimum",
        "or the maximum is to be changed, the actual value may be set to the",
        "new minmum or maximum if it falls outside the new bounds."
    };
    gmx_output_env_t                *oenv;
    t_filenm                         fnm[] = {
        { efDAT, "-f",  "pdin" , ffREAD  },
        { efDAT, "-f2", "pdin2", ffOPTRD },
        { efDAT, "-o", "pdout",  ffOPTWR }
    };
    static char *parameter  = (char *)"";
    static char *particle   = (char *)"";
    static char *mutability = (char *)"";
    real         pmin       = 0;
    real         pmax       = 0;
    real         pval       = 0;
    gmx_bool     force      = false;
    gmx_bool     analyze    = false;
    t_pargs                          pa[]     = 
    {
        { "-p",      FALSE, etSTR,  {&parameter},
          "Type of parameter to change, e.g. zeta." },
        { "-a",      FALSE, etSTR,  {&particle},
          "Particle type to change, if not set all particles of the correct type will be changed." },
        { "-min",    FALSE, etREAL, {&pmin},
          "Minimum value of parameter." },
        { "-val",    FALSE, etREAL, {&pval},
          "New value of parameter." },
        { "-max",    FALSE, etREAL,  {&pmax},
          "Maximum value of parameter." },
        { "-mut",    FALSE, etSTR,  {&mutability},
          "Set the mutability for the given parameters to the value. The following values are supported: Free, Fixed, Bounded" },
        { "-force",  FALSE, etBOOL, {&force},
          "Will change also non-mutable parameters. Use with care!" },
        { "-analyze", FALSE, etBOOL, {&analyze},
          "Analyze the parameters in a simple manner" }
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
                      opt2parg_bSet("-val", npargs, pa), pval,
                      opt2parg_bSet("-max", npargs, pa), pmax,
                      opt2parg_bSet("-mut", npargs, pa), mutability,
                      force);
        if (analyze)
        {
            analyzePoldata(&pd);
        }
        if (opt2bSet("-f2", NFILE, fnm))
        {
            alexandria::Poldata pd2;
            try 
            {
                alexandria::readPoldata(opt2fn("-f2", NFILE, fnm), &pd2);
            }
            GMX_CATCH_ALL_AND_EXIT_WITH_FATAL_ERROR;
            compare_pd(&pd, &pd2, parameter);
        }
    }
    if (opt2bSet("-o", NFILE, fnm))
    {
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
    }
    return 0;
}

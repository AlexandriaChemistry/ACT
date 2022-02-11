/*
 * This source file is part of the Alexandria Chemistry Toolkit.
 *
 * Copyright (C) 2014-2021 
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

#include <cstdlib>

#include <map>
#include <set>

#include "gromacs/commandline/pargs.h"
#include "gromacs/fileio/oenv.h"
#include "gromacs/utility/arraysize.h"

#include "act/poldata/act_checksum.h"
#include "alexandria/alex_modules.h"
#include "act/poldata/forcefieldparameter.h"
#include "act/basics/interactiontype.h"
#include "act/basics/mutability.h"
#include "act/poldata/poldata.h"
#include "act/poldata/poldata_xml.h"

namespace alexandria
{

static void setMinMaxMut(FILE *fp,
                         ForceFieldParameter *pp,
                         bool bSetMin, double pmin,
                         bool bSetVal, double pval,
                         bool bSetMax, double pmax,
                         bool bSetMut, const std::string &mutability,
                         bool bScale,  double scale,
                         bool stretch, const std::string &particleId)
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
    if (bSetMin)
    {
        pp->setMinimum(pmin);
        pp->setMaximum(std::max(pmin, pp->maximum()));
        pp->setValue(std::max(pmin, pp->value()));
        if (fp)
        {
            fprintf(fp, "Minimum set to %g for %s\n", pmin, particleId.c_str());
        }
    }
    if (bSetMax)
    {
        pp->setMaximum(pmax);
        pp->setValue(std::min(pmax, pp->value()));
        pp->setMinimum(std::min(pmax, pp->minimum()));
        if (fp)
        {
            fprintf(fp, "Maximum set to %g %s\n", pmax, particleId.c_str());
        }
    }
    if (bSetVal || bScale)
    {
        if (bScale)
        {
            printf("scaling by %g\n", scale);
            if (bSetVal && fp)
            {
                fprintf(fp, "Ignoring -val when using -scale\n");
            }
            pval = scale * pp->value();
        }
        if (pval < pp->minimum())
        {
            pp->setMinimum(pval);
        }
        if (pval > pp->maximum())
        {
            pp->setMaximum(pval);
        }
        pp->setValue(pval);
        if (fp)
        {
            fprintf(fp, "Value set to %g %s\n", pval, particleId.c_str());
        }
    }
    if (bSetMut)
    {
        Mutability mut;
        if (nameToMutability(mutability, &mut))
        {
            pp->setMutability(mut);
        }
    }
    if (stretch)
    {
        auto range = pp->maximum()-pp->minimum();
        if (std::fabs(pp->value() - pp->minimum()) < 0.01*range && range > 0)
        {
            if (!pp->setMinimum(pp->minimum()-0.2*range))
            {
                pp->setMinimum(0.8*pp->minimum());
            }
            if (fp)
            {
                fprintf(fp, "Minimum stretched to %g for %s\n",
                        pp->minimum(), particleId.c_str());
            }
        }
        if (std::fabs(pp->value() - pp->maximum()) < 0.01*range && range > 0)
        {
            pp->setMaximum(pp->maximum()+0.2*range);
            if (fp)
            {
                fprintf(fp, "Maximum stretched to %g for %s\n",
                        pp->maximum(), particleId.c_str());
            }
        }
    }
}

static void modifyParticle(const std::string &paramType,
                           ParticleTypeIterator particle,
                           bool bSetMin, double pmin,
                           bool bSetVal, double pval,
                           bool bSetMax, double pmax,
                           bool bSetMut, const std::string &mutability,
                           bool bScale,  double scale,
                           bool force, bool stretch)
{
    // Check particletypes instead
    if (particle->hasParameter(paramType))
    {
        auto ff = particle->parameter(paramType);
        if (force || ff->isMutable())
        {
            std::string particleId = 
                gmx::formatString("%s - %s",
                                  particle->id().id().c_str(),
                                  paramType.c_str());
            setMinMaxMut(debug, ff,
                         bSetMin, pmin, 
                         bSetVal, pval, bSetMax, pmax,
                         bSetMut, mutability, 
                         bScale,  scale,
                         stretch, particleId);
        }
    }
    else
    {
        printf("No interaction type corresponding to %s\n", paramType.c_str());
    }
}

static void modifyInteraction(Poldata *pd,
                              InteractionType itype,
                              const std::string &paramType,
                              const Identifier &pId,
                              bool bSetMin, double pmin,
                              bool bSetVal, double pval,
                              bool bSetMax, double pmax,
                              bool bSetMut, const std::string &mutability,
                              bool bScale,  double scale,
                              bool force, bool stretch)
{
    auto fs = pd->findForces(itype)->parameters();
    for(auto &ffs : *fs)
    {
        if (pId == ffs.first)
        {
            for(auto &pp : ffs.second)
            {
                if (pp.first == paramType)
                {
                    if (force || pp.second.isMutable())
                    {
                        std::string myId = 
                            gmx::formatString("%s - %s - %s",
                                              pId.id().c_str(),
                                              paramType.c_str(),
                                              interactionTypeToString(itype).c_str());
                        setMinMaxMut(debug, &pp.second,
                                     bSetMin, pmin, 
                                     bSetVal, pval, bSetMax, pmax,
                                     bSetMut, mutability, 
                                     bScale,  scale,
                                     stretch, myId);
                    }
                }
            }
        }
    }
}

static void modifyPoldata(Poldata *pd,
                          const std::string &paramType,
                          const std::string &particle,
                          bool bSetMin, double pmin,
                          bool bSetVal, double pval,
                          bool bSetMax, double pmax,
                          bool bSetMut, const std::string &mutability,
                          bool bScale,  double scale,
                          bool force, bool stretch)
{
    if (!(bSetVal || bSetMin || bSetMax || bSetMut || stretch || bScale))
    {
        printf("No parameter to change.\n");
        return;
    }
    if (paramType.empty())
    {
        printf("Empty parameter type. Nothing to modify.\n");
        return;
    }
    std::vector<ParticleTypeIterator> myParticles;
    if (!particle.empty())
    {
        for (const auto &s : gmx::splitString(particle))
        {
            if (pd->hasParticleType(s))
            {
                myParticles.push_back(pd->findParticleType(s));
            }
        }
    }
    else
    {
        auto particleTypes = pd->particleTypes();
        for (auto pt=particleTypes->begin(); pt < particleTypes->end(); ++pt)
        {
            myParticles.push_back(pt);
        }
    }
    if (myParticles.empty())
    {
        printf("Cannot find any particle %s\n", particle.c_str());
    }
    for (auto p : myParticles)
    {
        InteractionType itype;
        if (pd->typeToInteractionType(paramType, &itype))
        {
            if (!p->hasInteractionType(itype))
            {
                continue;
            }
            auto pId    = p->interactionTypeToIdentifier(itype);
            if (pId.id().empty())
            {
                continue;
            }
            auto natoms = interactionTypeToNatoms(itype);
            switch (natoms)
            {
            case 1:
                {
                    modifyInteraction(pd, itype, paramType, pId,
                                      bSetMin, pmin,
                                      bSetVal, pval,
                                      bSetMax, pmax,
                                      bSetMut, mutability,
                                      bScale,  scale,
                                      force, stretch);
                    break;
                }
            case 2:
                {
                    auto q1id = pId.id();
                    for (auto q2: myParticles)
                    {
                        auto q2id = q2->interactionTypeToIdentifier(itype).id();
                        const double bondorders[] = { 1, 1.5, 2, 3 };
                        const size_t nBondorder   = std::extent<decltype(bondorders)>::value;
                        for(size_t bb = 0; bb < nBondorder; bb++)
                        {
                            auto qId = Identifier({q1id, q2id}, { bondorders[bb] }, CanSwap::No);
                            modifyInteraction(pd, itype, paramType, qId,
                                              bSetMin, pmin,
                                              bSetVal, pval,
                                              bSetMax, pmax,
                                              bSetMut, mutability,
                                              bScale,  scale,
                                              force, stretch);
                        }
                    }
                    break;
                }
            default:
                fprintf(stderr, "Don't know how to handle interactions with %d atoms", natoms);
            }
        }
        else
        {
            modifyParticle(paramType, p,
                           bSetMin, pmin,
                           bSetVal, pval,
                           bSetMax, pmax,
                           bSetMut, mutability,
                           bScale,  scale,
                           force, stretch);
        }
    }
}

static const std::set<InteractionType> &findInteractionMap(const std::string &analyze,
                                                           bool              *found)
{
    static std::set<InteractionType> bonds = {
        InteractionType::BONDS,
        InteractionType::ANGLES,
        InteractionType::LINEAR_ANGLES,
        InteractionType::PROPER_DIHEDRALS,
        InteractionType::IMPROPER_DIHEDRALS };
    static std::set<InteractionType> other = {
        InteractionType::VDW,
        InteractionType::CONSTR,
        InteractionType::VSITE2,
        InteractionType::VSITE3FAD,
        InteractionType::VSITE3OUT };
    static std::set<InteractionType> eem = {
        InteractionType::POLARIZATION,
        InteractionType::CHARGEDISTRIBUTION,
        InteractionType::BONDCORRECTIONS,
        InteractionType::ELECTRONEGATIVITYEQUALIZATION };
    static std::map<const std::string, std::set<InteractionType> > mymap = {
        { "EEM", eem }, { "BONDS", bonds }, { "OTHER", other } 
    };
    if (mymap.find(analyze)== mymap.end())
    {
        printf("Don't know how to analyze this group %s\n", analyze.c_str());
        *found = false;
        return mymap["EEM"];
    }
    else
    {
        *found = true;
        return mymap[analyze];
    }
}

static void analyzePoldata(Poldata           *pd,
                           const std::string &analyze)
{
    bool found;
    
    auto myset = findInteractionMap(analyze, &found);
    if (!found)
    {
        return;
    }
    int    mindata   = 1;
    double tolerance = 0.001;
    for(auto &fc : pd->forcesConst())
    {
        auto itype = fc.first;
        if (myset.find(itype) == myset.end())
        {
            continue;
        }
        for(auto &fm : fc.second.parametersConst())
        {
            auto myid = fm.first;
            for(auto &param : fm.second)
            {
                auto ptype = param.first;
                auto ppp   = param.second;
                if (ppp.ntrain() >= mindata && ppp.isMutable())
                {
                    if (fabs(ppp.value() - ppp.minimum()) < tolerance)
                    {
                        printf("%s %s %s at minimum %g\n",
                               interactionTypeToString(itype).c_str(),
                               myid.id().c_str(),
                               ptype.c_str(), ppp.minimum());
                    }
                    else if (fabs(ppp.value() - ppp.maximum()) < tolerance)
                    {
                        printf("%s %s %s at maximum %g\n",
                               interactionTypeToString(itype).c_str(),
                               myid.id().c_str(),
                               ptype.c_str(), ppp.maximum());
                    }
                }
            }
        }
    }
}

static void dumpPoldata(Poldata           *pd,
                        const std::string &analyze,
                        const std::string &particle,
                        const std::string &filenm)
{
    bool found;
    
    auto myset = findInteractionMap(analyze, &found);
    if (!found)
    {
        return;
    }
    auto particles = gmx::splitString(particle);
    FILE *fp = fopen(filenm.c_str(), "w");
    int nparm = 0;
    for(auto &fc : pd->forcesConst())
    {
        auto itype = fc.first;
        if (myset.find(itype) == myset.end())
        {
            continue;
        }
        for(auto &fm : fc.second.parametersConst())
        {
            auto myid = fm.first;
            if (std::find(particles.begin(), particles.end(), myid.id()) != particles.end())
            {
                for(auto &param : fm.second)
                {
                    auto ptype = param.first;
                    auto ppp   = param.second;
                    fprintf(fp, "%s|%s|%g\n",
                            ptype.c_str(), myid.id().c_str(), ppp.value());
                    nparm++;
                }
            }
        }
    }
    fclose(fp);
    printf("Stored %d parameters in %s\n", nparm, filenm.c_str());
}

static void copy_missing(const Poldata     *pdref,
                         Poldata           *pdout,
                         const std::string &analyze,
                         bool               replace)
{
    bool found;
    auto myset = findInteractionMap(analyze, &found);
    if (!found)
    {
        return;
    }
    if (replace)
    {
        printf("Replacing %s interactions from %s to %s\n",
               analyze.c_str(),
               pdref->filename().c_str(), pdout->filename().c_str());
    }
    else
    {
        printf("Copying %s missing interactions from %s to %s\n",
               analyze.c_str(),
               pdref->filename().c_str(), pdout->filename().c_str());
    }
    for(auto &fc : pdref->forcesConst())
    {
        auto itype = fc.first;
        if (myset.find(itype) == myset.end())
        {
            continue;
        }
        if (!pdout->interactionPresent(itype))
        {
            printf("Cannot find interaction %s in %s, giving up\n",
                   interactionTypeToDescription(itype).c_str(), pdout->filename().c_str());
            return;
        }
        auto fcout = pdout->findForces(itype);
        if (replace)
        {
            fcout->clearParameters();
        }
        for(auto &fm : fc.second.parametersConst())
        {
            auto myid = fm.first;
            if (!fcout->parameterExists(myid))
            {
                printf("Parameter %s missing from %s for %s\n",
                       myid.id().c_str(), pdout->filename().c_str(),
                       interactionTypeToDescription(itype).c_str());
                for(auto &param : fm.second)
                {
                    auto ptype = param.first;
                    auto ppp   = param.second;
                    fcout->addParameter(myid, ptype, ppp);
                }
            }
        }
    }
}

static void implant_values(const Poldata     *pdref,
                           Poldata           *pdout,
                           const std::string &analyze)
{
    bool found;
    auto myset = findInteractionMap(analyze, &found);
    if (!found)
    {
        return;
    }
    printf("Implanting optimized values for %s interactions from %s to %s\n",
           analyze.c_str(),
           pdref->filename().c_str(), pdout->filename().c_str());
    
    for(auto &fc : pdref->forcesConst())
    {
        auto itype = fc.first;
        if (myset.find(itype) == myset.end())
        {
            continue;
        }
        if (!pdout->interactionPresent(itype))
        {
            printf("Cannot find interaction %s in %s, giving up\n",
                   interactionTypeToDescription(itype).c_str(), pdout->filename().c_str());
            return;
        }
        auto fcout = pdout->findForces(itype);
        for(auto &fm : fc.second.parametersConst())
        {
            auto myid = fm.first;
            if (fcout->parameterExists(myid))
            {
                for(auto &fmp : fm.second)
                {
                    auto pout = fcout->findParameterType(myid, fmp.first);
                
                    if (fmp.second.ntrain() > pout->ntrain())
                    {
                        printf("Parameter %s , interaction %s will be replaced in %s\n",
                               myid.id().c_str(), interactionTypeToDescription(itype).c_str(),
                               pdout->filename().c_str());
                        pout->copy(fmp.second);
                    }
                    else
                    {
                        printf("Parameter %s , interaction %s will not be replaced in %s, old one has ntrain %d new %d\n",
                               myid.id().c_str(), interactionTypeToDescription(itype).c_str(),
                               pdout->filename().c_str(), pout->ntrain(), fmp.second.ntrain());
                    }
                }
            }
        }
    }
}

static void compare_pd(Poldata *pd1,
                       Poldata *pd2,
                       const std::string   &ptype)
{
    printf("Comparing %s and %s\n",
           pd1->filename().c_str(), pd2->filename().c_str());
    InteractionType itype1, itype2;
    if (!pd1->typeToInteractionType(ptype, &itype1) ||
        !pd2->typeToInteractionType(ptype, &itype2))
    {
        return;
    }
    auto f1     = pd1->findForcesConst(itype1);
    auto f2     = pd2->findForcesConst(itype2);
    printf("%-17s %11s %9s %9s %9s %4s | %9s %9s %9s %4s | %9s\n",
           "Param", "Ptype", 
           "Minimum_1", "Value_1", "Maximum_1", "N_1",
           "Minimum_2", "Value_2", "Maximum_2", "N_2", "Differ.");
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
            if (p1.ntrain() > 0 && p2.ntrain() > 0 &&
                p1.value() != p2.value())
            {
                double d = p1.value()-p2.value();
                psum2   += d*d;
                np2     += 1;
                printf("%-17s %11s %9.4f %9.4f %9.4f %4d | %9.4f %9.4f %9g %4d | %9.4f\n",
                       ptype.c_str(), id1.first.id().c_str(),
                       p1.minimum(), p1.value(), p1.maximum(), p1.ntrain(),
                       p2.minimum(), p2.value(), p2.maximum(), p2.ntrain(),
                       d);
            }
        }
    }
    if (np2 > 1)
    {                                           
        printf("RMSD: %9g\n", std::sqrt(psum2/np2));
    }
}

int poldata_edit(int argc, char*argv[])
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
        "new minmum or maximum if it falls outside the new bounds.",
    };
    gmx_output_env_t                *oenv;
    t_filenm                         fnm[] = {
        { efXML, "-f",  "pdin" , ffREAD  },
        { efXML, "-f2", "pdin2", ffOPTRD },
        { efXML, "-o", "pdout",  ffOPTWR },
        { efDAT, "-dump", "params", ffOPTWR }
    };

    static char *parameter  = (char *)"";
    static char *particle   = (char *)"";
    static char *mutability = (char *)"";
    real         pmin       = 0;
    real         pmax       = 0;
    real         pval       = 0;
    real         scale      = 1;
    gmx_bool     force      = false;
    gmx_bool     stretch    = false;
    static char *missing    = (char *)"";
    static char *replace    = (char *)"";
    static char *implant    = (char *)"";
    static char *analyze    = (char *)"";
    t_pargs      pa[]       =
    {
        { "-p",      FALSE, etSTR,  {&parameter},
          "Type of parameter to change, e.g. zeta." },
        { "-a",      FALSE, etSTR,  {&particle},
          "Particle type to change, if not set all particles of the correct type will be changed. A quoted string may be passed to adjust multiple particle types at once, e.g. -a 'h1_z hc_z'" },
        { "-min",    FALSE, etREAL, {&pmin},
          "Minimum value of parameter." },
        { "-val",    FALSE, etREAL, {&pval},
          "New value of parameter." },
        { "-scale",  FALSE, etREAL, {&scale},
          "Scale the value of a parameter by this number, e.g. -scale 0.99." },
        { "-max",    FALSE, etREAL,  {&pmax},
          "Maximum value of parameter." },
        { "-mut",    FALSE, etSTR,  {&mutability},
          "Set the mutability for the given parameters to the value. The following values are supported: Free, Fixed, Bounded" },
        { "-force",  FALSE, etBOOL, {&force},
          "Will change also non-mutable parameters. Use with care!" },
        { "-stretch", FALSE, etBOOL, {&stretch},
          " Will automatically stretch boundaries for individual parameters" },
        { "-ana", FALSE, etSTR, {&analyze},
          "Analyze either the EEM, the BONDED or OTHER parameters in a simple manner" },
        { "-copy_missing", FALSE, etSTR, {&missing},
          "Copy either the EEM, the BONDED or OTHER parameters from file two [TT]-f2[tt] that are missing from file one [TT]-f[tt] to another [TT]-o[tt]." },
        { "-replace", FALSE, etSTR, {&replace},
          "Replace either the EEM, the BONDED or OTHER parameters in file one [TT]-f[ff] by those from file two [TT]-f2[tt] and store in another [TT]-o[tt]." },
        { "-implant", FALSE, etSTR, {&implant},
          "Implant (write over) either the EEM, the BONDED or OTHER parameters in file one [TT]-f[ff] by those from file two [TT]-f2[tt] and store in another [TT]-o[tt]." }
    };
    int                 npargs = asize(pa);
    int                 NFILE  = asize(fnm);
    alexandria::Poldata pd;
    
    CommunicationRecord cr;
    if (cr.isMaster())
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
        (void) pd.verifyCheckSum(stderr, poldataCheckSum(&pd));
        if (opt2bSet("-f2", NFILE, fnm))
        {
            alexandria::Poldata pd2;
            try
            {
                alexandria::readPoldata(opt2fn("-f2", NFILE, fnm), &pd2);
            }
            GMX_CATCH_ALL_AND_EXIT_WITH_FATAL_ERROR;
            (void) pd2.verifyCheckSum(stderr, poldataCheckSum(&pd2));
            
            if (strlen(missing) > 0)
            {
                copy_missing(&pd2, &pd, missing, false);
            }
            else if (strlen(replace) > 0)
            {
                copy_missing(&pd2, &pd, replace, true);
            }
            else if (strlen(implant) > 0)
            {
                implant_values(&pd2, &pd, implant);
            }
            else
            {
                compare_pd(&pd, &pd2, parameter);
            }
        }
        else if (strlen(analyze) > 0)
        {
            if (strlen(particle) > 0)
            {
                dumpPoldata(&pd, analyze, particle, opt2fn("-dump", NFILE, fnm));
            }
            else
            {
                analyzePoldata(&pd, analyze);
            }
        }
        else
        {
            modifyPoldata(&pd, parameter, particle,
                          opt2parg_bSet("-min", npargs, pa), pmin,
                          opt2parg_bSet("-val", npargs, pa), pval,
                          opt2parg_bSet("-max", npargs, pa), pmax,
                          opt2parg_bSet("-mut", npargs, pa), mutability,
                          opt2parg_bSet("-scale", npargs, pa), scale,
                          force, stretch);
        }
    }
    if (opt2bSet("-o", NFILE, fnm))
    {
        if (cr.isParallel())
        {
            if (cr.isMaster())
            {
                pd.Send(&cr, 1);
                std::string outfile(opt2fn("-o", NFILE, fnm));
                cr.send_str(1, &outfile);
            }
            else if (cr.rank() == 1)
            {
                pd.Receive(&cr, 0);
                std::string outfile;
                cr.recv_str(0, &outfile);
                alexandria::writePoldata(outfile, &pd, 0);
            }
        }
        else
        {
            std::string checkSum = poldataCheckSum(&pd);
            if (checkSum == pd.checkSum())
            {
                printf("No changes to poldata structure, not writing a new file.");
            }
            else
            {
                pd.updateTimeStamp();
                pd.setCheckSum(checkSum);
                writePoldata(opt2fn("-o", NFILE, fnm), &pd, 0);
            }
        }
    }
    return 0;
}

} // namespace alexandria


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

#include <cstdlib>

#include <algorithm>
#include <map>
#include <set>

#include "act/alexandria/alex_modules.h"
#include "act/basics/interactiontype.h"
#include "act/basics/mutability.h"
#include "act/forces/forcecomputer.h"
#include "act/forcefield/act_checksum.h"
#include "act/forcefield/combruleutil.h"
#include "act/forcefield/forcefield_parameter.h"
#include "act/forcefield/forcefield_parametername.h"
#include "act/forcefield/forcefield.h"
#include "act/forcefield/forcefield_xml.h"
#include "gromacs/commandline/pargs.h"
#include "gromacs/fileio/oenv.h"
#include "gromacs/topology/ifunc.h"
#include "gromacs/utility/arraysize.h"

namespace alexandria
{

static void setMinMaxMut(FILE *fp,
                         ForceFieldParameter *pp,
                         bool bSetMin, double pmin,
                         bool bSetVal, double pval,
                         bool bSetMax, double pmax,
                         bool bSetMut, const std::string &mutability,
                         bool bScale,  double scale,
                         bool stretch, const std::string &particleId,
                         bool bLimits, double factor)
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
        GMX_RELEASE_ASSERT(pmin <= pmax, "Please choose minimum to be at most the maximum");
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
        double range = std::fabs(pp->maximum()-pp->minimum());
        fprintf(stderr, "Trying to stretch parameter with unit %s range %g\n", pp->unit().c_str(), range);
        if (std::fabs(pp->value() - pp->minimum()) < 0.01*range && range > 0)
        {
            double newmin = pp->minimum()-0.2*range;
            if (newmin <= 0 && pp->nonNegative())
            {
                newmin = 0.8*pp->minimum();
            }
            if (pp->setMinimum(newmin))
            {
                if (fp)
                {
                    fprintf(fp, "Minimum stretched to %g for %s\n",
                            pp->minimum(), particleId.c_str());
                }
            }
            else if (fp)
            {
                fprintf(fp, "Could not change minimum from %g to %g\n", pp->minimum(), newmin);
            }
        }
        if (std::fabs(pp->value() - pp->maximum()) < 0.01*range && range > 0)
        {
            double newmax = pp->maximum()+0.2*range;
            pp->setMaximum(newmax);
            if (fp)
            {
                fprintf(fp, "Maximum stretched to %g for %s\n",
                        newmax, particleId.c_str());
            }
        }
    }
    if (bLimits)
    {
        double mm = pp->value()*factor;
        double mx = pp->value()/factor;
        pp->setMinimum(std::min(mx, mm));
        pp->setMaximum(std::max(mx, mm));
    }
}

static void modifyParticle(const std::string &paramType,
                           ParticleType      *particle,
                           bool bSetMin, double pmin,
                           bool bSetVal, double pval,
                           bool bSetMax, double pmax,
                           bool bSetMut, const std::string &mutability,
                           bool bScale,  double scale,
                           bool force,   bool stretch,
                           bool bLimit,  double factor)
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
                         stretch, particleId,
                         bLimit,  factor);
        }
    }
    else
    {
        printf("No interaction type corresponding to %s\n", paramType.c_str());
    }
}

static void modifyInteraction(ForceField *pd,
                              InteractionType itype,
                              const std::string &paramType,
                              const Identifier &pId,
                              bool bSetMin, double pmin,
                              bool bSetVal, double pval,
                              bool bSetMax, double pmax,
                              bool bSetMut, const std::string &mutability,
                              bool bScale,  double scale,
                              bool force,   bool stretch,
                              bool bLimit,  double factor)
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
                                     stretch, myId,
                                     bLimit,  factor);
                    }
                }
            }
        }
    }
}

static void modifyForceField(ForceField *pd,
                             const std::string &paramType,
                             const std::string &particle,
                             bool bSetMin, double pmin,
                             bool bSetVal, double pval,
                             bool bSetMax, double pmax,
                             bool bSetMut, const std::string &mutability,
                             bool bScale,  double scale,
                             bool force,   bool stretch,
                             bool bLimit,  double factor)
{
    if (!(bSetVal || bSetMin || bSetMax || bSetMut || stretch || bScale || bLimit))
    {
        printf("No parameter to change.\n");
        return;
    }
    if (paramType.empty())
    {
        printf("Empty parameter type. Nothing to modify.\n");
        return;
    }
    InteractionType itype;
    if (pd->typeToInteractionType(paramType, &itype))
    {
        if (pd->interactionPresent(itype))
        {
            // A true interaction
            auto fs = pd->findForces(itype);
            for (auto &ffp : fs->parametersConst())
            {
                auto id = ffp.first;
                if (particle.empty() || id.id() == particle)
                {
                    modifyInteraction(pd, itype, paramType, id,
                                      bSetMin, pmin,
                                      bSetVal, pval,
                                      bSetMax, pmax,
                                      bSetMut, mutability,
                                      bScale,  scale,
                                      force,   stretch,
                                      bLimit,  factor);
                }
            }
        }
        else
        {
            std::vector<ParticleType *> myParticles;
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
                for (auto &pt : *pd->particleTypes())
                {
                    myParticles.push_back(&pt.second);
                }
            }
            // A particle?
            for (auto &p : myParticles)
            {
                modifyParticle(paramType, p,
                               bSetMin, pmin,
                               bSetVal, pval,
                               bSetMax, pmax,
                               bSetMut, mutability,
                               bScale,  scale,
                               force,   stretch,
                               bLimit,  factor);
            }
        }
    }
}

static const std::set<InteractionType> &findInteractionMap(const std::string &analyze,
                                                           bool              *found)
{
    static std::set<InteractionType> bonds = {
        InteractionType::VDW,
        InteractionType::ELECTROSTATICS,
        InteractionType::BONDS,
        InteractionType::ANGLES,
        InteractionType::LINEAR_ANGLES,
        InteractionType::PROPER_DIHEDRALS,
        InteractionType::IMPROPER_DIHEDRALS };
    static std::set<InteractionType> other = {
        InteractionType::CONSTR,
        InteractionType::VSITE1,
        InteractionType::VSITE2,
        InteractionType::VSITE2FD,
        InteractionType::VSITE3FD,
        InteractionType::VSITE3FAD,
        InteractionType::VSITE3OUT,
        InteractionType::VSITE3OUTS };
    static std::set<InteractionType> eem = {
        InteractionType::POLARIZATION,
        InteractionType::ELECTROSTATICS,
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

static void analyzeForceField(ForceField *pd)
{
    int    mindata   = 1;
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

static void dumpForceField(ForceField        *pd,
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

static void plotInteractions(ForceField           *pd,
                             const std::string &analyze)
{
    bool found;
    auto myset = findInteractionMap(analyze, &found);
    if (!found)
    {
        return;
    }
    ForceComputer fc;
    for(auto &m : myset)
    {
        fc.plot(pd, m);
    }
}

static void copy_missing(const ForceField  *pdref,
                         ForceField        *pdout,
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

static void implant_values(const ForceField  *pdref,
                           ForceField        *pdout,
                           const std::string &analyze,
                           const std::string &particle,
                           const std::string &parameter,
                           bool               verbose)
{
    bool found;
    auto myset = findInteractionMap(analyze, &found);
    if (!found)
    {
        return;
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
            if (verbose)
            {
                printf("Cannot find interaction %s in %s\n",
                       interactionTypeToDescription(itype).c_str(), pdout->filename().c_str());
            }
            continue;
        }
        if (verbose)
        {
            printf("Will try to implant values for %s interactions from %s to %s\n",
                   interactionTypeToString(itype).c_str(),
                   pdref->filename().c_str(), pdout->filename().c_str());
        }
        auto fcout = pdout->findForces(itype);
        for(auto &fm : fc.second.parametersConst())
        {
            auto myid    = fm.first;
            auto myatoms = myid.atoms();
            // If we have a specified atom, it should be in the identifier, if not skip this one
            if (!particle.empty() && std::find(myatoms.begin(), myatoms.end(), particle) == myatoms.end())
            {
                continue;
            }
            if (fcout->parameterExists(myid))
            {
                for(auto &fmp : fm.second)
                {
                    auto pout = fcout->findParameterType(myid, fmp.first);
                    // If we have a specified parameter, it should be that one, if not skip it
                    if (!parameter.empty() && fmp.first != parameter)
                    {
                        continue;
                    }

                    if (fmp.second.ntrain() > 0)
                    {
                        printf("Parameter %s for particle %s will be replaced in %s\n",
                               fmp.first.c_str(), myid.id().c_str(),
                               interactionTypeToDescription(itype).c_str());
                        pout->copy(fmp.second);
                    }
                    else if (verbose)
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

static void compare_pd(ForceField *pd1,
                       ForceField *pd2,
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

static void copyDeToD0(ForceField *pd)
{
    auto fs = pd->findForces(InteractionType::BONDS);
    if (fs->potential() != Potential::MORSE_BONDS)
    {
        printf("Not using Morse in force field file %s\n", pd->filename().c_str());
        return;
    }
    auto ppp = fs->parameters();
    for (auto &p : (*ppp))
    {
        auto &param = p.second;
        if (param.find(morse_name[morseDE]) != param.end() && 
            param.find(morse_name[morseD0]) != param.end())
        {
            double De = param.find(morse_name[morseDE])->second.value();
            param.find(morse_name[morseD0])->second.setValue(-De);
        }
    }
}

static void addBondEnergy(ForceField *pd)
{
    auto fs = pd->findForces(InteractionType::BONDS);
    if (fs->potential() != Potential::HARMONIC_BONDS)
    {
        printf("Not using Bonds in force field file %s\n",
               pd->filename().c_str());
        return;
    }
    auto ppp = fs->parameters();
    for (auto &p : (*ppp))
    {
        auto &param = p.second;
        if (param.find(bond_name[bondLENGTH]) != param.end() && 
            param.find(bond_name[bondKB]) != param.end() &&
            param.find(bond_name[bondENERGY]) == param.end())
        {
            ForceFieldParameter be("kJ/mol", -200, 0, 1, -1000, -20, 
                                   Mutability::Bounded, false,
                                   false, true);
            param.insert({ bond_name[bondENERGY], be });
        }
    }
}

int edit_ff(int argc, char*argv[])
{
    std::vector<const char *> desc =
    {
        "edit_ff reads a force field file and can do a number of things.[PAR]",
        "If the flag [TT]-ana[TT] is set it will analyze the file and print",
        "if parameters are close to their bounds.[PAR]",
        "It may write a new file that can be compared to the input.",
        "If run in parallel, this utility will read on the first processor,",
        "then process the data and send it to the second processor. That one",
        "will then write the new file. This is for debugging force field",
        "communication.[PAR]",
        "Modifications of parameters can be made by specifying both",
        "input and output files and what parameters to change.",
        "If the value, the minimum",
        "or the maximum is to be changed, the actual value may be set to the",
        "new minmum or maximum if it falls outside the new bounds.",
    };
    CombRuleUtil crule;
    crule.addInfo(&desc);
    gmx_output_env_t                *oenv;
    t_filenm                         fnm[] = {
        { efXML, "-ff",   "aff_in" , ffREAD  },
        { efXML, "-ff2",  "aff_in2", ffOPTRD },
        { efXML, "-o",    "aff_out", ffOPTWR },
        { efDAT, "-dump", "params",  ffOPTWR }
    };

    static char *parameter  = (char *)"";
    static char *particle   = (char *)"";
    static char *mutability = (char *)"";
    real         pmin       = 0;
    real         pmax       = 0;
    real         pval       = 0;
    real         scale      = 1;
    real         limits     = 1;
    gmx_bool     force      = false;
    gmx_bool     stretch    = false;
    gmx_bool     plot       = false;
    gmx_bool     De2D0      = false;
    gmx_bool     bondenergy = false;
    gmx_bool     forceWrite = false;
    gmx_bool     bcast      = false;
    gmx_bool     bounds     = false;
    gmx_bool     verbose    = false;
    static char *missing    = (char *)"";
    static char *replace    = (char *)"";
    static char *implant    = (char *)"";
    static char *analyze    = (char *)"";
    std::vector<t_pargs> pa =
    {
        { "-p",      FALSE, etSTR,  {&parameter},
          "Type of parameter to change, e.g. zeta." },
        { "-a",      FALSE, etSTR,  {&particle},
          "Particle/interaction type to change, if not set all particles/interaction of the correct type will be changed. A quoted string may be passed to adjust multiple particle types at once, e.g. -a 'h1_z hc_z'" },
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
          "Will automatically stretch boundaries for individual parameters" },
        { "-limits",  FALSE, etREAL, {&limits},
          "Reset the limits for a parameter (class) to the current value of the parameter times this number (between 0 and 1) and one over the value. If you set e.g. -limits 0.8 the parameter min and max will be set to 0.8 respectively 1.25 times the present value." },
        { "-ana", FALSE, etSTR, {&analyze},
          "Analyze either the EEM, the BONDS or OTHER parameters in a simple manner" },
        { "-bounds", FALSE, etBOOL, {&bounds},
          "Check whether any parameter is hitting the wall, i.e. is at one of the boundaries." },
        { "-copy_missing", FALSE, etSTR, {&missing},
          "Copy either the EEM, the BONDS or OTHER parameters from file two [TT]-f2[tt] that are missing from file one [TT]-f[tt] to another [TT]-o[tt]." },
        { "-replace", FALSE, etSTR, {&replace},
          "Replace either the EEM, the BONDS or OTHER parameters in file one [TT]-f[ff] by those from file two [TT]-f2[tt] and store in another [TT]-o[tt]." },
        { "-implant", FALSE, etSTR, {&implant},
          "Implant (write over) either the EEM, the BONDS or OTHER parameters in file one [TT]-f[ff] by those from file two [TT]-f2[tt] and store in another [TT]-o[tt]. Only parameters with ntrain larger than zero will be copied. This can be used to merge training data from multiple runs." },
        { "-de2d0", FALSE, etBOOL, {&De2D0},
          "This is a hack to copy -De to D0 in the Morse potential" },
        { "-bondenergy", FALSE, etBOOL, {&bondenergy},
          "This is a hack to a bondenergy field to the BOND potential" },
        { "-plot",    FALSE, etBOOL, {&plot},
          "Plot many interactions as a function of distance or angle" },
        { "-bcast",   FALSE, etBOOL, {&bcast},
          "Use broadcast rather than send/receive to communicate force field" },
        { "-v", FALSE, etBOOL, {&verbose},
          "Print more stuff during processing" },
        { "-write",   FALSE, etBOOL, {&forceWrite},
          "Write out a force field file even if there were no changes" }
    };
    int NFILE  = asize(fnm);
    crule.addPargs(&pa);
    if (!parse_common_args(&argc, argv, 0, NFILE, fnm, pa.size(), pa.data(),
                           desc.size(), desc.data(), 0, nullptr, &oenv))
    {
        return 0;
    }
    CommunicationRecord cr;
    cr.init(cr.size());
    ForceField pd;
    if (cr.isMaster())
    {
        try 
        {
            alexandria::readForceField(opt2fn("-ff", NFILE, fnm), &pd);
        }
        GMX_CATCH_ALL_AND_EXIT_WITH_FATAL_ERROR;
        (void) pd.verifyCheckSum(stderr, forcefieldCheckSum(&pd));
        if (opt2bSet("-ff2", NFILE, fnm))
        {
            alexandria::ForceField pd2;
            try
            {
                alexandria::readForceField(opt2fn("-ff2", NFILE, fnm), &pd2);
            }
            GMX_CATCH_ALL_AND_EXIT_WITH_FATAL_ERROR;
            (void) pd2.verifyCheckSum(stderr, forcefieldCheckSum(&pd2));
            
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
                implant_values(&pd2, &pd, implant, particle, parameter, verbose);
            }
            else
            {
                compare_pd(&pd, &pd2, parameter);
            }
        }
        else if (bounds)
        {
            analyzeForceField(&pd);
        }
        else if (strlen(analyze) > 0)
        {
            auto dumpfn = opt2fn_null("-dump", NFILE, fnm);
            if (strlen(particle) > 0 && dumpfn)
            {
                dumpForceField(&pd, analyze, particle, dumpfn);
            }
            else if (plot)
            {
                plotInteractions(&pd, analyze);
            }
        }
        else if (De2D0)
        {
            copyDeToD0(&pd);
        }
        else if (bondenergy)
        {
            addBondEnergy(&pd);
        }
        else
        {
            bool bLimits = opt2parg_bSet("-limits", pa.size(), pa.data());
            if (limits <= 0)
            {
                fprintf(stderr, "%g is an inappropriate value for the limits option\n", limits);
                limits  = 1;
                bLimits = false;
            }
            else if (limits > 1)
            {
                limits = 1.0/limits;
            }
            modifyForceField(&pd, parameter, particle,
                             opt2parg_bSet("-min", pa.size(), pa.data()), pmin,
                             opt2parg_bSet("-val", pa.size(), pa.data()), pval,
                             opt2parg_bSet("-max", pa.size(), pa.data()), pmax,
                             opt2parg_bSet("-mut", pa.size(), pa.data()), mutability,
                             opt2parg_bSet("-scale", pa.size(), pa.data()), scale,
                             force, stretch, bLimits, limits);
        }
    }
    // Fetch new combination rules if necessary
    std::map<InteractionType, ForceFieldParameterList *> its =
        {
            { InteractionType::VDW,            nullptr },
            { InteractionType::CHARGETRANSFER, nullptr }
        };

    for(auto &fst : its)
    {
        if (pd.interactionPresent(fst.first))
        {
            fst.second = pd.findForces(fst.first);
        }
    }
    int nRuleChanged = crule.extract(its[InteractionType::VDW],
                                     its[InteractionType::CHARGETRANSFER]);
    if (nRuleChanged > 0)
    {
        printf("Inserted %d new style combination rules from command line.\n",
               nRuleChanged);
    }
    else
    {
        nRuleChanged = crule.convert(its[InteractionType::VDW]);
        printf("Converted old style comb rule to %d new style combination rules.\n", nRuleChanged);
    }
    its[InteractionType::VDW]->removeOption("combination_rule");

    if (opt2bSet("-o", NFILE, fnm))
    {
        if (cr.isParallel())
        {
            int root = 0;
            CommunicationStatus cs;
            if (cr.isMaster())
            {
            
                std::string outfile(opt2fn("-o", NFILE, fnm));
                printf("Will send force field to my helpers to write %s\n", outfile.c_str());
                if (bcast)
                {
                    cs = pd.BroadCast(&cr, root, cr.comm_world());
                    cr.bcast(&outfile, cr.comm_world());
                }
                else
                {
                    for(int dest = 1; dest < cr.size(); dest++)
                    {
                        cs = pd.Send(&cr, dest);
                        cr.send(dest, outfile);
                    }
                }
            }
            else
            {
                std::string outfile;
                if (bcast)
                {
                    cs = pd.BroadCast(&cr, root, cr.comm_world());
                    cr.bcast(&outfile, cr.comm_world());
                }
                else
                {
                    cs = pd.Receive(&cr, 0);
                    cr.recv(0, &outfile);
                }
                if (CommunicationStatus::OK == cs && cr.rank() == 2)
                {
                    alexandria::writeForceField(outfile, &pd, 0);
                }
            }
        }
        else
        {
            std::string checkSum = forcefieldCheckSum(&pd);
            if (checkSum == pd.checkSum() && !forceWrite && nRuleChanged == 0)
            {
                printf("No changes to forcefield structure, not writing a new file.");
            }
            else
            {
                pd.updateTimeStamp();
                pd.setCheckSum(checkSum);
                writeForceField(opt2fn("-o", NFILE, fnm), &pd, 0);
            }
        }
    }
    return 0;
}

} // namespace alexandria


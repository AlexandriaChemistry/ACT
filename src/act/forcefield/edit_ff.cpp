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

#include <algorithm>
#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <map>
#include <set>
#include <string>
#include <utility>
#include <vector>

#include "act/alexandria/alex_modules.h"
#include "act/basics/identifier.h"
#include "act/basics/interactiontype.h"
#include "act/basics/msg_handler.h"
#include "act/basics/mutability.h"
#include "act/forcefield/act_checksum.h"
#include "act/forcefield/combruleutil.h"
#include "act/forcefield/forcefield.h"
#include "act/forcefield/forcefield_json.h"
#include "act/forcefield/forcefield_parameter.h"
#include "act/forcefield/forcefield_parameterlist.h"
#include "act/forcefield/forcefield_validation.h"
#include "act/forcefield/forcefield_xml.h"
#include "act/forcefield/generate_dependent.h"
#include "act/forcefield/particletype.h"
#include "act/forces/forcecomputer.h"
#include "act/utility/jsontree.h"
#include "edit_ff.h"
#include "gromacs/commandline/filenm.h"
#include "gromacs/commandline/pargs.h"
#include "gromacs/fileio/filetypes.h"
#include "gromacs/fileio/oenv.h"
#include "gromacs/utility/basedefinitions.h"
#include "gromacs/utility/exceptions.h"
#include "gromacs/utility/fatalerror.h"
#include "gromacs/utility/gmxassert.h"
#include "gromacs/utility/real.h"
#include "gromacs/utility/stringutil.h"

namespace alexandria
{

//! \brief Set parameter limits
static void setMinMaxMut(FILE *fp,
                         ForceFieldParameter *pp,
                         bool bSetMin, double pmin,
                         bool bSetVal, double pval,
                         bool bSetMax, double pmax,
                         bool bSetMut, const std::string &mutability,
                         bool bScale,  double scale,
                         bool force,   bool stretch,
                         const std::string &particleId,
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
        double clampedVal = std::max(pmin, pp->value());
        if (force)
        {
            pp->forceSetValue(clampedVal);
        }
        else if (pp->isMutable())
        {
            pp->setValue(clampedVal);
        }
        if (fp)
        {
            fprintf(fp, "Minimum set to %g for %s\n", pmin, particleId.c_str());
        }
    }
    if (bSetMax)
    {
        // As for the minimum, the bounds are consistent before the value
        // is moved into them
        pp->setMaximum(pmax);
        pp->setMinimum(std::min(pmax, pp->minimum()));
        double clampedVal = std::min(pmax, pp->value());
        if (force)
        {
            pp->forceSetValue(clampedVal);
        }
        else if (pp->isMutable())
        {
            pp->setValue(clampedVal);
        }
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
        if (force)
        {
            pp->forceSetValue(pval);
        }
        else
        {
            pp->setValue(pval);
        }
        if (fp)
        {
            fprintf(fp, "Value set to %g %s\n", pval, particleId.c_str());
        }
    }
    if (bSetMut)
    {
        Mutability mut;
        if (!nameToMutability(mutability, &mut))
        {
            GMX_THROW(gmx::InvalidInputError(gmx::formatString("Unknown mutability '%s'", mutability.c_str()).c_str()));
        }
        pp->setMutability(mut);
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

//! \brief Change parameters for one particle
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
        // Allow min/max/stretch/limits changes even for Fixed parameters when
        // no value or scale change is requested, since setMinimum/setMaximum
        // are independent of mutability.
        bool hasValueOrScaleChange = bSetVal || bScale;
        bool hasLimitChange        = bSetMin || bSetMax || stretch || bLimit;
        if (force || ff->isMutable() || (hasLimitChange && !hasValueOrScaleChange))
        {
            const std::string particleId = 
                      gmx::formatString("%s - %s",
                                        particle->id().id().c_str(),
                                        paramType.c_str());
            setMinMaxMut(debug, ff,
                         bSetMin, pmin, 
                         bSetVal, pval, bSetMax, pmax,
                         bSetMut, mutability, 
                         bScale,  scale,
                         force,   stretch, particleId,
                         bLimit,  factor);
        }
    }
    else
    {
        printf("No interaction type corresponding to %s\n", paramType.c_str());
    }
}

//! \brief Modify parameters for one interaction
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
                        const std::string myId = 
                                  gmx::formatString("%s - %s - %s",
                                                    pId.id().c_str(),
                                                    paramType.c_str(),
                                                    interactionTypeToString(itype).c_str());
                        setMinMaxMut(debug, &pp.second,
                                     bSetMin, pmin, 
                                     bSetVal, pval, bSetMax, pmax,
                                     bSetMut, mutability, 
                                     bScale,  scale,
                                     force,   stretch, myId,
                                     bLimit,  factor);
                    }
                }
            }
        }
    }
}

//! \brief Remove particles from a FF
static void removeParticles(ForceField *pd,
                            const std::string &particle)
{
    for (const auto &s : gmx::splitString(particle))
    {
        // First remove all the forces involving s
        auto forces = pd->forces();
        for (auto f = forces->begin(); f != forces->end(); ++f)
        {
            auto params = f->second.parameters();
            for (auto fp = params->begin(); fp != params->end(); )
            {
                const auto myatoms = fp->first.atoms();
                bool found = false;
                for (const auto &mya : myatoms)
                {
                    found = found || mya == s;
                }
                if (found)
                {
                    fp = params->erase(fp);
                }
                else
                {
                    ++fp;
                }
            }
        }
        // Then remove s itself
        auto ptypes = pd->particleTypes();
        for (auto p = ptypes->begin(); p != ptypes->end(); )
        {
            if (s == p->first.id())
            {
                p = ptypes->erase(p);
            }
            else
            {
                ++p;
            }
        }
    }
}

/*! \brief Set or clear the attribute nonnegative of parameters
 * Parameters of interactions as for modifyForceField. Reading a force field
 * makes all parameters of particle types non-negative, except the charge,
 * so for those only setting it lasts.
 * \param[inout] pd          The force field
 * \param[in]    paramType   The type of parameter, such as bondenergy
 * \param[in]    particle    The identifiers, all if empty
 * \param[in]    nonNegative Whether the parameters may not be negative
 * \return the number of parameters changed
 */
static int setNonNegative(ForceField        *pd,
                          const std::string &paramType,
                          const std::string &particle,
                          bool               nonNegative)
{
    if (paramType.empty())
    {
        printf("Empty parameter type, use -p. Nothing to modify.\n");
        return 0;
    }
    auto ids     = gmx::splitString(particle);
    const auto wanted  = [&ids](const std::string &id)
          {
              return ids.empty() || std::find(ids.begin(), ids.end(), id) != ids.end();
          };
    int  changed = 0;
    const auto change  = [nonNegative, &changed](ForceFieldParameter *p)
          {
              if (p->nonNegative() == nonNegative)
              {
                  return;
              }
              if (nonNegative)
              {
                  p->setNonNegative();
              }
              else
              {
                  p->clearNonNegative();
              }
              changed += 1;
          };
    InteractionType itype;
    if (pd->typeToInteractionType(paramType, &itype) && pd->interactionPresent(itype))
    {
        auto fs = pd->findForces(itype);
        for (auto &entry : *fs->parameters())
        {
            const auto param = entry.second.find(paramType);
            if (wanted(entry.first.id()) && param != entry.second.end())
            {
                change(&param->second);
            }
        }
    }
    else
    {
        for (auto &pt : *pd->particleTypes())
        {
            const auto param = pt.second.parameters()->find(paramType);
            if (wanted(pt.first.id()) && param != pt.second.parameters()->end())
            {
                change(&param->second);
            }
        }
    }
    printf("Made %d parameter(s) of type %s %s.\n", changed, paramType.c_str(),
           nonNegative ? "non-negative" : "allowed to be negative");
    return changed;
}

//! \brief Modify everything that the user has requested
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
            // A true interaction, of which one or more identifiers may be given
            auto ids = gmx::splitString(particle);
            auto fs  = pd->findForces(itype);
            for (auto &ffp : fs->parametersConst())
            {
                const auto id = ffp.first;
                if (ids.empty() || std::find(ids.begin(), ids.end(), id.id()) != ids.end())
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

//! \return subset of interactions corresponding to user input
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
        InteractionType::VSITE3,
        InteractionType::VSITE3S,
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

//! \brief Investigate FF and print interesting observations
static void analyzeForceField(ForceField *pd)
{
    unsigned int mindata   = 1;
    double tolerance = 0.001;
    for(auto &fc : pd->forcesConst())
    {
        auto itype = fc.first;
        for(auto &fm : fc.second.parametersConst())
        {
            const auto myid = fm.first;
            for(auto &param : fm.second)
            {
                const auto ptype = param.first;
                const auto ppp   = param.second;
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

//! \brief Dump select information from a force field to a file
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
            const auto myid = fm.first;
            if (std::find(particles.begin(), particles.end(), myid.id()) != particles.end())
            {
                for(auto &param : fm.second)
                {
                    const auto ptype = param.first;
                    const auto ppp   = param.second;
                    fprintf(fp, "%s|%s|%g\n",
                            ptype.c_str(), myid.id().c_str(), ppp.value());
                    nparm++;
                }
            }
        }
    }
    // The widths of the charge distributions are those of the particle types
    if (myset.find(InteractionType::ELECTROSTATICS) != myset.end())
    {
        for(const auto &pt : pd->particleTypesConst())
        {
            if (pt.second.hasParameter(zetaName) &&
                std::find(particles.begin(), particles.end(), pt.first.id()) != particles.end())
            {
                fprintf(fp, "%s|%s|%g\n", zetaName, pt.first.id().c_str(),
                        pt.second.parameterConst(zetaName).value());
                nparm++;
            }
        }
    }
    fclose(fp);
    printf("Stored %d parameters in %s\n", nparm, filenm.c_str());
}

//! \brief Plot potential curves to many xvg files
static void plotInteractions(MsgHandler        *msghandler,
                             ForceField        *pd,
                             const std::string &analyze)
{
    bool found;
    const auto myset = findInteractionMap(analyze, &found);
    if (!found)
    {
        return;
    }
    const ForceComputer fc(1e-3, 100, 1.0);
    for(auto &m : myset)
    {
        fc.plot(msghandler, pd, m);
    }
    msghandler->msg(ACTStatus::Info,
                    gmx::formatString("ForceComputer: %zu evaluations, %zu shell iterations and %zu shell convergence failures",
                                      fc.numEvaluations(), fc.numShellIterations(), fc.numShellConvergenceFailed()));
}

//! \brief Copy missing information from a reference FF to an output FF
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
            const auto myid = fm.first;
            if (!fcout->parameterExists(myid))
            {
                printf("Parameter %s missing from %s for %s\n",
                       myid.id().c_str(), pdout->filename().c_str(),
                       interactionTypeToDescription(itype).c_str());
                for(auto &param : fm.second)
                {
                    const auto ptype = param.first;
                    const auto ppp   = param.second;
                    fcout->addParameter(myid, ptype, ppp);
                }
            }
        }
    }
}

//! \brief Implant specific values to a FF from a reference
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
            const auto myid    = fm.first;
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

//! \brief Compare to force field for a specific parameter
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
    const auto& f1     = pd1->findForcesConst(itype1);
    const auto& f2     = pd2->findForcesConst(itype2);
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
            auto fp1 = f1.findParameterMapConst(id1.first);
            auto fp2 = f2.findParameterMapConst(id1.first);
            const auto p1  = fp1[ptype];
            const auto p2  = fp2[ptype];
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

//! The result of validating a force field
struct Validation
{
    //! The number of problems, a wrong checksum included
    size_t nissues = 0;
    //! The number of problems that writing the force field removes
    size_t nfixed  = 0;
};

/*! \brief Print the problems with a force field
 * These are a checksum that does not match the contents and the problems
 * with the parameters, see validateForceField.
 * \param[in] pd              The force field
 * \param[in] fileName        Its file
 * \param[in] contentCheckSum The checksum of the contents as read
 * \param[in] writing         Whether the force field will be written
 * \return the number of problems and of those that writing removes
 */
static Validation printValidation(const ForceField  *pd,
                                  const std::string &fileName,
                                  const std::string &contentCheckSum,
                                  bool               writing)
{
    Validation result;
    if (pd->checkSum() == contentCheckSum)
    {
        printf("Checksum of %s: %s, correct\n", fileName.c_str(), contentCheckSum.c_str());
    }
    else
    {
        // Writing gives the file the checksum of its contents
        result.nissues += 1;
        result.nfixed  += 1;
        printf("Checksum of %s: the file has '%s', the contents as read give %s.\n"
               "This is corrected by writing the force field with -o, or with edit_ff -sign.\n",
               fileName.c_str(), pd->checkSum().c_str(), contentCheckSum.c_str());
    }
    const auto issues = validateForceField(*pd);
    const std::vector<std::pair<ForceFieldIssue::Kind, const char *>> kinds =
    {
        { ForceFieldIssue::Kind::ResetToBounds,       "Values outside their bounds, set to the nearest bound" },
        { ForceFieldIssue::Kind::Duplicate,           "Parameters defined more than once, the first definition is used" },
        { ForceFieldIssue::Kind::OutOfBounds,         "Values outside their bounds that may not be changed" },
        { ForceFieldIssue::Kind::MinimumAboveMaximum, "Minimum larger than maximum" },
        { ForceFieldIssue::Kind::Negative,            "Negative values of parameters that may not be negative" }
    };
    size_t nfixed = 0;
    printf("Validation of the parameters of %s: %zu problem(s)\n", fileName.c_str(), issues.size());
    for (const auto &kind : kinds)
    {
        size_t n = 0;
        for (const auto &issue : issues)
        {
            if (issue.kind != kind.first)
            {
                continue;
            }
            if (n++ == 0)
            {
                printf("\n%s:\n", kind.second);
            }
            printf("  %-28s %-24s %-12s %s\n", issue.where.c_str(), issue.identifier.c_str(),
                   issue.parameter.c_str(), issue.text.c_str());
            nfixed += issue.repaired ? 1 : 0;
        }
    }
    if (nfixed > 0)
    {
        printf("\n%zu problem(s) %s by writing the force field with -o.\n", nfixed,
               writing ? "are corrected" : "would be corrected");
    }
    if (issues.size() > nfixed)
    {
        printf("%zu problem(s) have to be corrected by hand, for instance with edit_ff -p -a -min -max -val -nonneg.\n",
               issues.size() - nfixed);
    }
    result.nissues += issues.size();
    result.nfixed  += nfixed;
    return result;
}

/*! \brief Tool to manipulate force field files
 * \param[in] argc Number of arguments on the cmd line
 * \param[in] argv The actual arguments
 * \return 0 if ok, 1 otherwise
 */
void EditFFOptions::check() const
{
    if (forceFieldFile.empty())
    {
        GMX_THROW(gmx::InvalidInputError("There is no force field file (-ff)"));
    }
    const std::vector<std::pair<const char *, std::optional<double>>> numbers = {
        { "-min", minimum }, { "-val", value }, { "-max", maximum }, { "-scale", scale }, { "-limits", limits }
    };
    for (const auto &n : numbers)
    {
        if (n.second && !std::isfinite(*n.second))
        {
            GMX_THROW(gmx::InvalidInputError(gmx::formatString("The value of %s is not a number", n.first)));
        }
    }
    if (mutability)
    {
        Mutability mut;
        if (!nameToMutability(*mutability, &mut))
        {
            GMX_THROW(gmx::InvalidInputError(gmx::formatString("Unknown mutability '%s' of -mut", mutability->c_str())));
        }
    }
    // Signing writes the force field as it was read, nothing else
    if (sign)
    {
        const std::vector<std::pair<const char *, bool>> others = {
            { "-ff2", !secondForceFieldFile.empty() }, { "-min", minimum.has_value() },
            { "-val", value.has_value() }, { "-max", maximum.has_value() },
            { "-scale", scale.has_value() }, { "-mut", mutability.has_value() },
            { "-nonneg", nonNegative.has_value() }, { "-limits", limits.has_value() },
            { "-del", remove }, { "-stretch", stretch }, { "-bounds", bounds },
            { "-validate", validate }, { "-ana", !analyze.empty() }
        };
        for (const auto &o : others)
        {
            if (o.second)
            {
                GMX_THROW(gmx::InvalidInputError(gmx::formatString("-sign cannot be combined with %s", o.first)));
            }
        }
    }
    // These need the second force field, without it they would do nothing
    const std::vector<std::pair<std::string, std::string>> needSecond = {
        { "-copy_missing", copyMissing }, { "-replace", replace }, { "-implant", implant }
    };
    for (const auto &kind : needSecond)
    {
        if (!kind.second.empty() && secondForceFieldFile.empty())
        {
            GMX_THROW(gmx::InvalidInputError(gmx::formatString("%s needs a second force field (-ff2)", kind.first.c_str())));
        }
    }
}

int runEditFF(MsgHandler *msghandler, const EditFFOptions &options)
{
    int status = 0;
    options.check();
    ForceField pd;
    // A force field that is repaired with parameters of a second one need
    // not be consistent before, but it has to be afterwards.
    const bool repair = (!options.secondForceFieldFile.empty() &&
                         (!options.copyMissing.empty() || !options.replace.empty() || !options.implant.empty()));
    alexandria::readForceField(options.forceFieldFile, &pd, msghandler, true, !repair);
    // The checksum of the force field as it was read, which differs from
    // that in the file if it was converted, see ForceField::convertCoulomb
    const std::string readCheckSum = forcefieldCheckSum(&pd);
    const bool checkSumOk = pd.verifyCheckSum(msghandler, readCheckSum);
    bool forceWrite = options.forceWrite;
    if (!options.secondForceFieldFile.empty())
    {
        alexandria::ForceField pd2;
        // The second force field is only a source of parameters
        alexandria::readForceField(options.secondForceFieldFile, &pd2, msghandler, true, false);
        (void) pd2.verifyCheckSum(msghandler, forcefieldCheckSum(&pd2));

        if (repair)
        {
            if (!options.copyMissing.empty())
            {
                copy_missing(&pd2, &pd, options.copyMissing, false);
            }
            else if (!options.replace.empty())
            {
                copy_missing(&pd2, &pd, options.replace, true);
            }
            else
            {
                implant_values(&pd2, &pd, options.implant, options.particles, options.parameter, false);
            }
            pd.checkConsistency(msghandler);
        }
        else
        {
            compare_pd(&pd, &pd2, options.parameter);
        }
    }
    else if (options.bounds)
    {
        analyzeForceField(&pd);
    }
    else if (options.validate)
    {
        const auto validation = printValidation(&pd, options.forceFieldFile, readCheckSum,
                                                !options.outputFile.empty());
        // The corrections were made when reading, write them
        forceWrite = forceWrite || validation.nfixed > 0;
        // Also if writing corrects the problems: they are in the file that was read
        status     = validation.nissues > 0 ? 1 : 0;
    }
    else if (!options.analyze.empty())
    {
        if (!options.particles.empty() && !options.dumpFile.empty())
        {
            dumpForceField(&pd, options.analyze, options.particles, options.dumpFile);
        }
        else if (options.plot)
        {
            generateDependentParameter(msghandler, &pd);
            plotInteractions(msghandler, &pd, options.analyze);
        }
    }
    else if (options.remove)
    {
        removeParticles(&pd, options.particles);
    }
    else if (options.sign)
    {
        // Nothing to change, only to write below
    }
    else
    {
        bool   bLimits = options.limits.has_value();
        double limits  = options.limits.value_or(1);
        if (bLimits && limits <= 0)
        {
            fprintf(stderr, "%g is an inappropriate value for the limits option\n", limits);
            limits  = 1;
            bLimits = false;
        }
        else if (limits > 1)
        {
            limits = 1.0/limits;
        }
        const bool bNonNeg = options.nonNegative.has_value();
        if (bNonNeg)
        {
            (void) setNonNegative(&pd, options.parameter, options.particles, *options.nonNegative);
        }
        const bool other = (options.minimum || options.value || options.maximum || options.mutability ||
                            options.scale || options.stretch || bLimits);
        if (other || !bNonNeg)
        {
            modifyForceField(&pd, options.parameter, options.particles,
                             options.minimum.has_value(), options.minimum.value_or(0),
                             options.value.has_value(), options.value.value_or(0),
                             options.maximum.has_value(), options.maximum.value_or(0),
                             options.mutability.has_value(), options.mutability.value_or(""),
                             options.scale.has_value(), options.scale.value_or(1),
                             options.force, options.stretch, bLimits, limits);
        }
    }
    // Fetch new combination rules if necessary
    std::map<InteractionType, ForceFieldParameterList *> its =
        {
            { InteractionType::VDW,                  nullptr },
            { InteractionType::VDWCORRECTION ,       nullptr },
            { InteractionType::INDUCTIONCORRECTION , nullptr }
        };

    for(auto &fst : its)
    {
        if (pd.interactionPresent(fst.first))
        {
            fst.second = pd.findForces(fst.first);
        }
    }
    CombRuleUtil crule;
    crule.setRules(options.combinationRules);
    int nRuleChanged = crule.extract(&pd);
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
    if (its[InteractionType::VDW])
    {
        its[InteractionType::VDW]->removeOption("combination_rule");
    }
    // Whether the force field differs from the one that was read
    const std::string checkSum = forcefieldCheckSum(&pd);
    const bool        changed  = (checkSum != readCheckSum || nRuleChanged > 0);
    bool              written  = false;
    // Signing without -o writes over the file that was read
    const std::string outputFile = (options.sign && options.outputFile.empty()) ? pd.filename()
                                                                               : options.outputFile;
    if (!outputFile.empty())
    {
        if (!changed && !forceWrite && !options.sign)
        {
            printf("No changes to forcefield structure, not writing a new file.");
        }
        else
        {
            pd.updateTimeStamp();
            pd.setCheckSum(checkSum);
            writeForceField(outputFile, &pd, false);
            written = true;
            if (options.sign)
            {
                printf("Wrote %s with checksum %s and time stamp %s.\n", outputFile.c_str(),
                       pd.checkSum().c_str(), pd.timeStamp().c_str());
            }
        }
    }
    // After writing, such that the checksum is that of the new file
    if (!options.jsonFile.empty())
    {
        // The force field as it is now, for programs that show or edit it,
        // and whether what was asked for changed anything
        JsonTree jtree("edit_ff");
        jtree.addBool("changed", changed);
        jtree.addBool("written", written);
        // Whether the checksum in the file is that of the contents as read
        jtree.addBool("checksum_ok", checkSumOk);
        jtree.addObject(forceFieldToJson(pd));
        jtree.addObject(forceFieldIssuesToJson(validateForceField(pd)));
        jtree.write(options.jsonFile, true);
    }
    return status;
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
        "new minimum or maximum if it falls outside the new bounds.[PAR]",
        "With [TT]-json[tt] the force field is written in json format to the",
        "given file, edit_ff.json by default, as it is after the changes that",
        "were asked for: the particle types and the interactions with their",
        "options, combination rules and parameters (value, uncertainty,",
        "minimum, maximum, mutability), the problems that [TT]-validate[tt]",
        "looks for, whether the force field differs from the one that was",
        "read ([TT]changed[tt]) and whether a file was written with [TT]-o[tt]",
        "([TT]written[tt]). A force field without changes is not written.[PAR]",
        "With [TT]-validate[tt] the parameters are checked: values have to be",
        "within their bounds, the minimum must not be larger than the maximum,",
        "parameters that may not be negative must not be, and each parameter",
        "may be defined only once. Reading a force field sets values of trained",
        "parameters outside their bounds to the nearest bound and uses the first",
        "definition of a parameter that is defined more than once, so a force",
        "field written with [TT]-o[tt] no longer has these problems. The others",
        "are listed to be corrected by hand. The checksum in the file is",
        "compared with that of the contents as well; writing the force field",
        "or [TT]-sign[tt] corrects it. If there is any problem, the exit status",
        "is 1, also when the force field written with [TT]-o[tt] no longer has it.[PAR]",
        "With [TT]-sign[tt] the force field is written with a new checksum and",
        "time stamp, such that programs read it without the warning that the",
        "checksum does not match, for instance after editing the file by hand.",
        "The file is written to [TT]-o[tt] or, without [TT]-o[tt], over the file",
        "that was read. Older files are converted to the present format first.",
        "Nothing else is done, so [TT]-sign[tt] cannot be combined with options",
        "that change or analyze parameters."
    };
    CombRuleUtil crule;
    crule.addInfo(&desc);
    gmx_output_env_t                *oenv;
    std::vector<t_filenm> fnm = {
        { efXML, "-ff",   "aff_in" , ffREAD  },
        { efXML, "-ff2",  "aff_in2", ffOPTRD },
        { efXML, "-o",    "aff_out", ffOPTWR },
        { efDAT, "-dump", "params",  ffOPTWR },
        { efJSON, "-json", "edit_ff", ffOPTWR }
    };

    char        *parameter  = (char *)"";
    char        *particle   = (char *)"";
    char        *mutability = (char *)"";
    real         pmin       = 0;
    real         pmax       = 0;
    real         pval       = 0;
    real         scale      = 1;
    real         limits     = 1;
    gmx_bool     force      = false;
    gmx_bool     stretch    = false;
    gmx_bool     plot       = false;
    gmx_bool     forceWrite = false;
    gmx_bool     bounds     = false;
    gmx_bool     validate   = false;
    gmx_bool     sign       = false;
    const char  *nonneg[]   = { nullptr, "unchanged", "yes", "no", nullptr };
    gmx_bool     remove     = false;
    char        *missing    = (char *)"";
    char        *replace    = (char *)"";
    char        *implant    = (char *)"";
    char        *analyze    = (char *)"";
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
        { "-nonneg", FALSE, etENUM, {nonneg},
          "Set (yes) or clear (no) the attribute nonnegative of the parameters given with -p and -a, which prevents setting a negative minimum. Parameters of particle types other than the charge are always non-negative." },
        { "-mut",    FALSE, etSTR,  {&mutability},
          "Set the mutability for the given parameters to the value. The following values are supported: Fixed, Bounded" },
        { "-del",    FALSE, etBOOL,  {&remove},
          "Together with the -a flag, will remove all interactions containing this particle" },
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
        { "-validate", FALSE, etBOOL, {&validate},
          "Check the checksum and that the parameters are within their bounds and defined only once, and write the corrected force field with [TT]-o[tt], see above. The exit status is 1 if there is a problem." },
        { "-copy_missing", FALSE, etSTR, {&missing},
          "Copy either the EEM, the BONDS or OTHER parameters from file two [TT]-f2[tt] that are missing from file one [TT]-f[tt] to another [TT]-o[tt]." },
        { "-replace", FALSE, etSTR, {&replace},
          "Replace either the EEM, the BONDS or OTHER parameters in file one [TT]-f[ff] by those from file two [TT]-f2[tt] and store in another [TT]-o[tt]." },
        { "-implant", FALSE, etSTR, {&implant},
          "Implant (write over) either the EEM, the BONDS or OTHER parameters in file one [TT]-f[ff] by those from file two [TT]-f2[tt] and store in another [TT]-o[tt]. Only parameters with ntrain larger than zero will be copied. This can be used to merge training data from multiple runs." },
        { "-plot",    FALSE, etBOOL, {&plot},
          "Plot many interactions as a function of distance or angle" },
        { "-write",   FALSE, etBOOL, {&forceWrite},
          "Write out a force field file even if there were no changes" },
        { "-sign",    FALSE, etBOOL, {&sign},
          "Write the force field with a new checksum and time stamp, to [TT]-o[tt] or over the input file, see above." }
    };
    crule.addPargs(&pa);
    MsgHandler msghandler;
    msghandler.addOptions(&pa, &fnm, "edit_ff");
    if (!parse_common_args(&argc, argv, 0, fnm.size(), fnm.data(), pa.size(), pa.data(),
                           desc.size(), desc.data(), 0, nullptr, &oenv))
    {
        return 0;
    }

    CommunicationRecord cr;
    cr.init(cr.size());
    msghandler.optionsFinished(fnm, &cr);
    // Only the master edits the force field
    if (!cr.isMaster())
    {
        return 0;
    }

    try
    {
        EditFFOptions options;
        options.forceFieldFile = opt2fn("-ff", fnm.size(), fnm.data());
        const auto setFile     = [&fnm](const char *option)
        {
            return opt2bSet(option, fnm.size(), fnm.data()) ? std::string(opt2fn(option, fnm.size(), fnm.data()))
                                                            : std::string();
        };
        options.secondForceFieldFile = setFile("-ff2");
        options.outputFile           = setFile("-o");
        const char *dump             = opt2fn_null("-dump", fnm.size(), fnm.data());
        options.dumpFile             = dump ? dump : "";
        const char *json             = opt2fn_null("-json", fnm.size(), fnm.data());
        options.jsonFile             = json ? json : "";
        options.parameter            = parameter;
        options.particles            = particle;
        const auto given             = [&pa](const char *option)
        {
            return opt2parg_bSet(option, pa.size(), pa.data());
        };
        if (given("-min"))
        {
            options.minimum = pmin;
        }
        if (given("-val"))
        {
            options.value = pval;
        }
        if (given("-max"))
        {
            options.maximum = pmax;
        }
        if (given("-scale"))
        {
            options.scale = scale;
        }
        if (given("-mut"))
        {
            options.mutability = mutability;
        }
        if (given("-limits"))
        {
            options.limits = limits;
        }
        if (std::string(nonneg[0]) != "unchanged")
        {
            options.nonNegative = std::string(nonneg[0]) == "yes";
        }
        options.remove           = remove;
        options.force            = force;
        options.stretch          = stretch;
        options.bounds           = bounds;
        options.validate         = validate;
        options.plot             = plot;
        options.forceWrite       = forceWrite;
        options.sign             = sign;
        options.analyze          = analyze;
        options.copyMissing      = missing;
        options.replace          = replace;
        options.implant          = implant;
        options.combinationRules = crule.rules();
        return runEditFF(&msghandler, options);
    }
    catch (const gmx::InvalidInputError &ex)
    {
        fprintf(stderr, "%s\n", ex.what());
        return 1;
    }
}

} // namespace alexandria

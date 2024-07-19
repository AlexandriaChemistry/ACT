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

#include "particletype.h"

#include <cstdlib>

#include <list>
#include <map>

#include "gromacs/utility/exceptions.h"
#include "act/utility/communicationrecord.h"
#include "act/basics/identifier.h"

namespace alexandria
{

//! Map to convert a std::string to an InteractionType
static std::map<std::string, InteractionType> stringToItype =
    {
        { "acmtype",     InteractionType::ELECTRONEGATIVITYEQUALIZATION },
        { "zetatype",    InteractionType::ELECTROSTATICS },
        { "poltype",     InteractionType::POLARIZATION },
        { "bondtype",    InteractionType::BONDS },
        { "vdwtype",     InteractionType::VDW },
        { "vdwcorrtype", InteractionType::VDWCORRECTION },
        { "qttype",      InteractionType::CHARGETRANSFER }
    };
    
//! List of potential options for a particle
static std::list<std::string> particleOptions = 
    { "element", "atomnumber", "row" };

const std::string &ParticleType::optionValue(const std::string &type) const
{
    auto pm = option_.find(type);
    if (pm == option_.end())
    {
        GMX_THROW(gmx::InternalError(gmx::formatString("No such option type %s", type.c_str()).c_str()));
    }
    return pm->second;
}

void ParticleType::setOption(const std::string &key,
                             const std::string &value)
{
    if (stringToItype.find(key) != stringToItype.end() ||
        std::find(particleOptions.begin(), particleOptions.end(), key) != particleOptions.end())
    {
        option_.insert({key, value});
    }
    else
    {
        GMX_THROW(gmx::InternalError(gmx::formatString("Will not add unknown option %s to particle type %s", key.c_str(), id_.id().c_str()).c_str()));
    }
}

static InteractionType remapInteractionType(InteractionType itype)
{
    if (itype == InteractionType::BONDCORRECTIONS)
    {
        itype = InteractionType::ELECTRONEGATIVITYEQUALIZATION;
    }
    else if (itype == InteractionType::ANGLES || 
             itype == InteractionType::LINEAR_ANGLES ||
             itype == InteractionType::IMPROPER_DIHEDRALS ||
             itype == InteractionType::PROPER_DIHEDRALS)
    {
        itype = InteractionType::BONDS;
    } 
    return itype;
}
    
bool ParticleType::hasInteractionType(InteractionType itype) const
{
    itype = remapInteractionType(itype);
    for(const auto &s2i : stringToItype)
    {
        if (s2i.second == itype)
        {
            return hasOption(s2i.first);
        }
    }
    return false;
}

Identifier ParticleType::interactionTypeToIdentifier(InteractionType itype) const
{
    auto itypeOrig = itype;
    itype = remapInteractionType(itype);
    for(const auto &s2i : stringToItype)
    {
        if (s2i.second == itype && hasOption(s2i.first))
        {
            return Identifier(optionValue(s2i.first));
        }
    }
    GMX_THROW(gmx::InvalidInputError(gmx::formatString("Interaction type %s not present in particle %s", interactionTypeToString(itypeOrig).c_str(), id_.id().c_str()).c_str()));
    // To make the compiler happy.
    Identifier dummy;
    return dummy;
}

double ParticleType::paramValue(const std::string &type) const
{
    auto pm = parameterMap_.find(type);
    if (pm == parameterMap_.end())
    {
        GMX_THROW(gmx::InternalError(gmx::formatString("No such parameter type %s", type.c_str()).c_str()));
    }
    return pm->second.value();
}    

const ForceFieldParameter &ParticleType::parameterConst(const std::string &type) const
{
    if (parameterMap_.find(type) == parameterMap_.end())
    {
        GMX_THROW(gmx::InternalError(gmx::formatString("No such parameter %s in particle %s", type.c_str(), id().id().c_str()).c_str()));
    }
    return parameterMap_.find(type)->second;
}

ForceFieldParameter *ParticleType::parameter(const std::string &type)
{
    if (parameterMap_.find(type) == parameterMap_.end())
    {
        GMX_THROW(gmx::InternalError(gmx::formatString("No such parameter %s in particle %s", type.c_str(), id().id().c_str()).c_str()));
    }
    return &parameterMap_.find(type)->second;
}

double ParticleType::mass() const
{
    std::string mass("mass");
    if (hasParameter(mass))
    {
        return paramValue(mass);
    }
    return 0.0;
}

double ParticleType::charge() const
{
    std::string charge("charge");
    if (hasParameter(charge))
    {
        return paramValue(charge);
    }
    return 0.0;
}

static int myatoi(const std::string &str)
{
    int d;
    
    if (sscanf(str.c_str(), "%d", &d) == 1)
    {
        return d;
    }
    return 0;
}

int ParticleType::atomnumber() const
{
    std::string anr("atomnumber");
    if (hasOption(anr))
    {
        return myatoi(optionValue(anr));
    }
    return 0;
}

int ParticleType::row() const
{
    std::string anr("row");
    if (hasOption(anr))
    {
        return myatoi(optionValue(anr));
    }
    return 0;
}

std::string ParticleType::element() const
{
    std::string elem("element");
    if (hasOption(elem))
    {
        return optionValue(elem);
    }
    std::string empty;
    return empty;
}

CommunicationStatus ParticleType::Send(const CommunicationRecord *cr, int dest)
{
    CommunicationStatus cs = CommunicationStatus::OK;
    id_.Send(cr, dest);
    cr->send(dest, desc_);
    cr->send(dest, actParticleToString(apType_));
    cr->send(dest, option_.size());
    for(const auto &opt : option_)
    {
        cr->send(dest, opt.first);
        cr->send(dest, opt.second);
    }
    cr->send(dest, parameterMap_.size());
    for(const auto &param : parameterMap_)
    {
        cr->send(dest, param.first);
        cs = param.second.Send(cr, dest);
    }
    return cs;
}
    
CommunicationStatus ParticleType::BroadCast(const CommunicationRecord *cr,
                                            int                        root,
                                            MPI_Comm                   comm)
{
    CommunicationStatus cs = cr->bcast_data(comm);
    if (CommunicationStatus::OK == cs)
    {
        cs = id_.BroadCast(cr, root, comm);
        cr->bcast(&desc_, comm);
        std::string aptype;
        if (cr->isMaster())
        {
            aptype = actParticleToString(apType_);
        }
        cr->bcast(&aptype, comm);
        if (!stringToActParticle(aptype, &apType_))
        {
            GMX_THROW(gmx::InternalError("Communicating ActParticle"));
        }
        size_t nopt = option_.size();
        cr->bcast(&nopt, comm);
        if (cr->rank() == root)
        {
            for(auto &opt : option_)
            {
                std::string key(opt.first), value(opt.second);
                cr->bcast(&key, comm);
                cr->bcast(&value, comm);
            }
        }
        else
        {
            option_.clear();
            for(size_t i = 0; i < nopt; i++)
            {
                std::string key, value;
                cr->bcast(&key, comm);
                cr->bcast(&value, comm);
                setOption(key, value);
            }
        }
    }
    size_t nparm = parameterMap_.size();
    cr->bcast(&nparm, comm);
    if (cr->rank() == root)
    {
        for(auto &pp : parameterMap_)
        {
            std::string type(pp.first);
            cr->bcast(&type, comm);
            pp.second.BroadCast(cr, root, comm);
        }
    }
    else
    {
        parameterMap_.clear();
        for(size_t i = 0; i < nparm; i++)
        {
            std::string type;
            cr->bcast(&type, comm);
            ForceFieldParameter ff;
            ff.BroadCast(cr, root, comm);
            parameterMap_.insert({type, ff});
        }
    }
    return cs;
}

CommunicationStatus ParticleType::Receive(const CommunicationRecord *cr, int src)
{
    CommunicationStatus cs = CommunicationStatus::OK;
    cs = id_.Receive(cr, src);
    cr->recv(src, &desc_);
    std::string aptype;
    cr->recv(src, &aptype);
    if (!stringToActParticle(aptype, &apType_))
    {
        GMX_THROW(gmx::InternalError("Problem receiving ActParticle"));
    }
    size_t nopt;
    cr->recv(src, &nopt);
    option_.clear();
    for(size_t i = 0; i < nopt; i++)
    {
        std::string key, value;
        cr->recv(src, &key);
        cr->recv(src, &value);
        setOption(key, value);
    }
    size_t nparm;
    cr->recv(src, &nparm);
    parameterMap_.clear();
    for(size_t i = 0; i < nparm; i++)
    {
        std::string type;
        cr->recv(src, &type);
        ForceFieldParameter ff;
        ff.Receive(cr, src);
        parameterMap_.insert({type, ff});
    }
    return cs;
}

} // namespace


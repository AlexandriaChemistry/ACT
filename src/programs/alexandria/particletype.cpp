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

#include "particletype.h"

#include <cstdlib>

#include <list>
#include <map>

#include "gromacs/utility/exceptions.h"
#include "gmx_simple_comm.h"
#include "identifier.h"

namespace alexandria
{

std::map<std::string, InteractionType> stringToItype =
    {
        { "acmtype",  InteractionType::ELECTRONEGATIVITYEQUALIZATION },
        { "zetatype", InteractionType::CHARGEDISTRIBUTION },
        { "poltype",  InteractionType::POLARIZATION },
        { "bondtype", InteractionType::BONDS },
        { "vdwtype",  InteractionType::VDW }
    };
    
std::list<std::string> particleOptions = 
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
    
bool ParticleType::hasInteractionType(InteractionType itype) const
{
    for(const auto &s2i : stringToItype)
    {
        if (s2i.second == itype)
        {
            return true;
        }
    }
    return false;
}

Identifier ParticleType::interactionTypeToIdentifier(InteractionType itype) const
{
    for(const auto &s2i : stringToItype)
    {
        if (s2i.second == itype && hasOption(s2i.first))
        {
            return Identifier({optionValue(s2i.first)}, CanSwap::No);
        }
    }
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

const ForceFieldParameter &ParticleType::parameter(const std::string &type) const
{
    if (parameterMap_.find(type) == parameterMap_.end())
    {
        GMX_THROW(gmx::InternalError(gmx::formatString("No such parameter %s in particle %s", type.c_str(), id().id().c_str()).c_str()));
    }
    return parameterMap_.find(type)->second;
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

double ParticleType::refEnthalpy() const
{ 
    std::string refenthalpy("ref_enthalpy");
    if (hasParameter(refenthalpy))
    {
        return paramValue(refenthalpy);
    }
    return 0.0;
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

CommunicationStatus ParticleType::Send(const t_commrec *cr, int dest)
{
    CommunicationStatus cs;
    id_.Send(cr, dest);
    gmx_send_str(cr, dest, &desc_);
    gmx_send_int(cr, dest, gmxParticleType_);
    gmx_send_int(cr, dest, option_.size());
    for(const auto &opt : option_)
    {
        gmx_send_str(cr, dest, &opt.first);
        gmx_send_str(cr, dest, &opt.second);
    }
    gmx_send_int(cr, dest, parameterMap_.size());
    for(const auto &param : parameterMap_)
    {
        gmx_send_str(cr, dest, &param.first);
        cs = param.second.Send(cr, dest);
    }
    return cs;
}
    
CommunicationStatus ParticleType::Receive(const t_commrec *cr, int src)
{
    CommunicationStatus cs;
    cs = id_.Receive(cr, src);
    gmx_recv_str(cr, src, &desc_);
    gmxParticleType_ = gmx_recv_int(cr, src);
    int nopt = gmx_recv_int(cr, src);
    for(int i = 0; i < nopt; i++)
    {
        std::string key, value;
        gmx_recv_str(cr, src, &key);
        gmx_recv_str(cr, src, &value);
        option_.insert({key, value});
    }
    int nparm = gmx_recv_int(cr, src);
    for(int i = 0; i < nparm; i++)
    {
        std::string type;
        gmx_recv_str(cr, src, &type);
        ForceFieldParameter ff;
        ff.Receive(cr, src);
        parameterMap_.insert({type, ff});
    }
    return cs;
}

} // namespace


/*
 * This source file is part of the Alexandria Chemistry Toolkit.
 *
 * Copyright (C) 2020 
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
 * \author David van der Spoel <david.vanderspoel@icm.uu.se>
 */

#include "forcefieldparameterlist.h"

#include <map>

#include "gromacs/topology/idef.h"
#include "gromacs/topology/ifunc.h"
#include "gromacs/utility/exceptions.h"
#include "gromacs/utility/fatalerror.h"
#include "gromacs/utility/stringutil.h"

namespace alexandria
{

ForceFieldParameterList::ForceFieldParameterList(const std::string &function,
                                                 CanSwap            canSwap) : function_(function), canSwap_(canSwap)
{
    if (function.empty())
    {
        fType_ = F_NRE;
    }
    else
    {
        size_t funcType;
        for (funcType = 0; funcType < F_NRE; funcType++)
        {
            if (strcasecmp(interaction_function[funcType].name, function_.c_str()) == 0)
            {
                break;
            }
        }
        if (funcType == F_NRE)
        {
            GMX_THROW(gmx::InvalidInputError(gmx::formatString("Force function '%s' does not exist in gromacs", function_.c_str()).c_str()));
        }
        fType_ = funcType;
    }
}

void ForceFieldParameterList::addParameter(const Identifier          &identifier,
                                           const std::string         &type,
                                           const ForceFieldParameter &param)
{
    auto params                  = parameters_.find(identifier);
    ForceFieldParameter newParam = param;    
    newParam.setIndex(counter_);
    if (params == parameters_.end())
    {
        // New parameter!
        std::map<std::string, ForceFieldParameter> entry;
        entry[type] = newParam;
        parameters_[identifier] = entry;
        // Increase counter since we have a new parameter
        counter_ += 1;
    }
    else
    {
        if (params->second.find(type) != params->second.end())
        {
            GMX_THROW(gmx::InvalidInputError(gmx::formatString("A parameter with type %s was defined for identifier %s already.", type.c_str(), identifier.id().c_str()).c_str()));
        }
        parameters_[identifier][type] = newParam;
    }
}

size_t ForceFieldParameterList::parameterId(const Identifier &identifier) const
{
    auto p = parameters_.find(identifier);
    if (p == parameters_.end())
    {
        GMX_THROW(gmx::InternalError(gmx::formatString("Cannot find parameter %s for %s",
                                                       identifier.id().c_str(), 
                                                       interaction_function[fType()].name).c_str()));
    }
    auto ffp = p->second.find("unit");
    if (ffp != p->second.end())
    {
        return ffp->second.index();
    }
    else
    {
        GMX_THROW(gmx::InternalError(gmx::formatString("Empty parameter list for %s in %s",
                                                       identifier.id().c_str(), 
                                                       interaction_function[fType()].name).c_str()));
    }
}

const std::map<std::string, ForceFieldParameter> &ForceFieldParameterList::findParametersConst(const Identifier &identifier) const
{
    auto params = parameters_.find(identifier);
    
    if (params == parameters_.end())
    {
        GMX_THROW(gmx::InvalidInputError(gmx::formatString("No such identifier %s in const parameter list with %d entries for function '%s'",
                                                           identifier.id().c_str(),
                                                           static_cast<int>(parameters_.size()),
                                                           function_.c_str()).c_str()));
    }
    
    return params->second;
}

const ForceFieldParameter &ForceFieldParameterList::findParameterTypeConst(const Identifier  &identifier,
                                                                           const std::string &type) const
{
    auto params = parameters_.find(identifier);
    
    if (params == parameters_.end())
    {
        GMX_THROW(gmx::InvalidInputError(gmx::formatString("No such identifier %s in parameter list for %s looking for type %s", identifier.id().c_str(), function_.c_str(), type.c_str()).c_str()));
    }
    auto ffparam = params->second.find(type);
    if (ffparam == params->second.end())
    {
        GMX_THROW(gmx::InvalidInputError(gmx::formatString("No such type %s in parameter list for %s", type.c_str(), function_.c_str()).c_str()));
    }
    return ffparam->second;
}

ForceFieldParameter *ForceFieldParameterList::findParameterType(const Identifier  &identifier,
                                                                const std::string &type)
{
    auto params = parameters_.find(identifier);
    
    if (params == parameters_.end())
    {
        GMX_THROW(gmx::InvalidInputError(gmx::formatString("No such identifier %s in parameter list for %s looking for type %s", identifier.id().c_str(), function_.c_str(), type.c_str()).c_str()));
    }
    auto ffparam = params->second.find(type);
    if (ffparam == params->second.end())
    {
        GMX_THROW(gmx::InvalidInputError(gmx::formatString("No such type %s in parameter list for %s", type.c_str(), function_.c_str()).c_str()));
    }
    return &ffparam->second;
}

ForceFieldParameterMap *ForceFieldParameterList::findParameters(const Identifier &identifier)
{
    auto params = parameters_.find(identifier);
    
    if (params == parameters_.end())
    {
        GMX_THROW(gmx::InvalidInputError(gmx::formatString("No such identifier %s in mutable parameter list for %s", identifier.id().c_str(), function_.c_str()).c_str()));
    }
    
    return &params->second;
}

CommunicationStatus ForceFieldParameterList::Send(const t_commrec *cr, int dest) const
{
    CommunicationStatus cs;
    cs = gmx_send_data(cr, dest);
    if (CS_OK == cs)
    {
        gmx_send_str(cr, dest, &function_);
        std::string canSwapString = canSwapToString(canSwap_);
        gmx_send_str(cr, dest, &canSwapString);
        gmx_send_int(cr, dest, fType_);
        gmx_send_int(cr, dest, options_.size());
        for(auto const &x : options_)
        {
            gmx_send_str(cr, dest, &x.first);
            gmx_send_str(cr, dest, &x.second);
        }
        gmx_send_int(cr, dest, parameters_.size());
        for(auto const &x : parameters_)
        {
            cs = x.first.Send(cr, dest);
            if (CS_OK == cs)
            {
                gmx_send_int(cr, dest, x.second.size());
                for(auto const &p : x.second)
                {
                    gmx_send_str(cr, dest, &p.first);
                    cs = p.second.Send(cr, dest);
                    if (CS_OK != cs)
                    {
                        break;
                    }
                }
            }
            else
            {
                break;
            }
        }
    }
    return cs;
}

CommunicationStatus ForceFieldParameterList::Receive(const t_commrec *cr, int src)
{
    CommunicationStatus cs;
    cs = gmx_recv_data(cr, src);
    if (CS_OK == cs)
    {
        gmx_recv_str(cr, src, &function_);
        std::string canSwapString;
        gmx_recv_str(cr, src, &canSwapString);
        canSwap_     = stringToCanSwap(canSwapString);
        fType_       = gmx_recv_int(cr, src);
        int noptions = gmx_recv_int(cr, src);
        options_.clear();
        for(int i = 0; i < noptions; i++)
        {
            std::string key, value;
            gmx_recv_str(cr, src, &key);
            gmx_recv_str(cr, src, &value);
            options_.insert({key, value});
        }
        if (debug)
        {
            fprintf(debug, "Done receiving options\n");
            fflush(debug);
        }
        int nparam =  gmx_recv_int(cr, src);
        parameters_.clear();
        for(int i = 0; i < nparam; i++)
        {
            Identifier key;
            cs = key.Receive(cr, src);
            if (debug)
            {
                fprintf(debug, "Done receiving key %s\n", key.id().c_str());
                fflush(debug);
            }
            if (CS_OK == cs)
            {
                int ntype = gmx_recv_int(cr, src);
                if (debug)
                {
                    fprintf(debug, "Done receiving ntype = %d\n", ntype);
                    fflush(debug);
                }
                for(int j = 0; j < ntype; j++)
                {
                    std::string type;
                    gmx_recv_str(cr, src, &type);
                    ForceFieldParameter p;
                    cs = p.Receive(cr, src);
                    if (CS_OK == cs)
                    {
                        parameters_[key].insert({type, p});
                        if (debug)
                        {
                            fprintf(debug, "Done receiving parameter\n");
                            fflush(debug);
                        }
                    }
                    else
                    {
                        break;
                    }
                }
                if (debug)
                {
                    fprintf(debug, "Done receiving %d parameters\n", nparam);
                    fflush(debug);
                }
            }
            else
            {
                break;
            }
        }
    }
    return cs;
}

} // namespace

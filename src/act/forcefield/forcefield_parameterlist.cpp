/*
 * This source file is part of the Alexandria Chemistry Toolkit.
 *
 * Copyright (C) 2020-2022
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
 * \author David van der Spoel <david.vanderspoel@icm.uu.se>
 */

#include "forcefield_parameterlist.h"

#include <map>

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
    auto params = parameters_.find(identifier);
    
    if (params == parameters_.end())
    {
        // New parameter!
        std::map<std::string, ForceFieldParameter> newParam = {
            { type, param }
        };
        newParam[type].setIndex(counter_);
        parameters_.insert({identifier, newParam});
        // Increase counter since we have a new parameter
        counter_ += 1;
    }
    else
    {
        if (params->second.find(type) == params->second.end())
        {
            params->second.insert({ type, param });
        }
        else
        {
            fprintf(stderr, "Ignoring parameter with type %s defined for identifier %s since it exists already.\n", type.c_str(), identifier.id().c_str());
        }
    }
}

size_t ForceFieldParameterList::parameterId(const Identifier &identifier) const
{
    auto p = parameters_.find(identifier);

    if (p == parameters_.end())
    {
        GMX_THROW(gmx::InternalError(gmx::formatString("Cannot find parameter %s for %s",
                                                       identifier.id().c_str(), 
                                                       interaction_function[gromacsType()].name).c_str()));
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
                                                       interaction_function[gromacsType()].name).c_str()));
    }
}

const std::map<std::string, ForceFieldParameter> &ForceFieldParameterList::findParametersConst(const Identifier &identifier) const
{
    auto iter = parameters_.find(identifier);
    if (parameters_.end() == iter)
    {
        GMX_THROW(gmx::InvalidInputError(gmx::formatString("1. No such identifier '%s' in const parameter list with %d entries for function '%s'",
                                                           identifier.id().c_str(),
                                                           static_cast<int>(parameters_.size()),
                                                           function_.c_str()).c_str()));
    }
    
    return iter->second;
}

const ForceFieldParameter &ForceFieldParameterList::findParameterTypeConst(const Identifier  &identifier,
                                                                           const std::string &type) const
{
    auto params = parameters_.find(identifier);
    
    if (params == parameters_.end())
    {
        GMX_THROW(gmx::InvalidInputError(gmx::formatString("2. No such identifier %s in parameter list for %s looking for type %s", identifier.id().c_str(), function_.c_str(), type.c_str()).c_str()));
    }
    auto ffparam = params->second.find(type);
    if (ffparam == params->second.end())
    {
        GMX_THROW(gmx::InvalidInputError(gmx::formatString("No such type '%s' in parameter list for %s", type.c_str(), function_.c_str()).c_str()));
    }
    return ffparam->second;
}

ForceFieldParameter *ForceFieldParameterList::findParameterType(const Identifier  &identifier,
                                                                const std::string &type)
{
    auto params = parameters_.find(identifier);
    
    if (params == parameters_.end())
    {
        GMX_THROW(gmx::InvalidInputError(gmx::formatString("3. No such identifier %s in parameter list for %s looking for type %s", identifier.id().c_str(), function_.c_str(), type.c_str()).c_str()));
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
        GMX_THROW(gmx::InvalidInputError(gmx::formatString("4. No such identifier %s in mutable parameter list for %s", identifier.id().c_str(), function_.c_str()).c_str()));
    }
    
    return &params->second;
}

bool ForceFieldParameterList::parameterExists(const Identifier &identifier) const
{
    return parameters_.end() != parameters_.find(identifier);
}

void ForceFieldParameterList::dump(FILE *fp) const
{
    if (!fp)
    {
        return;
    }
    fprintf(fp, "Function: %s\n", function_.c_str());
    fprintf(fp, "CanSwap: %s\n", canSwapToString(canSwap_).c_str());
    fprintf(fp, "Ftype: %d\n", fType_);
    for(const auto &opt : options_)
    {
        fprintf(fp, "Option: type='%s' value='%s'\n", opt.first.c_str(),
                opt.second.c_str());
    }
    for(const auto &param : parameters_)
    {
        fprintf(fp, "Identifier: %s\n", param.first.id().c_str());
        for(const auto &key: param.second)
        {
            fprintf(fp, "  Type: %s Value %g\n",
                    key.first.c_str(),
                    key.second.value());
        }
    }
}

CommunicationStatus ForceFieldParameterList::Send(const CommunicationRecord *cr, int dest) const
{
    CommunicationStatus cs = CommunicationStatus::OK;
    if (CommunicationStatus::SEND_DATA == cr->send_data(dest))
    {
        cr->send_str(dest, &function_);
        std::string canSwapString = canSwapToString(canSwap_);
        if (canSwapString.empty())
        {
            GMX_THROW(gmx::InternalError("Empty canSwapString"));
        }
        cr->send_str(dest, &canSwapString);
        cr->send_int(dest, fType_);
        cr->send_int(dest, options_.size());
        for(auto const &x : options_)
        {
            cr->send_str(dest, &x.first);
            cr->send_str(dest, &x.second);
        }
        cr->send_int(dest, parameters_.size());
        for(auto const &x : parameters_)
        {
            cs = x.first.Send(cr, dest);
            if (CommunicationStatus::OK == cs)
            {
                cr->send_int(dest, x.second.size());
                for(auto const &p : x.second)
                {
                    cr->send_str(dest, &p.first);
                    cs = p.second.Send(cr, dest);
                    if (CommunicationStatus::OK != cs)
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
        cr->send_int(dest, counter_);
    }
    return cs;
}

CommunicationStatus ForceFieldParameterList::BroadCast(const CommunicationRecord *cr,
                                                       int                        root,
                                                       MPI_Comm                   comm)
{
    CommunicationStatus cs = cr->bcast_data(comm);
    if (CommunicationStatus::OK == cs)
    {
        cr->bcast(&function_, comm);
        std::string canSwapString;
        if (cr->rank() == root)
        {
            canSwapString.assign(canSwapToString(canSwap_));
        }
        cr->bcast(&canSwapString, comm);
        canSwap_ = stringToCanSwap(canSwapString);
        int ftype = fType_;
        cr->bcast(&ftype, comm);
        fType_ = ftype;
        int noptions = options_.size();
        cr->bcast(&noptions, comm);
        if (cr->rank() == root)
        {
            for(const auto &opt : options_)
            {
                std::string key   = opt.first;
                std::string value = opt.second;
                cr->bcast(&key, comm);
                cr->bcast(&value, comm);
            }
        }
        else
        {
            options_.clear();
            for(int i = 0; i < noptions; i++)
            {
                std::string key, value;
                cr->bcast(&key, comm);
                cr->bcast(&value, comm);
                options_.insert({key, value});
            }
        }
        if (debug)
        {
            fprintf(debug, "Done receiving options\n");
            fflush(debug);
        }
        int nparam = parameters_.size();
        cr->bcast(&nparam, comm);
        if (cr->rank() == root)
        {
            for(auto &p : parameters_)
            {
                Identifier pp = p.first;
                pp.BroadCast(cr, root, comm);
                int ntype = p.second.size();
                cr->bcast(&ntype, comm);
                for(auto &q : p.second)
                {
                    std::string type = q.first;
                    cr->bcast(&type, comm);
                    q.second.BroadCast(cr, root, comm);
                }
            }
        }
        else
        {
            parameters_.clear();
            for(int i = 0; i < nparam; i++)
            {
                Identifier key;
                cs = key.BroadCast(cr, root, comm);
                if (debug)
                {
                    fprintf(debug, "Done broadcasting key %s\n", key.id().c_str());
                    fflush(debug);
                }
                if (CommunicationStatus::OK == cs)
                {
                    int ntype;
                    cr->bcast(&ntype, comm);
                    for(int j = 0; j < ntype; j++)
                    {
                        std::string type;
                        cr->bcast(&type, comm);
                        ForceFieldParameter p;
                        cs = p.BroadCast(cr, root, comm);
                        if (CommunicationStatus::OK == cs)
                        {
                            parameters_[key].insert({type, p});
                            if (debug)
                            {
                                fprintf(debug, "Done receiving parameter\n");
                                fflush(debug);
                            }
                        }
                    }
                }
                else
                {
                    break;
                }
                if (debug)
                {
                    fprintf(debug, "Done receiving %d parameters\n", nparam);
                    fflush(debug);
                }
            }
        }
        int ctr = counter_;
        cr->bcast(&ctr, comm);
        counter_ = ctr;
    }
    return cs;
}

CommunicationStatus ForceFieldParameterList::Receive(const CommunicationRecord *cr, int src)
{
    CommunicationStatus cs = CommunicationStatus::OK;
    if (CommunicationStatus::RECV_DATA == cr->recv_data(src))
    {
        cr->recv_str(src, &function_);
        std::string canSwapString;
        cr->recv_str(src, &canSwapString);
        canSwap_     = stringToCanSwap(canSwapString);
        fType_       = cr->recv_int(src);
        int noptions = cr->recv_int(src);
        options_.clear();
        for(int i = 0; i < noptions; i++)
        {
            std::string key, value;
            cr->recv_str(src, &key);
            cr->recv_str(src, &value);
            options_.insert({key, value});
        }
        if (debug)
        {
            fprintf(debug, "Done receiving options\n");
            fflush(debug);
        }
        int nparam =  cr->recv_int(src);
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
            if (CommunicationStatus::OK == cs)
            {
                int ntype = cr->recv_int(src);
                if (debug)
                {
                    fprintf(debug, "Done receiving ntype = %d\n", ntype);
                    fflush(debug);
                }
                for(int j = 0; j < ntype; j++)
                {
                    std::string type;
                    cr->recv_str(src, &type);
                    ForceFieldParameter p;
                    cs = p.Receive(cr, src);
                    if (CommunicationStatus::OK == cs)
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
        counter_ = cr->recv_int(src);
    }
    return cs;
}

} // namespace

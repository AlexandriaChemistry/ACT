/*
 * This source file is part of the Alexandria Chemistry Toolkit.
 *
 * Copyright (C) 2020-2025
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

#include "act/basics/chargemodel.h"
#include "act/forcefield/potential.h"
#include "gromacs/topology/ifunc.h"
#include "gromacs/utility/exceptions.h"
#include "gromacs/utility/fatalerror.h"
#include "gromacs/utility/stringutil.h"

namespace alexandria
{

ForceFieldParameterList::ForceFieldParameterList(const std::string &function,
                                                 CanSwap            canSwap) : canSwap_(canSwap)
{
    setFunction(function);
}

void ForceFieldParameterList::setFunction(const std::string &function)
{
    if (!stringToPotential(function, &pot_))
    {
        GMX_THROW(gmx::InvalidInputError(gmx::formatString("Potential function '%s' unknown to ACT", function.c_str()).c_str()));
    }
}

void ForceFieldParameterList::addOption(const std::string &option, const std::string &value)
{
    if ((Potential::COULOMB_GAUSSIAN == pot_ ||
         Potential::COULOMB_POINT    == pot_ ||
         Potential::COULOMB_SLATER   == pot_ ) &&
        option.compare("chargetype") == 0)
    {
        auto qdist = name2ChargeDistributionType(value);
        pot_ =  chargeDistributionTypeToPotential(qdist);
    }
    else
    {
        options_.insert({option, value});
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
                                                       potentialToString(pot_).c_str()).c_str()));
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
                                                       potentialToString(pot_).c_str()).c_str()));
    }
}

//! \brief Placeholder for an empty ForceFieldParameterMap
static ForceFieldParameterMap EmptyForceFieldParameterMap{};

const ForceFieldParameterMap &ForceFieldParameterList::findParameterMapConst(const Identifier &identifier) const
{
    auto iter = parameters_.find(identifier);
    if (parameters_.end() == iter)
    {
        return EmptyForceFieldParameterMap;
    }
    
    return iter->second;
}

const ForceFieldParameterMap &ForceFieldParameterList::findParametersConst(const Identifier &identifier) const
{
    auto iter = parameters_.find(identifier);
    if (parameters_.end() == iter)
    {
        GMX_THROW(gmx::InvalidInputError(gmx::formatString("1. No such identifier '%s' in const parameter list with %zu entries for function '%s'",
                                                           identifier.id().c_str(),
                                                           parameters_.size(),
                                                           potentialToString(pot_).c_str()).c_str()));
    }
    
    return iter->second;
}

const ForceFieldParameter &ForceFieldParameterList::findParameterTypeConst(const Identifier  &identifier,
                                                                           const std::string &type) const
{
    auto params = parameters_.find(identifier);
    
    if (params == parameters_.end())
    {
        GMX_THROW(gmx::InvalidInputError(gmx::formatString("2. No such identifier %s in parameter list for %s looking for type %s", identifier.id().c_str(), potentialToString(pot_).c_str(), type.c_str()).c_str()));
    }
    auto ffparam = params->second.find(type);
    if (ffparam == params->second.end())
    {
        GMX_THROW(gmx::InvalidInputError(gmx::formatString("No such type '%s' in parameter list for %s", type.c_str(), potentialToString(pot_).c_str()).c_str()));
    }
    return ffparam->second;
}

ForceFieldParameter *ForceFieldParameterList::findParameterType(const Identifier  &identifier,
                                                                const std::string &type)
{
    auto params = parameters_.find(identifier);
    
    if (params == parameters_.end())
    {
        GMX_THROW(gmx::InvalidInputError(gmx::formatString("3. No such identifier %s in parameter list for %s looking for type %s", identifier.id().c_str(), potentialToString(pot_).c_str(), type.c_str()).c_str()));
    }
    auto ffparam = params->second.find(type);
    if (ffparam == params->second.end())
    {
        GMX_THROW(gmx::InvalidInputError(gmx::formatString("No such type %s in parameter list for %s", type.c_str(), potentialToString(pot_).c_str()).c_str()));
    }
    return &ffparam->second;
}

ForceFieldParameterMap *ForceFieldParameterList::findParameters(const Identifier &identifier)
{
    auto params = parameters_.find(identifier);
    
    if (params == parameters_.end())
    {
        GMX_THROW(gmx::InvalidInputError(gmx::formatString("4. No such identifier %s in mutable parameter list for %s", identifier.id().c_str(), potentialToString(pot_).c_str()).c_str()));
    }
    
    return &params->second;
}

const ForceFieldParameterMap *ForceFieldParameterList::findParametersPtrConst(const Identifier &identifier) const
{
    auto params = parameters_.find(identifier);
    
    if (params == parameters_.end())
    {
        GMX_THROW(gmx::InvalidInputError(gmx::formatString("4. No such identifier %s in mutable parameter list for %s", identifier.id().c_str(), potentialToString(pot_).c_str()).c_str()));
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
    fprintf(fp, "Function: %s\n", potentialToString(pot_).c_str());
    fprintf(fp, "CanSwap: %s\n", canSwapToString(canSwap_).c_str());
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
        std::string ppp = potentialToString(pot_);
        cr->send(dest, ppp);
        std::string canSwapString = canSwapToString(canSwap_);
        if (canSwapString.empty())
        {
            GMX_THROW(gmx::InternalError("Empty canSwapString"));
        }
        cr->send(dest, canSwapString);
        cr->send(dest, options_.size());
        for(auto const &x : options_)
        {
            cr->send(dest, x.first);
            cr->send(dest, x.second);
        }
        cr->send(dest, combrules_.size());
        for(auto const &x : combrules_)
        {
            cr->send(dest, x.first);
            cr->send(dest, x.second);
        }
        cr->send(dest, parameters_.size());
        for(auto const &x : parameters_)
        {
            cs = x.first.Send(cr, dest);
            if (CommunicationStatus::OK == cs)
            {
                cr->send(dest, x.second.size());
                for(auto const &p : x.second)
                {
                    cr->send(dest, p.first);
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
        cr->send(dest, counter_);
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
        std::string ppp = potentialToString(pot_);
        cr->bcast(&ppp, comm);
        if (!stringToPotential(ppp, &pot_))
        {
            GMX_THROW(gmx::InternalError(gmx::formatString("Invalid potential function '%s'", ppp.c_str()).c_str()));
        }
        std::string canSwapString;
        if (cr->rank() == root)
        {
            canSwapString.assign(canSwapToString(canSwap_));
        }
        cr->bcast(&canSwapString, comm);
        canSwap_ = stringToCanSwap(canSwapString);
        size_t noptions = options_.size();
        cr->bcast(&noptions, comm);
        size_t ncrule = combrules_.size();
        cr->bcast(&ncrule, comm);
        if (cr->rank() == root)
        {
            for(const auto &opt : options_)
            {
                std::string key   = opt.first;
                std::string value = opt.second;
                cr->bcast(&key, comm);
                cr->bcast(&value, comm);
            }
            for(const auto &opt : combrules_)
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
            for(size_t i = 0; i < noptions; i++)
            {
                std::string key, value;
                cr->bcast(&key, comm);
                cr->bcast(&value, comm);
                options_.insert({key, value});
            }
            combrules_.clear();
            for(size_t i = 0; i < ncrule; i++)
            {
                std::string key, value;
                cr->bcast(&key, comm);
                cr->bcast(&value, comm);
                combrules_.insert({key, value});
            }
        }
        if (debug)
        {
            fprintf(debug, "Done receiving options\n");
            fflush(debug);
        }
        size_t nparam = parameters_.size();
        cr->bcast(&nparam, comm);
        if (cr->rank() == root)
        {
            for(auto &p : parameters_)
            {
                Identifier pp = p.first;
                pp.BroadCast(cr, root, comm);
                size_t ntype = p.second.size();
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
            for(size_t i = 0; i < nparam; i++)
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
                    size_t ntype;
                    cr->bcast(&ntype, comm);
                    for(size_t j = 0; j < ntype; j++)
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
                    fprintf(debug, "Done receiving %zu parameters\n", nparam);
                    fflush(debug);
                }
            }
        }
        size_t ctr = counter_;
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
        std::string ppp = potentialToString(pot_);
        cr->recv(src, &ppp);
        if (!stringToPotential(ppp, &pot_))
        {
            GMX_THROW(gmx::InternalError(gmx::formatString("Invalid potential function '%s'", ppp.c_str()).c_str()));
        }
        std::string canSwapString;
        cr->recv(src, &canSwapString);
        canSwap_     = stringToCanSwap(canSwapString);
        size_t noptions;
        cr->recv(src, &noptions);
        options_.clear();
        size_t ncrule;
        for(size_t i = 0; i < noptions; i++)
        {
            std::string key, value;
            cr->recv(src, &key);
            cr->recv(src, &value);
            options_.insert({key, value});
        }
        cr->recv(src, &ncrule);
        combrules_.clear();
        for(size_t i = 0; i < ncrule; i++)
        {
            std::string key, value;
            cr->recv(src, &key);
            cr->recv(src, &value);
            combrules_.insert({key, value});
        }
        if (debug)
        {
            fprintf(debug, "Done receiving options\n");
            fflush(debug);
        }
        size_t nparam;
        cr->recv(src, &nparam);
        parameters_.clear();
        for(size_t i = 0; i < nparam; i++)
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
                size_t ntype;
                cr->recv(src, &ntype);
                if (debug)
                {
                    fprintf(debug, "Done receiving ntype = %zu\n", ntype);
                    fflush(debug);
                }
                for(size_t j = 0; j < ntype; j++)
                {
                    std::string type;
                    cr->recv(src, &type);
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
                    fprintf(debug, "Done receiving %zu parameters\n", nparam);
                    fflush(debug);
                }
            }
            else
            {
                break;
            }
        }
        cr->recv(src, &counter_);
    }
    return cs;
}

} // namespace

/*
 * This source file is part of the Alexandria Chemistry Toolkit.
 *
 * Copyright (C) 2020-2026
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

#include "forcefield_parameter.h"

#include <map>

#include <inttypes.h>

#include "gromacs/utility/exceptions.h"
#include "gromacs/utility/fatalerror.h"
#include "gromacs/utility/stringutil.h"

bool stringToBoolean(const std::string &str)
{
    return (str.compare("yes") == 0 || str.compare("Yes") == 0 ||
            str.compare("YES") == 0);
}

namespace alexandria
{

void ForceFieldParameter::forceSetValue(double value)
{
    value_ = value;
    // Need to update the intervalue since this is what the code will extract
    calculateInternalValue();
}

void ForceFieldParameter::setValue(double value)
{
    switch (mutability_)
    {
    case Mutability::Free:
        value_ = value;
        break;
    case Mutability::ACM:
        // Setting this value will not affect anything but it is needed for edit_ff
        value_ = value;
        break;
    case Mutability::Bounded:
        if (value >= minimum_ && value <= maximum_)
        {
            value_ = value;
        }
        else
        {
            double newval = std::min(maximum_, std::max(minimum_, value));
            if (strict_)
            {
                auto buf = gmx::formatString("Can not modify value to %g as it would be outside its bounds of %g-%g. Setting it to %g %s",
                                             value, minimum_, maximum_, newval, unit_.c_str());
                GMX_THROW(gmx::InvalidInputError(buf));
            }
            value_ = newval;
        }
        break;
    case Mutability::Dependent:
    case Mutability::Fixed:
        if (strict_)
        {
            auto buf = gmx::formatString("Cannot set parameter value to %g %s since it has mutability %s",
                                         value, unit_.c_str(), mutabilityName(mutability_).c_str());
            GMX_THROW(gmx::InvalidInputError(buf));
        }
        break;
    }
    calculateInternalValue();
}

void ForceFieldParameter::setUncertainty(double uncertainty)
{ 
    if (mutability_ == Mutability::Free || mutability_ == Mutability::Bounded)
    {
        uncertainty_ = std::max(0.0, uncertainty);
    }
    else if (strict_)
    {
        auto buf = gmx::formatString("Cannot modify uncertainty since the parameter is fixed/dependent");
        GMX_THROW(gmx::InternalError(buf));
    }
}

void ForceFieldParameter::setNtrain(unsigned int ntrain)
{ 
    if (mutability_ != Mutability::Fixed)
    {
        ntrain_ = ntrain;
    }
    else if (strict_)
    {
        auto buf = gmx::formatString("Cannot modify ntrain since the parameter is fixed/dependent");
        GMX_THROW(gmx::InternalError(buf));
    }
}

void ForceFieldParameter::copy(const ForceFieldParameter &src)
{
    unit_                = src.unit();
    value_               = src.value();
    originalValue_       = src.originalValue();
    uncertainty_         = src.uncertainty();
    originalUncertainty_ = src.uncertainty();
    ntrain_              = src.ntrain();
    originalNtrain_      = src.ntrain();
    minimum_             = src.minimum();
    maximum_             = src.maximum();
    mutability_          = src.mutability();
    strict_              = src.strict();
}

CommunicationStatus ForceFieldParameter::Send(const CommunicationRecord *cr, int dest) const
{
    CommunicationStatus cs = CommunicationStatus::OK;
    if (CommunicationStatus::SEND_DATA == cr->send_data(dest))
    {
        cr->send(dest, unit_);
        cr->send(dest, value_);
        cr->send(dest, mutabilityName(mutability_));
        cr->send(dest, originalValue_);
        cr->send(dest, uncertainty_);
        cr->send(dest, originalUncertainty_);
        cr->send(dest, ntrain_);
        cr->send(dest, originalNtrain_);
        cr->send(dest, minimum_);
        cr->send(dest, maximum_);
        cr->send(dest, nonNegative_);
        cr->send(dest, strict_ ? 1 : 0);
        if (debug)
        {
            fprintf(debug, "Sent most of a parameter\n");
            fflush(debug);
        }
        cr->send(dest, mutabilityName(mutability_));

        if (nullptr != debug)
        {
            fprintf(debug, "Sent ForceFieldParameter %g %g %s ntrain %d, status %s\n",
                    value_, uncertainty_, unit_.c_str(), ntrain_, cs_name(cs).c_str());
            fflush(debug);
        }
    }
    return cs;
}

CommunicationStatus ForceFieldParameter::BroadCast(const CommunicationRecord *cr,
                                                   int                        root,
                                                   MPI_Comm                   comm)
{
    CommunicationStatus cs = cr->bcast_data(comm);
    if (CommunicationStatus::OK == cs)
    {
        cr->bcast(&unit_, comm);
        cr->bcast(&value_, comm);
        std::string mutstr;
        if (cr->rank() == root)
        {
            mutstr = mutabilityName(mutability_);
        }
        cr->bcast(&mutstr, comm);
        Mutability mut;
        if (!nameToMutability(mutstr, &mut))
        {
            GMX_THROW(gmx::InternalError(gmx::formatString("Invalid mutability %s", mutstr.c_str()).c_str()));
        }
        mutability_          = mut;
        cr->bcast(&originalValue_, comm);
        cr->bcast(&uncertainty_, comm);
        cr->bcast(&originalUncertainty_, comm);
        cr->bcast(&ntrain_, comm);
        cr->bcast(&originalNtrain_, comm);
        cr->bcast(&minimum_, comm);
        cr->bcast(&maximum_, comm);
        cr->bcast(&nonNegative_, comm);
        int strict = strict_;
        cr->bcast(&strict, comm);
        strict_ = strict;
        if (debug)
        {
            fprintf(debug, "Received most of ff param\n");
            fflush(debug);
        }
        if (nullptr != debug)
        {
            fprintf(debug, "Received ForceFieldParameter %g %g %s ntrain %d, status %s\n",
                    value_, uncertainty_, unit_.c_str(), ntrain_, cs_name(cs).c_str());
            fflush(debug);
        }
        calculateInternalValue();
    }
    return cs;
}

CommunicationStatus ForceFieldParameter::Receive(const CommunicationRecord *cr, int src)
{
    CommunicationStatus cs = CommunicationStatus::OK;
    if (CommunicationStatus::RECV_DATA == cr->recv_data(src))
    {
        cr->recv(src, &unit_);
        cr->recv(src, &value_);
        std::string mutstr;
        cr->recv(src, &mutstr);
        Mutability mut;
        if (!nameToMutability(mutstr, &mut))
        {
            GMX_THROW(gmx::InternalError(gmx::formatString("Invalid mutability %s", mutstr.c_str()).c_str()));
        }
        mutability_          = mut;
        cr->recv(src, &originalValue_ );
        cr->recv(src, &uncertainty_);
        cr->recv(src, &originalUncertainty_);
        cr->recv(src, &ntrain_);
        cr->recv(src, &originalNtrain_);
        cr->recv(src, &minimum_);
        cr->recv(src, &maximum_);
        cr->recv(src, &nonNegative_);
        cr->recv(src, &strict_);
        if (debug)
        {
            fprintf(debug, "Received most of ff param\n");
            fflush(debug);
        }
        std::string mname;
        cr->recv(src, &mname);
        if (!nameToMutability(mname, &mutability_))
        {
            GMX_THROW(gmx::InternalError("Communicating mutability"));
        }
        if (nullptr != debug)
        {
            fprintf(debug, "Received ForceFieldParameter %g %g %s ntrain %d, status %s\n",
                    value_, uncertainty_, unit_.c_str(), ntrain_, cs_name(cs).c_str());
            fflush(debug);
        }
        calculateInternalValue();
    }
    return cs;
}

} // namespace alexandria

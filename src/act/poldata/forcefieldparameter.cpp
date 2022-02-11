/*
 * This source file is part of the Alexandria Chemistry Toolkit.
 *
 * Copyright (C) 2020-2022
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

#include "forcefieldparameter.h"

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

void ForceFieldParameter::setValue(double value)
{
    switch (mutability_)
    {
    case Mutability::Free:
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
                auto buf = gmx::formatString("Can not modify value outside its bounds of %g-%g. Setting it to %g %s",
                                             minimum_, maximum_, newval, unit_.c_str());
                GMX_THROW(gmx::InvalidInputError(buf));
            }
            value_ = newval;
        }
        break;
    case Mutability::ACM:
    case Mutability::Dependent:
    case Mutability::Fixed:
        if (strict_)
        {
            auto buf = gmx::formatString("Cannot modify parameter since it is fixed/dependent");
            GMX_THROW(gmx::InvalidInputError(buf));
        }
        break;
    }
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

void ForceFieldParameter::setNtrain(int ntrain)
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
        cr->send_str(dest, &unit_);
        cr->send_double(dest, value_);
        cr->send_str(dest, &mutabilityName(mutability_));
        cr->send_double(dest, originalValue_);
        cr->send_double(dest, uncertainty_);
        cr->send_double(dest, originalUncertainty_);
        cr->send_int(dest, ntrain_);
        cr->send_int(dest, originalNtrain_);
        cr->send_double(dest, minimum_);
        cr->send_double(dest, maximum_);
        cr->send_int(dest, strict_ ? 1 : 0);
        if (debug)
        {
            fprintf(debug, "Sent most of a parameter\n");
            fflush(debug);
        }
        std::string mname = mutabilityName(mutability_);
        cr->send_str(dest, &mname);

        if (nullptr != debug)
        {
            fprintf(debug, "Sent ForceFieldParameter %g %g %s %d, status %s\n",
                    value_, uncertainty_, unit_.c_str(), ntrain_, cs_name(cs));
            fflush(debug);
        }
    }
    return cs;
}

CommunicationStatus ForceFieldParameter::Receive(const CommunicationRecord *cr, int src)
{
    CommunicationStatus cs = CommunicationStatus::OK;
    if (CommunicationStatus::RECV_DATA == cr->recv_data(src))
    {
        cr->recv_str(src, &unit_);
        value_               = cr->recv_double(src);
        std::string mutstr;
        cr->recv_str(src, &mutstr);
        Mutability mut;
        if (!nameToMutability(mutstr, &mut))
        {
            GMX_THROW(gmx::InternalError(gmx::formatString("Invalid mutability %s", mutstr.c_str()).c_str()));
        }
        mutability_          = mut;
        originalValue_       = cr->recv_double(src);
        uncertainty_         = cr->recv_double(src);
        originalUncertainty_ = cr->recv_double(src);
        ntrain_              = cr->recv_int(src);
        originalNtrain_      = cr->recv_int(src);
        minimum_             = cr->recv_double(src);
        maximum_             = cr->recv_double(src);
        strict_              = cr->recv_int(src);
        if (debug)
        {
            fprintf(debug, "Received most of ff param\n");
            fflush(debug);
        }
        std::string mname;
        cr->recv_str(src, &mname);
        if (!nameToMutability(mname, &mutability_))
        {
            GMX_THROW(gmx::InternalError("Communicating mutability"));
        }
        if (nullptr != debug)
        {
            fprintf(debug, "Received ForceFieldParameter %g %g %s %d, status %s\n",
                    value_, uncertainty_, unit_.c_str(), ntrain_, cs_name(cs));
            fflush(debug);
        }
    }
    return cs;
}

} // namespace alexandria

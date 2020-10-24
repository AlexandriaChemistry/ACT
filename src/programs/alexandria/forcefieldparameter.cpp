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

#include "forcefieldparameter.h"

#include <map>

#include <inttypes.h>

#include "gromacs/utility/exceptions.h"
#include "gromacs/utility/fatalerror.h"
#include "gromacs/utility/stringutil.h"

namespace alexandria
{

static std::map<Mutability, const std::string> mut2string =
    {
        { Mutability::Fixed,     "Fixed"     },
        { Mutability::Dependent, "Dependent" },
        { Mutability::Bounded,   "Bounded"   },
        { Mutability::Free,      "Free"      }
    };

static std::map<const std::string, Mutability> string2mut;
    
const std::string &mutabilityName(Mutability mutability)
{
    auto m2s = mut2string.find(mutability);
    
    return m2s->second;
}

bool nameToMutability(const std::string &name, Mutability *mutability)
{
    if (string2mut.empty())
    {
        for (auto iter = mut2string.begin(); iter != mut2string.end(); ++iter)
        {
            string2mut.insert({iter->second, iter->first});
        }
    }
    auto s2m = string2mut.find(name);
    if (s2m != string2mut.end())
    {
        *mutability = s2m->second;
        return true;
    }
    else
    {
        return false;
    }
}

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
                auto buf = gmx::formatString("Can not modify value outside its bounds of %g-%g. Setting it to %g.",
                                             minimum_, maximum_, newval);
                GMX_THROW(gmx::InvalidInputError(buf));
            }
            value_ = newval;
        }
        break;
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

void ForceFieldParameter::setNtrain(uint64_t ntrain)
{ 
    if (mutability_ == Mutability::Free || mutability_ == Mutability::Bounded)
    {
        ntrain_ = ntrain;
    }
    else if (strict_)
    {
        auto buf = gmx::formatString("Cannot modify ntrain since the parameter is fixed/dependent");
        GMX_THROW(gmx::InternalError(buf));
    }
}

CommunicationStatus ForceFieldParameter::Send(const t_commrec *cr, int dest) const
{
    CommunicationStatus cs;
    cs = gmx_send_data(cr, dest);
    if (CS_OK == cs)
    {
        gmx_send_str(cr, dest, &unit_);
        gmx_send_double(cr, dest, value_);
        gmx_send_double(cr, dest, originalValue_);
        gmx_send_double(cr, dest, uncertainty_);
        gmx_send_double(cr, dest, originalUncertainty_);
        gmx_send_int(cr, dest, static_cast<int>(ntrain_));
        gmx_send_int(cr, dest, static_cast<int>(originalNtrain_));
        gmx_send_double(cr, dest, minimum_);
        gmx_send_double(cr, dest, maximum_);
        gmx_send_int(cr, dest, strict_ ? 1 : 0);
        if (debug)
        {
            fprintf(debug, "Sent most of a parameter\n");
            fflush(debug);
        }
        switch (mutability_)
        {
        case Mutability::Free:
            gmx_send_int(cr, dest, 0);
            break;
        case Mutability::Bounded:
            gmx_send_int(cr, dest, 1);
            break;
        case Mutability::Dependent:
            gmx_send_int(cr, dest, 2);
            break;
        case Mutability::Fixed:
            gmx_send_int(cr, dest, 3);
            break;
        }

        if (nullptr != debug)
        {
            fprintf(debug, "Sent ForceFieldParameter %g %g %s %" PRIu64 ", status %s\n",
                    value_, uncertainty_, unit_.c_str(), ntrain_, cs_name(cs));
            fflush(debug);
        }
    }
    return cs;
}

CommunicationStatus ForceFieldParameter::Receive(const t_commrec *cr, int src)
{
    CommunicationStatus cs;
    cs = gmx_recv_data(cr, src);
    if (CS_OK == cs)
    {
        gmx_recv_str(cr, src, &unit_);
        value_               = gmx_recv_double(cr, src);
        originalValue_       = gmx_recv_double(cr, src);
        uncertainty_         = gmx_recv_double(cr, src);
        originalUncertainty_ = gmx_recv_double(cr, src);
        ntrain_              = static_cast<uint64_t>(gmx_recv_int(cr, src));
        originalNtrain_      = static_cast<uint64_t>(gmx_recv_int(cr, src));
        minimum_             = gmx_recv_double(cr, src);
        maximum_             = gmx_recv_double(cr, src);
        strict_              = gmx_recv_int(cr, src);
        if (debug)
        {
            fprintf(debug, "Received most of ff param\n");
            fflush(debug);
        }
        int mut              = gmx_recv_int(cr,src);
        switch (mut)
        {
        case 0:
            mutability_ = Mutability::Free;
            break;
        case 1:
            mutability_ = Mutability::Bounded;
            break;
        case 2:
            mutability_ = Mutability::Dependent;
            break;
        case 3:
            mutability_ = Mutability::Fixed;
            break;
        default:
            // Should not happen.
            GMX_THROW(gmx::InternalError(gmx::formatString("Unexpected mutability %d", mut))); 
        }
        if (nullptr != debug)
        {
            fprintf(debug, "Received ForceFieldParameter %g %g %s %" PRIu64 ", status %s\n",
                    value_, uncertainty_, unit_.c_str(), ntrain_, cs_name(cs));
            fflush(debug);
        }
    }
    return cs;
}

} // namespace alexandria

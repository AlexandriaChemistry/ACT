/*
 * This source file is part of the Alexandria Chemistry Toolkit.
 *
 * Copyright (C) 2021
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

#include "vsite.h"

#include <map>
#include <string>

#include "gromacs/utility/exceptions.h"
#include "gromacs/utility/fatalerror.h"

namespace alexandria
{

static std::map<VsiteType, std::string> evt_names = {
    { VsiteType::LINEAR,       "linear"       },
    { VsiteType::PLANAR,       "planar"       },
    { VsiteType::RING_PLANAR,  "ring_planar"  },
    { VsiteType::IN_PLANE,     "in_plane"     },
    { VsiteType::OUT_OF_PLANE, "out_of_plane" },
    { VsiteType::ALL,          "all"          }
};

const char *vsiteType2string(VsiteType vType)
{
    return evt_names[vType].c_str();
}

VsiteType string2vsiteType(const char *mystring)
{
    for (auto &i : evt_names)
    {
        if (i.second.compare(mystring) == 0)
        {
            return i.first;
        }
    }
    GMX_THROW(gmx::InternalError(gmx::formatString("Invalid vsite type %s", mystring).c_str()));
}

Vsite::Vsite(const std::string &atype,
             const std::string &type,
             int                number,
             double             distance,
             double             angle,
             int                ncontrolatoms)
    :
      atype_(atype),
      type_(string2vsiteType(type.c_str())),
      number_(number),
      distance_(distance),
      angle_(angle),
      ncontrolatoms_(ncontrolatoms)
{}

CommunicationStatus Vsite::Send(const CommunicationRecord *cr, int dest)
{
    CommunicationStatus cs = CommunicationStatus::OK;
    if (CommunicationStatus::SEND_DATA == cr->send_data(dest))
    {
        std::string vtype;
        vtype.assign(vsiteType2string(type_));
        cr->send_str(dest, &atype_);
        cr->send_str(dest, &vtype);
        cr->send_int(dest, number_);
        cr->send_double(dest, distance_);
        cr->send_double(dest, angle_);
        cr->send_int(dest, ncontrolatoms_);
        if (nullptr != debug)
        {
            fprintf(debug, "Sent Vsite %s %s %d %g %g %d, status %s\n",
                    atype_.c_str(), vsiteType2string(type_), number_,
                    distance_, angle_, ncontrolatoms_, cs_name(cs));
            fflush(debug);
        }
    }
    return cs;
}

CommunicationStatus Vsite::Receive(const CommunicationRecord *cr, int src)
{
    CommunicationStatus cs = CommunicationStatus::OK;
    if (CommunicationStatus::RECV_DATA == cr->recv_data(src))
    {
        std::string type;
        cr->recv_str(src, &atype_);
        cr->recv_str(src, &type);
        type_          = string2vsiteType(type.c_str());
        number_        = cr->recv_int(src);
        distance_      = cr->recv_double(src);
        angle_         = cr->recv_double(src);
        ncontrolatoms_ = cr->recv_int(src);
        if (nullptr != debug)
        {
            fprintf(debug, "Received Vsite %s %s %d %g %g %d, status %s\n",
                    atype_.c_str(), vsiteType2string(type_), number_,
                    distance_, angle_, ncontrolatoms_, cs_name(cs));
            fflush(debug);
        }
    }
    return cs;
}

} // namespace alexandria

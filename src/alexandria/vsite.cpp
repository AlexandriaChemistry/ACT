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

#include "gmx_simple_comm.h"

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

CommunicationStatus Vsite::Send(const t_commrec *cr, int dest)
{
    CommunicationStatus cs;
    cs = gmx_send_data(cr, dest);
    if (CS_OK == cs)
    {
        std::string vtype;
        vtype.assign(vsiteType2string(type_));
        gmx_send_str(cr, dest, &atype_);
        gmx_send_str(cr, dest, &vtype);
        gmx_send_int(cr, dest, number_);
        gmx_send_double(cr, dest, distance_);
        gmx_send_double(cr, dest, angle_);
        gmx_send_int(cr, dest, ncontrolatoms_);
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

CommunicationStatus Vsite::Receive(const t_commrec *cr, int src)
{
    CommunicationStatus cs;
    cs = gmx_recv_data(cr, src);
    if (CS_OK == cs)
    {
        std::string type;
        gmx_recv_str(cr, src, &atype_);
        gmx_recv_str(cr, src, &type);
        type_          = string2vsiteType(type.c_str());
        number_        = gmx_recv_int(cr, src);
        distance_      = gmx_recv_double(cr, src);
        angle_         = gmx_recv_double(cr, src);
        ncontrolatoms_ = gmx_recv_int(cr, src);
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
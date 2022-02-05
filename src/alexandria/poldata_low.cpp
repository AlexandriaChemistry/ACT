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

#include "poldata_low.h"

#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <cstring>

#include <algorithm>
#include <map>
#include <vector>

#include "communicationrecord.h"
#include "utility/stringutil.h"

namespace alexandria
{

Bosque::Bosque(const std::string &bosque, double polarizability)
    :
      bosque_(bosque),
      polarizability_(polarizability)
{}

CommunicationStatus Bosque::Send(const CommunicationRecord *cr, int dest)
{
    CommunicationStatus cs = CommunicationStatus::OK;
    if (CommunicationStatus::SEND_DATA == cr->send_data(dest))
    {
        cr->send_str(dest, &bosque_);
        cr->send_double(dest, polarizability_);
        if (nullptr != debug)
        {
            fprintf(debug, "Sent Bosque %s %g, status %s\n",
                    bosque_.c_str(), polarizability_, cs_name(cs));
            fflush(debug);
        }
    }
    return cs;
}

CommunicationStatus Bosque::Receive(const CommunicationRecord *cr, int src)
{
    CommunicationStatus cs = CommunicationStatus::OK;
    if (CommunicationStatus::RECV_DATA == cr->recv_data(src))
    {
        cr->recv_str(src, &bosque_);
        polarizability_ = cr->recv_double(src);
        if (nullptr != debug)
        {
            fprintf(debug, "Received Bosque %s %g, status %s\n",
                    bosque_.c_str(), polarizability_, cs_name(cs));
            fflush(debug);
        }
    }
    return cs;
}

Miller::Miller(const std::string &miller,
               int                atomnumber,
               double             tauAhc,
               double             alphaAhp,
               const std::string &alexandria_equiv)
    :
      miller_(miller),
      atomnumber_(atomnumber),
      tauAhc_(tauAhc),
      alphaAhp_(alphaAhp),
      alexandria_equiv_(alexandria_equiv)
{}

CommunicationStatus Miller::Send(const CommunicationRecord *cr, int dest)
{
    CommunicationStatus cs = CommunicationStatus::OK;
    if (CommunicationStatus::SEND_DATA == cr->send_data(dest))
    {
        cr->send_str(dest, &miller_);
        cr->send_int(dest, atomnumber_);
        cr->send_double(dest, tauAhc_);
        cr->send_double(dest, alphaAhp_);
        cr->send_str(dest, &alexandria_equiv_);

        if (nullptr != debug)
        {
            fprintf(debug, "Sent Miller %s %d %g %g %s, status %s\n",
                    miller_.c_str(), atomnumber_, tauAhc_,
                    alphaAhp_, alexandria_equiv_.c_str(), cs_name(cs));
            fflush(debug);
        }
    }
    return cs;
}

CommunicationStatus Miller::Receive(const CommunicationRecord *cr, int src)
{
    CommunicationStatus cs = CommunicationStatus::OK;
    if (CommunicationStatus::RECV_DATA == cr->recv_data(src))
    {
        cr->recv_str(src, &miller_);
        atomnumber_ = cr->recv_int(src);
        tauAhc_     = cr->recv_double(src);
        alphaAhp_   = cr->recv_double(src);
        cr->recv_str(src, &alexandria_equiv_);
        if (nullptr != debug)
        {
            fprintf(debug, "Received Miller %s %d %g %g %s, status %s\n",
                    miller_.c_str(), atomnumber_, tauAhc_,
                    alphaAhp_, alexandria_equiv_.c_str(), cs_name(cs));
            fflush(debug);
        }
    }
    return cs;
}

Symcharges::Symcharges(const std::string &central,
                       const std::string &attached,
                       int                numattach)
    :
      central_(central),
      attached_(attached),
      numattach_(numattach)
{}

CommunicationStatus Symcharges::Send(const CommunicationRecord *cr, int dest)
{
    CommunicationStatus cs = CommunicationStatus::OK;
    if (CommunicationStatus::SEND_DATA == cr->send_data(dest))
    {
        cr->send_str(dest, &central_);
        cr->send_str(dest, &attached_);
        cr->send_int(dest, numattach_);

        if (nullptr != debug)
        {
            fprintf(debug, "Sent Symcharges %s %s %d, status %s\n",
                    central_.c_str(), attached_.c_str(), numattach_, cs_name(cs));
            fflush(debug);
        }
    }
    return cs;
}

CommunicationStatus Symcharges::Receive(const CommunicationRecord *cr, int src)
{
    CommunicationStatus cs = CommunicationStatus::OK;
    if (CommunicationStatus::RECV_DATA == cr->recv_data(src))
    {
        cr->recv_str(src, &central_);
        cr->recv_str(src, &attached_);
        numattach_ = cr->recv_int(src);
        if (nullptr != debug)
        {
            fprintf(debug, "Received Symcharges %s %s %d, status %s\n",
                    central_.c_str(), attached_.c_str(), numattach_, cs_name(cs));
            fflush(debug);
        }
    }
    return cs;
}

} // namespace alexandria

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

//#include "gromacs/topology/ifunc.h"
//#include "gromacs/utility/cstringutil.h"
//#include "gromacs/utility/stringutil.h"

#include "gmx_simple_comm.h"
//#include "identifier.h"
//#include "plistwrapper.h"
#include "stringutil.h"

namespace alexandria
{


Bosque::Bosque(const std::string &bosque, double polarizability)
    :
      bosque_(bosque),
      polarizability_(polarizability)
{}

CommunicationStatus Bosque::Send(const t_commrec *cr, int dest)
{
    CommunicationStatus cs;
    cs = gmx_send_data(cr, dest);
    if (CS_OK == cs)
    {
        gmx_send_str(cr, dest, &bosque_);
        gmx_send_double(cr, dest, polarizability_);
        if (nullptr != debug)
        {
            fprintf(debug, "Sent Bosque %s %g, status %s\n",
                    bosque_.c_str(), polarizability_, cs_name(cs));
            fflush(debug);
        }
    }
    return cs;
}

CommunicationStatus Bosque::Receive(const t_commrec *cr, int src)
{
    CommunicationStatus cs;
    cs = gmx_recv_data(cr, src);
    if (CS_OK == cs)
    {
        gmx_recv_str(cr, src, &bosque_);
        polarizability_ = gmx_recv_double(cr, src);
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

CommunicationStatus Miller::Send(const t_commrec *cr, int dest)
{
    CommunicationStatus cs;
    cs = gmx_send_data(cr, dest);
    if (CS_OK == cs)
    {
        gmx_send_str(cr, dest, &miller_);
        gmx_send_int(cr, dest, atomnumber_);
        gmx_send_double(cr, dest, tauAhc_);
        gmx_send_double(cr, dest, alphaAhp_);
        gmx_send_str(cr, dest, &alexandria_equiv_);

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

CommunicationStatus Miller::Receive(const t_commrec *cr, int src)
{
    CommunicationStatus cs;
    cs = gmx_recv_data(cr, src);
    if (CS_OK == cs)
    {
        gmx_recv_str(cr, src, &miller_);
        atomnumber_ = gmx_recv_int(cr, src);
        tauAhc_     = gmx_recv_double(cr, src);
        alphaAhp_   = gmx_recv_double(cr, src);
        gmx_recv_str(cr, src, &alexandria_equiv_);
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

CommunicationStatus Symcharges::Send(const t_commrec *cr, int dest)
{
    CommunicationStatus cs;
    cs = gmx_send_data(cr, dest);
    if (CS_OK == cs)
    {
        gmx_send_str(cr, dest, &central_);
        gmx_send_str(cr, dest, &attached_);
        gmx_send_int(cr, dest, numattach_);

        if (nullptr != debug)
        {
            fprintf(debug, "Sent Symcharges %s %s %d, status %s\n",
                    central_.c_str(), attached_.c_str(), numattach_, cs_name(cs));
            fflush(debug);
        }
    }
    return cs;
}

CommunicationStatus Symcharges::Receive(const t_commrec *cr, int src)
{
    CommunicationStatus cs;
    cs = gmx_recv_data(cr, src);
    if (CS_OK == cs)
    {
        gmx_recv_str(cr, src, &central_);
        gmx_recv_str(cr, src, &attached_);
        numattach_ = gmx_recv_int(cr, src);
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

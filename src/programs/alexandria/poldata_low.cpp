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

#include "gromacs/topology/ifunc.h"
#include "gromacs/utility/cstringutil.h"
#include "gromacs/utility/stringutil.h"

#include "gmx_simple_comm.h"
#include "identifier.h"
#include "plistwrapper.h"
#include "stringutil.h"

namespace alexandria
{

Ffatype::Ffatype(const std::string &desc,
                 const std::string &type,
                 const std::string &ptype,
                 const std::string &btype,
                 const std::string &ztype,
                 const std::string &acmtype,
                 const std::string &elem,
                 double             mass,
                 int                atomnumber,
                 double             charge,
                 int                row,
                 Mutability         mutability,
                 const std::string &refEnthalpy) :
    desc_(desc), type_(type),
    elem_(elem),
    refEnthalpy_(refEnthalpy),
    mass_(mass),
    atomnumber_(atomnumber),
    charge_(charge),
    row_(row),
    mutability_(mutability)
{
    subType_.insert({InteractionType::VDW, type});
    subType_.insert({InteractionType::POLARIZATION, ptype});
    for(auto &it : {InteractionType::BONDS, InteractionType::ANGLES, InteractionType::LINEAR_ANGLES,
                    InteractionType::PROPER_DIHEDRALS,  InteractionType::IMPROPER_DIHEDRALS })
    {
        subType_.insert({it, btype});
    }
    subType_.insert({InteractionType::CHARGEDISTRIBUTION, ztype});
    subType_.insert({InteractionType::ELECTRONEGATIVITYEQUALIZATION, acmtype});
    subType_.insert({InteractionType::BONDCORRECTIONS, acmtype});
    if (type_.empty())
    {
        GMX_THROW(gmx::InternalError(gmx::formatString("Trying to a type with no name. ptype = %s btype = %s ztype = %s.",
                                                       ptype.c_str(), btype.c_str(),
                                                       ztype.c_str()).c_str()));
    }
}

Identifier Ffatype::id(InteractionType iType) const
{
    auto ss = subType_.find(iType);
    if (ss == subType_.end())
    {
        GMX_THROW(gmx::InternalError(gmx::formatString("No such subType %s for atomtype '%s'",
                                                       interactionTypeToString(iType).c_str(),
                                                       type_.c_str()).c_str()));
    }
    else
    {
        return Identifier({ss->second}, CanSwap::No);
    }
}

static const char * evt_names[evtNR] = {
    "linear",
    "planar",
    "ring_planar",
    "in_plane",
    "out_of_plane",
    "all"
};

const char *vsiteType2string(VsiteType vType)
{
    if (vType < evtNR)
    {
        return evt_names[vType];
    }
    return nullptr;
}

VsiteType string2vsiteType(const char *string)
{
    int i;
    for (i = 0; i < evtNR; i++)
    {
        if (gmx_strcasecmp(string, evt_names[i]) == 0)
        {
            return static_cast<VsiteType>(i);
        }
    }
    return evtNR;
}

CommunicationStatus Ffatype::Send(const t_commrec *cr, int dest)
{
    CommunicationStatus cs;
    cs = gmx_send_data(cr, dest);
    if (CS_OK == cs)
    {
        gmx_send_str(cr, dest, &desc_);
        gmx_send_str(cr, dest, &type_);
        gmx_send_int(cr, dest, subType_.size());
        for(auto &it : subType_)
        {
            std::string key = interactionTypeToString(it.first);
            gmx_send_str(cr, dest, &key);
            gmx_send_str(cr, dest, &it.second);
        }
        gmx_send_str(cr, dest, &elem_);
        gmx_send_double(cr, dest, mass_);
        gmx_send_int(cr, dest, atomnumber_);
        gmx_send_double(cr, dest, charge_);
        gmx_send_int(cr, dest, row_);
        std::string mut = mutabilityName(mutability_);
        gmx_send_str(cr, dest, &mut);
        gmx_send_str(cr, dest, &refEnthalpy_);
        if (nullptr != debug)
        {
            fprintf(debug, "Sent Fftype %s %s %s %s %s %g %g %d %s %s, status %s\n",
                    desc_.c_str(), type_.c_str(), 
                    subType_[InteractionType::BONDS].c_str(),
                    subType_[InteractionType::POLARIZATION].c_str(),
                    elem_.c_str(), mass_, charge_, row_, mut.c_str(),
                    refEnthalpy_.c_str(), cs_name(cs));
            fflush(debug);
        }
    }
    return cs;
}

CommunicationStatus Ffatype::Receive(const t_commrec *cr, int src)
{
    CommunicationStatus cs;
    cs = gmx_recv_data(cr, src);
    if (CS_OK == cs)
    {
        gmx_recv_str(cr, src, &desc_);
        gmx_recv_str(cr, src, &type_);
        int nsub = gmx_recv_int(cr, src);
        subType_.clear();
        for(int i = 0; i < nsub; i++)
        {
            std::string key, value;
            gmx_recv_str(cr, src, &key);
            gmx_recv_str(cr, src, &value);
            InteractionType iType = stringToInteractionType(key);
            subType_.insert({iType, value});
        }
        gmx_recv_str(cr, src, &elem_);
        mass_ = gmx_recv_double(cr, src);
        atomnumber_ = gmx_recv_int(cr, src);
        charge_ = gmx_recv_double(cr, src);
        row_ = gmx_recv_int(cr, src);
        std::string mut;
        gmx_recv_str(cr, src, &mut);
        if (!nameToMutability(mut, &mutability_))
        {
            GMX_THROW(gmx::InternalError(gmx::formatString("Should not happen but did. mut = %s", mut.c_str()).c_str()));
        }
        gmx_recv_str(cr, src, &refEnthalpy_);

        if (nullptr != debug)
        {
            fprintf(debug, "Received Fftype %s %s %s %s %s %s %g %g %d %s, status %s\n",
                    desc_.c_str(), type_.c_str(), 
                    subType_[InteractionType::BONDS].c_str(),
                    subType_[InteractionType::POLARIZATION].c_str(),
                    elem_.c_str(), mass_, charge_, row_, mut.c_str(),
                    refEnthalpy_.c_str(), cs_name(cs));
            fflush(debug);
        }
    }
    return cs;
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
    std::string         central, attached;
    cs = gmx_recv_data(cr, src);
    if (CS_OK == cs)
    {
        gmx_recv_str(cr, src, &central);
        gmx_recv_str(cr, src, &attached);
        const_cast<std::string &>(central_)   = central;
        const_cast<std::string &>(attached_)  = attached;
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

/*
 * This source file is part of the Alexandria Chemistry Toolkit.
 *
 * Copyright (C) 2014-2021
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

#include "composition.h"

#include <cstdio>

#include <string>
#include <vector>

#include "communication.h"
#include "gmx_simple_comm.h"
#include "units.h"

namespace alexandria
{

void CalcAtom::SetUnit(const std::string &unit)
{
    if ((unit_.size() == 0) && (unit.size() > 0))
    {
        unit_ = unit;
    }
    else
    {
        if (unit_.size() == 0)
        {
            fprintf(stderr, "Replacing CalcAtom unit '%s' by '%s'\n", unit_.c_str(), unit.c_str());
        }
    }
}

bool CalcAtom::Equal(CalcAtom ca)
{
    return !((name_.compare(ca.getName()) != 0) ||
             (obType_.compare(ca.getObtype()) != 0) ||
             (x_ != ca.getX()) ||
             (y_ != ca.getY()) ||
             (z_ != ca.getZ()) ||
             (atomID_ != ca.getAtomid()));
}

CommunicationStatus CalcAtom::Receive(t_commrec *cr, int src)
{
    CommunicationStatus cs;
    int                 Ncharge;

    cs = gmx_recv_data(cr, src);
    if (CS_OK == cs)
    {
        gmx_recv_str(cr, src, &name_);
        gmx_recv_str(cr, src, &obType_);
        gmx_recv_str(cr, src, &residueName_);
        residueNumber_ = gmx_recv_int(cr, src);
        atomID_ = gmx_recv_int(cr, src);
        gmx_recv_str(cr, src, &unit_);
        x_      = gmx_recv_double(cr, src);
        y_      = gmx_recv_double(cr, src);
        z_      = gmx_recv_double(cr, src);
        Ncharge = gmx_recv_int(cr, src);

        for (int n = 0; (CS_OK == cs) && (n < Ncharge); n++)
        {
            std::string type;
            gmx_recv_str(cr, src, &type);
            double q = gmx_recv_double(cr, src);
            AddCharge(stringToQtype(type), q);
        }
    }
    if (nullptr != debug)
    {
        fprintf(debug, "Received CalcAtom, status %s\n", cs_name(cs));
     }
    return cs;
}

CommunicationStatus CalcAtom::Send(t_commrec *cr, int dest) const
{
    CommunicationStatus  cs;

    cs = gmx_send_data(cr, dest);
    if (CS_OK == cs)
    {
        gmx_send_str(cr, dest, &name_);
        gmx_send_str(cr, dest, &obType_);
        gmx_send_str(cr, dest, &residueName_);
        gmx_send_int(cr, dest, residueNumber_);
        gmx_send_int(cr, dest, atomID_);
        gmx_send_str(cr, dest, &unit_);
        gmx_send_double(cr, dest, x_);
        gmx_send_double(cr, dest, y_);
        gmx_send_double(cr, dest, z_);
        gmx_send_int(cr, dest, q_.size());

        for (const auto &qi : q_)
        {
            std::string type = qTypeName(qi.first);
            gmx_send_str(cr, dest, &type);
            gmx_send_double(cr, dest, qi.second);
        }
    }
    if (nullptr != debug)
    {
        fprintf(debug, "Sent CalcAtom, status %s\n", cs_name(cs));
    }
    return cs;
}

CommunicationStatus AtomNum::Send(t_commrec *cr, int dest) const
{
    CommunicationStatus cs;

    cs = gmx_send_data(cr, dest);
    if (CS_OK == cs)
    {
        gmx_send_str(cr, dest, &catom_);
        gmx_send_int(cr, dest, cnumber_);
        if (nullptr != debug)
        {
            fprintf(debug, "Sent AtomNum %s %d, status %s\n",
                    catom_.c_str(), cnumber_, cs_name(cs));
            fflush(debug);
        }
    }
    return cs;
}

CommunicationStatus AtomNum::Receive(t_commrec *cr, int src)
{
    CommunicationStatus cs;

    cs = gmx_recv_data(cr, src);
    if (CS_OK == cs)
    {
        gmx_recv_str(cr, src, &catom_);
        cnumber_ = gmx_recv_int(cr, src);
        if (nullptr != debug)
        {
            fprintf(debug, "Received AtomNum %s %d, status %s\n",
                    catom_.c_str(), cnumber_, cs_name(cs));
            fflush(debug);
        }
    }
    return cs;
}

void MolecularComposition::AddAtom(AtomNum an)
{
    AtomNumIterator mci = searchAtom(an.getAtom());
    if (mci == atomnum_.end())
    {
        atomnum_.push_back(an);
    }
    else
    {
        mci->SetNumber(mci->getNumber()+an.getNumber());
    }
}

void MolecularComposition::DeleteAtom(const std::string &catom)
{
    AtomNumIterator ani;

    if ((ani = searchAtom(catom)) != atomnum_.end())
    {
        atomnum_.erase(ani);
    }
}

AtomNumConstIterator MolecularComposition::searchAtomConst(const std::string &an) const
{
    for (auto ani = atomnum_.begin(); ani < atomnum_.end(); ++ani)
    {
        if (an.compare(ani->getAtom()) == 0)
        {
            return ani;
        }
    }
    return atomnum_.end();
}

AtomNumIterator MolecularComposition::searchAtom(const std::string &an)
{
    for (auto ani = atomnum_.begin(); ani < atomnum_.end(); ++ani)
    {
        if (an.compare(ani->getAtom()) == 0)
        {
            return ani;
        }
    }
    return atomnum_.end();
}

void MolecularComposition::ReplaceAtom(const std::string &oldatom,
                                       const std::string &newatom)
{

    for (auto &i : atomnum_)
    {
        if (oldatom.compare(i.getAtom()) == 0)
        {
            i.SetAtom(newatom);
            break;
        }
    }
}

int MolecularComposition::CountAtoms(const std::string &atom) const
{
    for (auto &i : atomnum_)
    {
        if (atom.compare(i.getAtom()) == 0)
        {
            return i.getNumber();
        }
    }
    return 0;
}

int MolecularComposition::CountAtoms() const
{
    int             nat = 0;

    for (auto &i : atomnum_)
    {
        nat += i.getNumber();
    }
    return nat;
}

CommunicationStatus MolecularComposition::Send(t_commrec *cr, int dest) const
{
    CommunicationStatus cs = gmx_send_data(cr, dest);
    if (CS_OK == cs)
    {
        gmx_send_int(cr, dest, atomnum_.size());
        gmx_send_str(cr, dest, &compname_);
        for (auto &ani : atomnum_)
        {
            cs = ani.Send(cr, dest);
            if (CS_OK != cs)
            {
                break;
            }
        }
        if (nullptr != debug)
        {
            fprintf(debug, "Sent MolecularComposition %s, status %s\n",
                    compname_.c_str(), cs_name(cs));
            fflush(debug);
        }
    }
    return cs;
}

CommunicationStatus MolecularComposition::Receive(t_commrec *cr, int src)
{
    int                 Natomnum;
    CommunicationStatus cs = gmx_recv_data(cr, src);
    if (CS_OK == cs)
    {
        Natomnum = gmx_recv_int(cr, src);
        gmx_recv_str(cr, src, &compname_);
        CommunicationStatus cs2;
        for (int n = 0; n < Natomnum; n++)
        {
            AtomNum an;
            cs2 = an.Receive(cr, src);
            if (CS_OK == cs2)
            {
                AddAtom(an);
            }
        }
        if (nullptr != debug)
        {
            fprintf(debug, "Received MolecularComposition %s, status %s\n",
                    compname_.c_str(), cs_name(cs));
            fflush(debug);
        }
    }
    return cs;
}


} // namespace alexandria

/*
 * This source file is part of the Alexandria Chemistry Toolkit.
 *
 * Copyright (C) 2014-2022
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

#include "communicationrecord.h"
#include "utility/units.h"

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

CommunicationStatus CalcAtom::Receive(const CommunicationRecord *cr, int src)
{
    CommunicationStatus cs       = CommunicationStatus::OK;
    int                 Ncharge;

    if (CommunicationStatus::RECV_DATA == cr->recv_data(src))
    {
        cr->recv_str(src, &name_);
        cr->recv_str(src, &obType_);
        cr->recv_str(src, &residueName_);
        residueNumber_ = cr->recv_int(src);
        atomID_ = cr->recv_int(src);
        cr->recv_str(src, &unit_);
        x_      = cr->recv_double(src);
        y_      = cr->recv_double(src);
        z_      = cr->recv_double(src);
        Ncharge = cr->recv_int(src);

        for (int n = 0; (CommunicationStatus::OK == cs) && (n < Ncharge); n++)
        {
            std::string type;
            cr->recv_str(src, &type);
            double q = cr->recv_double(src);
            AddCharge(stringToQtype(type), q);
        }
    }
    if (nullptr != debug)
    {
        fprintf(debug, "Received CalcAtom, status %s\n", cs_name(cs));
     }
    return cs;
}

CommunicationStatus CalcAtom::Send(const CommunicationRecord *cr, int dest) const
{
    CommunicationStatus  cs = CommunicationStatus::OK;

    if (CommunicationStatus::SEND_DATA == cr->send_data(dest))
    {
        cr->send_str(dest, &name_);
        cr->send_str(dest, &obType_);
        cr->send_str(dest, &residueName_);
        cr->send_int(dest, residueNumber_);
        cr->send_int(dest, atomID_);
        cr->send_str(dest, &unit_);
        cr->send_double(dest, x_);
        cr->send_double(dest, y_);
        cr->send_double(dest, z_);
        cr->send_int(dest, q_.size());

        for (const auto &qi : q_)
        {
            std::string type = qTypeName(qi.first);
            cr->send_str(dest, &type);
            cr->send_double(dest, qi.second);
        }
    }
    if (nullptr != debug)
    {
        fprintf(debug, "Sent CalcAtom, status %s\n", cs_name(cs));
    }
    return cs;
}

CommunicationStatus AtomNum::Send(const CommunicationRecord *cr, int dest) const
{
    CommunicationStatus cs = CommunicationStatus::OK;

    if (CommunicationStatus::SEND_DATA == cr->send_data(dest))
    {
        cr->send_str(dest, &catom_);
        cr->send_int(dest, cnumber_);
        if (nullptr != debug)
        {
            fprintf(debug, "Sent AtomNum %s %d, status %s\n",
                    catom_.c_str(), cnumber_, cs_name(cs));
            fflush(debug);
        }
    }
    return cs;
}

CommunicationStatus AtomNum::Receive(const CommunicationRecord *cr, int src)
{
    CommunicationStatus cs = CommunicationStatus::OK;

    if (CommunicationStatus::RECV_DATA == cr->recv_data(src))
    {
        cr->recv_str(src, &catom_);
        cnumber_ = cr->recv_int(src);
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

CommunicationStatus MolecularComposition::Send(const CommunicationRecord *cr, int dest) const
{
    CommunicationStatus cs = CommunicationStatus::OK;
    if (CommunicationStatus::SEND_DATA == cr->send_data(dest))
    {
        cr->send_int(dest, atomnum_.size());
        cr->send_str(dest, &compname_);
        for (auto &ani : atomnum_)
        {
            cs = ani.Send(cr, dest);
            if (CommunicationStatus::OK != cs)
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

CommunicationStatus MolecularComposition::Receive(const CommunicationRecord *cr, int src)
{
    int                 Natomnum;
    CommunicationStatus cs = CommunicationStatus::OK;
    if (CommunicationStatus::RECV_DATA == cr->recv_data(src))
    {
        Natomnum = cr->recv_int(src);
        cr->recv_str(src, &compname_);
        CommunicationStatus cs2;
        for (int n = 0; n < Natomnum; n++)
        {
            AtomNum an;
            cs2 = an.Receive(cr, src);
            if (CommunicationStatus::OK == cs2)
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

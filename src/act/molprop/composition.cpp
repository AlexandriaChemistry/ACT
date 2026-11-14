/*
 * This source file is part of the Alexandria Chemistry Toolkit.
 *
 * Copyright (C) 2014-2026
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
 * \author Mohammad Mehdi Ghahremanpour <mohammad.ghahremanpour@icm.uu.se>
 * \author David van der Spoel <david.vanderspoel@icm.uu.se>
 */

#include "composition.h"

#include <cstdio>

#include <string>
#include <vector>

#include "act/utility/communicationrecord.h"
#include "act/utility/units.h"
#include "gromacs/utility/stringutil.h"

namespace alexandria
{

void CalcAtom::setCoordUnit(const std::string &unit)
{
    if ((coord_unit_.empty()) && (!unit.empty()))
    {
        coord_unit_ = unit;
    }
    else
    {
        if (coord_unit_.empty())
        {
            fprintf(stderr, "Trying to replace CalcAtom unit '%s' by '%s'\n", coord_unit_.c_str(), unit.c_str());
        }
    }
}

void CalcAtom::setForceUnit(const std::string &unit)
{
    force_unit_ = unit;
}

bool CalcAtom::Equal(const CalcAtom &ca)
{
    return !((name_.compare(ca.getName()) != 0) ||
             (obType_.compare(ca.getObtype()) != 0) ||
             (x_ != ca.getX()) ||
             (y_ != ca.getY()) ||
             (z_ != ca.getZ()) ||
             (atomID_ != ca.getAtomid()));
}

CommunicationStatus CalcAtom::BroadCast(const CommunicationRecord *cr,
                                        int                        root,
                                        MPI_Comm                   comm)
{
    CommunicationStatus cs = cr->bcast_data(comm);

    if (CommunicationStatus::OK == cs)
    {
        cr->bcast(&name_, comm);
        cr->bcast(&obType_, comm);
        cr->bcast(&residueName_, comm);
        cr->bcast(&residueNumber_, comm);
        cr->bcast(&atomID_, comm);
        cr->bcast(&coord_unit_, comm);
        cr->bcast(&x_, comm);
        cr->bcast(&y_, comm);
        cr->bcast(&z_, comm);
        cr->bcast(&force_unit_, comm);
        cr->bcast(&fx_, comm);
        cr->bcast(&fy_, comm);
        cr->bcast(&fz_, comm);
        size_t Ncharge = q_.size();
        cr->bcast(&Ncharge, comm);
        if (cr->rank() == root)
        {
            for (auto &qi : q_)
            {
                std::string type(qi.first);
                cr->bcast(&type, comm);
                cr->bcast(&qi.second, comm);
            }
        }
        else
        {
            for (size_t n = 0; (CommunicationStatus::OK == cs) && (n < Ncharge); n++)
            {
                std::string type;
                cr->bcast(&type, comm);
                double q;
                cr->bcast(&q, comm);
                AddCharge(type, q);
            }
        }
    }
    if (cr->mh())
    {
        cr->mh()->writeDebug(gmx::formatString("Received CalcAtom, status %s\n", cs_name(cs).c_str()));
    }
    return cs;
}

void MolecularComposition::AddAtom(const AtomNum &an)
{
    const AtomNumIterator mci = searchAtom(an.getAtom());
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

} // namespace alexandria

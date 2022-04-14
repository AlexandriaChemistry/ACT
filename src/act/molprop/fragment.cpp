/*
 * This source file is part of the Alexandria Chemistry Toolkit.
 *
 * Copyright (C) 2022
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

#include "fragment.h"

#include "gromacs/utility/stringutil.h"

namespace alexandria
{

void Fragment::makeAtomString()
{
    atomString_.clear();
    for(auto &a : atoms_)
    {
        atomString_ += gmx::formatString(" %d", a+1);
    }
}

void Fragment::makeTexFormula()
{
    texform_.clear();
    for(const auto c : formula_)
    {
        if (isdigit(c))
        {
            texform_ += gmx::formatString("$_%c$", c);
        }
        else
        {
            texform_ += c;
        }
    }
}

void Fragment::dump(FILE *fp) const
{
    if (!fp)
    {
        return;
    }
    fprintf(fp, "Fragment %s mass %g formula %s charge %d multiplicity %d atoms %s\n",
            id_.c_str(), mass_, formula_.c_str(), charge_, multiplicity_,
            atomString_.c_str());
}
    
CommunicationStatus Fragment::Receive(const CommunicationRecord *cr, int src)
{
    CommunicationStatus cs = CommunicationStatus::OK;
    mass_         = cr->recv_double(src);
    charge_       = cr->recv_int(src);
    multiplicity_ = cr->recv_int(src);
    cr->recv_str(src, &formula_);
    cr->recv_str(src, &texform_);
    cr->recv_str(src, &id_);
    int natom     = cr->recv_int(src);
    atoms_.clear();
    for(int i = 0; i < natom; i++)
    {
        atoms_.push_back(cr->recv_int(src));
    }
    return cs;
}

CommunicationStatus Fragment::Send(const CommunicationRecord *cr, int dest) const
{
    CommunicationStatus cs = CommunicationStatus::OK;
    cr->send_double(dest, mass_);
    cr->send_int(dest, charge_);
    cr->send_int(dest, multiplicity_);
    cr->send_str(dest, &formula_);
    cr->send_str(dest, &texform_);
    cr->send_str(dest, &id_);
    cr->send_int(dest, atoms_.size());
    for(auto &a : atoms_)
    {
        cr->send_int(dest, a);
    }
    return cs;
}

} // namespace alexandria

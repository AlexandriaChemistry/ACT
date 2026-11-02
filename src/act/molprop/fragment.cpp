/*
 * This source file is part of the Alexandria Chemistry Toolkit.
 *
 * Copyright (C) 2022-2024
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
 * \author David van der Spoel <david.vanderspoel@icm.uu.se>
 */

#include "fragment.h"

#include "gromacs/utility/stringutil.h"
#include "gromacs/utility/textwriter.h"

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

void Fragment::dump(gmx::TextWriter *tw) const
{
    if (tw)
    {
        tw->writeStringFormatted("Fragment %s mass %g formula %s charge %d multiplicity %d atoms %s\n",
                                 inchi_.c_str(), mass_, formula_.c_str(), charge_, multiplicity_,
                                 atomString_.c_str());
    }
}
    
CommunicationStatus Fragment::BroadCast(const CommunicationRecord *cr,
                                        int                        root,
                                        MPI_Comm                   comm)
{
    CommunicationStatus cs = cr->bcast_data(comm);
    if (CommunicationStatus::OK == cs)
    {
        cr->bcast(&mass_, comm);
        cr->bcast(&charge_, comm);
        cr->bcast(&multiplicity_, comm);
        cr->bcast(&symmetryNumber_, comm);
        cr->bcast(&formula_, comm);
        cr->bcast(&texform_, comm);
        cr->bcast(&inchi_, comm);
        cr->bcast(&iupac_, comm);
        int natom     = atoms_.size();
        cr->bcast(&natom, comm);
        if (cr->rank() != root)
        {
            atoms_.resize(natom);
        }
        for(int i = 0; i < natom; i++)
        {
            cr->bcast(&atoms_[i], comm);
        }
    }
    if (cr->rank() != root)
    {
        makeAtomString();
    }
    return cs;
}

} // namespace alexandria

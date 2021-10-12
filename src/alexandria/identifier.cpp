/*
 * This source file is part of the Alexandria Chemistry Toolkit.
 *
 * Copyright (C) 2020 
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

#include "identifier.h"

#include <map>
#include <string>
#include <vector>

#include "gromacs/utility/exceptions.h"
#include "gromacs/utility/stringutil.h"

#include "communication.h"
#include "gmx_simple_comm.h"
#include "stringutil.h"

namespace alexandria
{

std::map<CanSwap, std::string> cs2string = 
    {
        { CanSwap::No,       "false"    },
        { CanSwap::Yes,      "true"     }
    };
    
CanSwap stringToCanSwap(const std::string &str)
{
    for(auto &cs : cs2string)
    {
        if (cs.second == str)
        {
            return cs.first;
        }
    }
    GMX_THROW(gmx::InvalidInputError(gmx::formatString("The value '%s' cannot be converted to a CanSwap value", str.c_str()).c_str()));
}

const std::string &canSwapToString(CanSwap canSwap)
{
    return cs2string[canSwap];
}

static const char *IdDelimeter = "#";

Identifier::Identifier(const std::vector<std::string> &atoms,
                       CanSwap                         canSwap)
{
    for(size_t k = 0; k < atoms.size(); k++)
    {
        id_.append(atoms[k]);
        if (k < atoms.size()-1)
        {
            id_.append(IdDelimeter);
        }
    }
    if (canSwap == CanSwap::Yes)
    {
        for(size_t k = 0; k < atoms.size(); k++)
        {
            swappedId_.append(atoms[atoms.size()-k-1]);
            if (k < atoms.size()-1)
            {
                swappedId_.append(IdDelimeter);
            }
        }
    }
    orderAtoms();
    atoms_ = atoms;
}

void Identifier::orderAtoms()
{
    if (swappedId_.size() > 0 && swappedId_ < id_)
    {
        std::string tmp = id_;
        id_             = swappedId_;
        swappedId_      = tmp;
    }
}

Identifier::Identifier(InteractionType    iType,
                       const std::string &id,
                       CanSwap            canSwap)
{
    id_        = id;
    auto ids   = split(id, IdDelimeter[0]);
    auto idEnd = ids.end();
    if (iType == InteractionType::BONDS ||
        iType == InteractionType::BONDCORRECTIONS)
    {
        idEnd--;
        bondOrder_ = atof(ids[ids.size()-1].c_str());
    }
    else
    {
        bondOrder_ = 0;
    }
    atoms_.clear();
    for (auto i = ids.begin(); i < idEnd; ++i)
    {
        atoms_.push_back(*i);
    }
    if (canSwap == CanSwap::Yes)
    {
        if (iType == InteractionType::BONDS ||
            iType == InteractionType::BONDCORRECTIONS)
        {
            swappedId_ = gmx::formatString("%s%s%s%s%g", ids[1].c_str(), IdDelimeter,
                                           ids[0].c_str(), IdDelimeter, bondOrder_);
        }
        else
        {
            swappedId_.clear();
            for(auto i = ids.rbegin(); i < ids.rend(); ++i)
            {
                swappedId_.append(*i);
                if (i < ids.rend()-1)
                {
                    swappedId_.append(IdDelimeter);
                }
            }
        }
    }
    orderAtoms();
}

Identifier::Identifier(const std::vector<std::string> &atoms,
                       double                          bondOrder,
                       CanSwap                         canSwap)
{
    GMX_RELEASE_ASSERT(atoms.size() == 2,
                       "Bonds should have two atoms exactly");
    id_.assign(gmx::formatString("%s%s%s%s%g",
                                 atoms[0].c_str(), IdDelimeter,
                                 atoms[1].c_str(), IdDelimeter,
                                 bondOrder));
    if (canSwap == CanSwap::Yes)
    {
        swappedId_.assign(gmx::formatString("%s%s%s%s%g",
                                            atoms[1].c_str(), IdDelimeter,
                                            atoms[0].c_str(), IdDelimeter,
                                            bondOrder));
    }
    atoms_     = atoms;
    if (bondOrder < 1 || bondOrder > 3)
    {
        GMX_THROW(gmx::InvalidInputError(gmx::formatString("Bond order should be 1, 2 or 3, not %g", bondOrder).c_str()));
    }
    bondOrder_ = bondOrder;
    orderAtoms();
}

CommunicationStatus Identifier::Send(const t_commrec *cr, int dest) const
{
    gmx_send_str(cr, dest, &id_);
    gmx_send_str(cr, dest, &swappedId_);
    gmx_send_double(cr, dest, bondOrder_);
    gmx_send_int(cr, dest, static_cast<int>(atoms_.size()));
    for(auto &a : atoms_)
    {
        gmx_send_str(cr, dest, &a);
    }
    return CS_OK;
}

CommunicationStatus Identifier::Receive(const t_commrec *cr, int src)
{
    gmx_recv_str(cr, src, &id_);
    gmx_recv_str(cr, src, &swappedId_);
    bondOrder_ = gmx_recv_double(cr, src);
    int natoms = gmx_recv_int(cr, src);
    atoms_.clear();
    for(int i = 0; i < natoms; i++)
    {
        std::string a;
        gmx_recv_str(cr, src, &a);
        atoms_.push_back(a);
    }
    return CS_OK;
}

bool operator==(const Identifier &a, const Identifier &b)
{
#ifdef DEBUG
    if (a.atoms().size() != b.atoms().size())
    {
        GMX_THROW(gmx::InternalError(gmx::formatString("Comparing identifiers of different type with %d resp. %d components", static_cast<int>(a.atoms().size()), static_cast<int>(b.atoms().size())).c_str()));
    }
#endif
    return (a.id() == b.id());
}

bool operator<(const Identifier &a, const Identifier &b)
{
    if (a == b)
    {
        return false;
    }
    else
    {
        return (a.id() < b.id());
    }
}

} // namespace alexandria

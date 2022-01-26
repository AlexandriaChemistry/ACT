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
#include "gromacs/utility/fatalerror.h"
#include "gromacs/utility/stringutil.h"

#include "communicationrecord.h"
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

/*! To distinguish bond orders in identifiers
 * These are the SMARTS Bond Primitives defined here
 * https://www.daylight.com/dayhtml/doc/theory/theory.smarts.html
 * 1.5 refers to aromatic bonds but also resonant systems, such
 * as a carboxylate.
 */
static std::map<double, char> BondOrderDelimeter =
    { { 1, '~' }, { 1.5, ':' }, { 2, '=' }, { 3, '#' } };

/*! Reverse map of the above.
 * Will be created on first use
 */
static std::map<char, double> ReverseBondOrderDelimeter;

static void createReverseBondOrderDelimeter()
{
    if (ReverseBondOrderDelimeter.empty())
    {
        for (auto &bod : BondOrderDelimeter)
        {
            ReverseBondOrderDelimeter.insert(std::pair<char, double>(bod.second, bod.first));
        }
    }
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


void Identifier::createSwapped(CanSwap canSwap)
{
    if (canSwap == CanSwap::Yes)
    {
        swappedId_ = atoms_[atoms_.size()-1];
        for(size_t i = atoms_.size()-1; i > 0; i--)
        {
            swappedId_ += BondOrderDelimeter[bondOrders_[i-1]];
            swappedId_ += atoms_[i-1];
        }
    }
}

Identifier::Identifier(const std::string &atom)
{
    atoms_.push_back(atom);
    id_ = atom;
}

Identifier::Identifier(InteractionType    iType,
                       const std::string &id,
                       CanSwap            canSwap)
{
    if (debug)
    {
        fprintf(debug, "Trying to create %s identifier from '%s'\n", interactionTypeToString(iType).c_str(), id.c_str());
    }
    id_       = id;
    createReverseBondOrderDelimeter();
    size_t c0 = 0;
    for (size_t c = 0; c < id.size(); c++)
    {
        auto ptr = ReverseBondOrderDelimeter.find(id[c]);
        if (ptr != ReverseBondOrderDelimeter.end())
        {
            atoms_.push_back(id.substr(c0,c-c0));
            bondOrders_.push_back(ptr->second);
            c0 = c+1;
        }
    }
    if (c0 < id.size())
    {
        atoms_.push_back(id.substr(c0, id.size()-c0));
    }
    int  natoms = interactionTypeToNatoms(iType);
    // For natoms atoms we should have natoms-1 bond orders
    // which means 2*natoms - 1 items in ids
    if (static_cast<int>(atoms_.size()) != natoms)
    {
        GMX_THROW(gmx::InvalidInputError(gmx::formatString("Expected %d atoms but found %d. Id = %s", natoms, static_cast<int>(atoms_.size()), id.c_str()).c_str()));
    }
    if (static_cast<int>(bondOrders_.size()) != natoms-1)
    {
        GMX_THROW(gmx::InvalidInputError(gmx::formatString("Expected %d bondOrders but found %d. Id = %s", natoms-1, static_cast<int>(bondOrders_.size()), id.c_str()).c_str()));
    }
    createSwapped(canSwap);
    orderAtoms();
}

Identifier::Identifier(const std::vector<std::string> &atoms,
                       const std::vector<double>      &bondOrders,
                       CanSwap                         canSwap)
{
    if (bondOrders.size()+1 != atoms.size())
    {
        GMX_THROW(gmx::InvalidInputError(gmx::formatString("Expecting %d bond orders for %d atoms, but got %d", static_cast<int>(atoms.size()-1), static_cast<int>(atoms.size()), static_cast<int>(bondOrders.size())).c_str()));
    }
    id_ = atoms[0];
    for(size_t i = 1; i < atoms.size(); i++)
    {
        if (BondOrderDelimeter.find(bondOrders[i-1]) == BondOrderDelimeter.end())
        {
            gmx_fatal(FARGS, "Don't understand a bond order of %g", bondOrders[i-1]);
        }
        if (!atoms[i].empty())
        {
            id_ += BondOrderDelimeter[bondOrders[i-1]];
            id_ += atoms[i];
        }
    }
    atoms_      = atoms;
    bondOrders_ = bondOrders;
    createSwapped(canSwap);
    orderAtoms();
}

CommunicationStatus Identifier::Send(const CommunicationRecord *cr, int dest) const
{
    cr->send_str(dest, &id_);
    cr->send_str(dest, &swappedId_);
    cr->send_int(dest, static_cast<int>(atoms_.size()));
    for(auto &a : atoms_)
    {
        cr->send_str(dest, &a);
    }
    cr->send_int(dest, static_cast<int>(bondOrders_.size()));
    for(auto &b : bondOrders_)
    {
        cr->send_double(dest, b);
    }
    return CS_OK;
}

CommunicationStatus Identifier::Receive(const CommunicationRecord *cr, int src)
{
    cr->recv_str(src, &id_);
    cr->recv_str(src, &swappedId_);
    int natoms = cr->recv_int(src);
    atoms_.clear();
    for(int i = 0; i < natoms; i++)
    {
        std::string a;
        cr->recv_str(src, &a);
        atoms_.push_back(a);
    }
    int nbo = cr->recv_int(src);
    bondOrders_.clear();
    for(int i = 0; i < nbo; i++)
    {
        bondOrders_.push_back(cr->recv_double(src));
    }
    return CS_OK;
}

bool operator==(const Identifier &a, const Identifier &b)
{
    return (a.id() == b.id());
}

bool operator<(const Identifier &a, const Identifier &b)
{
    return (a.id() < b.id());
}

} // namespace alexandria

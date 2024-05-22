/*
 * This source file is part of the Alexandria Chemistry Toolkit.
 *
 * Copyright (C) 2020-2023
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

#include "identifier.h"

#include <map>
#include <string>
#include <vector>

#include "gromacs/utility/exceptions.h"
#include "gromacs/utility/fatalerror.h"
#include "gromacs/utility/stringutil.h"

#include "act/utility/communicationrecord.h"
#include "act/utility/stringutil.h"

namespace alexandria
{

std::map<CanSwap, std::string> cs2string =
    {
        { CanSwap::No,       "false"    },
        { CanSwap::Yes,      "true"     },
        { CanSwap::Idih,     "idih"     },
        { CanSwap::Linear,   "linear"   },
        { CanSwap::Vsite2,   "vsite2"   },
        { CanSwap::Vsite3,   "vsite3"   }
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
    { { 1, '~' }, { 1.5, ':' }, { 2, '=' }, { 3, '#' }, { 9, '!' } };

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

const std::string &Identifier::id() const
{
    if (ids_.empty())
    {
        GMX_THROW(gmx::InternalError("Identifier without ids"));
    }
    else
    {
        return ids_[0];
    }
}

const std::string &Identifier::swapped() const
{
    if (ids_.size() == 2)
    {
        return ids_[1];
    }
    else
    {
        GMX_THROW(gmx::InternalError(gmx::formatString("Identifier with insufficient #ids (%zu)", ids_.size()).c_str()));
    }
}

void Identifier::orderAtoms()
{
    if ((canSwap_ == CanSwap::Yes && atoms_.size() > 0) ||
        (canSwap_ == CanSwap::Linear && atoms_.size() == 3))
    {
        std::string swapped = atoms_[atoms_.size()-1];
        for(size_t i = atoms_.size()-1; i > 0; i--)
        {
            swapped += BondOrderDelimeter[bondOrders_[i-1]];
            swapped += atoms_[i-1];
        }
        if (swapped >= ids_[0])
        {
            ids_.push_back(swapped);
        }
        else
        {
            // Reverse the order on things
            auto tmpid = ids_[0];
            ids_[0] = swapped;
            ids_.push_back(tmpid);
            auto tmpat = atoms_;
            atoms_.clear();
            for(size_t i = tmpat.size(); i > 0; i--)
            {
                atoms_.push_back(tmpat[i-1]);
            }
            auto tmpbo = bondOrders_;
            bondOrders_.clear();
            for(size_t i = tmpbo.size(); i > 0; i--)
            {
                bondOrders_.push_back(tmpbo[i-1]);
            }
        }
    }
    else if (canSwap_ == CanSwap::Vsite2 && atoms_.size() == 3)
    {
        std::string swapped = atoms_[1] + BondOrderDelimeter[bondOrders_[0]] + atoms_[0] + BondOrderDelimeter[bondOrders_[1]] + atoms_[2];
        ids_.push_back(swapped);
    }
    else if (canSwap_ == CanSwap::Idih && atoms_.size() == 4)
    {
        // For this variant the first (central) atom remains the same, but the rest rotates.
        std::vector<std::vector<int>> swap = {
            { 0, 1, 2, 3 }, { 0, 2, 3, 1 }, { 0, 2, 1, 3 },
            { 0, 3, 1, 2 }, { 0, 3, 2, 1 }, { 0, 1, 3, 2 }
        };
        size_t smallest = 0;
        for(size_t s = 0; s < swap.size(); s++)
        {
            const auto &ss = swap[s];
            std::string swapped = atoms_[0];
            for(size_t i = 1; i < ss.size(); i++)
            {
                swapped += BondOrderDelimeter[bondOrders_[ss[i]-1]];
                swapped += atoms_[ss[i]];
            }
            if (swapped < ids_[smallest])
            {
                smallest = s;
            }
            // TODO Only insert ones we do not already have
            if (s > 0)
            {
                ids_.push_back(swapped);
            }
        }
        if (smallest > 0)
        {
            // Reverse the order on things
            auto tmpid = ids_[0];
            ids_[0] = ids_[smallest];
            ids_[smallest] = tmpid;
            auto tmpat = atoms_;
            atoms_.clear();
            for(size_t i = 0; i < tmpat.size(); i++)
            {
                atoms_.push_back(tmpat[swap[smallest][i]]);
            }
            auto tmpbo = bondOrders_;
            bondOrders_.clear();
            for(size_t i = 0; i < tmpbo.size(); i++)
            {
                auto boIndex = swap[smallest][i+1]-1;
                bondOrders_.push_back(tmpbo[boIndex]);
            }
        }
    }
}

Identifier::Identifier(const std::string &atom)
{
    atoms_.push_back(atom);
    ids_.push_back(atom);
    canSwap_ = CanSwap::Yes;
    update();
}

Identifier::Identifier(InteractionType    iType,
                       const std::string &id,
                       CanSwap            canSwap)
{
    if (debug)
    {
        fprintf(debug, "Trying to create %s identifier from '%s'\n", interactionTypeToString(iType).c_str(), id.c_str());
    }
    ids_.push_back(id);
    canSwap_  = canSwap;
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
    if (InteractionType::VDW != iType && InteractionType::ELECTROSTATICS != iType &&
        InteractionType::VSITE3OUT != iType && InteractionType::VSITE3OUTS != iType)
    {
        // Those InteractionTypes that have a fixed number of atom types
        // are tested here. (VDW and COULOMB both store single atoms and
        // atom pairs in the force field file).
        size_t natoms = interactionTypeToNatoms(iType);
        // For natoms atoms we should have natoms-1 bond orders
        // which means 2*natoms - 1 items in ids
        if (atoms_.size() != natoms)
        {
            GMX_THROW(gmx::InvalidInputError(gmx::formatString("Expected %zu atoms but found %zu. Id = %s", natoms, atoms_.size(), id.c_str()).c_str()));
        }
        if (bondOrders_.size() != natoms-1)
        {
            GMX_THROW(gmx::InvalidInputError(gmx::formatString("Expected %zu bondOrders but found %zu. Id = %s", natoms-1, bondOrders_.size(), id.c_str()).c_str()));
        }
    }
    orderAtoms();
}

void Identifier::update()
{
    std::string newid(atoms_[0]);
    for(size_t i = 1; i < atoms_.size(); i++)
    {
        if (BondOrderDelimeter.find(bondOrders_[i-1]) == BondOrderDelimeter.end())
        {
            gmx_fatal(FARGS, "Don't understand a bond order of %g", bondOrders_[i-1]);
        }
        if (!atoms_[i].empty())
        {
            newid += BondOrderDelimeter[bondOrders_[i-1]];
            newid += atoms_[i];
        }
    }
    if (ids_.empty())
    {
        ids_.push_back(newid);
    }
    else
    {
        ids_[0] = newid;
    }
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
    atoms_      = atoms;
    bondOrders_ = bondOrders;
    canSwap_    = canSwap;
    update();
}

CommunicationStatus Identifier::Send(const CommunicationRecord *cr, int dest) const
{
    cr->send_int(dest, ids_.size());
    for(const auto &ii : ids_)
    {
        cr->send_str(dest, &ii);
    }
    auto tmp = canSwapToString(canSwap_);
    cr->send_str(dest, &tmp);
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
    return CommunicationStatus::OK;
}

CommunicationStatus Identifier::BroadCast(const CommunicationRecord *cr,
                                          int                        root,
                                          MPI_Comm                   comm)
{
    int nids = ids_.size();
    cr->bcast(&nids, comm);
    std::string tmp;
    for(int i = 0; i < nids; i++)
    {
        if (cr->rank() == root)
        {
            tmp.assign(ids_[i]);
        }
        cr->bcast(&tmp, comm);
        if (cr->rank() != root)
        {
            ids_.push_back(tmp);
        }
    }

    if (cr->rank() == root)
    {
        tmp.assign(canSwapToString(canSwap_));
    }
    cr->bcast(&tmp, comm);
    if (cr->rank() != root)
    {
        canSwap_ = stringToCanSwap(tmp);
    }
    int natoms = atoms_.size();
    cr->bcast(&natoms, comm);
    if (cr->rank() == root)
    {
        for(auto a = atoms_.begin(); a < atoms_.end(); ++a)
        {
            std::string aa(*a);
            cr->bcast(&aa, comm);
        }
    }
    else
    {
        atoms_.clear();
        for(int i = 0; i < natoms; i++)
        {
            std::string a;
            cr->bcast(&a, comm);
            atoms_.push_back(a);
        }
    }
    int nbo = bondOrders_.size();
    cr->bcast(&nbo, comm);
    if (cr->rank() != root)
    {
        bondOrders_.resize(nbo);
    }
    cr->bcast(&bondOrders_, comm);

    return CommunicationStatus::OK;
}

CommunicationStatus Identifier::Receive(const CommunicationRecord *cr, int src)
{
    int nids = cr->recv_int(src);
    std::string tmp;
    for(int i = 0; i < nids; i++)
    {
        cr->recv_str(src, &tmp);
        ids_.push_back(tmp);
    }

    cr->recv_str(src, &tmp);
    canSwap_ = stringToCanSwap(tmp);
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
    return CommunicationStatus::OK;
}

} // namespace alexandria

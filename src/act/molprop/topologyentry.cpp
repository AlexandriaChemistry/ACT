/*
 * This source file is part of the Alexandria Chemistry Toolkit.
 *
 * Copyright (C) 2021-2026
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

#include "topologyentry.h"

#include <mutex>
#include <set>
#include <tuple>
#include <vector>

#include "act/basics/interactiontype.h"
#include "act/basics/msg_handler.h"
#include "act/forcefield/forcefield_parametername.h"

namespace alexandria
{

void TopologyEntry::check(size_t nAtom) const
{
    if (indices_.size() != nAtom)
    {
        GMX_THROW(gmx::InternalError(gmx::formatString("Expected %lu atom indices, found %lu",
                                                       nAtom, indices_.size()).c_str()));
    }
}

void TopologyEntry::renumberAtoms(const std::vector<int> &renumber)
{
    for(int &a : indices_)
    {
        if (a < static_cast<int>(renumber.size()))
        {
            a = renumber[a];
        }
        else
        {
            printf("Warning atom %d not in renumber array\n", a);
        }
    }
}

void TopologyEntry::setBondOrder(size_t ai, double bo)
{
    if (ai >= bondOrder_.size())
    {
        GMX_THROW(gmx::InvalidInputError("Incorrect index"));
    }
    bondOrder_[ai] = bo;
}

AtomPair AtomPair::swap() const
{
    AtomPair te;
    auto &indices = atomIndices();
    for (size_t i = 0; i < indices.size(); i++)
    {
        te.addAtom(indices[indices.size()-1-i]);
    }
    return te;
}

Bond Bond::swap() const
{
    Bond te;
    auto &indices = atomIndices();
    for (size_t i = 0; i < indices.size(); i++)
    {
        te.addAtom(indices[indices.size()-1-i]);
    }
    auto &bondOrder = bondOrders();
    for (size_t i = 0; i < bondOrder.size(); i++)
    {
        te.addBondOrder(bondOrder[bondOrder.size()-1-i]);
    }
    te.setGromacsType(gromacsType());
    te.setId(id());
    return te;
}

void Angle::setBonds(const Bond &bij, const Bond &bjk)
{
    if (bij.aJ() != bjk.aI())
    {
        GMX_THROW(gmx::InvalidInputError("Incorrect bonds passed"));
    }
    addAtom(bij.aI());
    addAtom(bij.aJ());
    addAtom(bjk.aJ());
    addBondOrder(bij.bondOrder());
    addBondOrder(bjk.bondOrder());
}

const Identifier *sharedIdentifier(const Identifier *id)
{
    static const Identifier empty;
    if (nullptr == id)
    {
        return &empty;
    }
    // All of the identifier, unlike its comparison operator, which takes
    // an identifier and its swapped one to be the same
    struct Exact
    {
        bool operator()(const Identifier &a, const Identifier &b) const
        {
            return std::forward_as_tuple(a.canSwap(), a.ids(), a.atoms(), a.bondOrders()) <
                std::forward_as_tuple(b.canSwap(), b.ids(), b.atoms(), b.bondOrders());
        }
    };
    static std::set<Identifier, Exact> shared;
    static std::mutex                   mutex;
    const std::lock_guard<std::mutex>   lock(mutex);
    return &*shared.insert(*id).first;
}

CommunicationStatus TopologyEntry::BroadCast(const CommunicationRecord *cr,
                                             int                        root,
                                             MPI_Comm                   comm)
{
    CommunicationStatus cs = cr->bcast_data(comm);

    if (CommunicationStatus::OK == cs)
    {
        int nai = indices_.size();
        cr->bcast(&nai, comm);
        if (cr->rank() == root)
        {
            for (int i= 0; i < nai; i++)
            {
                cr->bcast(&indices_[i], comm);
            }
        }
        else
        {
            indices_.resize(nai, 0);
            for (int i = 0; i < nai; i++)
            {
                cr->bcast(&indices_[i], comm);
            }
        }
        int nbo = bondOrder_.size();
        cr->bcast(&nbo, comm);
        if (cr->rank() != root)
        {
            bondOrder_.resize(nbo);
        }
        cr->bcast(&bondOrder_, comm);
    }
    else if (cr->mh())
    {
        cr->mh()->msg(ACTStatus::Error,
                      gmx::formatString("Failed receiving TopologyEntry, status %s\n", cs_name(cs).c_str()));
    }
    return cs;
}

bool AtomPair::operator==(const AtomPair &other) const
{
    return ((aI() == other.aI() && aJ() == other.aJ()) ||
            (aJ() == other.aI() && aI() == other.aJ()));
}

bool AtomPair::operator<(const AtomPair &other) const
{
    return (aI() < other.aI() ||
            (aI() == other.aI() && aJ() < other.aJ()) ||
            (aI() == other.aJ() && aJ() < other.aI()));
}

bool SingleAtom::operator==(const SingleAtom &other) const
{
    return aI() == other.aI();
}

bool SingleAtom::operator<(const SingleAtom &other) const
{
    return aI() < other.aI();
}

bool Bond::operator==(const Bond &other) const
{
    return ((aI() == other.aI() && aJ() == other.aJ()) ||
            (aJ() == other.aI() && aI() == other.aJ()));
}

void AtomPair::get(int *ai, int *aj) const
{
    check(2);
    *ai        = atomIndex(0);
    *aj        = atomIndex(1);
}

void Bond::get(int *ai, int *aj, double *bondorder) const
{
    check(2);
    *ai        = atomIndex(0);
    *aj        = atomIndex(1);
    *bondorder = bondOrders()[0];
}

void Vsite1::get(int *ai, int *vs) const
{
    check(2);
    *ai = atomIndex(0);
    *vs = atomIndex(1);
}

void Vsite2::get(int *ai, int *aj, int *vs) const
{
    check(3);
    *ai = atomIndex(0);
    *aj = atomIndex(1);
    *vs = atomIndex(2);
}

void Vsite3::get(int *ai, int *aj, int *ak, int *vs) const
{
    check(4);
    *ai = atomIndex(0);
    *aj = atomIndex(1);
    *ak = atomIndex(2);
    *vs = atomIndex(3);
}


void Vsite3OUT::get(int *ai, int *aj, int *ak, int *vs, int *sign) const
{
    check(4);
    *ai = atomIndex(0);
    *aj = atomIndex(1);
    *ak = atomIndex(2);
    *vs = atomIndex(3);
    *sign = sign_;

}



Angle::Angle(const Bond bij, const Bond bjk)
{
    setBonds(bij, bjk);
}

Improper::Improper(const Bond &bij, const Bond &bik, const Bond &bil)
{
    if (bij.aI() != bik.aI())
    {
        GMX_THROW(gmx::InvalidInputError("Incorrect second bond passed"));
    }
    if (bij.aI() != bil.aI())
    {
        GMX_THROW(gmx::InvalidInputError("Incorrect third bond passed"));
    }
    addAtom(bij.aI());
    addAtom(bij.aJ());
    addAtom(bik.aJ());
    addAtom(bil.aJ());
    GMX_RELEASE_ASSERT(atomIndices().size() == 4, "Something weird with impropers");
    addBondOrder(bij.bondOrder());
    addBondOrder(bik.bondOrder());
    addBondOrder(bil.bondOrder());
}

Proper::Proper(const Bond &bij, const Bond &bjk, const Bond &bkl)
{
    if (bij.aJ() != bjk.aI())
    {
        GMX_THROW(gmx::InvalidInputError("Incorrect second bond passed"));
    }
    if (bjk.aJ() != bkl.aI())
    {
        GMX_THROW(gmx::InvalidInputError("Incorrect third bond passed"));
    }
    addAtom(bij.aI());
    addAtom(bij.aJ());
    addAtom(bjk.aJ());
    addAtom(bkl.aJ());
    addBondOrder(bij.bondOrder());
    addBondOrder(bjk.bondOrder());
    addBondOrder(bkl.bondOrder());
}

} // namespace

/*
 * This source file is part of the Alexandria Chemistry Toolkit.
 *
 * Copyright (C) 2021-2022
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

#include "topology.h"

#include <algorithm>
#include <set>
#include <vector>

#include "gromacs/math/vectypes.h"
#include "gromacs/utility/fatalerror.h"

#include "mymol_low.h"

namespace alexandria
{

void TopologyEntry::check(size_t nAtom) const
{
    if (indices_.size() != nAtom)
    {
        GMX_THROW(gmx::InternalError(gmx::formatString("Expected %d atom indices, found %d",
                static_cast<int>(nAtom), static_cast<int>(indices_.size())).c_str()));
    }
}

void TopologyEntry::renumberAtoms(const std::vector<int> &renumber)
{
    for(auto &a : indices_)
    {
        a = renumber[a];
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
    return te;
}

CommunicationStatus TopologyEntry::Send(const CommunicationRecord *cr, int dest) const
{
    CommunicationStatus cs = CommunicationStatus::OK;

    if (CommunicationStatus::SEND_DATA == cr->send_data(dest))
    {
        cr->send_int(dest, indices_.size());
        for(auto &ai : indices_)
        {
            cr->send_int(dest, ai);
        }
        cr->send_int(dest, bondOrder_.size());
        for(auto &bo : bondOrder_)
        {
            cr->send_double(dest, bo);
        }
    }
    else if (nullptr != debug)
    {
        fprintf(debug, "Trying to send TopologyEntry, status %s\n", cs_name(cs));
        fflush(debug);
    }
    return cs;
}

CommunicationStatus TopologyEntry::Receive(const CommunicationRecord *cr, int src)
{
    CommunicationStatus cs = CommunicationStatus::OK;

    if (CommunicationStatus::RECV_DATA == cr->recv_data(src))
    {
        int nai = cr->recv_int(src);
        for (int i=0; i < nai; i++)
        {
            addAtom(cr->recv_int(src));
        }
        int nbo = cr->recv_int(src);
        for (int i=0; i < nbo; i++)
        {
            addBondOrder(cr->recv_double(src));
        }
    }
    else if (nullptr != debug)
    {
        fprintf(debug, "Trying to receive TopologyEntry, status %s\n", cs_name(cs));
        fflush(debug);
    }
    return cs;
}

bool AtomPair::operator==(const AtomPair &other) const
{
    return ((aI() == other.aI() && aJ() == other.aJ()) ||
            (aJ() == other.aI() && aI() == other.aJ()));
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

Angle::Angle(Bond bij, Bond bjk)
{
    b_[0] = bij;
    b_[1] = bjk;
    if (b_[0].aJ() != b_[1].aI())
    {
        GMX_THROW(gmx::InvalidInputError("Incorrect bonds passed"));
    }
    addAtom(b_[0].aI());
    addAtom(b_[0].aJ());
    addAtom(b_[1].aJ());
    addBondOrder(bij.bondOrder());
    addBondOrder(bjk.bondOrder());
}

void Angle::renumberAtoms(const std::vector<int> &renumber)
{
    b_[0].renumberAtoms(renumber);
    b_[1].renumberAtoms(renumber);
}

Improper::Improper(Bond bij, Bond bik, Bond bil)
{
    b_[0] = bij;
    b_[1] = bik;
    b_[2] = bil;
    if (b_[0].aI() != b_[1].aI())
    {
        GMX_THROW(gmx::InvalidInputError("Incorrect second bond passed"));
    }
    if (b_[0].aI() != b_[2].aI())
    {
        GMX_THROW(gmx::InvalidInputError("Incorrect third bond passed"));
    }
    addAtom(b_[0].aI());
    addAtom(b_[0].aJ());
    addAtom(b_[1].aJ());
    addAtom(b_[2].aJ());
    addBondOrder(bij.bondOrder());
    addBondOrder(bik.bondOrder());
    addBondOrder(bil.bondOrder());
}

void Improper::renumberAtoms(const std::vector<int> &renumber)
{
    b_[0].renumberAtoms(renumber);
    b_[1].renumberAtoms(renumber);
    b_[2].renumberAtoms(renumber);
}

Proper::Proper(Bond bij, Bond bjk, Bond bkl)
{
    b_[0] = bij;
    b_[1] = bjk;
    b_[2] = bkl;
    if (b_[0].aJ() != b_[1].aI())
    {
        GMX_THROW(gmx::InvalidInputError("Incorrect second bond passed"));
    }
    if (b_[1].aJ() != b_[2].aI())
    {
        GMX_THROW(gmx::InvalidInputError("Incorrect third bond passed"));
    }
    addAtom(b_[0].aI());
    addAtom(b_[0].aJ());
    addAtom(b_[1].aJ());
    addAtom(b_[2].aJ());
    addBondOrder(bij.bondOrder());
    addBondOrder(bjk.bondOrder());
    addBondOrder(bkl.bondOrder());
}

void Proper::renumberAtoms(const std::vector<int> &renumber)
{
    b_[0].renumberAtoms(renumber);
    b_[1].renumberAtoms(renumber);
    b_[2].renumberAtoms(renumber);
}

Topology::Topology(const std::vector<Bond> &bonds)
{
    entries_.insert(std::pair<InteractionType, std::vector<TopologyEntry *> > (InteractionType::BONDS, {}));
    for(auto &b : bonds)
    {
        entries_[InteractionType::BONDS].push_back(new Bond(b));
    }
}

void Topology::setAtoms(const t_atoms *atoms)
{
    atoms_.clear();
    for(int i = 0; i < atoms->nr; i++)
    {
        atoms_.push_back(ActAtom(*atoms->atomname[i],
                                 *atoms->atomtype[i],
                                 atoms->atom[i].ptype,
                                 atoms->atom[i].m,
                                 atoms->atom[i].q));
    }
}

const Bond *Topology::findBond(int ai, int aj) const
{
    Bond b(ai, aj, 1.0);
    auto bonds = entries_.find(InteractionType::BONDS)->second;
    auto bptr  = std::find_if(bonds.begin(), bonds.end(),
                              [b](const TopologyEntry *bb)
    { return (bb->atomIndex(0) == b.atomIndex(0) && bb->atomIndex(1) == b.atomIndex(1)) ||
             (bb->atomIndex(0) == b.atomIndex(1) && bb->atomIndex(1) == b.atomIndex(0)); });
    if (bptr == bonds.end())
    {
        GMX_THROW(gmx::InternalError("Cannot find bond"));
    }
    auto myptr = static_cast<Bond *>(*bptr);
    if (myptr->aI() == ai)
    {
        return myptr;
    }
    else
    {
        return findBond(aj, ai);
    }
}

const TopologyEntry *Topology::findTopologyEntry(InteractionType            itype,
                                                 const std::vector<int>    &aindex,
                                                 const std::vector<double> &bondOrder,
                                                 CanSwap                    cs) const
{
    if (!hasEntry(itype))
    {
        return nullptr;
    }
    auto entry = entries_.find(itype)->second;
    TopologyEntry te;
    for (auto &a : aindex)
    {
        te.addAtom(a);
    }
    for (auto &b : bondOrder)
    {
        te.addBondOrder(b);
    }
    auto bptr  = std::find_if(entry.begin(), entry.end(),
                              [te, cs](const TopologyEntry *bb)
    {
        bool same = te.atomIndices().size() == bb->atomIndices().size();
        if (same)
        {
            for(size_t i = 0; i < te.atomIndices().size() and same; i++)
            {
                same = same && (te.atomIndices()[i] == bb->atomIndices()[i]);
            }
            for(size_t i = 0; i < te.bondOrders().size() and same; i++)
            {
                same = same && (te.bondOrders()[i] == bb->bondOrders()[i]);
            }
        }
        if (!same && CanSwap::Yes == cs)
        {
            same = true;
            for(size_t i = 0; i < te.atomIndices().size() and same; i++)
            {
                same = same && (te.atomIndices()[i] == bb->atomIndices()[bb->atomIndices().size()-1-i]);
            }
            for(size_t i = 0; i < te.bondOrders().size() and same; i++)
            {
                same = same && (te.bondOrders()[i] == bb->bondOrders()[bb->bondOrders().size()-1-i]);
            }
        }
        return same;
    });
    if (bptr == entry.end())
    {
        return nullptr;
    }
    else
    {
        return *bptr;
    }
}
 
void Topology::makeAngles(const gmx::HostVector<gmx::RVec> &x,
                          double                            LinearAngleMin)
{
    auto bonds = entry(InteractionType::BONDS);
    entries_.insert(std::pair<InteractionType, std::vector<TopologyEntry *>>(InteractionType::ANGLES, {}));
    entries_.insert(std::pair<InteractionType, std::vector<TopologyEntry *>>(InteractionType::LINEAR_ANGLES, {}));
    auto &angles    = entries_.find(InteractionType::ANGLES)->second;
    auto &linangles = entries_.find(InteractionType::LINEAR_ANGLES)->second;
    for(size_t i = 0; i < bonds.size(); i++)
    {
        auto  b1  = static_cast<Bond *>(bonds[i]);
        int   ai1 = b1->atomIndex(0);
        int   aj1 = b1->atomIndex(1);
        for(size_t j = 0; j < bonds.size(); j++)
        {
            if (i == j)
            {
                continue;
            }
            auto  b2     = static_cast<Bond *>(bonds[j]);
            int   ai2    = b2->atomIndex(0);
            int   aj2    = b2->atomIndex(1);
            Angle *angle = nullptr;
            if (aj1 == ai2)
            {
                angle = new Angle(*b1, *b2);
            }
            else if (aj1 == aj2)
            {
                angle = new Angle(*b1, b2->swap());
            }
            else if (ai1 == ai2)
            {
                angle = new Angle(b1->swap(), *b2);
            }
            else if (ai1 == aj2)
            {
                angle = new Angle(b1->swap(), b2->swap());
            }
            if (angle)
            {
                
                if (is_linear(x[angle->atomIndex(0)], x[angle->atomIndex(1)], x[angle->atomIndex(2)],
                              nullptr, LinearAngleMin))
                {
                    angle->setLinear(true);
                    if (nullptr == findTopologyEntry(InteractionType::LINEAR_ANGLES, angle->atomIndices(),
                                          angle->bondOrders(), CanSwap::Yes))
                    {
                        linangles.push_back(angle);
                    }
                }
                else
                {
                    if (nullptr == findTopologyEntry(InteractionType::ANGLES, angle->atomIndices(),
                                          angle->bondOrders(), CanSwap::Yes))
                    {
                        angles.push_back(angle);
                    }
                }
            }
        }
    }
}

void Topology::makeImpropers(const gmx::HostVector<gmx::RVec> &x,
                             double                            PlanarAngleMax)
{
    size_t nAtoms = x.size();
    // Premature optimization is the root...
    // Store the bonds in a bidrectional data structure
    // using a vector of sets.
    std::vector<std::set<int> > mybonds;
    mybonds.resize(nAtoms);
    auto &bonds = entry(InteractionType::BONDS);
    for(auto &b : bonds)
    {
        mybonds[b->atomIndex(0)].insert(b->atomIndex(1));
        mybonds[b->atomIndex(1)].insert(b->atomIndex(0));
    }
    entries_.insert(std::pair<InteractionType, std::vector<TopologyEntry *>>(InteractionType::IMPROPER_DIHEDRALS, {}));
    auto &impropers = entries_.find(InteractionType::IMPROPER_DIHEDRALS)->second;
    // Now loop over the atoms
    const rvec *myx = as_rvec_array(x.data());
    for(size_t i = 0; i < nAtoms; i++)
    {
        if (mybonds[i].size() == 3)
        {
            std::vector<int> jkl;
            for (auto &s : mybonds[i])
            {
                jkl.push_back(s);
            }
            if (is_planar(myx[i], myx[jkl[0]], myx[jkl[1]], myx[jkl[2]], nullptr, PlanarAngleMax))
            {
                Bond bjkl[3];
                for(int m = 0; m < 3; m++)
                {
                    bjkl[m] = *findBond(i, jkl[m]);
                    if (bjkl[m].atomIndex(0) != static_cast<int>(i))
                    {
                        bjkl[m] = bjkl[m].swap();
                    }
                }
                impropers.push_back(new Improper(bjkl[0], bjkl[1], bjkl[2]));
            }
        }
    }
}

void Topology::makePairs(int natoms)
{
    for(const auto &itype : { InteractionType::VDW,
                              InteractionType::COULOMB })
    {
        entries_.insert({itype, {} });
        auto &pairs = entries_.find(itype)->second;
        for(int i = 0; i < natoms; i++)
        {
            // Check for exclusions is done later.
            for(int j = i+1; j < natoms; j++)
            {
                pairs.push_back(new AtomPair(i, j));
            }
        }
    }
}

void  Topology::addShellPairs()
{
    if (hasEntry(InteractionType::POLARIZATION))
    {
        auto pol = entries_.find(InteractionType::POLARIZATION)->second;
        for (const auto &itype : { InteractionType::VDW,
                                   InteractionType::COULOMB } )
        {
            if (!hasEntry(itype))
            {
                continue;
            }
            auto &ffpl = entries_.find(itype)->second;
            std::vector<AtomPair *> newf;
            for(const auto &p_i : pol)
            {
                auto core_i  = p_i->atomIndices()[0];
                auto shell_i = p_i->atomIndices()[1];
                for(const auto &p_j : pol)
                {
                    auto core_j  = p_j->atomIndices()[0];
                    auto shell_j = p_j->atomIndices()[1];
                    for(auto &f : ffpl)
                    {
                        auto indices = f->atomIndices();
                        if (!(core_i == indices[0] && core_j == indices[1]) && 
                            !(core_i == indices[1] && core_j == indices[0]))
                        {
                            // This pair of cores interacts, now add the shells.
                            newf.push_back(new AtomPair(core_i, shell_j));
                            newf.push_back(new AtomPair(core_j, shell_i));
                            newf.push_back(new AtomPair(shell_i, shell_j));
                        }
                    }
                }
            }
            for(const auto &n : newf)
            {
                ffpl.push_back(std::move(n));
            }
        }
    }   
}

void Topology::makePropers()
{
    entries_.insert(std::pair<InteractionType, std::vector<TopologyEntry *>>(InteractionType::PROPER_DIHEDRALS, {}));
    auto &propers = entries_.find(InteractionType::PROPER_DIHEDRALS)->second;
    auto &angles  = entries_.find(InteractionType::ANGLES)->second;
    for(size_t i = 0; i < angles.size(); i++)
    {
        auto a1 = static_cast<Angle *>(angles[i]);
        if (a1->isLinear())
        {
            continue;
        }
        auto bij1 = a1->bij();
        auto bjk1 = a1->bjk();
        for(size_t j = i+1; j < angles.size(); j++)
        {
            auto a2 = static_cast<Angle *>(angles[j]);
            if (a2->isLinear())
            {
                continue;
            }
            auto bij2 = a2->bij();
            auto bjk2 = a2->bjk();
            // We are looking for linear dihedrals with order
            // i-j-k-l, and there are four possibilities:
            // 1) the first angle is i-j-k and the second j-k-l
            if (bjk1.aI() == bij2.aI() && bjk1.aJ() == bij2.aJ())
            {
                propers.push_back(new Proper(bij1, bjk1, bjk2));
            }
            // 2) the first angle is k-j-i and the second j-k-l
            else if (bij1.aI() == bij2.aJ() && bij1.aJ() == bij2.aI())
            {
                propers.push_back(new Proper(bjk1.swap(), bij2, bjk2));
            }
            // 2) the first angle is k-j-i and the second l-k-j
            else if (bij1.aI() == bij2.aJ() && bij1.aJ() == bjk2.aJ() )
            {
                propers.push_back(new Proper(bjk1.swap(), bij1.swap(), bij2.swap()));
            }
            // 4) the first angle is i-j-k and the second l-k-j
            else if (bij1.aJ() == bjk2.aJ() && bjk1.aJ() == bij2.aJ())
            {
                propers.push_back(new Proper(bij1, bjk1, bij2.swap()));
            }
        }
    }
}

void Topology::makeVsite2s(const ForceFieldParameterList &vsite2)
{
    gmx_fatal(FARGS, "Vsite2 not implemented yet");
}

void Topology::renumberAtoms(const std::vector<int> &renumber)
{
    for(auto &myEntry: entries_)
    {
        for(auto &b : myEntry.second)
        {
            b->renumberAtoms(renumber);
        }
    }
}

const std::vector<TopologyEntry *> &Topology::entry(InteractionType itype) const
{
    return entries_.find(itype)->second;
}

void Topology::addEntry(InteractionType                     itype,
                        const std::vector<TopologyEntry *> &entry)
{
    if (hasEntry(itype))
    {
        GMX_THROW(gmx::InternalError(gmx::formatString("InteractionType %s already present",
                                     interactionTypeToString(itype).c_str()).c_str()));
    }
    entries_.insert(std::pair<InteractionType, std::vector<TopologyEntry *>>(itype, entry));
}

void Topology::generateExclusions(int nrexcl,
                                  int nratoms)
{
    exclusions_.resize(nratoms);
    
    for(auto &myEntry: entries_)
    {
        switch (myEntry.first)
        {
        case InteractionType::BONDS:
            {
                if (nrexcl > 0)
                {
                    for(auto &b : myEntry.second)
                    {
                        auto a = b->atomIndices();
                        exclusions_[a[0]].push_back(a[1]);
                        exclusions_[a[1]].push_back(a[0]);
                    }
                }
            }
            break;
        case InteractionType::ANGLES:
        case InteractionType::LINEAR_ANGLES:
            {
                if (nrexcl > 1)
                {
                    for(auto &b : myEntry.second)
                    {
                        auto a = b->atomIndices();
                        exclusions_[a[0]].push_back(a[2]);
                        exclusions_[a[2]].push_back(a[0]);
                    }
                }
            } 
            break;
        default:
            break;
        }
    }
    // Update our own VDW and COULOMB pairs if present
    for(const auto &itype : { InteractionType::VDW,
                              InteractionType::COULOMB })
    {
        if (hasEntry(itype))
        {
            auto itt = &entries_.find(itype)->second;
            for(size_t ai = 0; ai < exclusions_.size(); ai++)
            {
                for(size_t j = 0; j < exclusions_[ai].size(); j++)
                {
                    AtomPair ap(ai, exclusions_[ai][j]);
                    for(auto pp = itt->begin(); pp < itt->end(); ++pp)
                    {
                    if (ap == *static_cast<AtomPair *>(*pp))
                    {
                        itt->erase(pp);
                        break;
                    }
                    }
                }
            }
        }
    }
}

t_excls *Topology::gromacsExclusions()
{
    t_excls *gmx;
    snew(gmx, exclusions_.size());
    for(size_t i = 0; i < exclusions_.size(); i++)
    {
        snew(gmx[i].e, exclusions_[i].size());
        for (size_t j = 0; j < exclusions_[i].size(); j++)
        {
            gmx[i].e[j] = exclusions_[i][j];
        }
        gmx[i].nr = exclusions_[i].size();
    }
    return gmx;
}

} // namespace alexandria

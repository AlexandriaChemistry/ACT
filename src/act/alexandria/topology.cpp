/*
 * This source file is part of the Alexandria Chemistry Toolkit.
 *
 * Copyright (C) 2021-2023
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

#include "topology.h"

#include <algorithm>
#include <map>
#include <set>
#include <vector>

#include "act/basics/interactiontype.h"
#include "act/forcefield/forcefield_parametername.h"
#include "act/molprop/experiment.h"
#include "gromacs/gmxpreprocess/grompp-impl.h"
#include "gromacs/math/vectypes.h"
#include "gromacs/utility/fatalerror.h"

#include "actmol_low.h"

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
    else if (nullptr != debug)
    {
        fprintf(debug, "Trying to receive TopologyEntry, status %s\n", cs_name(cs));
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

bool AtomPair::operator<(const AtomPair &other) const
{
    return (aI() < other.aI() ||
            (aI() == other.aI() && aJ() < other.aJ()) ||
            (aI() == other.aJ() && aJ() < other.aI()));
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

void Vsite2::get(int *ai, int *aj, int *vs) const
{
    check(3);
    *ai = atomIndex(0);
    *aj = atomIndex(1);
    *vs = atomIndex(2);
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
    GMX_RELEASE_ASSERT(atomIndices().size() == 4, "Something weird with impropers");
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
    if (!bonds.empty())
    {
        entries_.insert({ InteractionType::BONDS, {} });
        for(auto &b : bonds)
        {
            if (b.aI() < 0 || b.aJ() < 0)
            {
                GMX_THROW(gmx::InternalError("Negative atom indices when adding bonds"));
            }
            entries_[InteractionType::BONDS].push_back(new Bond(b));
        }
    }
}

static void dump_entry(FILE *fp, 
                       const std::vector<TopologyEntry *> &entries,
                       const std::string &label)
{
    if (nullptr == fp)
    {
        return;
    }
    for(auto &entry : entries)
    {
        fprintf(fp, "%s", label.c_str());
        for (auto &ai : entry->atomIndices())
        {
            fprintf(fp, " %d", ai);
        }
        fprintf(fp, "\n");
    }
}

void Topology::addShells(FILE                   *fp,
                         const ForceField       *pd,
                         std::vector<gmx::RVec> *x)
{
    /* Calculate the total number of Atom and Vsite particles and
     * generate the renumbering array.
     */
    std::vector<int>             shellRenumber;
    std::map<int, int>           originalAtomIndex;

    int                          nshell = 0;
    shellRenumber.resize(atoms_.size(), 0);
    for (size_t i = 0; i < atoms_.size(); i++)
    {
        auto atype       = pd->findParticleType(atoms_[i].ffType());
        shellRenumber[i] = i + nshell;
        originalAtomIndex.insert({shellRenumber[i], i});
        if (atype->hasInteractionType(InteractionType::POLARIZATION))
        {
            // TODO: Update if particles can have more than one shell
            nshell++;
        }
        if (debug)
        {
            fprintf(debug, "Shell renumber %s %zu is %d\n", atoms_[i].name().c_str(), i, shellRenumber[i]);
        }
    }
    if (fp)
    {
        fprintf(fp, "Found %d shells to be added\n", nshell);
    }
    
    /* Add Polarization to the plist. */
    std::vector<TopologyEntry *> pols;
    auto &fs  = pd->findForcesConst(InteractionType::POLARIZATION);
    
    // Atoms first
    size_t i = 0;
    auto xiter = x->begin();
    for (auto iter = atoms_.begin(); iter < atoms_.end() && xiter < x->end(); ++iter, ++xiter)
    {
        if (iter->pType() == eptAtom || iter->pType() == eptVSite)
        {
            std::string atomtype(iter->ffType());
            if (pd->hasParticleType(atomtype))
            {
                // TODO: Allow adding multiple shells
                auto fa = pd->findParticleType(atomtype);
                if (fa->hasInteractionType(InteractionType::POLARIZATION))
                {
                    auto ptype = fa->interactionTypeToIdentifier(InteractionType::POLARIZATION);
                    auto param  = fs.findParameterTypeConst(ptype, pol_name[polALPHA]);
                    auto fp     = pd->findParticleType(ptype);
                    auto charge = fp->charge();
                    auto pol    = convertToGromacs(param.value(), param.unit());
                    if (pol <= 0)
                    {
                        GMX_THROW(gmx::InvalidInputError(gmx::formatString("Polarizability should be positive for %s", fa->id().id().c_str()).c_str()));
                    }
                    // TODO Multiple shell support
                    auto pp = new TopologyEntry();
                    auto core = shellRenumber[i];
                    auto shell = core+1;
                    pp->addAtom(core);
                    pp->addAtom(shell);
                    pp->addBondOrder(1.0);
                    pols.push_back(pp);

                    // Insert shell atom
                    auto shellName = iter->name() + "s";
                    ActAtom newshell(shellName, "EP", ptype.id(), eptShell,
                                     0, 0, charge);
                    newshell.setResidueNumber(iter->residueNumber());
                    iter = atoms_.insert(iter+1, newshell);
                    // Insert new coordinate
                    xiter = x->insert(xiter+1, (*x)[shellRenumber[i]]); 
                    // Create link between shell and atom
                    atoms_[core].addShell(shell);
                    atoms_[shell].setCore(core);
                    if (debug)
                    {
                        fprintf(debug, "Just made an atom-shell pair %s %d %s %d\n",
                                atoms_[core].name().c_str(), core,
                                atoms_[shell].name().c_str(), shell);
                    }
                }
            }
            else
            {
                GMX_THROW(gmx::InvalidInputError(gmx::formatString("Cannot find atomtype %s in forcefield\n", 
                                                                   atomtype.c_str())));
            }
            // Increase counter for atoms and vsites only.
            i += 1;
        }
    }
    // Check length of arrays.
    if (atoms_.size() != x->size())
    {
        GMX_THROW(gmx::InternalError(gmx::formatString("Inconsistency: %zu atoms but %zu coordinates.", atoms_.size(),
                                                       x->size()).c_str()));
    }
    if (!pols.empty())
    {
        addEntry(InteractionType::POLARIZATION, pols);
        dump_entry(debug, pols, "Before renumber");
        renumberAtoms(shellRenumber);
        dump_entry(debug, pols, "After renumber");
    }
}


void Topology::addBond(const Bond &bond)
{
    if (bond.aI() < 0 || bond.aJ() < 0)
    {
        GMX_THROW(gmx::InternalError("Negative atom indices when adding bonds"));
    }
    entries_.insert({ InteractionType::BONDS, {} });
    entries_[InteractionType::BONDS].push_back(new Bond(bond));
}

int Topology::makeVsites(const ForceField       *pd,
                         std::vector<gmx::RVec> *x)
{
    int vsitesAdded = 0;
    // Check whether we have virtual sites in the force field.
    if (pd->interactionPresent(InteractionType::VSITE2))
    {
        vsitesAdded += makeVsite2s(pd);
    }
    if (vsitesAdded > 0)
    {
        x->resize(x->size()+vsitesAdded);
    }
    return vsitesAdded;
}

double Topology::mass() const
{
    double mtot = 0;
    for(size_t i = 0; i < atoms_.size(); i++)
    {
        mtot += atoms_[i].mass();
    }
    return mtot;
}

const Bond *Topology::findBond(int ai, int aj) const
{
    Bond b(ai, aj, 1.0);
    std::vector<TopologyEntry *>::const_iterator bptr;
    auto bondsptr = entries_.find(InteractionType::BONDS);
    if (entries_.end() != bondsptr)
    {
        const std::vector<TopologyEntry *> &bonds = bondsptr->second;
        bptr  = std::find_if(bonds.begin(), bonds.end(),
                             [b](const TopologyEntry *bb)
                             { return (bb->atomIndex(0) == b.atomIndex(0) && bb->atomIndex(1) == b.atomIndex(1)) ||
                               (bb->atomIndex(0) == b.atomIndex(1) && bb->atomIndex(1) == b.atomIndex(0)); });
        if (bptr == bonds.end())
        {
            GMX_THROW(gmx::InternalError(gmx::formatString("Cannot find bond between %d and %d", ai, aj).c_str()));
        }
    }
    else
    {
        GMX_THROW(gmx::InternalError("There are no bonds at all."));
    }
    auto myptr = static_cast<const Bond *>(*bptr);
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
    auto eptr  = entries_.find(itype);
    if (entries_.end() == eptr)
    {
        GMX_THROW(gmx::InternalError(gmx::formatString("No such entry %s in topology", interactionTypeToString(itype).c_str()).c_str()));
    }
    auto entry = eptr->second;
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
 
void Topology::makeAngles(const std::vector<gmx::RVec> &x,
                          double                        LinearAngleMin)
{
    auto ib = InteractionType::BONDS;
    if (entries_.find(ib) == entries_.end())
    {
        return;
    }
    auto &bonds = entry(ib);
    entries_.insert({ InteractionType::ANGLES, {} });
    entries_.insert({ InteractionType::LINEAR_ANGLES, {} });
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

void Topology::makeImpropers(const std::vector<gmx::RVec> &x,
                             double                        PlanarAngleMax)
{
    auto ib = InteractionType::BONDS;
    if (entries_.find(ib) == entries_.end())
    {
        return;
    }
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
    entries_.insert({ InteractionType::IMPROPER_DIHEDRALS, {} });
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

void Topology::fixShellExclusions()
{
    for (const auto &itype : { InteractionType::VDW,
                               InteractionType::COULOMB } )
    {
        if (!hasEntry(itype) || exclusions_.find(itype) == exclusions_.end())
        {
            continue;
        }
        auto &ffpl = entries_.find(itype)->second;
        // Loop over all exclusions
        for(size_t ai = 0; ai < exclusions_[itype].size(); ++ai)
        {
            for(size_t jj = 0; jj < exclusions_[itype][ai].size(); ++jj)
            {
                size_t aj = exclusions_[itype][ai][jj];
                // Check whether these particles have shells
                for(size_t si : atoms_[ai].shells())
                {
                    for(size_t sj : atoms_[aj].shells())
                    {
                        // See whether this interaction exists
                        auto it = ffpl.begin();
                        while (ffpl.end() != it)
                        {
                            size_t aai = (*it)->atomIndices()[0];
                            size_t aaj = (*it)->atomIndices()[1];
                            if (((aai == si || aai == ai) && (aaj == sj || aaj == aj)) || 
                                ((aaj == si || aaj == ai) && (aai == sj || aai == aj)))
                            {
                                it = ffpl.erase(it);
                            }
                            else
                            {
                                ++it;
                            }
                        }
                    }
                }                
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
            std::set<AtomPair> newf;
            for(auto &f : ffpl)
            {
                auto indices = f->atomIndices();
                auto core_i = indices[0];
                auto core_j = indices[1];
                int shell_i = -1;
                int shell_j = -1;
                for(const auto &p_i : pol)
                {
                    if (core_i == p_i->atomIndices()[0])
                    {
                        shell_i = p_i->atomIndices()[1];
                        break;
                    }
                }
                for(const auto &p_j : pol)
                {
                    if (core_j == p_j->atomIndices()[0])
                    {
                        shell_j = p_j->atomIndices()[1];
                        break;
                    }
                }
                if (shell_i != -1)
                {
                    newf.emplace(AtomPair(shell_i, core_j));
                }
                if (shell_j != -1)
                {
                    newf.emplace(AtomPair(core_i, shell_j));
                }
                if (shell_i != -1 && shell_j != -1)
                {
                    newf.emplace(AtomPair(shell_i, shell_j));
                }
            }
            for(const auto &n : newf)
            {
                auto ap = new AtomPair(n);
                ffpl.push_back(std::move(ap));
            }
        }
    }
}

void Topology::makePropers()
{
    auto ia = InteractionType::ANGLES;
    if (entries_.find(ia) == entries_.end())
    {
        return;
    }
    entries_.insert({ InteractionType::PROPER_DIHEDRALS, {} });
    auto &propers = entries_.find(InteractionType::PROPER_DIHEDRALS)->second;
    auto &angles  = entries_.find(ia)->second;
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

int Topology::makeVsite2s(const ForceField *pd)
{
    auto itype_vs2 = InteractionType::VSITE2;
    auto &ffvs     = pd->findForcesConst(itype_vs2);
    if (ffvs.empty())
    {
        return 0;
    }
    // Loop over all the bonds to make vsites.
    auto itype_b = InteractionType::BONDS;
    if (entries_.find(itype_b) == entries_.end())
    {
        return 0;
    }
    auto &bonds     = entry(itype_b);
    entries_.insert({ itype_vs2, {} });
    auto &vsite2    = entries_.find(itype_vs2)->second;
    int vsitesAdded = 0;
    for(size_t i = 0; i < bonds.size(); i++)
    {
        auto bid    = static_cast<Bond *>(bonds[i])->id();
        auto batoms = bid.atoms();
        auto border = bid.bondOrders();
        GMX_RELEASE_ASSERT(2 == batoms.size(), "Not enough atoms");
        printf("Found bond %s %s\n", batoms[0].c_str(), batoms[1].c_str());
        
        for(auto fvs : ffvs.parametersConst())
        {
            printf("Checking vsite %s\n", fvs.first.id().c_str());
            auto vsatoms = fvs.first.atoms();
            auto vsbo    = fvs.first.bondOrders();
            bool equal   = batoms.size()+1 == vsatoms.size();
            for(size_t i = 0; i < batoms.size() && equal; ++i)
            {
                equal = equal && (batoms[i] == vsatoms[i]);
            }
            for(size_t i = 0; i < border.size() && equal; ++i)
            {
                equal = equal && (border[i] == vsbo[i]);
            }
            if (equal)
            {
                auto vsname = vsatoms[vsatoms.size()-1];
                if (!pd->hasParticleType(vsname))
                {
                    printf("No such particle type %s as found in vsite %s\n", vsname.c_str(), fvs.first.id().c_str());
                }
                else
                {
                    auto ptype = pd->findParticleType(vsname);
                    
                    ActAtom newatom(ptype->id().id(), ptype->element(), ptype->id().id(),
                                    ptype->gmxParticleType(), 0, ptype->mass(), ptype->charge());
                    atoms_.push_back(newatom);
                    vsite2.push_back(new Vsite2(bonds[i]->atomIndex(0), bonds[i]->atomIndex(1), atoms_.size()-1));
                    vsitesAdded += 1;
                }
            }
        }    
    }
    return vsitesAdded;
}

void Topology::renumberAtoms(const std::vector<int> &renumber)
{
    for(auto &myEntry: entries_)
    {
        if (InteractionType::POLARIZATION != myEntry.first)
        {
            for(auto &b : myEntry.second)
            {
                b->renumberAtoms(renumber);
            }
        }
    }
}

void Topology::dumpPairlist(FILE *fp, InteractionType itype) const
{
    if (nullptr == fp)
    {
        return;
    }
    auto plist = entry(itype);
    for(auto pl = plist.begin(); pl < plist.end(); ++pl)
    {
        fprintf(fp, "PAIRLIST %s %d %d\n", interactionTypeToString(itype).c_str(),
                (*pl)->atomIndices()[0], (*pl)->atomIndices()[1]);
    }
}

immStatus Topology::GenerateAtoms(const ForceField       *pd,
                                  const MolProp          *mol,
                                  std::vector<gmx::RVec> *x)
{
    immStatus imm   = immStatus::OK;
    
    auto ci = mol->findExperimentConst(JobType::OPT);
    if (!ci)
    {
        ci = mol->findExperimentConst(JobType::TOPOLOGY);
    }
    if (ci)
    {
        if (ci->NAtom() == 0)
        {
            return immStatus::NoAtoms;
        }
        int res0 = -666;
        int nres = 0;
        for (auto &cai : ci->calcAtomConst())
        {
            auto myunit = cai.coordUnit();
            double    xx, yy, zz;
            cai.coords(&xx, &yy, &zz);
            int resnr = cai.residueNumber();
            if (resnr != res0)
            {
                res0  = resnr;
                addResidue(nres, cai.residueName());
                nres += 1;
            }
            gmx::RVec xxx = { convertToGromacs(xx, myunit),
                              convertToGromacs(yy, myunit),
                              convertToGromacs(zz, myunit) };
            x->push_back(xxx);

            if (pd->hasParticleType(cai.getObtype()))
            {
                auto atype = pd->findParticleType(cai.getObtype());
                
                ActAtom newatom(cai.getName(), atype->element(), cai.getObtype(), eptAtom,
                                atype->atomnumber(), atype->mass(), atype->charge());
                newatom.setResidueNumber(nres-1);
                addAtom(newatom);
            }
            else
            {
                fprintf(stderr, "Cannot find atomtype %s (atom %zu) in forcefield, there are %d atomtypes.\n", 
                        cai.getObtype().c_str(), atoms_.size(), pd->nParticleTypes());
                
                return immStatus::AtomTypes;
            }
        }
    }
    else
    {
        imm = immStatus::Topology;
    }
    if (nullptr != debug)
    {
        fprintf(debug, "Tried to convert %s to gromacs. LOT is %s/%s. Natoms is %zu\n",
                mol->getMolname().c_str(),
                ci->getMethod().c_str(),
                ci->getBasisset().c_str(), atoms_.size());
    }

    return imm;
}

void Topology::build(const ForceField             *pd,
                     std::vector<gmx::RVec>       *x,
                     double                        LinearAngleMin,
                     double                        PlanarAngleMax,
                     missingParameters             missing)
{
    makeAngles(*x, LinearAngleMin);
    makeImpropers(*x, PlanarAngleMax);
    // Check whether we have dihedrals in the force field.
    if (pd->interactionPresent(InteractionType::PROPER_DIHEDRALS))
    {
        // Store temporary address of variable for performance
        auto &fs = pd->findForcesConst(InteractionType::PROPER_DIHEDRALS);
        if (!fs.empty() || missingParameters::Generate == missing)
        {
            makePropers();
        }
    }
    // Check whether we have virtual sites in the force field.
    if (pd->interactionPresent(InteractionType::VSITE2))
    {
        auto &fs = pd->findForcesConst(InteractionType::VSITE2);
        if (!fs.empty())
        {
            makeVsite2s(pd);
        }
    }
    makePairs(x->size());
    int nexclvdw;
    if (!ffOption(*pd, InteractionType::VDW, 
                  "nexcl", &nexclvdw))
    {
        nexclvdw = 0;
    }
    int nexclqq;
    if (!ffOption(*pd, InteractionType::COULOMB, 
                  "nexcl", &nexclqq))
    {
        nexclqq = 0;
    }
    generateExclusions(nexclvdw, nexclqq, x->size());
    if (pd->polarizable())
    {
        addShells(debug, pd, x);
        addShellPairs();
        fixShellExclusions();
    }
    setIdentifiers(pd);
    if (missing != missingParameters::Generate)
    {
        fillParameters(pd);
    }
    if (debug)
    {
        dumpPairlist(debug, InteractionType::COULOMB);
        dumpPairlist(debug, InteractionType::VDW);
    }
}

const std::vector<TopologyEntry *> &Topology::entry(InteractionType itype) const
{
    return entries_.find(itype)->second;
}

void Topology::addResidue(int                residueNumber, 
                          const std::string &residueName)
{
    if (residueNumber >= static_cast<int>(residueNames_.size()))
    {
        residueNames_.resize(residueNumber+1);
        residueNames_[residueNumber] = residueName;
    }
    else if (residueNames_[residueNumber] != residueName)
    {
        GMX_THROW(gmx::InternalError(gmx::formatString("Residue name mismatch. Have residues %d = %s, but new name %s",
                                                       residueNumber, residueNames_[residueNumber].c_str(),
                                                       residueName.c_str()).c_str()));
    }
}
    
void Topology::addEntry(InteractionType                     itype,
                        const std::vector<TopologyEntry *> &entry)
{
    if (hasEntry(itype))
    {
        GMX_THROW(gmx::InternalError(gmx::formatString("InteractionType %s already present",
                                                       interactionTypeToString(itype).c_str()).c_str()));
    }
    if (entry.size() > 0)
    {
        entries_.insert({ itype, entry });
    }
}

void Topology::generateExclusions(int nrexclvdw,
                                  int nrexclqq,
                                  int nratoms)
{
    // Update our own VDW and COULOMB pairs if present
    for(const auto &itype : { InteractionType::VDW,
                              InteractionType::COULOMB })
    {
        if (!hasEntry(itype))
        {
            continue;
        }
        exclusions_.insert({itype, {}});
        exclusions_[itype].resize(nratoms);
        int nrexcl = nrexclqq;
        if (InteractionType::VDW == itype)
        {
            nrexcl = nrexclvdw;
        }
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
                            exclusions_[itype][a[0]].push_back(a[1]);
                            exclusions_[itype][a[1]].push_back(a[0]);
                        }
                    }
                }
                break;
            case InteractionType::POLARIZATION:
                {
                    for(auto &b : myEntry.second)
                    {
                        auto a = b->atomIndices();
                        exclusions_[itype][a[0]].push_back(a[1]);
                        exclusions_[itype][a[1]].push_back(a[0]);
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
                            exclusions_[itype][a[0]].push_back(a[2]);
                            exclusions_[itype][a[2]].push_back(a[0]);
                        }
                    }
                } 
                break;
            default:
                break;
            }
        }
        auto itt = &entries_.find(itype)->second;
        for(size_t ai = 0; ai < exclusions_[itype].size(); ai++)
        {
            for(size_t j = 0; j < exclusions_[itype][ai].size(); j++)
            {
                AtomPair ap(ai, exclusions_[itype][ai][j]);
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

t_excls *Topology::gromacsExclusions()
{
    t_excls *gmx;
    auto     itype = InteractionType::VDW;
    snew(gmx, exclusions_[itype].size());
    for(size_t i = 0; i < exclusions_[itype].size(); i++)
    {
        snew(gmx[i].e, exclusions_[itype][i].size());
        for (size_t j = 0; j < exclusions_[itype][i].size(); j++)
        {
            gmx[i].e[j] = exclusions_[itype][i][j];
        }
        gmx[i].nr = exclusions_[itype][i].size();
    }
    return gmx;
}

void Topology::dump(FILE *fp) const
{
    if (nullptr == fp)
    {
        return;
    }
    for(auto &myEntry: entries_)
    {
        fprintf(fp, "%s\n", interactionTypeToString(myEntry.first).c_str());
        for (auto &tt : myEntry.second)
        {
            for(auto &aa : tt->atomIndices())
            {
                fprintf(fp, " %d", aa+1);
            }
            fprintf(fp, "\n");
        }
    }
}

static void fillParams(const ForceFieldParameterList &fs,
                       const Identifier              &btype,
                       int                            nr,
                       const char                    *param_names[],
                       std::vector<double>           *param)
{
    for (int i = 0; i < nr; i++)
    {
        if (!fs.parameterExists(btype))
        {
            GMX_THROW(gmx::InternalError(gmx::formatString("Cannot find parameters for %s and topid %s", fs.function().c_str(),
                                                           btype.id().c_str()).c_str()));
        }
        else
        {
            // If a particular parameter is missing, we set it to zero.
            // This is needed in particular to deal with combination rules, where
            // some force fields have the atomic sigma, epsilon, etc. while others
            // store the combined values, but not atomic ones (because that would
            // not make sense).
            auto ff      = fs.findParametersConst(btype);
            auto fp      = ff.find(param_names[i]);
            double value = 0;
            if (ff.end() != fp)
            {
                value = fp->second.internalValue();
            }
            param->push_back(value);
        }
    }
}

void Topology::fillParameters(const ForceField *pd)
{
    for(auto &entry : entries_)
    {
        auto &fs = pd->findForcesConst(entry.first);
        for(auto &topentry : entry.second)
        {
            const auto &topID = topentry->id();

            std::vector<double> param;
            switch (fs.gromacsType())
            {
            case F_LJ:
                fillParams(fs, topID, ljNR, lj_name, &param);
                break;
            case F_BHAM:
                fillParams(fs, topID, wbhNR, wbh_name, &param);
                break;
            case F_GBHAM:
                fillParams(fs, topID, gbhNR, gbh_name, &param);
                break;
            case F_COUL_SR:
                fillParams(fs, topID, coulNR, coul_name, &param);
                break;
            case F_MORSE:
                fillParams(fs, topID, morseNR, morse_name, &param);
                break;
            case F_CUBICBONDS:
                fillParams(fs, topID, cubicNR, cubic_name, &param);
                break;
            case F_BONDS:
                fillParams(fs, topID, bondNR, bond_name, &param);
                break;
            case F_ANGLES:
                fillParams(fs, topID, angleNR, angle_name, &param);
                break;
            case F_UREY_BRADLEY:
                fillParams(fs, topID, ubNR, ub_name, &param);
                break;
            case F_LINEAR_ANGLES:
                fillParams(fs, topID, linangNR, linang_name, &param);
                break;
            case F_IDIHS:
                fillParams(fs, topID, idihNR, idih_name, &param);
                break;
            case F_FOURDIHS:
                fillParams(fs, topID, fdihNR, fdih_name, &param);
                break;
            case F_POLARIZATION:
                fillParams(fs, topID, polNR, pol_name, &param);
                break;
            case F_PDIHS:
                fillParams(fs, topID, pdihNR, pdih_name, &param);
                break;
            default:
                GMX_THROW(gmx::InternalError(gmx::formatString("Missing case %s when filling the topology structure.", interaction_function[fs.gromacsType()].name).c_str()));
            }
            topentry->setParams(param);
        }
    }
}

void Topology::setIdentifiers(const ForceField *pd)
{
    for(auto &entry : entries_)
    {
        auto &fs = pd->findForcesConst(entry.first);
        for(auto &topentry : entry.second)
        {
            std::string              atypes;
            std::vector<std::string> btype;
            for(const size_t &jj : topentry->atomIndices())
            {
                if (jj >= atoms_.size())
                {
                    GMX_THROW(gmx::InternalError(gmx::formatString("Atom index %zu should be less than %zu for %s",
                                                                   jj, atoms_.size(),
                                                                   interactionTypeToString(entry.first).c_str()).c_str()
                                                 ));
                }
                // Check whether we do not e.g. have a shell here
                //if (!pd->hasParticleType(atoms_[jj].ffType()))
                //{
                //  continue;
                //}
                auto atype = pd->findParticleType(atoms_[jj].ffType());
                if (!atypes.empty())
                {
                    atypes += " ";
                }
                atypes += atype->id().id();
                switch (entry.first)
                {
                case InteractionType::VDW:
                case InteractionType::VSITE2:
                    {
                        btype.push_back(atoms_[jj].ffType());
                        break;
                    }
                case InteractionType::POLARIZATION:
                    {
                        // For COULOMB there are two particles,
                        // but for polarization just one.
                        if (atype->hasInteractionType(entry.first))
                        {
                            btype.push_back(atype->interactionTypeToIdentifier(entry.first).id());
                        }
                        break;
                    }
                case InteractionType::COULOMB:
                    {
                        // For COULOMB there are two particles,
                        // but for polarization just one.
                        if (atype->hasInteractionType(entry.first))
                        {
                            btype.push_back(atype->interactionTypeToIdentifier(entry.first).id());
                        }
                        break;
                    }
                default:
                    {
                        auto itype = InteractionType::BONDS;
                        if (atype->hasInteractionType(itype))
                        {
                            btype.push_back(atype->interactionTypeToIdentifier(itype).id());
                        }
                        break;
                    }
                }
            }
            if (!btype.empty())
            {
                if (btype.size() == 1)
                {
                    topentry->setId(Identifier(btype[0]));
                }
                else
                {
                    topentry->setId({ Identifier(btype, topentry->bondOrders(),
                                                 fs.canSwap()) });
                    if ((btype.size() == 3 && btype[2].compare("h_b") == 0 && topentry->bondOrders()[1] != 1) || 
                        (topentry->id().id().compare("c2_b:c2_b:h_b~c2_b") == 0))
                    {
                        fprintf(stderr, "Found a funny one: %s\n", topentry->id().id().c_str());
                    }
                }
            }
            else if (debug)
            {
                fprintf(debug, "Could not find identifier for %s for atomtype '%s'", interactionTypeToString(entry.first).c_str(), atypes.c_str());
            }
        }
    }
}

} // namespace alexandria

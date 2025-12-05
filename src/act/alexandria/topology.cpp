/*
 * This source file is part of the Alexandria Chemistry Toolkit.
 *
 * Copyright (C) 2021-2025
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
#include <list>
#include <map>
#include <set>
#include <string>
#include <vector>

#include "act/basics/interactiontype.h"
#include "act/basics/msg_handler.h"
#include "act/forcefield/forcefield_parametername.h"
#include "act/molprop/experiment.h"
#include "act/molprop/topologyentry.h"
#include "gromacs/gmxpreprocess/grompp-impl.h"
#include "gromacs/math/vectypes.h"
#include "gromacs/utility/fatalerror.h"

namespace alexandria
{

class ActAtomListItem
{
private:
    //! The atom information
    ActAtom   atom_;
    //! The original index
    size_t    index_ = 0;
    //! Atomic coordinates
    gmx::RVec x_     = { 0, 0, 0 };
public:
    //! constructor
    ActAtomListItem(const ActAtom   atom,
                    size_t          index,
                    const gmx::RVec x) : atom_(atom), index_(index), x_(x) {}

    //! \return the atom
    const ActAtom atom() const { return atom_; }

    //! \return pointer to the atom for editing
    ActAtom *atomPtr() { return &atom_; }

    //! \return the index
    size_t index() const { return index_; }

    //! \return the coordinates
    const gmx::RVec x() const { return x_; }

    //! Equal index
    bool operator==(size_t index) const { return index == index_; }
};

void Topology::init(const std::vector<Bond> &bonds)
{
    for(const auto &b : bonds)
    {
        addBond(b);
    }
}

//! Print to msghandler
static void dump_entry(MsgHandler                *msghandler,
                       const TopologyEntryVector &entries,
                       const std::string         &label)
{
    if (msghandler->printLevel() != ACTStatus::Debug)
    {
        return;
    }
    for(auto &entry : entries)
    {
        std::string str = label;
        for (auto &ai : entry->atomIndices())
        {
            str.append(" " + std::to_string(ai));
        }
        msghandler->writeDebug(str);
    }
}

bool Topology::hasVsites() const
{
    return (hasEntry(InteractionType::VSITE1) ||
            hasEntry(InteractionType::VSITE2) ||
            hasEntry(InteractionType::VSITE2FD) ||
            hasEntry(InteractionType::VSITE3) ||
            hasEntry(InteractionType::VSITE3S) ||
            hasEntry(InteractionType::VSITE3FD) ||
            hasEntry(InteractionType::VSITE3FAD) ||
            hasEntry(InteractionType::VSITE3OUT) ||
            hasEntry(InteractionType::VSITE3OUTS));
}

void Topology::addShells(MsgHandler       *msghandler,
                         const ForceField *pd,
                         AtomList         *atomList)
{
    if (!pd->polarizable())
    {
        return;
    }
    /* Add Polarization to the plist. */
    TopologyEntryVector pols;
    auto &fs  = pd->findForcesConst(InteractionType::POLARIZATION);
    auto pol_name = potentialToParameterName(fs.potential());
    if (pol_name.empty())
    {
        msghandler->msg(ACTStatus::Fatal, "Polarization function not defined in FF file");
        return;
    }
    // Loop through the atomList.
    for (auto iter = atomList->begin(); iter != atomList->end(); iter = std::next(iter))
    {
        if (iter->atom().pType() == ActParticle::Atom || iter->atom().pType() == ActParticle::Vsite)
        {
            std::string atomtype(iter->atom().ffType());
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
                        continue;
                    }
                    // TODO Multiple shell support
                    TopologyEntry pp;
                    auto core  = iter->index();
                    auto shell = atomList->size();
                    pp.addAtom(core);
                    pp.addAtom(shell);
                    pp.addBondOrder(1.0);
                    pols.push_back(std::move(pp));

                    // Insert shell atom
                    auto shellName = iter->atom().name() + "s";
                    ActAtom newshell(shellName, "EP", ptype.id(), ActParticle::Shell,
                                     0, 0, charge, fa->row());
                    newshell.setResidueNumber(iter->atom().residueNumber());
                    auto olditer = iter;
                    iter = atomList->insert(std::next(iter), ActAtomListItem(newshell, shell, olditer->x()));
                    // Create link between shell and atom
                    olditer->atomPtr()->addShell(shell);
                    iter->atomPtr()->addCore(core);
                }
            }
            else
            {
                msghandler->msg(ACTStatus::Warning,
                                ACTMessage::AtomTypes,
                                gmx::formatString("Cannot find atomtype %s in forcefield\n",
                                                  atomtype.c_str()));
            }
        }
    }

    if (!pols.empty())
    {
        addEntry(InteractionType::POLARIZATION, pols);
        dump_entry(msghandler, pols, "The pols");
    }
}

void Topology::addBond(const Bond &bond)
{
    if (bond.aI() < 0 || bond.aJ() < 0)
    {
        GMX_THROW(gmx::InternalError("Negative atom indices when adding bonds"));
    }
    auto itb = InteractionType::BONDS;
    if (entries_.find(itb) == entries_.end())
    {
        entries_.insert({ itb, TopologyEntryVector{} });
    }
    entries_[itb].push_back(std::any_cast<Bond>(std::move(bond)));
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

Bond Topology::findBond(int ai, int aj) const
{
    Bond b(ai, aj, 1.0);
    auto bondsptr = entries_.find(InteractionType::BONDS);
    TopologyEntryVector::const_iterator bptr;
    if (entries_.end() == bondsptr)
    {
        GMX_THROW(gmx::InternalError("There are no bonds at all."));
    }
    for(bptr = bondsptr->second.begin(); bptr < bondsptr->second.end(); ++bptr)
    {
        if (((*bptr)->atomIndex(0) == ai && (*bptr)->atomIndex(1) == aj) ||
            ((*bptr)->atomIndex(0) == aj && (*bptr)->atomIndex(1) == ai))
        {
            break;
        }
    }
    if (bptr == bondsptr->second.end())
    {
        GMX_THROW(gmx::InternalError(gmx::formatString("Cannot find bond between %d and %d", ai, aj).c_str()));
    }

    auto c = (*bptr)->self();
    if ((*bptr)->atomIndex(0) == ai)
    {
        return *static_cast<const Bond *>(c);
    }
    else
    {
        return static_cast<const Bond *>(c)->swap();
    }
}

const TopologyEntry *Topology::findTopologyEntry(const TopologyEntryVector &entries,
                                                 const std::vector<int>    &aindex,
                                                 const std::vector<double> &bondOrder,
                                                 CanSwap                    cs) const
{
    TopologyEntry te;
    for (auto &a : aindex)
    {
        te.addAtom(a);
    }
    for (auto &b : bondOrder)
    {
        te.addBondOrder(b);
    }
    TopologyEntryVector::const_iterator bptr;
    for(bptr = entries.begin(); bptr < entries.end(); ++bptr)
    {
        bool same = te.atomIndices().size() == (*bptr)->atomIndices().size();
        if (same)
        {
            for(size_t i = 0; i < te.atomIndices().size() and same; i++)
            {
                same = same && (te.atomIndices()[i] == (*bptr)->atomIndices()[i]);
            }
            for(size_t i = 0; i < te.bondOrders().size() and same; i++)
            {
                same = same && (te.bondOrders()[i] == (*bptr)->bondOrders()[i]);
            }
        }
        if (!same && CanSwap::Yes == cs)
        {
            same = true;
            for(size_t i = 0; i < te.atomIndices().size() and same; i++)
            {
                same = same && (te.atomIndex(i) == (*bptr)->atomIndex((*bptr)->atomIndices().size()-1-i));
            }
            for(size_t i = 0; i < te.bondOrders().size() and same; i++)
            {
                same = same && (te.bondOrder(i) == (*bptr)->bondOrder((*bptr)->bondOrders().size()-1-i));
            }
        }
        if (same)
        {
            break;
        }
    }
    if (bptr == entries.end())
    {
        return nullptr;
    }
    else
    {
        return std::any_cast<const TopologyEntry *>((*bptr)->self());
    }
}

void Topology::makeAngles(MsgHandler                   *msghandler,
                          const ForceField             *pd,
                          const std::vector<gmx::RVec> &x,
                          double                        LinearAngleMin)
{
    auto ib = InteractionType::BONDS;
    if (entries_.find(ib) == entries_.end())
    {
        return;
    }
    auto &bonds = entry(ib);
    TopologyEntryVector angles{};
    TopologyEntryVector linangles{};
    for(size_t i = 0; i < bonds.size(); i++)
    {
        auto b1  = static_cast<const Bond &>(*bonds[i]->self());
        int  ai1 = bonds[i]->atomIndex(0);
        int  aj1 = bonds[i]->atomIndex(1);
        for(size_t j = 0; j < bonds.size(); j++)
        {
            if (i == j)
            {
                continue;
            }
            auto b2  = static_cast<const Bond &>(*bonds[j]->self());
            int  ai2 = bonds[j]->atomIndex(0);
            int  aj2 = bonds[j]->atomIndex(1);
            Angle angle{};
            // Check if the bonds share atoms
            if (aj1 == ai2)
            {
                angle.setBonds(b1, b2);
            }
            else if (aj1 == aj2)
            {
                angle.setBonds(b1, b2.swap());
            }
            else if (ai1 == ai2)
            {
                angle.setBonds(b1.swap(), b2);
            }
            else if (ai1 == aj2)
            {
                angle.setBonds(b1.swap(), b2.swap());
            }
            if (angle.atomIndices().size() == 3)
            {
                if (is_linear(x[angle.atomIndex(0)], x[angle.atomIndex(1)], x[angle.atomIndex(2)],
                              nullptr, LinearAngleMin))
                {
                    if (nullptr == findTopologyEntry(linangles, angle.atomIndices(),
                                                     angle.bondOrders(), CanSwap::Yes))
                    {
                        angle.setLinear(true);
                        linangles.push_back(std::move(angle));
                    }
                }
                else
                {
                    if (nullptr == findTopologyEntry(angles, angle.atomIndices(),
                                                     angle.bondOrders(), CanSwap::Yes))
                    {
                        angles.push_back(std::move(angle));
                    }
                }
            }
        }
    }
    if (!angles.empty())
    {
        entries_.insert({ InteractionType::ANGLES, std::move(angles) });
        if (pd)
        {
            setEntryIdentifiers(msghandler, pd, InteractionType::ANGLES);
        }
    }
    if (!linangles.empty())
    {
        entries_.insert({ InteractionType::LINEAR_ANGLES, std::move(linangles) });
        if (pd)
        {
            setEntryIdentifiers(msghandler, pd, InteractionType::LINEAR_ANGLES);
        }
    }
}

void Topology::makeImpropers(MsgHandler                   *msghandler,
                             const ForceField             *pd,
                             const std::vector<gmx::RVec> &x,
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
    TopologyEntryVector impropers{};
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
                    auto bbb = findBond(i, jkl[m]);
                    if (bbb.atomIndex(0) != static_cast<int>(i))
                    {
                        bjkl[m] = bbb.swap();
                    }
                    else
                    {
                        bjkl[m] = bbb;
                    }
                }
                impropers.push_back(Improper(bjkl[0], bjkl[1], bjkl[2]));
            }
        }
    }
    if (!impropers.empty())
    {
        entries_.insert({ InteractionType::IMPROPER_DIHEDRALS, std::move(impropers) });
        if (pd)
        {
            setEntryIdentifiers(msghandler, pd, InteractionType::IMPROPER_DIHEDRALS);
        }
    }
}

void Topology::makePairs(MsgHandler                          *msghandler,
                         const ForceField                    *pd,
                         InteractionType                      itype,
                         const std::vector<std::set<size_t>> &exclusions)
{
    TopologyEntryVector pairs{};
    for(size_t i = 0; i < atoms_.size(); i++)
    {
        for(size_t j = i+1; j < atoms_.size(); j++)
        {
            if (exclusions[i].find(j) == exclusions[i].end())
            {
                pairs.push_back(AtomPair(i, j));
            }
        }
    }
    // Finally insert the identifiers
    if (!pairs.empty())
    {
        entries_.insert({itype, std::move(pairs) });
        if (pd)
        {
            setEntryIdentifiers(msghandler, pd, itype);
        }
    }
}

void  Topology::addShellPairs()
{
    if (hasEntry(InteractionType::POLARIZATION))
    {
        auto &pol = entries_.find(InteractionType::POLARIZATION)->second;
        for (const auto &itype : { InteractionType::VDW,
                                   InteractionType::ELECTROSTATICS } )
        {
            if (!hasEntry(itype))
            {
                continue;
            }
            auto &ffpl = entries_.find(itype)->second;
            std::set<AtomPair> newf;
            for(auto &f : ffpl)
            {
                auto core_i = f->atomIndex(0);
                auto core_j = f->atomIndex(1);
                int shell_i = -1;
                int shell_j = -1;
                for(const auto &p_i : pol)

                {
                    if (core_i == p_i->atomIndex(0))
                    {
                        shell_i = p_i->atomIndex(1);
                        break;
                    }
                }
                for(const auto &p_j : pol)
                {
                    if (core_j == p_j->atomIndex(0))
                    {
                        shell_j = p_j->atomIndex(1);
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
                AtomPair ap(n);
                ffpl.push_back(std::move(ap));
            }
        }
    }
}

void Topology::makePropers(MsgHandler       *msghandler,
                           const ForceField *pd)
{
    auto ia = InteractionType::ANGLES;
    if (entries_.find(ia) == entries_.end())
    {
        return;
    }
    TopologyEntryVector propers{};
    auto &angles  = entries_.find(ia)->second;
    for(size_t i = 0; i < angles.size(); i++)
    {
        auto a1 = static_cast<const Angle *>(angles[i]->self());
        if (a1->isLinear())
        {
            continue;
        }
        auto bij1 = a1->bij();
        auto bjk1 = a1->bjk();
        for(size_t j = i+1; j < angles.size(); j++)
        {
            auto a2 = static_cast<const Angle *>(angles[j]->self());
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
                propers.push_back(Proper(bij1, bjk1, bjk2));
            }
            // 2) the first angle is k-j-i and the second j-k-l
            else if (bij1.aI() == bij2.aJ() && bij1.aJ() == bij2.aI())
            {
                propers.push_back(Proper(bjk1.swap(), bij2, bjk2));
            }
            // 2) the first angle is k-j-i and the second l-k-j
            else if (bij1.aI() == bij2.aJ() && bij1.aJ() == bjk2.aJ() )
            {
                propers.push_back(Proper(bjk1.swap(), bij1.swap(), bij2.swap()));
            }
            // 4) the first angle is i-j-k and the second l-k-j
            else if (bij1.aJ() == bjk2.aJ() && bjk1.aJ() == bij2.aJ())
            {
                propers.push_back(Proper(bij1, bjk1, bij2.swap()));
            }
        }
    }
    if (!propers.empty())
    {
        entries_.insert({ InteractionType::PROPER_DIHEDRALS, std::move(propers) });
        if (pd)
        {
            setEntryIdentifiers(msghandler, pd, InteractionType::PROPER_DIHEDRALS);
        }
    }
}

std::map<InteractionType, size_t> Topology::makeVsite1s(MsgHandler       *msghandler,
                                                        const ForceField *pd,
                                                        AtomList         *atomList)
{
    if (!pd)
    {
        GMX_THROW(gmx::InternalError("Why did you call makeVsite1s without a force field?"));
        return {};
    }
    auto itype = InteractionType::VSITE1;
    if (!pd->interactionPresent(itype))
    {
        return {};
    }
    const auto &fs = pd->findForcesConst(itype);
    if (fs.empty())
    {
        return {};
    }
    TopologyEntryVector v1top;
    for(auto atom = atomList->begin(); atom != atomList->end(); atom++)
    {
        auto aa = atom->atom();
        // First find the bond type for this atom
        const auto ptype1 = pd->findParticleType(aa.ffType());
        for(const auto &mm : fs.parametersConst())
        {
            // Use the type of the atom to compare to the force field
            auto faa = mm.first.atoms();
            if (ptype1->id().id() == faa[0])
            {
                // The vsite however, should have the same bond as atom type.
                const auto ptype = pd->findParticleType(faa[1]);
                std::string vstype;
                ActAtom newatom(ptype->id().id(), vstype, ptype->id().id(),
                                ptype->apType(),
                                0, ptype->mass(), ptype->charge(),
                                ptype->row());
                // Put virtual site straight after the last atom.
                int vs1 = atomList->size();
                newatom.addCore(atom->index());
                // Residue number
                newatom.setResidueNumber(aa.residueNumber());
                gmx::RVec vzero = { 0, 0, 0 };
                atomList->insert(std::next(atom), ActAtomListItem(newatom, vs1, vzero));
                // Create new topology entry
                Vsite1 vsnew(atom->index(), vs1);
                msghandler->writeDebug(gmx::formatString("Adding vs1 %s-%lu %d\n",
                                                         aa.element().c_str(), atom->index(), vs1));
                // Special bond order for vsites
                vsnew.addBondOrder(9);
                v1top.push_back(std::any_cast<Vsite1>(std::move(vsnew)));
            }
        }
    }
    // If we did find any vsite1 instances, add the whole vector to the topology.
    // A very subtle programming issue arises here:
    // after the std::move operation, the vector is empty
    // and therefore vsite1.size() == 0. Hence we have to store the size in a variable.
    std::map<InteractionType, size_t> num_v1;
    if (!v1top.empty())
    {
        num_v1.insert({ itype, v1top.size() });
        entries_.insert({ itype, std::move(v1top) });
    }
    return num_v1;
}

std::map<InteractionType, size_t> Topology::makeVsite2s(MsgHandler       *msghandler,
                                                        const ForceField *pd,
                                                        AtomList         *atomList)
{
    if (!pd)
    {
        GMX_THROW(gmx::InternalError("Why did you call makeVsites2 without a force field?"));
        return {};
    }
    std::vector<InteractionType> v2s = { InteractionType::VSITE2,
                                         InteractionType::VSITE2FD };
    std::map<InteractionType, ForceFieldParameterList> ffvs;
    for (const auto itype : v2s)
    {
        if (pd->interactionPresent(itype))
        {
            const auto &fs = pd->findForcesConst(itype);
            if (!fs.empty())
            {
                ffvs.insert({itype, fs });
            }
        }
    }
    if (ffvs.empty())
    {
        msghandler->msg(ACTStatus::Verbose,
                        "Force field does not contain any two particle virtual sites.");
        return {};
    }
    else
    {
        msghandler->msg(ACTStatus::Verbose,
                        gmx::formatString("There are %zu non-empty two particle vsite entries in the force field.",
                                          ffvs.size()));
    }
    auto itype_bonds = InteractionType::BONDS;
    if (entries_.find(itype_bonds) == entries_.end())
    {
        msghandler->writeDebug("There are no bonds to generate vsites3 from.");
        return {};
    }
    auto &bonds     = entry(itype_bonds);
    // Count the number of interactions of each type that we find
    std::map<InteractionType, size_t> num_v2;
    for (const auto &myffvs : ffvs)
    {
        TopologyEntryVector v2top;
        for (const auto &fvs : myffvs.second.parametersConst())
        {
            // Force field info for vsites2 type
            auto ff_vs2  = pd->findForcesConst(myffvs.first);
            auto vsatoms = fvs.first.atoms();
            auto vsbo    = fvs.first.bondOrders();

            // Loop over bonds to see if there are any pairs of atoms
            // that should be augmented with a virtual site.
            // Since we insert atoms and coordinates, the bonds array needs to be
            // interpreted with care.
            std::string bondtype("bondtype");
            for(size_t i = 0; i < bonds.size(); i++)
            {
                auto mybond = static_cast<const Bond *>(bonds[i]->self());
                mybond->check(2);
                auto bid    = mybond->id();
                auto border = bid.bondOrders();
                int  ai     = mybond->aI();
                int  aj     = mybond->aJ();
                auto pti    = pd->findParticleType(atoms_[ai].ffType());
                auto ptj    = pd->findParticleType(atoms_[aj].ffType());
                if (!pti->hasOption(bondtype))
                {
                    GMX_THROW(gmx::InternalError(gmx::formatString("Particle type %s has no bondtype option but is bonded to %s",
                                                                   atoms_[ai].ffType().c_str(),
                                                                   atoms_[aj].ffType().c_str()).c_str()));
                }
                else if (!ptj->hasOption(bondtype))
                {
                    GMX_THROW(gmx::InternalError(gmx::formatString("Particle type %s has no bondtype option but is bonded to %s",
                                                                   atoms_[aj].ffType().c_str(),
                                                                   atoms_[ai].ffType().c_str()).c_str()));
                }
                auto bai    = pti->optionValue("bondtype");
                auto baj    = ptj->optionValue("bondtype");
                
                msghandler->writeDebug(gmx::formatString("Found bond %s %s", bai.c_str(), baj.c_str()));

                bool found   = false;
                if (border[0] == vsbo[0])
                {
                    if (bai == vsatoms[0] && baj == vsatoms[1])
                    {
                        found = true;
                    }
                    else if (baj == vsatoms[0] && bai == vsatoms[1])
                    {
                        found   = true;
                        int tmp = ai; ai = aj; aj = tmp;
                    }
                }
                // Make dummmy bond identifier using just the two atoms and the real bond order.
                if (found)
                {
                    auto vsname = vsatoms[vsatoms.size()-1];
                    if (!pd->hasParticleType(vsname))
                    {
                        printf("No such particle type %s as found in vsite %s\n",
                               vsname.c_str(), fvs.first.id().c_str());
                    }
                    else
                    {
                        auto ptype = pd->findParticleType(vsname);
                        if (!ptype->hasOption(bondtype))
                        {
                            GMX_THROW(gmx::InternalError(gmx::formatString("particle type %s has no bondtype option", vsname.c_str()).c_str()));
                            
                        }
                        std::string vstype = ptype->optionValue(bondtype);
                        ActAtom newatom(ptype->id().id(), vstype, ptype->id().id(),
                                        ptype->apType(),
                                        0, ptype->mass(), ptype->charge(),
                                        ptype->row());
                        // Put virtual site straight after the last atom.
                        int vs2 = atomList->size();
                        newatom.addCore(ai);
                        newatom.addCore(aj);
                        // Residue number
                        newatom.setResidueNumber(atoms_[ai].residueNumber());
                        gmx::RVec vzero = { 0, 0, 0 };
                        size_t after = std::max(ai, aj);
                        auto iter = std::find(atomList->begin(), atomList->end(), after);
                        atomList->insert(std::next(iter), ActAtomListItem(newatom, vs2, vzero));
                        // Create new topology entry
                        Vsite2 vsnew(ai, aj, vs2);
                        msghandler->writeDebug(gmx::formatString("Adding vs2 %s-%d %s-%d %d",
                                                                 atoms_[ai].element().c_str(), ai,
                                                                 atoms_[aj].element().c_str(), aj, vs2));

                        // Add bond orders, copied from the bond.
                        for (auto b : border)
                        {
                            vsnew.addBondOrder(b);
                        }
                        // Special bond order for vsites
                        vsnew.addBondOrder(9);
                        v2top.push_back(std::any_cast<Vsite2>(std::move(vsnew)));
                    }
                }
            }
        }
        // If we did find any vsite2 instances, add the whole vector to the topology.
        // A very subtle programming issue arises here:
        // after the std::move operation, the vector is empty
        // and therefore vsite2.size() == 0. Hence we have to store the size in a variable.
        if (!v2top.empty())
        {
            num_v2.insert({ myffvs.first, v2top.size() });
            entries_.insert({ myffvs.first, std::move(v2top) });
        }
    }
    return num_v2;
}

std::map<InteractionType, size_t> Topology::makeVsite3s(MsgHandler       *msghandler,
                                                        const ForceField *pd,
                                                        AtomList         *atomList)
{
    if (!pd)
    {
        GMX_THROW(gmx::InternalError("Why did you call makeVsites3 without a force field?"));
        return {};
    }
    std::vector<InteractionType> v3s = { InteractionType::VSITE3,
                                         InteractionType::VSITE3S,
                                         InteractionType::VSITE3FD,
                                         InteractionType::VSITE3OUT,
                                         InteractionType::VSITE3OUTS };
    std::map<InteractionType, ForceFieldParameterList> ffvs;
    for (const auto itype : v3s)
    {
        if (pd->interactionPresent(itype))
        {
            const auto &fs = pd->findForcesConst(itype);
            if (!fs.empty())
            {
                ffvs.insert({itype, fs });
            }
        }
    }
    if (ffvs.empty())
    {
        msghandler->writeDebug("Force field does not contain any three particle virtual sites.");
        return {};
    }
    else
    {
        msghandler->writeDebug(gmx::formatString("There are %zu non-empty three particle vsite entries in the force field.", ffvs.size()));
    }
    auto itype_angles = InteractionType::ANGLES;
    if (entries_.find(itype_angles) == entries_.end())
    {
        msghandler->writeDebug("There are no angles to generate vsites3 from.");
        return {};
    }
    // Count the number of interactions of each type that we find
    std::map<InteractionType, size_t> num_v3;
    for (const auto &myffvs : ffvs)
    {
        TopologyEntryVector v3top;
        for (const auto &fvs : myffvs.second.parametersConst())
        {
            msghandler->writeDebug(gmx::formatString("Checking vsite %s", fvs.first.id().c_str()));

            auto vsatoms = fvs.first.atoms();
            auto vsbo    = fvs.first.bondOrders();
            auto &angles     = entry(itype_angles);
            for (size_t i = 0; i < angles.size(); i++)
            {
                auto mybond = static_cast<const Angle *>(angles[i]->self());

                mybond->check(3);
                auto bid    = mybond->id();
                auto border = bid.bondOrders();
                int  ai     = mybond->atomIndex(0);
                int  aj     = mybond->atomIndex(1);
                int  ak     = mybond->atomIndex(2);

                auto bai    = pd->findParticleType(atoms_[ai].ffType())->optionValue("bondtype");
                auto baj    = pd->findParticleType(atoms_[aj].ffType())->optionValue("bondtype");
                auto bak    = pd->findParticleType(atoms_[ak].ffType())->optionValue("bondtype");
                msghandler->writeDebug(gmx::formatString("Found angle %s %s %s", bai.c_str(), baj.c_str(), bak.c_str()));


                bool found   = false;
                if (border[0] == vsbo[0] && border[1] == vsbo[1] &&
                    bai == vsatoms[0] && baj == vsatoms[1] && bak == vsatoms[2])
                {
                    found = true;
                }
                else if (baj == vsatoms[1] && bak == vsatoms[0] &&
                         bai == vsatoms[2] && border[1] == vsbo[0] && border[0] == vsbo[1] )
                {
                    found   = true;
                    int tmp = ak; ak = ai; ai = tmp;
                }
                if (!found)
                {
                    continue;
                }

                auto vsname = vsatoms[vsatoms.size() - 1];
                if (!pd->hasParticleType(vsname))
                {
                    printf("No such particle type %s as found in vsite %s\n",
                           vsname.c_str(), fvs.first.id().c_str());
                }
                else
                {
                    auto ptype = pd->findParticleType(vsname);
                    std::string vstype = ptype->optionValue("bondtype");
                    // Determine how many particles to add
                    int maxpid = 1;
                    if (InteractionType::VSITE3OUT == myffvs.first ||
                        InteractionType::VSITE3OUTS == myffvs.first)
                    {
                        maxpid = 2;
                    }
                    for (int pid=0; pid<maxpid; pid++)
                    {
                        ActAtom newatom(ptype->id().id(), vstype, ptype->id().id(),
                                        ptype->apType(),
                                        0, ptype->mass(), ptype->charge(), ptype->row());

                        int vs3 = atomList->size();
                        newatom.addCore(ai);
                        newatom.addCore(aj);
                        newatom.addCore(ak);
                        newatom.setResidueNumber(atoms_[ai].residueNumber());
                        msghandler->writeDebug(gmx::formatString("Adding %s %s%d %s%d %s%d %d",
                                                                 interactionTypeToString(myffvs.first).c_str(),
                                                                 atoms_[ai].element().c_str(), ai,
                                                                 atoms_[aj].element().c_str(), aj,
                                                                 atoms_[ak].element().c_str(), ak, vs3));

                        gmx::RVec vzero = {0, 0, 0};
                        size_t after = std::max({ai, aj, ak});
                        auto iter = std::find(atomList->begin(), atomList->end(), after);
                        atomList->insert(std::next(iter), ActAtomListItem(newatom, vs3, vzero));

                        // Create new topology entry
                        switch(myffvs.first)
                        {
                        case InteractionType::VSITE3:
                        case InteractionType::VSITE3S:
                        case InteractionType::VSITE3FD:
                            {
                                Vsite3 vsnew(ai, aj, ak, vs3);
                                // Add bond orders, cp from the angle.
                                for (auto b : border)
                                {
                                    vsnew.addBondOrder(b);
                                }

                                // Special bond order for vsites
                                vsnew.addBondOrder(9);
                                v3top.push_back(std::any_cast<Vsite3>(std::move(vsnew)));
                            }
                            break;
                        case InteractionType::VSITE3OUT:
                        case InteractionType::VSITE3OUTS:
                            {
                                // We are creating two of these, with different sign on the c parameter.
                                int       sign = 2*pid - 1;
                                Vsite3OUT vsnew(ai, aj, ak, vs3, sign);
                                // Add bond orders, cp from the angle.
                                for (auto b : border)
                                {
                                    vsnew.addBondOrder(b);
                                }

                                // Special bond order for vsites
                                vsnew.addBondOrder(9);
                                v3top.push_back(std::any_cast<Vsite3OUT>(std::move(vsnew)));
                            }
                            break;
                        default: // This is ok?
                            break;
                        }
                    }
                }
            }
        }

        if (!v3top.empty())
        {
            num_v3.insert({ myffvs.first, v3top.size() });
            entries_.insert({ myffvs.first, std::move(v3top) });
        }
    }
    return num_v3;
}

std::map<InteractionType, size_t> Topology::makeVsite4s(MsgHandler       *msghandler,
                                                        const ForceField *pd,
                                                        AtomList         *atomList)
{
    if (!pd)
    {
        GMX_THROW(gmx::InternalError("Why did you call makeVsites4 without a force field?"));
        return {};
    }
    std::vector<InteractionType> v4s = { InteractionType::VSITE4,
                                         InteractionType::VSITE4S,
                                         InteractionType::VSITE4S3 };
    std::map<InteractionType, ForceFieldParameterList> ffvs;
    for (const auto itype : v4s)
    {
        if (pd->interactionPresent(itype))
        {
            const auto &fs = pd->findForcesConst(itype);
            if (!fs.empty())
            {
                ffvs.insert({itype, fs });
            }
        }
    }
    if (ffvs.empty())
    {
        msghandler->writeDebug("Force field does not contain any four particle virtual sites.");
        return {};
    }
    else
    {
        msghandler->writeDebug(gmx::formatString("There are %zu non-empty four particle vsite entries in the force field.", ffvs.size()));
    }
    auto itype_bonds = InteractionType::BONDS;
    if (entries_.find(itype_bonds) == entries_.end())
    {
        msghandler->writeDebug("There are no bonds to generate vsites4 from.");
        return {};
    }
    auto &bonds  = entry(itype_bonds);
    // Count the number of interactions of each type that we find
    std::map<InteractionType, size_t> num_v4;
    for (const auto &myffvs : ffvs)
    {
        TopologyEntryVector v4top;
        for (const auto &fvs : myffvs.second.parametersConst())
        {
            msghandler->writeDebug(gmx::formatString("Checking for vsite %s",
                                                     fvs.first.id().c_str()));
            auto vsatoms = fvs.first.atoms();
            auto vsbo    = fvs.first.bondOrders();
            // Loop over atoms
            for(auto aaa = atomList->begin(); aaa != atomList->end(); ++aaa)
            {
                // Shortcut
                auto atom = aaa->atom();
                if (atom.pType() != ActParticle::Atom)
                {
                    continue;
                }
                std::string btype("bondtype");
                auto ptype = pd->findParticleType(atom.ffType());
                if (!ptype->hasOption(btype))
                {
                    continue;
                }
                auto atom_b = ptype->optionValue(btype);
                if (atom_b != vsatoms[0])
                {
                    continue;
                }
                // Found a matching first atom bondtype, let's check whether its bonds match the vsite
                std::vector<int> vsbond(vsatoms.size(), -1);
                vsbond[0] = aaa->index();
                for(const auto &b : bonds)
                {
                    int aj = -1;
                    if (static_cast<size_t>(b->atomIndex(0)) == aaa->index())
                    {
                        aj = b->atomIndex(1);
                    }
                    else if (static_cast<size_t>(b->atomIndex(1)) == aaa->index())
                    {
                        aj = b->atomIndex(1);
                    }
                    if (aj == -1)
                    {
                        continue;
                    }
                    // Lookup aj
                    for(auto a2 = atomList->begin(); a2 != atomList->end(); ++a2)
                    {
                        // Complicated way of comparing size_t and int
                        if (a2->index() - aj == 0)
                        {
                            for (size_t i = 1; i < vsatoms.size(); i++)
                            {
                                if (a2->atom().id().id() == vsatoms[i])
                                {
                                    if (vsbond[i] == -1)
                                    {
                                        vsbond[i] = aj;
                                    }
                                }
                            }
                        }
                    }
                }
                bool foundAllBonds = true;
                for(const auto nb : vsbond)
                {
                    foundAllBonds = foundAllBonds && (nb >= 0);
                }
                if (foundAllBonds)
                {
                    // Yay! We found a four particle virtual site
                    msghandler->writeDebug(gmx::formatString("Found vsite4 %d %d %d %d",
                                                             vsbond[0], vsbond[1],
                                                             vsbond[2], vsbond[3]));
                    auto vsname = vsatoms[vsatoms.size() - 1];
                    if (!pd->hasParticleType(vsname))
                    {
                        printf("No such particle type %s as found in vsite %s\n",
                               vsname.c_str(), fvs.first.id().c_str());
                    }
                    else
                    {
                        auto ptype = pd->findParticleType(vsname);
                        std::string vstype = ptype->optionValue("bondtype");
                        // Add one particle
                        ActAtom newatom(ptype->id().id(), vstype, ptype->id().id(),
                                        ptype->apType(),
                                        0, ptype->mass(), ptype->charge(), ptype->row());

                        int vs4 = atomList->size();
                        
                        newatom.addCore(vsbond[0]);
                        newatom.addCore(vsbond[1]);
                        newatom.addCore(vsbond[2]);
                        newatom.addCore(vsbond[3]);
                        newatom.setResidueNumber(atoms_[vsbond[0]].residueNumber());
                        msghandler->writeDebug(gmx::formatString("Adding %s %s%d %s%d %s%d %s%d %d",
                                                                 interactionTypeToString(myffvs.first).c_str(),
                                                                 atoms_[vsbond[0]].element().c_str(), vsbond[0],
                                                                 atoms_[vsbond[1]].element().c_str(), vsbond[1],
                                                                 atoms_[vsbond[2]].element().c_str(), vsbond[2], 
                                                                 atoms_[vsbond[3]].element().c_str(), vsbond[3],
                                                                 vs4));

                        gmx::RVec vzero = {0, 0, 0};
                        size_t after = std::max({vsbond[0], vsbond[1], vsbond[2], vsbond[3]});
                        auto iter = std::find(atomList->begin(), atomList->end(), after);
                        atomList->insert(std::next(iter), ActAtomListItem(newatom, vs4, vzero));

                        // Create new topology entry
                        {
                            Vsite4 vsnew(vsbond[0], vsbond[1], vsbond[2], vsbond[3], vs4);
                            // Add bond orders, cp from the FF
                            for (auto b : vsbo)
                            {
                                vsnew.addBondOrder(b);
                            }
                            
                            // Special bond order for vsites
                            vsnew.addBondOrder(9);
                            v4top.push_back(std::any_cast<Vsite4>(std::move(vsnew)));
                        }
                    }
                }
            }
        }

        if (!v4top.empty())
        {
            num_v4.insert({ myffvs.first, v4top.size() });
            entries_.insert({ myffvs.first, std::move(v4top) });
        }
    }
    return num_v4;
}

void Topology::addVsitesToCores()
{
    for(size_t i = 0; i < atoms_.size(); i++)
    {
        if (ActParticle::Vsite == atoms_[i].pType())
        {
            for(const auto cj : atoms_[i].cores())
            {
                atoms_[cj].addVsite(i);
            }
        }
    }
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

void Topology::dumpPairlist(gmx::TextWriter *tw, InteractionType itype) const
{
    if (nullptr == tw)
    {
        return;
    }
    if (!hasEntry(itype))
    {
        return;
    }
    for(const auto &pl : entry(itype))
    {
        tw->writeLine(gmx::formatString("PAIRLIST %s %d %d",
                                        interactionTypeToString(itype).c_str(),
                                        pl->atomIndex(0), pl->atomIndex(1)));
    }
}

void Topology::GenerateAtoms(MsgHandler             *msghandler,
                             const ForceField       *pd,
                             const MolProp          *mol,
                             std::vector<gmx::RVec> *x)
{
    ACTMessage imm   = ACTMessage::OK;

    auto ci = mol->findExperimentConst(JobType::OPT);
    if (!ci)
    {
        ci = mol->findExperimentConst(JobType::TOPOLOGY);
    }
    if (!ci)
    {
        ci = mol->findExperimentConst(JobType::SP);
    }
    if (!ci)
    {
        msghandler->msg(ACTStatus::Error, ACTMessage::Topology,
                        gmx::formatString("No molecule %s", mol->getMolname().c_str()).c_str());
        return;
    }
    if (ci)
    {
        if (ci->NAtom() == 0)
        {
            msghandler->msg(ACTStatus::Error, ACTMessage::NoAtoms, "No atoms");
            return;
        }
        int res0  = -666;
        int nres  = 0;
        int natom = 0;
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

                ActAtom newatom(cai.getName(), atype->element(), cai.getObtype(), ActParticle::Atom,
                                atype->atomnumber(), atype->mass(), atype->charge(), atype->row());
                newatom.setResidueNumber(nres-1);
                addAtom(newatom);
                realAtoms_.push_back(natom++);
            }
            else
            {
                msghandler->msg(ACTStatus::Error,
                                ACTMessage::AtomTypes,
                                gmx::formatString("Cannot find atomtype %s (atom %zu) in forcefield, there are %d atomtypes.",
                                                  cai.getObtype().c_str(), atoms_.size(), pd->nParticleTypes()).c_str());

                return;
            }
        }
    }
    else
    {
        msghandler->msg(ACTStatus::Error, ACTMessage::Topology,
                        "No structure to make topology for.");
    }
    msghandler->writeDebug(gmx::formatString("Tried to convert '%s' to ACT. LOT is '%s/%s'. Natoms is %zu. Result: %s.",
                                             mol->getMolname().c_str(),
                                             ci->getMethod().c_str(),
                                             ci->getBasisset().c_str(), atoms_.size(),
                                             actMessage(imm)).c_str());
}

void Topology::build(MsgHandler             *msghandler,
                     const ForceField       *pd,
                     std::vector<gmx::RVec> *x,
                     double                  LinearAngleMin,
                     double                  PlanarAngleMax,
                     missingParameters       missing)
{
    // We use a list structure to keep track of the atoms and their coordinates.
    // This is needed since we may insert both shells and virtual sites. Although
    // one can insert items in std::vector structures it may lead to poor (N^2)
    // scaling to insert for instance one shell per particle. At the end of this
    // algorithm we copy the contents of the list back to std::vectors.
    AtomList atomList;
    auto     nRealAtoms = atoms_.size();
    for(size_t index = 0; index < atoms_.size(); index++)
    {
        ActAtomListItem item(atoms_[index], index, (*x)[index]);
        atomList.push_back(std::move(item));
    }

    // Before we can do anything we need to "identify" the bonds
    setEntryIdentifiers(msghandler, pd, InteractionType::BONDS);

    // Check whether we have virtual sites in the force field.
    auto nv1 = makeVsite1s(msghandler, pd, &atomList);
    auto nv2 = makeVsite2s(msghandler, pd, &atomList);

    // Before we can make three-particle vsites, we need to create
    // angles, but only temporarily.
    makeAngles(msghandler, pd, *x, LinearAngleMin);
    auto nv3 = makeVsite3s(msghandler, pd, &atomList);
    auto nv4 = makeVsite4s(msghandler, pd, &atomList);
    if (msghandler->debug() && (nv2.size() > 0 || nv3.size() > 0 || nv4.size() > 0) )
    {
        std::string str = "Added";
        for(const auto &mm : { nv1, nv2, nv3, nv4 })
        {
            for(const auto &nn : mm)
            {
                str += " " + std::to_string(nn.second) + " " + interactionTypeToString(nn.first);
            }
        }
        msghandler->writeDebug(str);
    }
    // Now throw away the angles again since we need to reorder the
    // lists of atoms and vsites.
    auto itype_a = InteractionType::ANGLES;
    if (entries_.end() != entries_.find(itype_a))
    {
        entries_.erase(itype_a);
    }

    // Now time for shells.
    addShells(msghandler, pd, &atomList);
    if (!msghandler->ok())
    {
        return;
    }
    if (msghandler->debug())
    {
        auto nshell = atomList.size()-nRealAtoms;
        for(const auto &mm : { nv1, nv2, nv3 })
        {
            for(const auto &nn : mm)
            {
                nshell -= nn.second;
            }
        }
        msghandler->writeDebug(gmx::formatString("Added %zu shells", nshell));
    }
    // Now there will be no more changes to the atomList and we can copy it back.
    std::vector<int> renumber(atomList.size(), 0);
    size_t atom = 0;
    atoms_.clear();
    x->clear();
    for(auto iter = atomList.begin(); iter != atomList.end(); iter = std::next(iter))
    {
        atoms_.push_back(iter->atom());
        (*x).push_back(iter->x());
        renumber[iter->index()] = atom++;
    }
    // Renumber the list of real atoms
    for (size_t i = 0; i < realAtoms_.size(); i++)
    {
        realAtoms_[i] = renumber[realAtoms_[i]];
    }
    // Now renumber the atoms' internal links to cores and shells
    for (auto &atom: atoms_)
    {
        auto ccc = atom.coresPtr();
        for(auto c = ccc->begin(); c < ccc->end(); ++c)
        {
            *c = renumber[*c];
        }
        auto sss = atom.shellsPtr();
        for(auto s = sss->begin(); s < sss->end(); ++s)
        {
            *s = renumber[*s];
        }
    }
    // Renumber the atoms in the TopologyEntries that have been created so far.
    renumberAtoms(renumber);
    setEntryIdentifiers(msghandler, pd, InteractionType::POLARIZATION);
    for(const auto nn : nv1)
    {
        setEntryIdentifiers(msghandler, pd, nn.first);
    }
    for(const auto nn : nv2)
    {
        setEntryIdentifiers(msghandler, pd, nn.first);
    }
    for(const auto nn : nv3)
    {
        setEntryIdentifiers(msghandler, pd, nn.first);
    }
    // Add vsite ids to cores
    addVsitesToCores();

    // Now make angles etc.
    makeAngles(msghandler, pd, *x, LinearAngleMin);
    makeImpropers(msghandler, pd, *x, PlanarAngleMax);
    // Check whether we have dihedrals in the force field.
    if (pd->interactionPresent(InteractionType::PROPER_DIHEDRALS))
    {
        // Store temporary address of variable for performance
        auto &fs = pd->findForcesConst(InteractionType::PROPER_DIHEDRALS);
        if (!fs.empty() || missingParameters::Generate == missing)
        {
            makePropers(msghandler, pd);
        }
    }
    int nexcl;
    auto myItype = InteractionType::ELECTROSTATICS;
    if (!ffOption(*pd, myItype, "nexcl", &nexcl))
    {
        msghandler->fatal(gmx::formatString("The number of exclusions is not specified for %s in %s", interactionTypeToString(myItype).c_str(), pd->filename().c_str()));
    }
    auto exclusions = generateExclusions(msghandler, nexcl);
    makePairs(msghandler, pd, InteractionType::VDW, exclusions);
    makePairs(msghandler, pd, InteractionType::ELECTROSTATICS, exclusions);
    auto itqt = InteractionType::VDWCORRECTION;
    // If the interaction has no parameters even though it is present, ignore
    if (pd->interactionPresent(itqt) && !pd->findForcesConst(itqt).empty())
    {
        makePairs(msghandler, pd, itqt, exclusions);
    }
    auto itic = InteractionType::INDUCTIONCORRECTION;
    // If the interaction has no parameters even though it is present, ignore
    if (pd->interactionPresent(itic) && !pd->findForcesConst(itic).empty())
    {
        makePairs(msghandler, pd, itic, exclusions);
    }
    if (missing != missingParameters::Generate)
    {
        fillParameters(msghandler, pd, missing);
    }
    if (msghandler->debug())
    {
        dumpPairlist(msghandler->twDebug(), InteractionType::ELECTROSTATICS);
        dumpPairlist(msghandler->twDebug(), InteractionType::VDW);
        dumpPairlist(msghandler->twDebug(), InteractionType::VDWCORRECTION);
        dumpPairlist(msghandler->twDebug(), InteractionType::INDUCTIONCORRECTION);
    }
}

const TopologyEntryVector &Topology::entry(InteractionType itype) const
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

void Topology::addEntry(InteractionType            itype,
                        const TopologyEntryVector &entry)
{
    if (hasEntry(itype))
    {
        GMX_THROW(gmx::InternalError(gmx::formatString("InteractionType %s already present",
                                                       interactionTypeToString(itype).c_str()).c_str()));
    }
    if (entry.size() > 0)
    {
        // Make a copy first.
        TopologyEntryVector newv;
        for(const auto &ee : entry)
        {
            newv.push_back(ee);
        }
        entries_.insert({ itype, std::move(newv) });
    }
}

std::vector<std::set<size_t>> Topology::generateExclusions(MsgHandler *msghandler,
                                                           int         nrexcl)
{
    msghandler->msg(ACTStatus::Debug,
                    gmx::formatString("Will generate %d exclusions", nrexcl));
    std::vector<std::set<size_t>> exclusions;
    exclusions.resize(atoms_.size());
    // First directly connected shells and vsites
    for (size_t i = 0; i < atoms_.size(); i++)
    {
        if (atoms_[i].pType() == ActParticle::Atom)
        {
            for (const auto s : atoms_[i].shells())
            {
                // Shells typically have one core but vsites can have more than one
                exclusions[i].insert(s);
                exclusions[s].insert(i);
            }
            for (const auto v : atoms_[i].vsites())
            {
                // Shells typically have one core but vsites can have more than one
                exclusions[i].insert(v);
                exclusions[v].insert(i);
            }
        }    
    }
    auto bonds = entries_.find(InteractionType::BONDS);
    if (bonds != entries_.end())
    {
        if (nrexcl > 0)
        {
            // Add first neighbors, assuming these are atoms only
            for(auto &b : bonds->second)
            {
                auto a = b->atomIndices();
                exclusions[a[0]].insert(a[1]);
                exclusions[a[1]].insert(a[0]);
            }
        }
        if (nrexcl > 1)
        {
            for(const auto itype : { InteractionType::ANGLES, InteractionType::LINEAR_ANGLES })
            {
                auto angles = entries_.find(itype);
                if (angles != entries_.end())
                {
                    // Add first neighbors, assuming these are atoms only
                    for(auto &b : angles->second)
                    {
                        auto a = b->atomIndices();
                        exclusions[a[0]].insert(a[2]);
                        exclusions[a[2]].insert(a[0]);
                    }
                }
            }
        }
        if (nrexcl > 2)
        {
            auto dihs = entries_.find(InteractionType::PROPER_DIHEDRALS);
            if (dihs != entries_.end())
            {
                // Add first neighbors, assuming these are atoms only
                for(auto &b : dihs->second)
                {
                    auto a = b->atomIndices();
                    exclusions[a[0]].insert(a[3]);
                    exclusions[a[3]].insert(a[0]);
                }
            }
        }
        if (nrexcl > 3)
        {
            msghandler->msg(ACTStatus::Warning, "Cannot handle more than 3 exclusions.");
        }
    }
    // Now replicate the exclusions of atoms to their vsites and shells
    for (size_t i = 0; i < atoms_.size(); i++)
    {
        if (atoms_[i].pType() == ActParticle::Atom)
        {
            continue;
        }
        for (const auto cc : atoms_[i].cores())
        {
            for (const auto c : exclusions[cc])
            {
                // Shells typically have one core but vsites can have more than one
                exclusions[i].insert(c);
                exclusions[c].insert(i);
                // Now add the shells of c
                for (const auto s : atoms_[c].shells())
                {
                    exclusions[i].insert(s);
                    exclusions[s].insert(i);
                }
                // Now add the vsites of c
                for (const auto v : atoms_[c].vsites())
                {
                    exclusions[i].insert(v);
                    exclusions[v].insert(i);
                }
            }
        }
    }
    if (msghandler->printLevel() == ACTStatus::Debug)
    {
        auto tw = msghandler->twDebug();
        for(size_t i = 0; i < atoms_.size(); i++)
        {
            std::string str = gmx::formatString("%-10s  %-10s %2zu exclusions:",
                                                atoms_[i].name().c_str(),
                                                actParticleToString(atoms_[i].pType()).c_str(), i);
            for (const auto j : exclusions[i])
            {
                str += gmx::formatString("  %2zu", j);
            }
            tw->writeLine(str);
        }
    }
    return exclusions;
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

//! \brief Low-level parameter filling
static void fillParams(MsgHandler                      *msghandler,
                       const ForceFieldParameterList   &fs,
                       const Identifier                &btype,
                       int                              independent,
                       const std::vector<const char *> &param_names,
                       std::vector<double>             *param)
{
    auto ff = fs.findParameterMapConst(btype);
    std::string msg = potentialToString(fs.potential()) + " " + btype.id();
    if (ff.empty())
    {
        msghandler->msg(ACTStatus::Error, ACTMessage::MissingFFParameter, msg);
        return;
    }
    if (param->empty())
    {
        param->resize(param_names.size(), 0);
    }
    int found = 0;
    for (size_t i = 0; i < param_names.size(); i++)
    {
        auto fp      = ff.find(param_names[i]);
        double value = 0;
        if (ff.end() != fp)
        {
            value  = fp->second.internalValue();
            found += 1;
        }
        (*param)[i] = value;
    }
    if (found != independent)
    {
        msg += gmx::formatString(" found %d, expected %d independent parameters", found, independent);
        msghandler->msg(ACTStatus::Warning, ACTMessage::MissingFFParameter, msg);
    }
}

void Topology::fillParameters(MsgHandler        *msghandler,
                              const ForceField  *pd,
                              missingParameters  missing)
{
    for(auto &entry : entries_)
    {
        if (!pd->interactionPresent(entry.first))
        {
            continue;
        }
        auto &fs = pd->findForcesConst(entry.first);

        // Loop over entries, increment of tp handled at end.
        for(auto tp = entry.second.begin(); tp < entry.second.end(); )
        {
            auto &topentry = *tp;

            const auto          &topID = topentry->id();
            std::vector<double>  param;

            auto pp = potprops.find(fs.potential());
            if (potprops.end() != pp)
            {
                fillParams(msghandler, fs, topID, 
                           pp->second.param.size(), pp->second.param, &param);
            }
            else
            {
                GMX_THROW(gmx::InternalError(gmx::formatString("Missing case %s when filling the topology structure.",
                                                               potentialToString(fs.potential()).c_str()).c_str()));
            }
            if (param.empty())
            {
                tp = entry.second.erase(tp);
            }
            else
            {
                if (!msghandler->ok())
                {
                    auto msg = gmx::formatString("Force field only contains %zu/%zu %s parameters for %s, using zero's.", 
                                                 param.size(),
                                                 pp->second.param.size(),
                                                 potentialToString(fs.potential()).c_str(), topID.id().c_str());
                    if (missing == missingParameters::Error)
                    {
                        msghandler->msg(ACTStatus::Error, ACTMessage::MissingFFParameter, msg);
                        return;
                    }
                    else
                    {
                        msghandler->msg(ACTStatus::Warning, ACTMessage::MissingFFParameter, msg);
                    }
                    tp = entry.second.erase(tp);
                }
                else
                {
                    topentry->setParams(param);
                    ++tp;
                }
            }
        }
    }
}

void Topology::setEntryIdentifiers(MsgHandler       *msghandler,
                                   const ForceField *pd,
                                   InteractionType   itype)
{
    if (!hasEntry(itype) || !pd->interactionPresent(itype))
    {
        return;
    }
    auto fs     = pd->findForcesConst(itype);
    auto iEntry = entries_.find(itype);
    if (entries_.end() == iEntry)
    {
        GMX_THROW(gmx::InternalError(gmx::formatString("Cannot find entry %s in topology",
                                                       interactionTypeToString(itype).c_str()).c_str()));
    }
    for(auto topentry = iEntry->second.begin(); topentry < iEntry->second.end(); ++topentry)
    {
        std::string              atypes;
        std::vector<std::string> btype;
        for(const size_t jj : (*topentry)->atomIndices())
        {
            if (jj >= atoms_.size())
            {
                GMX_THROW(gmx::InternalError(gmx::formatString("Atom index %zu should be less than %zu for %s",
                                                               jj, atoms_.size(),
                                                               interactionTypeToString(itype).c_str()).c_str()
                                             ));
            }
            auto atype = pd->findParticleType(atoms_[jj].ffType());
            if (!atypes.empty())
            {
                atypes += " ";
            }
            atypes += atype->id().id();
            switch (itype)
            {
            case InteractionType::VDW:
            case InteractionType::INDUCTIONCORRECTION:
            case InteractionType::VDWCORRECTION:
            case InteractionType::POLARIZATION:
            case InteractionType::ELECTROSTATICS:
                {
                    // For COULOMB there are two particles,
                    // but for POLARIZATION or VDWCORRECTION just one.
                    if (atype->hasInteractionType(itype))
                    {
                        btype.push_back(atype->interactionTypeToIdentifier(itype).id());
                    }
                    break;
                }
            case InteractionType::VSITE1:
                // Special treatment for vsite1, see https://github.com/AlexandriaChemistry/ACT/issues/962
                {
                    btype.push_back(atype->id().id());
                }
                break;
            default: // does something
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
                (*topentry)->setId(Identifier(btype[0]));
            }
            else
            {
                (*topentry)->setId({ Identifier(btype, (*topentry)->bondOrders(),
                                                fs.canSwap()) });
            }
        }
        else if (msghandler->debug())
        {
            msghandler->writeDebug(gmx::formatString("Could not find identifier for %s for atomtype '%s'",
                                                     interactionTypeToString(itype).c_str(), atypes.c_str()));
        }
    }
}

const std::vector<std::vector<int>> &Topology::exclusions(InteractionType itype) const
{
    auto exclptr = exclusions_.find(itype);
    if (exclusions_.end() != exclptr)
    {
        return exclptr->second;
    }
    GMX_THROW(gmx::InternalError(gmx::formatString("No exclusions for interaction type %s", interactionTypeToString(itype).c_str()).c_str()));
}

} // namespace alexandria

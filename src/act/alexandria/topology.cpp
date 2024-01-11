/*
 * This source file is part of the Alexandria Chemistry Toolkit.
 *
 * Copyright (C) 2021-2024
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
#include <vector>

#include "act/basics/interactiontype.h"
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
    size_t    index_;
    //! Atomic coordinates
    gmx::RVec x_;
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

Topology::Topology(const std::vector<Bond> &bonds)
{
    for(const auto &b : bonds)
    {
        addBond(b);
    }
}

static void dump_entry(FILE                      *fp,
                       const TopologyEntryVector &entries,
                       const std::string         &label)
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

void Topology::addShells(const ForceField *pd,
                         AtomList         *atomList)
{
    if (!pd->polarizable())
    {
        return;
    }
    /* Add Polarization to the plist. */
    TopologyEntryVector pols;
    auto &fs  = pd->findForcesConst(InteractionType::POLARIZATION);

    // Loop through the atomList.
    for (auto iter = atomList->begin(); iter != atomList->end(); iter = std::next(iter))
    {
        if (iter->atom().pType() == eptAtom || iter->atom().pType() == eptVSite)
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
                        GMX_THROW(gmx::InvalidInputError(gmx::formatString("Polarizability should be positive for %s", fa->id().id().c_str()).c_str()));
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
                    ActAtom newshell(shellName, "EP", ptype.id(), eptShell,
                                     0, 0, charge);
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
                GMX_THROW(gmx::InvalidInputError(gmx::formatString("Cannot find atomtype %s in forcefield\n",
                                                                   atomtype.c_str())));
            }
        }
    }

    if (!pols.empty())
    {
        addEntry(InteractionType::POLARIZATION, pols);
        dump_entry(debug, pols, "The pols");
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

void Topology::makeAngles(const ForceField             *pd,
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
            setEntryIdentifiers(pd, InteractionType::ANGLES);
        }
    }
    if (!linangles.empty())
    {
        entries_.insert({ InteractionType::LINEAR_ANGLES, std::move(linangles) });
        if (pd)
        {
            setEntryIdentifiers(pd, InteractionType::LINEAR_ANGLES);
        }
    }
}

void Topology::makeImpropers(const ForceField             *pd,
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
            setEntryIdentifiers(pd, InteractionType::IMPROPER_DIHEDRALS);
        }
    }
}

void Topology::makePairs(const ForceField *pd,
                         InteractionType   itype)
{
    TopologyEntryVector pairs{};
    for(size_t i = 0; i < atoms_.size(); i++)
    {
        // Check for exclusions is done later.
        for(size_t j = i+1; j < atoms_.size(); j++)
        {
            pairs.push_back(AtomPair(i, j));
        }
    }
    if (!pairs.empty())
    {
        // Now time for exclusions
        int nexcl;
        if (!ffOption(*pd, itype, "nexcl", &nexcl))
        {
            nexcl = 0;
        }
        //! Non bonded exclusions, array is length of number of atoms
        auto exclusions = generateExclusions(&pairs, nexcl);
        fixExclusions(&pairs, exclusions);
        // Finally insert the identifiers
        if (!pairs.empty())
        {
            entries_.insert({itype, std::move(pairs) });
            if (pd)
            {
                setEntryIdentifiers(pd, itype);
            }
        }
    }
}

void Topology::fixExclusions(TopologyEntryVector                 *pairs,
                             const std::vector<std::vector<int>> &exclusions)
{
    // Loop over all exclusions
    for(size_t ai = 0; ai < exclusions.size(); ++ai)
    {
        if (eptAtom != atoms_[ai].pType())
        {
            continue;
        }
        for(size_t jj = 0; jj < exclusions[ai].size(); ++jj)
        {
            size_t aj = exclusions[ai][jj];
            if (eptAtom != atoms_[aj].pType())
            {
                continue;
            }
            // Check whether these particles have shells
            for(size_t si : atoms_[ai].shells())
            {
                for(size_t sj : atoms_[aj].shells())
                {
                    // See whether this interaction exists
                    auto it = pairs->begin();
                    while (pairs->end() != it)
                    {
                        size_t aai = (*it)->atomIndex(0);
                        size_t aaj = (*it)->atomIndex(1);
                        if (((aai == si || aai == ai) && (aaj == sj || aaj == aj)) ||
                            ((aaj == si || aaj == ai) && (aai == sj || aai == aj)))
                        {
                            it = pairs->erase(it);
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
    // Now remove interactions with vsites that are excluded
    // from the constructing atoms.
    for(size_t i = 0; i < atoms_.size(); i++)
    {
        if (eptVSite == atoms_[i].pType())
        {
            // Each vsite has two or more cores
            for (size_t core : atoms_[i].cores())
            {
                std::vector<int> core_and_shells = atoms_[core].shells();
                core_and_shells.push_back(core);
                for (size_t cas : core_and_shells)
                {
                    // Loop over the exclusions for this itype and core or its shells
                    for (size_t jj = 0; jj < exclusions[cas].size(); ++jj)
                    {
                        size_t aj = exclusions[cas][jj];
                        // Now check the pair list
                        auto it   = pairs->begin();
                        while (pairs->end() != it)
                        {
                            size_t aai = (*it)->atomIndex(0);
                            size_t aaj = (*it)->atomIndex(1);
                            // If we find an interaction with the vsite, remove it.
                            if ((aai == i && aaj == aj) || (aaj == i && aai == aj))
                            {
                                it = pairs->erase(it);
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
        auto &pol = entries_.find(InteractionType::POLARIZATION)->second;
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

void Topology::makePropers(const ForceField *pd)
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
            setEntryIdentifiers(pd, InteractionType::PROPER_DIHEDRALS);
        }
    }
}

int Topology::makeVsite2s(const ForceField *pd,
                          AtomList         *atomList)
{
    if (!pd)
    {
        GMX_THROW(gmx::InternalError("Why did you call makeVsites2 without a force field?"));
        return 0;
    }
    auto itype_vs2 = InteractionType::VSITE2;
    if (!pd->interactionPresent(itype_vs2))
    {
        if (debug)
        {
            fprintf(debug, "Force field does not contain vsites.\n");
        }
        return 0;
    }
    auto &ffvs     = pd->findForcesConst(itype_vs2);
    if (ffvs.empty())
    {
        if (debug)
        {
            fprintf(debug, "Force field has zero vsite2 interactions!\n");
        }
        return 0;
    }
    // Loop over all the bonds to make vsites.
    auto itype_b = InteractionType::BONDS;
    if (entries_.find(itype_b) == entries_.end())
    {
        if (debug)
        {
            fprintf(debug, "There are no bonds to generate vsites2 from.\n");
        }
        return 0;
    }
    auto &bonds     = entry(itype_b);
    // Add the virtual sites entry
    TopologyEntryVector vsite2;
    // Force field info for vsites2
    auto ff_vs2 = pd->findForcesConst(itype_vs2);

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

        if (debug)
        {
            fprintf(debug, "Found bond %s %s\n", bai.c_str(), baj.c_str());
        }
        for(auto &fvs : ffvs.parametersConst())
        {
            if (debug)
            {
                fprintf(debug, "Checking vsite %s\n", fvs.first.id().c_str());
            }
            auto vsatoms = fvs.first.atoms();
            auto vsbo    = fvs.first.bondOrders();
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
                                    ptype->gmxParticleType(),
                                    0, ptype->mass(), ptype->charge());
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
                    if (debug)
                    {
                        fprintf(debug, "Adding vs2 %s%d %s%d %d\n",
                                atoms_[ai].element().c_str(), ai,
                                atoms_[aj].element().c_str(), aj, vs2);
                    }
                    // Add bond orders, copied from the bond.
                    for (auto b : border)
                    {
                        vsnew.addBondOrder(b);
                    }
                    // Special bond order for vsites
                    vsnew.addBondOrder(9);
                    vsite2.push_back(std::any_cast<Vsite2>(std::move(vsnew)));
                }
            }
        }
    }
    // If we did find any vsite2 instances, add the whole vector to the topology.
    // A very subtle programming issue arises here:
    // after the std::move operation, the vector is empty
    // and therefore vsite2.size() == 0. Hence we have to store the size in a variable.
    int nvsites2 = vsite2.size();
    if (nvsites2 > 0)
    {
        entries_.insert({ itype_vs2, std::move(vsite2) });
    }
    return nvsites2;
}



int Topology::makeVsite3s(const ForceField *pd,
                          AtomList         *atomList)
{
    if (!pd)
    {
        GMX_THROW(gmx::InternalError("Why did you call makeVsites3 without a force field?"));
        return 0;
    }
    auto itype_vs3 = InteractionType::VSITE3;
    if (!pd->interactionPresent(itype_vs3))
    {
        if (debug)
        {
            fprintf(debug, "Force field does not contain vsites3.\n");
        }
        return 0;
    }
    auto &ffvs     = pd->findForcesConst(itype_vs3);
    if (ffvs.empty())
    {
        if (debug)
        {
            fprintf(debug, "Force field has zero vsite3 interactions!\n");
        }
        return 0;
    }
    auto itype_b = InteractionType::ANGLES;
    if (entries_.find(itype_b) == entries_.end())
    {
        if (debug)
        {
            fprintf(debug, "There are no bonds to generate vsites3 from.\n");
        }
        return 0;
    }
    auto &angles     = entry(itype_b);

    TopologyEntryVector vsite3;

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
        fprintf(stderr, "%s %d %s %d %s %d\n", bai.c_str(), ai, baj.c_str(),aj, bak.c_str(), ak );
        if (debug)
        {
            fprintf(debug, "Found angle %s %s %s\n", bai.c_str(), baj.c_str(), bak.c_str());
        }


        for (auto &fvs : ffvs.parametersConst())
        {
            if (debug)
            {
                fprintf(debug, "Checking vsite %s\n", fvs.first.id().c_str());
            }

            auto vsatoms = fvs.first.atoms();
            auto vsbo    = fvs.first.bondOrders();
            bool found   = false;

            if (border[0] == vsbo[0] && border[1] == vsbo[1] && bai == vsatoms[0] && baj == vsatoms[1] && bak == vsatoms[2])
            {
                    found = true;
            }
                else if (baj == vsatoms[1] && bak == vsatoms[0] && bai == vsatoms[2] && border[1] == vsbo[0] && border[0] == vsbo[1] )

                {
                    found   = true;
                    int tmp = ak; ak = ai; ai = tmp;
                  //  int tmp=aj; aj=ai; ai=tmp;
                }


            if (found)
            {
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
                    for (int pid=0; pid<1; pid++)

             	    {
                          ActAtom newatom(ptype->id().id(), vstype, ptype->id().id(),
                                    ptype->gmxParticleType(),
                                    0, ptype->mass(), ptype->charge());

                    int vs3 = atomList->size();
                    newatom.addCore(ai);
                    newatom.addCore(aj);
                    newatom.addCore(ak);
                    newatom.setResidueNumber(atoms_[ai].residueNumber());


                    gmx::RVec vzero = {0, 0, 0};
                    size_t after = std::max({ai, aj, ak});
                    auto iter = std::find(atomList->begin(), atomList->end(), after);
                    atomList->insert(std::next(iter), ActAtomListItem(newatom, vs3, vzero));


                    // Create new top
                     Vsite3 vsnew(ai, aj, ak, vs3);
                    //  Vsite3 vsnew(ai, aj, ak, vsIds);
                    if (debug)
                    {
                        fprintf(debug, "Adding vs3 %s%d %s%d %s%d %d\n",
                                atoms_[ai].element().c_str(), ai,
                                atoms_[aj].element().c_str(), aj,
                                atoms_[ak].element().c_str(), ak, vs3);
                    }

                    // Add bond orders, cp from the bond.
                    for (auto b : border)
                    {
                        vsnew.addBondOrder(b);
                    }

                    // Special bond order for vsites
                    vsnew.addBondOrder(9);

                    vsite3.push_back(std::any_cast<Vsite3>(std::move(vsnew)));
                  }
                }
            }
        }
    }

    // If we found any vsite3 instances, add the vec to the top
    int nvsites3 = vsite3.size();
    if (nvsites3 > 0)
    {
        entries_.insert({itype_vs3, std::move(vsite3)});
    }

    return nvsites3;
}



int Topology::makeVsite3OUTs(const ForceField *pd,
                          AtomList         *atomList)
{
    if (!pd)
    {
        GMX_THROW(gmx::InternalError("Why did you call makeVsites3out without a force field?"));
    }
    auto itype_vs3out = InteractionType::VSITE3OUT;
    if (!pd->interactionPresent(itype_vs3out))
    {
        if (debug)
        {
            fprintf(debug, "Force field does not contain vsites3out.\n");
        }
        return 0;
    }
    auto &ffvs     = pd->findForcesConst(itype_vs3out);
    if (ffvs.empty())
    {
        if (debug)
        {
            fprintf(debug, "Force field has zero vsite3out interactions!\n");
        }
        return 0;
    }

    auto itype_b = InteractionType::ANGLES;
    if (entries_.find(itype_b) == entries_.end())
    {
        if (debug)
        {
            fprintf(debug, "There are no angles to generate vsites3out from.\n");
        }
        return 0;
    }
    auto &angles     = entry(itype_b);

    TopologyEntryVector vsite3out;

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

        if (debug)
        {
            fprintf(debug, "Found angle %s %s %s\n", bai.c_str(), baj.c_str(), bak.c_str());
        }

        for (auto &fvs : ffvs.parametersConst())
        {
            if (debug)
            {
                fprintf(debug, "Checking vsite %s\n", fvs.first.id().c_str());
            }

            auto vsatoms = fvs.first.atoms();
            auto vsbo    = fvs.first.bondOrders();
            bool found   = false;

            if (border[0] == vsbo[0] && border[1] == vsbo[1] && bai == vsatoms[0] && baj == vsatoms[1] && bak == vsatoms[2])
            {
                    found = true;
            }
                else if (baj == vsatoms[1] && bak == vsatoms[0] && bai == vsatoms[2] && border[1] == vsbo[0] && border[0] == vsbo[1] )

                {
                    found   = true;
                    int tmp = ak; ak = ai; ai = tmp;
                }


            if (found)
            {
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
                    for (int pid=0; pid<2; pid++)

             	    {
                          ActAtom newatom(ptype->id().id(), vstype, ptype->id().id(),
                                          ptype->gmxParticleType(),
                                          0, ptype->mass(), ptype->charge());

                          int vs3out = atomList->size();
                          newatom.addCore(ai);
                          newatom.addCore(aj);
                          newatom.addCore(ak);

                          newatom.setResidueNumber(atoms_[ai].residueNumber());

                          gmx::RVec vzero = {0, 0, 0};
                          size_t after = std::max({ai, aj, ak});
                          auto iter = std::find(atomList->begin(), atomList->end(), after);
                          atomList->insert(std::next(iter), ActAtomListItem(newatom, vs3out, vzero));

                          // Create new topology entry
                          Vsite3OUT vsnew(ai, aj, ak, vs3out);
                          if (debug)
                          {
                              fprintf(debug, "Adding vs3out %s%d %s%d %s%d %d\n",
                                      atoms_[ai].element().c_str(), ai,
                                      atoms_[aj].element().c_str(), aj,
                                      atoms_[ak].element().c_str(), ak, vs3out);
                          }
                          
                          // Add bond orders, cp from the bond.
                          for (auto b : border)
                          {
                              vsnew.addBondOrder(b);
                          }

                          // Special bond order for vsites
                          vsnew.addBondOrder(9);
                          vsite3out.push_back(std::any_cast<Vsite3OUT>(std::move(vsnew)));
                    }
                }
            }
        }
    }

    // If we found any vsite3 instances, add the vec to the top
    int nvsites3out = vsite3out.size();
    if (nvsites3out > 0)
    {
        entries_.insert({itype_vs3out, std::move(vsite3out)});
    }

    return nvsites3out;
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

void Topology::dumpPairlist(FILE *fp, InteractionType itype) const
{
    if (nullptr == fp)
    {
        return;
    }
    if (!hasEntry(itype))
    {
        return;
    }
    for(const auto &pl : entry(itype))
    {
        fprintf(fp, "PAIRLIST %s %d %d\n", interactionTypeToString(itype).c_str(),
                pl->atomIndex(0), pl->atomIndex(1));
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
    if (!ci)
    {
        return immStatus::Topology;
    }
    if (ci)
    {
        if (ci->NAtom() == 0)
        {
            return immStatus::NoAtoms;
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

                ActAtom newatom(cai.getName(), atype->element(), cai.getObtype(), eptAtom,
                                atype->atomnumber(), atype->mass(), atype->charge());
                newatom.setResidueNumber(nres-1);
                addAtom(newatom);
                realAtoms_.push_back(natom++);
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
        fprintf(debug, "Tried to convert '%s' to ACT. LOT is '%s/%s'. Natoms is %zu\n",
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
    // We use a list structure to keep track of the atoms and their coordinates.
    // This is needed since we may insert both shells and virtual sites. Although
    // one can insert items in std::vector structures it may lead to poor (N^2)
    // scaling to insert for instance one shell per particle. At the end of this
    // algorithm we copy the contents of the list back to std::vectors.
    AtomList atomList;
    for(size_t index = 0; index < atoms_.size(); index++)
    {
        atomList.push_back(ActAtomListItem(atoms_[index], index, (*x)[index]));
    }

    // Before we can do anything we need to "identify" the bonds
    setEntryIdentifiers(pd, InteractionType::BONDS);

    // Check whether we have virtual sites in the force field.
    int nvsites2 = makeVsite2s(pd, &atomList);

    // Before we can make three-particle vsites, we need to create
    // angles, but only temporarily.
    makeAngles(pd, *x, LinearAngleMin);
    int nvsites3    = makeVsite3s(pd, &atomList);
    int nvsites3out = makeVsite3OUTs(pd, &atomList);
    if (debug)
    {
        fprintf(debug, "Added %d vsite2 %d vsite3 %d vsite3out\n", nvsites2, nvsites3, nvsites3out);
    }
    // Now throw away the angles again since we need to reorder the
    // lists of atoms and vsites.
    auto itype_a = InteractionType::ANGLES;
    if (entries_.end() != entries_.find(itype_a))
    {
        entries_.erase(itype_a);
    }

    // Now time for shells.
    addShells(pd, &atomList);
    if (debug)
    {
        fprintf(debug, "Added %zu shells\n", atomList.size()-nvsites2-nvsites3-nvsites3out-atoms_.size());
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
    setEntryIdentifiers(pd, InteractionType::POLARIZATION);
    if (nvsites2 > 0)
    {
        setEntryIdentifiers(pd, InteractionType::VSITE2);
    }

    if (nvsites3 > 0)
    {
        setEntryIdentifiers(pd, InteractionType::VSITE3);
    }

    if (nvsites3out > 0)
    {
        setEntryIdentifiers(pd, InteractionType::VSITE3OUT);
    }

    // Now make angles etc.
    makeAngles(pd, *x, LinearAngleMin);
    makeImpropers(pd, *x, PlanarAngleMax);
    // Check whether we have dihedrals in the force field.
    if (pd->interactionPresent(InteractionType::PROPER_DIHEDRALS))
    {
        // Store temporary address of variable for performance
        auto &fs = pd->findForcesConst(InteractionType::PROPER_DIHEDRALS);
        if (!fs.empty() || missingParameters::Generate == missing)
        {
            makePropers(pd);
        }
    }
    makePairs(pd, InteractionType::VDW);
    makePairs(pd, InteractionType::COULOMB);

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

std::vector<std::vector<int>> Topology::generateExclusions(TopologyEntryVector *pairs,
                                                           int                  nrexcl)
{
    std::vector<std::vector<int>> exclusions;
    exclusions.resize(atoms_.size());
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
                        exclusions[a[0]].push_back(a[1]);
                        exclusions[a[1]].push_back(a[0]);
                    }
                }
            }
            break;
        case InteractionType::POLARIZATION:
            {
                for(auto &b : myEntry.second)
                {
                    auto a = b->atomIndices();
                    exclusions[a[0]].push_back(a[1]);
                    exclusions[a[1]].push_back(a[0]);
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
                        exclusions[a[0]].push_back(a[2]);
                        exclusions[a[2]].push_back(a[0]);
                    }
                }
            }
            break;
        default:
            break;
        }
    }
    for(size_t ai = 0; ai < exclusions.size(); ai++)
    {
        for(size_t j = 0; j < exclusions[ai].size(); j++)
        {
            size_t aj = exclusions[ai][j];
            for(auto pp = pairs->begin(); pp < pairs->end(); ++pp)
            {
                size_t a0 = (*pp)->atomIndex(0);
                size_t a1 = (*pp)->atomIndex(1);
                if ((ai == a0 && aj == a1) || (ai == a1 && aj == a0))
                {
                    pairs->erase(pp);
                    break;
                }
            }
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

static void fillParams(const ForceFieldParameterList &fs,
                       const Identifier              &btype,
                       int                            nr,
                       const char                    *param_names[],
                       std::vector<double>           *param)
{
    if (param->empty())
    {
        param->resize(nr, 0);
    }
    auto ff = fs.findParameterMapConst(btype);
    if (!ff.empty())
    {
        for (int i = 0; i < nr; i++)
        {
            auto fp      = ff.find(param_names[i]);
            double value = 0;
            if (ff.end() != fp)
            {
                value = fp->second.internalValue();
            }
            (*param)[i] = value;
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
                fillParams(fs, topID, lj12_6NR, lj12_6_name, &param);
                break;
            case F_LJ8_6:
                fillParams(fs, topID, lj8_6NR, lj8_6_name, &param);
                break;
            case F_LJ14_7:
                fillParams(fs, topID, lj14_7NR, lj14_7_name, &param);
                break;
            case F_WBHAM:
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
            case F_VSITE2:
                fillParams(fs, topID, vsite2NR, vsite2_name, &param);
                break;
            case F_VSITE3:
                fillParams(fs, topID, vsite3NR, vsite3_name, &param);
                break;
            case F_VSITE3OUT:
                fillParams(fs, topID, vsite3outNR, vsite3out_name, &param);
                break;
            case F_VSITE3FAD:
                fillParams(fs, topID, vsite3fadNR, vsite3fad_name, &param);
                break;
            default:
                GMX_THROW(gmx::InternalError(gmx::formatString("Missing case %s when filling the topology structure.", interaction_function[fs.gromacsType()].name).c_str()));
            }
            topentry->setParams(param);
        }
    }
}

void Topology::setEntryIdentifiers(const ForceField *pd,
                                   InteractionType   itype)
{
    if (!hasEntry(itype))
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
                {
                    btype.push_back(atoms_[jj].ffType());
                    break;
                }
            case InteractionType::POLARIZATION:
                {
                    // For COULOMB there are two particles,
                    // but for polarization just one.
                    if (atype->hasInteractionType(itype))
                    {
                        btype.push_back(atype->interactionTypeToIdentifier(itype).id());
                    }
                    break;
                }
            case InteractionType::COULOMB:
                {
                    // For COULOMB there are two particles,
                    // but for polarization just one.
                    if (atype->hasInteractionType(itype))
                    {
                        btype.push_back(atype->interactionTypeToIdentifier(itype).id());
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
                (*topentry)->setId(Identifier(btype[0]));
            }
            else
            {
                (*topentry)->setId({ Identifier(btype, (*topentry)->bondOrders(),
                                                fs.canSwap()) });
            }
        }
        else if (debug)
        {
            fprintf(debug, "Could not find identifier for %s for atomtype '%s'",
                    interactionTypeToString(itype).c_str(), atypes.c_str());
        }
    }
}

void Topology::setIdentifiers(const ForceField *pd)
{
    for(auto &entry : entries_)
    {
        // TODO this if statement should not be needed.
        if (InteractionType::VSITE2 != entry.first && InteractionType::VSITE3 != entry.first && InteractionType::VSITE3OUT != entry.first)
        {
            setEntryIdentifiers(pd, entry.first);
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

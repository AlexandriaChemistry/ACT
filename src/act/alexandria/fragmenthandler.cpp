/*
 * This source file is part of the Alexandria Chemistry Toolkit.
 *
 * Copyright (C) 2014-2022
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
#include "fragmenthandler.h"

#include <vector>

#include "act/qgen/qgen_acm.h"
#include "act/molprop/fragment.h"
#include "act/utility/stringutil.h"
#include "gromacs/topology/atoms.h"

namespace alexandria
{    

FragmentHandler::FragmentHandler(const ForceField               *pd,
                                 const std::vector<gmx::RVec>   &coordinates,
                                 const std::vector<std::string> &residueNames,
                                 const std::vector<ActAtom>     &atoms,
                                 const std::vector<Bond>        &bonds,
                                 const std::vector<Fragment>    *fragments,
                                 const std::vector<int>         &shellRenumber,
                                 missingParameters               missing)

{
    GMX_RELEASE_ASSERT(fragments != nullptr,
                       "Empty fragments passed. Wazzuppwitdat?");
    GMX_RELEASE_ASSERT(fragments->size() > 0, "No fragments. Huh?");
    if (fragments->size() != residueNames.size())
    {
        GMX_THROW(gmx::InternalError(gmx::formatString("Number of fragments (%zu) should match number of residue names (%zu)", fragments->size(), residueNames.size()).c_str()));
    }
    topologies_.resize(fragments->size());

    bonds_.resize(fragments->size());
    natoms_           = 0;
    size_t  ff        = 0;
    for(auto f = fragments->begin(); f < fragments->end(); ++f)
    {
        int minatom = atoms.size();
        for(auto &i : f->atoms())
        {
            minatom = std::min(minatom, i);
        }
        atomStart_.push_back(minatom);
        qtotal_.push_back(f->charge());
    }
    std::vector<bool> atomFound(atoms.size(), false);
    for(auto f = fragments->begin(); f < fragments->end(); ++f)
    {
        if (f->atoms().size() == 0)
        {
            GMX_THROW(gmx::InternalError(gmx::formatString("No atoms in fragment %zu with formula %s", ff, f->formula().c_str()).c_str()));
        }
        // ID
        ids_.push_back(f->id());
        // If polarizable we will need to add these
        std::vector<TopologyEntry *> pols;
        int offset     = atomStart_[ff];
        int pol_offset = offset;
        if (!shellRenumber.empty())
        {
            pol_offset = shellRenumber[offset];
            atomStart_[ff] = pol_offset;
        }
        // Count the number of atoms
        std::vector<int> toAdd;
        for(auto &a : f->atoms())
        {
            size_t anew = a;
            // Check whether this atom is present already
            if (anew >= atoms.size()) 
            {
                GMX_THROW(gmx::InvalidInputError(gmx::formatString("Atom number %zu in fragment %zu too large (max %zu)",
                                                                   anew, fragments->size(), atoms.size()).c_str()));
            }
            else if (atomFound[anew])
            {
                GMX_THROW(gmx::InvalidInputError(gmx::formatString("Atom %zu occurs more than once in fragment description", anew).c_str()));
            }
            else
            {
                atomFound[anew] = true;
            }
            if (!shellRenumber.empty())
            {
                GMX_RELEASE_ASSERT(anew < shellRenumber.size(), "Atom number out of range");
                anew   = shellRenumber[anew];
            }
            
            // We add the new atom index, but relative to the first atom
            // in the compound.
            toAdd.push_back(anew-pol_offset);
            if (pd->polarizable())
            {
                auto fa = pd->findParticleType(atoms[anew].ffType());
                if (fa->hasInteractionType(InteractionType::POLARIZATION))
                {
                    for(const auto &ss : atoms[anew].shells())
                    {
                        toAdd.push_back(ss-pol_offset);
                        auto pp = new TopologyEntry();
                        pp->addAtom(anew-pol_offset);
                        pp->addAtom(ss-pol_offset);
                        pp->addBondOrder(1.0);
                        pols.push_back(pp);
                    }
                }
            }
        }
        if (pd->polarizable())
        {
            topologies_[ff].addEntry(InteractionType::POLARIZATION, pols);
        }

        int j = 0;
        for(auto &i : toAdd)
        {
            auto fa      = pd->findParticleType(atoms[i+pol_offset].ffType());
            int  anumber = my_atoi(fa->optionValue("atomnumber").c_str(), "atomic number");
            ActAtom newat(atoms[i+pol_offset].name(),
                          fa->optionValue("element"),
                          atoms[i+pol_offset].ffType(),
                          atoms[i+pol_offset].pType(),
                          anumber,
                          atoms[i+pol_offset].mass(),
                          atoms[i+pol_offset].charge());
            int resnr = std::max(0, atoms[i+pol_offset].residueNumber()-1);
            newat.setResidueNumber(resnr);
            topologies_[ff].addResidue(resnr, residueNames[resnr+1]);
            topologies_[ff].addAtom(newat);
            j++;
        }
        // Make sure the residue names and numbers are consistent
        std::vector<Fragment> copyFragment;
        copyFragment.push_back((*fragments)[ff]);
        topologies_[ff].shellsToAtoms();
        QgenAcm_.push_back(QgenAcm(pd, topologies_[ff].atoms(), f->charge()));
        for(const auto &b : bonds)
        {
            int ai = b.aI();
            int aj = b.aJ();
            if (std::find(f->atoms().begin(), f->atoms().end(),
                          ai) != f->atoms().end() &&
                std::find(f->atoms().begin(), f->atoms().end(),
                          aj) != f->atoms().end())
            {
                // Bonds should be numbered from the start of the atom.
                // Now the shells should not be taken into account, since
                // the ACM code will do it.
                Bond bb(ai - offset, aj - offset, b.bondOrder());
                bonds_[ff].push_back(bb);
                if (shellRenumber.empty())
                {
                    topologies_[ff].addBond(bb);
                }
                else
                {
                    // For the internal topologies however, the atoms are renumbered
                    // already if there are shell, and the bonds should be too.
                    ai = shellRenumber[ai] - pol_offset;
                    aj = shellRenumber[aj] - pol_offset;
                
                    Bond bb2(ai, aj, b.bondOrder());
                    topologies_[ff].addBond(bb2);
                }
            }
        }
        natoms_ += toAdd.size();
        ff      += 1;
    }
    for(size_t ff = 0; ff < topologies_.size(); ff++)
    {
        int natom = topologies_[ff].atoms().size();
        std::vector<gmx::RVec> myx(natom);
        int j = 0;
        for (size_t i = atomStart_[ff]; i < atomStart_[ff]+natom; i++)
        {
            copy_rvec(coordinates[i], myx[j++]);
        }

        topologies_[ff].build(pd, myx, 175.0, 5.0, missingParameters::Error);
        topologies_[ff].setIdentifiers(pd);
        if (missing != missingParameters::Generate)
        {
            topologies_[ff].fillParameters(pd);
        }
    }
    if (debug)
    {
        fprintf(debug, "FragmentHandler: atoms.size() %lu natoms %zu nbonds %zu nfragments %zu\n",
                atoms.size(), natoms_, bonds.size(), topologies_.size());
    }
}

void FragmentHandler::fetchCharges(std::vector<double> *qq)
{
    qq->resize(natoms_, 0);
    for(size_t ff = 0; ff < topologies_.size(); ++ff)
    {
        for (size_t a = 0; a < topologies_[ff].atoms().size(); a++)
        {
            // TODO: Check whether this works for polarizable models
            (*qq)[atomStart_[ff] + a] = QgenAcm_[ff].getQ(a);
        }
    }
}

eQgen FragmentHandler::generateCharges(FILE                         *fp,
                                       const std::string            &molname,
                                       const std::vector<gmx::RVec> &x,
                                       const ForceField             *pd,
                                       std::vector<ActAtom>         *atoms)
{
    auto   eqgen = eQgen::OK;
    for(size_t ff = 0; ff < topologies_.size(); ++ff)
    {
        // TODO only copy the coordinates if there is more than one fragment.
        std::vector<gmx::RVec> xx;
        xx.resize(topologies_[ff].atoms().size());
        for(size_t a = 0; a < topologies_[ff].atoms().size(); a++)
        {
            copy_rvec(x[atomStart_[ff]+a], xx[a]);
        }
        QgenAcm_[ff].setQtotal(qtotal_[ff]);
        eqgen = QgenAcm_[ff].generateCharges(fp, molname, pd, 
                                             topologies_[ff].atomsPtr(),
                                             xx, bonds_[ff]);
        if (eQgen::OK != eqgen)
        {
            break;
        }
        for(size_t a = 0; a < topologies_[ff].atoms().size(); a++)
        {
            (*atoms)[atomStart_[ff]+a].setCharge(topologies_[ff].atoms()[a].charge());
        }
    }
    return eqgen;
}

} // namespace alexandria

/*
 * This source file is part of the Alexandria Chemistry Toolkit.
 *
 * Copyright (C) 2022-2025
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

#include "act/alexandria/topology.h"
#include "act/alexandria/symmetrize_charges.h"
#include "act/qgen/qgen_acm.h"
#include "act/molprop/fragment.h"
#include "act/utility/stringutil.h"

namespace alexandria
{    

FragmentHandler::FragmentHandler(ForceField                   *pd,
                                 const std::vector<gmx::RVec> &coordinates,
                                 const std::vector<ActAtom>   &atoms,
                                 const std::vector<Bond>      &bonds,
                                 const std::vector<Fragment>  *fragments,
                                 missingParameters             missing)

{
    GMX_RELEASE_ASSERT(fragments != nullptr,
                       "Empty fragments passed. Wazzuppwitdat?");
    GMX_RELEASE_ASSERT(fragments->size() > 0, "No fragments. Huh?");
    if (atoms.size() != coordinates.size())
    {
        fprintf(stderr, "Received %zu atoms and %zu coordinates in fragmenthandler. Giving up.\n",
                atoms.size(), coordinates.size());
        return;
    }
    bonds_.resize(fragments->size());
    std::vector<bool> atomFound(coordinates.size(), false);
    // Total number of atoms
    natoms_ = 0;
    bool allWell = true;
    size_t  ff = 0;
    // Copy the coordinates
    std::vector<gmx::RVec> x = coordinates;
    for(auto f = fragments->begin(); allWell && f < fragments->end(); ++f)
    {
        // First reality check
        if (f->atoms().size() == 0)
        {
            fprintf(stderr, "No atoms in fragment %zu with formula %s", ff, f->formula().c_str());
            allWell = false;
            break;
        }
        // Determine start of each molecule
        int offset = x.size()+1;
        for(auto &i : f->atoms())
        {
            offset = std::min(offset, i);
            if (atomFound[i])
            {
                fprintf(stderr, "Atom %d (%s) occurs multiple times, most recently in fragment %s\n", i, atoms[i].name().c_str(), f->formula().c_str());
                allWell = false;
                break;
            }
            atomFound[i] = true;
        }
        if (!allWell)
        {
            break;
        }
        // Split bonds
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
            }
        }
        // Create new topology
        auto top = new Topology(bonds_[ff]);

        // Split coordinate array
        size_t natom = f->atoms().size();
        std::vector<gmx::RVec> xfrag(natom);
        // Copy the atoms from the global topology and make new coordinate array
        int j = 0;
        for(size_t i = offset; i < offset+natom; i++)
        {
            top->addAtom(atoms[i]);
            copy_rvec(x[i], xfrag[j++]);
        }
        // Now build the rest of the topology
        bool allParametersFound = top->build(pd, &xfrag, 175.0, 5.0, missing);
        allWell = allWell && allParametersFound;
        if (allWell)
        {
            // Array of total charges
            qtotal_.push_back(f->charge());
            // ID copied from Fragment ID.
            ids_.push_back(f->inchi());
            // Structure for charge generation
            QgenAcm_.push_back(QgenAcm(pd, top->atoms(), bonds_[ff], f->charge()));
            // Total number of atoms
            natoms_ += top->atoms().size();
            // Extend topologies_ array
            topologies_.push_back(std::move(top));
        }

        // Increase counter
        ff += 1;
    }
    if (allWell)
    {
        // Finaly determine molecule boundaries
        atomStart_.resize(fragments->size(), 0);
        for(size_t ff = 0; ff < topologies_.size(); ff++)
        {
            if (0 == ff)
            {
                atomStart_[ff] = 0;
            }
            else
            {
                atomStart_[ff] = atomStart_[ff-1] + topologies_[ff-1]->atoms().size();
            }
        }
    }
    else
    {
        topologies_.clear();
        ids_.clear();
        qtotal_.clear();
    }
}

void FragmentHandler::fetchCharges(std::vector<double> *qq)
{
    qq->resize(natoms_, 0);
    for(size_t ff = 0; ff < topologies_.size(); ++ff)
    {
        auto myatoms = topologies_[ff]->atoms();
        for (size_t a = 0; a < myatoms.size(); a++)
        {
            // TODO: Check whether this works for polarizable models
            size_t index = atomStart_[ff] + a;
            GMX_RELEASE_ASSERT(index < natoms_, 
                               gmx::formatString("Index %ld out of range %ld", index, natoms_).c_str());
            (*qq)[index] = myatoms[a].charge();
            if (debug)
            {
                fprintf(debug, "Charge %ld = %g\n", index,
                        (*qq)[index]);
            }
        }
    }
}

eQgen FragmentHandler::generateCharges(FILE                         *fp,
                                       const std::string            &molname,
                                       const std::vector<gmx::RVec> &x,
                                       const ForceField             *pd,
                                       std::vector<ActAtom>         *atoms,
                                       const std::vector<int>       &symmetric_charges)
{
    auto   eqgen = eQgen::OK;
    switch (algorithm_)
    {
    case ChargeGenerationAlgorithm::EEM:
    case ChargeGenerationAlgorithm::SQE:
        {
            for(size_t ff = 0; ff < topologies_.size(); ++ff)
            {
                // TODO only copy the coordinates if there is more than one fragment.
                std::vector<gmx::RVec> xx;
                xx.resize(topologies_[ff]->atoms().size());
                for(size_t a = 0; a < topologies_[ff]->atoms().size(); a++)
                {
                    copy_rvec(x[atomStart_[ff]+a], xx[a]);
                }
                QgenAcm_[ff].setQtotal(qtotal_[ff]);
                eqgen = QgenAcm_[ff].generateCharges(fp, molname, pd, 
                                                     topologies_[ff]->atomsPtr(),
                                                     xx, bonds_[ff]);
                if (eQgen::OK != eqgen)
                {
                    fprintf(stderr, "Failed to generate charges for %s: '%s'\n",
                            molname.c_str(), QgenAcm_[ff].status());
                    break;
                }
                // Fetch charges to one vector, then symmetrize them
                std::vector<double> qnew;
                fetchCharges(&qnew);
                apply_symmetrized_charges(&qnew, symmetric_charges);
                setCharges(qnew);
                for(size_t a = 0; a < topologies_[ff]->atoms().size(); a++)
                {
                    (*atoms)[atomStart_[ff]+a].setCharge(topologies_[ff]->atoms()[a].charge());
                }
            }
        }
        break;
    case ChargeGenerationAlgorithm::Read:
        break;
    default: // throws
        GMX_THROW(gmx::InvalidInputError(gmx::formatString("No support for %s algorithm for fragments", 
                                                           chargeGenerationAlgorithmName(algorithm_).c_str()).c_str()));
    }
    return eqgen;
}

bool FragmentHandler::fetchCharges(std::vector<ActAtom> *atoms)
{
    if (atoms->size() != natoms_)
    {
        return false;
    }
    size_t k = 0;
    for(size_t i = 0; i < topologies_.size(); i++)
    {
        auto ats = topologies_[i]->atoms();
        for(size_t j = 0; j < ats.size(); j++)
        {
            (*atoms)[k++].setCharge(ats[j].charge());
        }
    }
    return true;
}

bool FragmentHandler::setCharges(const chargeMap &qmap)
{
    bool success = true;
    for(size_t i = 0; i < ids_.size() && success; i++)
    {
        auto qptr = qmap.find(ids_[i]);
        if (qmap.end() != qptr)
        {
            auto aptr = topologies_[i]->atomsPtr();
            if (aptr->size() == qptr->second.size())
            {
                for(size_t a = 0; success && a < aptr->size(); a++)
                {
                    if (!((*aptr)[a].id() == qptr->second[a].first))
                    {
                        fprintf(stderr, "Atom mismatch when reading charges from chargemap for %s. Expected %s but found %s. Make sure the atoms in your chargemap match those in your molprop file.\n",
                                ids_[i].c_str(), (*aptr)[a].id().id().c_str(),
                                qptr->second[a].first.id().c_str());
                        success = false;
                    }
                    (*aptr)[a].setCharge(qptr->second[a].second);
                    if (debug)
                    {
                        fprintf(debug, "qmap Charge %zu = %g\n", a,
                                (*aptr)[a].charge());
                    }
                }
            }
            else
            {
                fprintf(stderr, "Compound '%s' (%zu particles) from the input does not match the same compound in the charge map (%zu particles).\n",
                        ids_[i].c_str(), aptr->size(), qptr->second.size());
                success = false;
            }
        }
        else
        {
            fprintf(stderr, "Cannot find compound '%s' in the charge map\n",
                    ids_[i].c_str());
            success = false;
        }
    }
    fixedQ_ = success;
    // Set algorithm to reading charges for the future
    if (debug)
    {
        fprintf(debug, "Copied charges from chargemap to fragments\n");
    }
    if (success)
    {
        algorithm_ = ChargeGenerationAlgorithm::Read;
    }
    return success;
}

void FragmentHandler::setCharges(const std::vector<ActAtom> &atoms)
{
    for(size_t ff = 0; ff < topologies_.size(); ++ff)
    {
        auto aptr = topologies_[ff]->atomsPtr();
        for(size_t a = 0; a < aptr->size(); a++)
        {
            (*aptr)[a].setCharge(atoms[atomStart_[ff]+a].charge());
        }
    }
}

void FragmentHandler::setCharges(const std::vector<double> &q)
{
    int j = 0;
    for(size_t ff = 0; ff < topologies_.size(); ++ff)
    {
        auto aptr = topologies_[ff]->atomsPtr();
        for(size_t a = 0; a < aptr->size(); a++)
        {
            (*aptr)[a].setCharge(q[j]);
            j += 1;
        }
    }
}

} // namespace alexandria

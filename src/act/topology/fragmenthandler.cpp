/*
 * This source file is part of the Alexandria Chemistry Toolkit.
 *
 * Copyright (C) 2022-2026
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

#include "act/molprop/fragment.h"
#include "act/qgen/qgen_acm.h"
#include "act/topology/symmetrize_charges.h"
#include "act/topology/topology.h"
#include "act/utility/stringutil.h"
#include "gromacs/math/vec.h"

namespace alexandria
{    

void FragmentHandler::init(MsgHandler                   *msghandler,
                           ForceField                   *pd,
                           const std::vector<gmx::RVec> &coordinates,
                           const std::vector<ActAtom>   &atoms,
                           const std::vector<Bond>      &bonds,
                           const std::vector<Fragment>  *fragments,
                           missingParameters             missing)
{
    GMX_RELEASE_ASSERT(fragments != nullptr,
                       "Empty fragments passed. Wazzuppwitdat?");
    GMX_RELEASE_ASSERT(!fragments->empty(), "No fragments. Huh?");
    if (atoms.size() != coordinates.size())
    {
        msghandler->msg(ACTStatus::Error,
                        gmx::formatString("Received %zu atoms and %zu coordinates in fragmenthandler. Giving up.",
                                          atoms.size(), coordinates.size()));
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
        if (f->atoms().empty())
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
                const Bond bb(ai - offset, aj - offset, b.bondOrder());
                bonds_[ff].push_back(bb);
            }
        }
        // Create new topology
        Topology top;
        top.init(bonds_[ff]);

        // Split coordinate array
        size_t natom = f->atoms().size();
        std::vector<gmx::RVec> xfrag(natom);
        // Copy the atoms from the global topology and make new coordinate array
        int j = 0;
        for(size_t i = offset; i < offset+natom; i++)
        {
            top.addAtom(atoms[i]);
            copy_rvec(x[i], xfrag[j++]);
        }
        // Now build the rest of the topology
        top.setLeaveOutUnchargedPairs(leaveOutUncharged_);
        top.build(msghandler, pd, &xfrag, 175.0, 5.0, missing);
        allWell = allWell && msghandler->ok();
        if (allWell)
        {
            // Array of total charges
            qtotal_.push_back(f->charge());
            // ID copied from Fragment ID.
            ids_.push_back(f->inchi());
            // Structure for charge generation
            QgenAcm_.emplace_back(msghandler, pd, top.atoms(), bonds_[ff], f->charge());
            // Total number of atoms
            natoms_ += top.atoms().size();
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
        polinfo_.resize(fragments->size());
        for(size_t ff = 0; ff < topologies_.size(); ff++)
        {
            if (0 == ff)
            {
                atomStart_[ff] = 0;
            }
            else
            {
                atomStart_[ff] = atomStart_[ff-1] + topologies_[ff-1].atoms().size();
            }
            polinfo_[ff].init(pd, &topologies_[ff]);
        }
    }
    else
    {
        topologies_.clear();
        ids_.clear();
        qtotal_.clear();
    }
}

void FragmentHandler::fetchCharges(MsgHandler *msg_handler, std::vector<double> *qq)
{
    qq->resize(natoms_, 0);
    for(size_t ff = 0; ff < topologies_.size(); ++ff)
    {
        const auto &myatoms = topologies_[ff].atoms();
        for (size_t a = 0; a < myatoms.size(); a++)
        {
            //! \todo: Check whether this works for polarizable models
            size_t index = atomStart_[ff] + a;
            GMX_RELEASE_ASSERT(index < natoms_, 
                               gmx::formatString("Index %ld out of range %ld", index, natoms_).c_str());
            (*qq)[index] = myatoms[a].charge();
            if (msg_handler->debug())
            {
                msg_handler->writeDebug(gmx::formatString("Charge %ld = %g\n", index,
                                                          (*qq)[index]));
            }
        }
    }
}

bool FragmentHandler::sqeOnGpu(const ForceField *pd) const
{
    if (fixedQ_ || algorithm_ != ChargeGenerationAlgorithm::SQE)
    {
        return fixedQ_;
    }
    for (size_t ff = 0; ff < QgenAcm_.size(); ++ff)
    {
        if (!QgenAcm_[ff].sqeOnGpu(pd, bonds_[ff]))
        {
            return false;
        }
    }
    return true;
}

const std::vector<int> &FragmentHandler::acmParticles(size_t ff) const
{
    static const std::vector<int> none;
    return fixedQ_ ? none : QgenAcm_[ff].acmParticles();
}

const std::vector<const ParticleType *> &FragmentHandler::particleTypes(const ForceField *pd, size_t ff)
{
    if (ptypesVersion_ != pd->particleTypeVersion() || ptypes_.size() != topologies_.size())
    {
        ptypes_.assign(topologies_.size(), {});
        for (size_t f = 0; f < topologies_.size(); f++)
        {
            for (const auto &a : topologies_[f].atoms())
            {
                ptypes_[f].push_back(pd->hasParticleType(a.ffType()) ? pd->findParticleType(a.ffType()) : nullptr);
            }
        }
        ptypesVersion_ = pd->particleTypeVersion();
    }
    return ptypes_[ff];
}

void FragmentHandler::updateFixedCharges(const ForceField *pd, std::vector<ActAtom> *atoms)
{
    for (size_t ff = 0; ff < topologies_.size(); ++ff)
    {
        auto ta = topologies_[ff].atomsPtr();
        if (!fixedQ_)
        {
            QgenAcm_[ff].updateFixedCharges(ta);
        }
        const auto &ptypes = particleTypes(pd, ff);
        for (size_t a = 0; a < ta->size(); a++)
        {
            auto &fa = (*ta)[a];
            auto &ca = (*atoms)[atomStart_[ff]+a];
            ca.setCharge(fa.charge());
            ca.setChargeVS(fa.chargeVS(), fa.zetaVS(), fa.rowVS());
            if (ptypes[a])
            {
                double zeta = ptypes[a]->zeta();
                fa.setZeta(zeta);
                ca.setZeta(zeta);
            }
        }
    }
}

eQgen FragmentHandler::generateCharges(MsgHandler                   *msg_handler,
                                       const std::string            &molname,
                                       const std::vector<gmx::RVec> &x,
                                       ForceField                   *pd,
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
                //! \todo only copy the coordinates if there is more than one fragment.
                std::vector<gmx::RVec> xx;
                xx.resize(topologies_[ff].atoms().size());
                for(size_t a = 0; a < topologies_[ff].atoms().size(); a++)
                {
                    copy_rvec(x[atomStart_[ff]+a], xx[a]);
                }
                QgenAcm_[ff].setQtotal(qtotal_[ff]);
                eqgen = QgenAcm_[ff].generateCharges(msg_handler, molname, pd, 
                                                     topologies_[ff].atomsPtr(),
                                                     xx, bonds_[ff]);
                if (eQgen::OK != eqgen)
                {
                    msg_handler->msg(ACTStatus::Warning,
                                     gmx::formatString("Failed to generate charges for %s: '%s'\n",
                                                       molname.c_str(), QgenAcm_[ff].status()));
                    break;
                }
                // Fetch charges to one vector, then symmetrize them
                std::vector<double> qnew;
                fetchCharges(msg_handler, &qnew);
                apply_symmetrized_charges(&qnew, symmetric_charges);
                setCharges(qnew);
                // Copy back charges for this fragment to combined topology atoms
                for(size_t a = 0; a < topologies_[ff].atoms().size(); a++)
                {
                    const auto &fa = topologies_[ff].atoms()[a];
                    auto       &ca = (*atoms)[atomStart_[ff]+a];
                    ca.setCharge(fa.charge());
                    ca.setChargeVS(fa.chargeVS(), fa.zetaVS(), fa.rowVS());
                }
            }
        }
        break;
    case ChargeGenerationAlgorithm::NONE:
        // Just copy the charges from the force field, in case they are updated.
        // We need to change both the fragments and the total.
        for(size_t ff = 0; ff < topologies_.size(); ++ff)
        {
            auto   ta   = topologies_[ff].atomsPtr();
            double qtot = 0;
            for(size_t a = 0; a < ta->size(); a++)
            {
                if (!pd->hasParticleType((*ta)[a].ffType()))
                {
                    msg_handler->msg(ACTStatus::Error,
                                     gmx::formatString("No such particletype %s in force field",
                                                       (*ta)[a].ffType().c_str()));
                }
                auto ptype = pd->findParticleType((*ta)[a].ffType());
                double qnew = ptype->charge();
                if (msg_handler->debug())
                {
                    msg_handler->writeDebug(gmx::formatString("Will update charge for %s from %g to %g",
                                                              (*ta)[a].ffType().c_str(), (*ta)[a].charge(), qnew));
                }
                (*ta)[a].setCharge(qnew);
                (*atoms)[atomStart_[ff]+a].setCharge(qnew);
                // The second charge, if any, is part of the total too
                (*ta)[a].setChargeVS(ptype->chargeVS(), ptype->zetaVS(), ptype->rowVS());
                (*atoms)[atomStart_[ff]+a].setChargeVS(ptype->chargeVS(), ptype->zetaVS(), ptype->rowVS());
                qtot += qnew + ptype->chargeVS();
            }
            // Make sure that the total charge is integer.
            double dq = (qtotal_[ff]-qtot)/ta->size();
            for(size_t a = 0; a < ta->size(); a++)
            {
                double qnew = (*ta)[a].charge() + dq;
                (*ta)[a].setCharge(qnew);
                (*atoms)[atomStart_[ff]+a].setCharge(qnew);
            }
        }
        break;
    case ChargeGenerationAlgorithm::Read:
        break;
    default: // throws
        GMX_THROW(gmx::InvalidInputError(gmx::formatString("No support for %s algorithm for fragments", 
                                                           chargeGenerationAlgorithmName(algorithm_).c_str()).c_str()));
    }
    // The widths of the charge distributions are those of the particle
    // types, which may have been trained, whatever the charges are
    for(size_t ff = 0; ff < topologies_.size(); ++ff)
    {
        auto        ta     = topologies_[ff].atomsPtr();
        const auto &ptypes = particleTypes(pd, ff);
        for(size_t a = 0; a < ta->size(); a++)
        {
            if (ptypes[a])
            {
                double zeta = ptypes[a]->zeta();
                (*ta)[a].setZeta(zeta);
                (*atoms)[atomStart_[ff]+a].setZeta(zeta);
            }
        }
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
        const auto &ats = topologies_[i].atoms();
        for(size_t j = 0; j < ats.size(); j++)
        {
            (*atoms)[k].setCharge(ats[j].charge());
            (*atoms)[k].setZeta(ats[j].zeta());
            (*atoms)[k++].setChargeVS(ats[j].chargeVS(), ats[j].zetaVS(), ats[j].rowVS());
        }
    }
    return true;
}

void FragmentHandler::setCharges(MsgHandler      *msghandler,
                                 const ChargeMap &qmap)
{
    for(size_t i = 0; i < ids_.size() && msghandler->ok(); i++)
    {
        const auto qptr = qmap.find(ids_[i]);
        if (qmap.end() != qptr)
        {
            auto aptr = topologies_[i].atomsPtr();
            // Complicated loop since model in topology may contain shells or vsites
            // but at least the atoms should match.
            size_t q = 0;
            for(size_t a = 0; msghandler->ok() && a < aptr->size(); a++)
            {
                if ((*aptr)[a].pType() == ActParticle::Atom)
                {
                    if (!((*aptr)[a].id() == qptr->second[q].first))
                    {
                        msghandler->msg(ACTStatus::Warning,
                                        gmx::formatString("Atom mismatch when reading charges from chargemap for %s. Expected %s but found %s. Make sure the atoms in your chargemap match those in your molprop file.\n",
                                                          ids_[i].c_str(), (*aptr)[a].id().id().c_str(),
                                                          qptr->second[a].first.id().c_str()));
                    }
                    if (msghandler->ok())
                    {
                        (*aptr)[a].setCharge(qptr->second[q].second);
                        msghandler->msg(ACTStatus::Debug,
                                        gmx::formatString("qmap Charge %zu = %g\n", a,
                                                          (*aptr)[a].charge()));
                    }
                }
                q += 1;
            }
        }
        else
        {
            msghandler->msg(ACTStatus::Warning,
                            gmx::formatString("Cannot find compound '%s' in the charge map\n",
                                              ids_[i].c_str()));
        }
    }
    fixedQ_ = msghandler->ok();
    // Set algorithm to reading charges for the future
    msghandler->msg(ACTStatus::Debug,
                    "Copied charges from chargemap to fragments\n");
    if (msghandler->ok())
    {
        algorithm_ = ChargeGenerationAlgorithm::Read;
    }
}

void FragmentHandler::setCharges(const std::vector<ActAtom> &atoms)
{
    for(size_t ff = 0; ff < topologies_.size(); ++ff)
    {
        auto aptr = topologies_[ff].atomsPtr();
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
        auto aptr = topologies_[ff].atomsPtr();
        for(size_t a = 0; a < aptr->size(); a++)
        {
            (*aptr)[a].setCharge(q[j]);
            j += 1;
        }
    }
}

} // namespace alexandria

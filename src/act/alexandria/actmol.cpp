/*
 * This source file is part of the Alexandria Chemistry Toolkit.
 *
 * Copyright (C) 2014-2024
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

#include "actmol.h"

#include <cstdio>

#include <algorithm>
#include <map>
#include <random>
#include <string>

#include "act/alexandria/actmol_low.h"
#include "act/alexandria/gromacs_top.h"
#include "act/alexandria/pdbwriter.h"
#include "act/alexandria/symmetrize_charges.h"
#include "act/forcefield/forcefield_parameter.h"
#include "act/forcefield/forcefield_parametername.h"
#include "act/molprop/experiment.h"
#include "act/molprop/molprop_util.h"
#include "act/molprop/multipole_names.h"
#include "act/utility/regression.h"
#include "act/utility/units.h"
#include "gromacs/gmxpreprocess/grompp-impl.h"
#include "gromacs/math/vec.h"
#include "gromacs/math/vectypes.h"
#include "gromacs/pbcutil/pbc.h"
#include "gromacs/utility/fatalerror.h"
#include "gromacs/utility/futil.h"

namespace alexandria
{

ACTQprop::ACTQprop(const std::vector<ActAtom>   &atoms,
                   const std::vector<gmx::RVec> &x)
{
    qPqm_.setQandX(atoms, x);
    qPact_.setQandX(atoms, x);
    qPact_.initializeMoments();
    for(const auto &aa : atoms)
    {
        bool b = (eptAtom == aa.pType());
        isAtom_.push_back(b);
    }
}

void ACTQprop::copyRespQ()
{
    if (QgenResp_.nEsp() > 0)
    {
        std::vector<double> q = qPqm_.charge();
        int j = 0;
        for(size_t i = 0; i < isAtom_.size(); i++)
        {
            if (isAtom_[i])
            {
                GMX_RELEASE_ASSERT(j < QgenResp_.natoms(), "Atom index out of range copying charges from RESP to Qprops");
                q[i] = QgenResp_.getCharge(j++);
            }
        }
        qPqm_.setQ(q);
        qPact_.setQ(q);
    }
}

ACTMol::ACTMol()
{
    clear_mat(box_);
    for (int m = 0; m < DIM; m++)
    {
        box_[m][m] = 10.0;
    }
}

void ACTMol::findInPlaneAtoms(int ca, std::vector<int> *atoms) const
{
    int bca = 0;
    /*First try to find the atom bound to the central atom (ca).*/
    for (auto &bi : bondsConst())
    {
        if (ca == bi.aJ() ||
            ca == bi.aI())
        {
            if (ca == bi.aI())
            {
                bca = bi.aJ();
                atoms->push_back(bca);
            }
            else
            {
                bca = bi.aI();
                atoms->push_back(bca);
            }
        }
    }
    /* Now try to find atoms bound to bca, except ca. */
    for (auto bi : bondsConst())
    {
        if ((ca != bi.aJ()   &&
             ca != bi.aI())  &&
            (bca == bi.aJ()  ||
             bca == bi.aI()))
        {
            if (bca == bi.aI())
            {
                atoms->push_back(bi.aJ());
            }
            else
            {
                atoms->push_back(bi.aI());
            }
        }
    }
}

void ACTMol::findOutPlaneAtoms(int ca, std::vector<int> *atoms) const
{
    for (auto &bi : bondsConst())
    {
        if (bi.bondOrder() == 1  &&
            (ca == bi.aJ() ||
             ca == bi.aI()))
        {
            if (ca == bi.aI())
            {
                atoms->push_back(bi.aJ());
            }
            else
            {
                atoms->push_back(bi.aI());
            }
        }
    }
}

bool ACTMol::IsVsiteNeeded(std::string       atype,
                           const ForceField *pd) const
{
    auto vsite = pd->findVsite(atype);
    return vsite != pd->getVsiteEnd();
}

immStatus ACTMol::checkAtoms(const ForceField *pd)
{
    int nmissing        = 0;
    int atomnumberTotal = 0;
    auto myatoms = topology_->atoms();
    for (size_t i = 0; i < myatoms.size(); i++)
    {
        auto atype = myatoms[i].ffType();
        if (!pd->hasParticleType(atype))
        {
            printf("Could not find a force field entry for atomtype %s atom %zu in compound '%s'\n",
                   atype.c_str(), i+1, getMolname().c_str());
            nmissing++;
        }
        else
        {
            atomnumberTotal += pd->findParticleType(atype)->atomnumber();
        }
    }
    if (nmissing > 0)
    {
        return immStatus::AtomTypes;
    }
    // Check multiplicity
    int multOk = atomnumberTotal + totalMultiplicity() + totalCharge();
    if (multOk % 2 == 0)
    {
        fprintf(stderr, "WARNING: atomnumberTotal %d, totalMultiplicity %d, totalCharge %d for %s\n",
                atomnumberTotal, totalMultiplicity(), totalCharge(),
                getMolname().c_str());
        return immStatus::Multiplicity;
    }
    return immStatus::OK;
}

static std::vector<gmx::RVec> experCoords(const std::vector<gmx::RVec> &xxx,
                                          const Topology               *topology)
{
    auto myatoms = topology->atoms();
    gmx::RVec fzero = { 0, 0, 0 };
    std::vector<gmx::RVec> coords(myatoms.size(), fzero);
    int j = 0;
    for(size_t i = 0; i < myatoms.size(); i++)
    {
        if (myatoms[i].pType() == eptAtom)
        {
            // Read coords from the experimental structure without shells or vsites
            copy_rvec(xxx[j], coords[i]);
            j += 1;
        }
    }
    // We need to split the loop to do all the atoms first
    // since we cannot be sure about the order of atoms and shells.
    for(size_t i = 0; i < myatoms.size(); i++)
    {
        switch (myatoms[i].pType())
        {
        case eptAtom:
            // Do nothing
            break;
        case eptShell:
        case eptVSite:
            {
                auto cores = myatoms[i].cores();
                GMX_RELEASE_ASSERT(!cores.empty(), "Shell or vsite without core");
                copy_rvec(coords[cores[0]], coords[i]);
            }
            break;
        default:
            {
                GMX_THROW(gmx::InternalError(gmx::formatString("Don't know how to handle %s particle type", ptype_str[myatoms[i].pType()]).c_str()));
            }
        }
    }
    ForceComputer fcomp;
    fcomp.generateVsites(topology, &coords);
    return coords;
}

std::vector<gmx::RVec> ACTMol::xOriginal() const
{
    auto exper = findExperimentConst(JobType::OPT);
    if (nullptr == exper)
    {
        exper = findExperimentConst(JobType::TOPOLOGY);
    }
    if (nullptr == exper)
    {
        exper = findExperimentConst(JobType::SP);
        if (nullptr == exper)
        {
            GMX_THROW(gmx::InternalError(gmx::formatString("No structure at all for %s", getMolname().c_str()).c_str()));
        }
        else
        {
            fprintf(stderr, "Warning: No calculation for %s with jobtype %s or %s, using first %s calc.\n",
                    getMolname().c_str(), jobType2string(JobType::OPT), jobType2string(JobType::TOPOLOGY),
                    jobType2string(JobType::SP));
        }
    }
    return experCoords(exper->getCoordinates(), topology_);
}

void ACTMol::forceEnergyMaps(const ForceField                                                    *pd,
                             const ForceComputer                                                 *forceComp,
                             std::vector<std::vector<std::pair<double, double> > >               *forceMap,
                             std::vector<ACTEnergy>                                              *energyMap,
                             std::vector<std::pair<double, std::map<InteractionType, double> > > *interactionEnergyMap,
                             std::vector<std::pair<double, std::map<InteractionType, double> > > *energyComponentMap) const
{
    auto       myatoms = topology_->atoms();
    forceMap->clear();
    energyMap->clear();
    interactionEnergyMap->clear();
    energyComponentMap->clear();
    gmx::RVec fzero = { 0, 0, 0 };
    for (auto &ei : experimentConst())
    {
        auto coords = experCoords(ei.getCoordinates(), topology_);
        // We compute either interaction energies or normal energies for one experiment
        if (ei.hasProperty(MolPropObservable::INTERACTIONENERGY))
        {
            auto eprops = ei.propertyConst(MolPropObservable::INTERACTIONENERGY);
            if (eprops.size() > 1)
            {
                gmx_fatal(FARGS, "Multiple interaction energies for this experiment");
            }
            else if (eprops.size() == 1)
            {
                if (forceComp)
                {
                    std::vector<gmx::RVec> interactionForces(myatoms.size(), fzero);
                    std::vector<gmx::RVec> mycoords(myatoms.size(), fzero);
                    auto exp_coords = ei.getCoordinates();
                    int eindex = 0;
                    for(size_t i = 0; i < myatoms.size(); i++)
                    {
                        if (myatoms[i].pType() == eptAtom)
                        {
                            copy_rvec(exp_coords[eindex], mycoords[i]);
                            eindex += 1;
                        }
                        else if (myatoms[i].pType() == eptShell && i > 0)
                        {
                            // TODO fix this hack
                            copy_rvec(mycoords[i-1], mycoords[i]);
                        }
                    }
                    std::map<InteractionType, double> einter;
                    calculateInteractionEnergy(pd, forceComp, &einter, &interactionForces, &mycoords);
                    interactionEnergyMap->push_back({ eprops[0]->getValue(), einter});
                    // TODO Store the interaction forces
                }
            }
        }
        else if (ei.hasProperty(MolPropObservable::DELTAE0))
        {
            std::map<InteractionType, double> energies;
            std::vector<gmx::RVec> forces(myatoms.size(), fzero);
            (void) forceComp->compute(pd, topology_, &coords, &forces, &energies);
            auto eprops = ei.propertyConst(MolPropObservable::DELTAE0);
            if (eprops.size() > 1)
            {
                gmx_fatal(FARGS, "Multiple energies for this experiment");
            }
            else if (eprops.size() == 1)
            {
                energyMap->push_back(ACTEnergy(ei.id(), eprops[0]->getValue(), energies[InteractionType::EPOT]));
                energyComponentMap->push_back({ eprops[0]->getValue(), std::move(energies) });
            }
        
            const std::vector<gmx::RVec> &fff = ei.getForces();
            if (!fff.empty())
            {
                size_t ifff = 0;
                std::vector<std::pair<double, double> > thisForce;
                for (size_t i = 0; i < myatoms.size(); i++)
                {
                    if (myatoms[i].pType() == eptAtom)
                    {
                        if (ifff >= fff.size())
                        {
                            GMX_THROW(gmx::InternalError(gmx::formatString("Inconsistency: there are %lu atoms and shells, but only %zu forces", myatoms.size(), fff.size())));
                        }
                        for(int m = 0; m < DIM; m++)
                        {
                            thisForce.push_back({ fff[ifff][m], forces[i][m] });
                        }
                        ifff += 1;
                    }
                }
                forceMap->push_back(std::move(thisForce));
            }
        }
    }
}

static bool isLinearMolecule(const std::vector<ActAtom>   &myatoms,
                             const std::vector<gmx::RVec> &coords)
{
    std::vector<int> core;
    for(size_t i = 0; i < myatoms.size(); i++)
    {
        if (myatoms[i].pType() == eptAtom)
        {
            core.push_back(i);
        }
    }
    if (core.size() <= 2)
    {
        return true;
    }
    t_pbc *pbc    = nullptr;
    real th_toler = 175; // Degrees
    bool linear   = true;
    for(size_t c = 2; c < core.size(); c++)
    {
        linear = linear && is_linear(coords[core[c-2]],
                                     coords[core[c-1]], 
                                     coords[core[c]], pbc, th_toler);
        if (!linear)
        {
            break;
        }
    }
    return linear;
}

immStatus ACTMol::GenerateTopology(gmx_unused FILE   *fp,
                                   const ForceField  *pd,
                                   missingParameters  missing)
{
    immStatus   imm = immStatus::OK;
    std::string btype1, btype2;

    if (nullptr != debug)
    {
        fprintf(debug, "Generating topology for %s\n", getMolname().c_str());
    }
    generateComposition();
    if (NAtom() <= 0)
    {
        imm = immStatus::AtomTypes;
    }
    /* Store bonds in harmonic potential list first, update type later */
    if (immStatus::OK == imm)
    {
        topology_ = new Topology(*bonds());
    }
    // Get the coordinates.
    std::vector<gmx::RVec> coords = xOriginal();

    if (immStatus::OK == imm)
    {
        imm = topology_->GenerateAtoms(pd, this,  &coords);
    }
    if (immStatus::OK == imm)
    {
        imm = checkAtoms(pd);
    }
    // Create fragments before adding shells!
    if (immStatus::OK == imm)
    {
        fraghandler_ = new FragmentHandler(pd, coords, topology_->atoms(),
                                           bondsConst(), fragmentPtr(), missing);
        // Finally, extract frequencies etc.
        getHarmonics();
    }
    
    if (immStatus::OK == imm)
    {
        topology_->build(pd, &coords, 175.0, 5.0, missing);
    }
    auto myatoms = topology()->atoms();
    if (immStatus::OK == imm)
    {
        realAtoms_.clear();
        for(size_t i = 0; i < myatoms.size(); i++)
        {
            if (myatoms[i].pType() == eptAtom)
            {
                realAtoms_.push_back(i);
            }
        }
    }
    if (immStatus::OK == imm && missing != missingParameters::Generate)
    {
        std::vector<InteractionType> itUpdate;
        for(auto &entry : *(topology_->entries()))
        {
            itUpdate.push_back(entry.first);
        }
    }
    if (immStatus::OK == imm)
    {
        imm = checkAtoms(pd);
    }
    if (immStatus::OK != imm && debug)
    {
        for(const auto &emsg : error_messages_)
        {
            fprintf(debug, "%s\n", emsg.c_str());
        }
    }
    if (immStatus::OK == imm)
    {
        isLinear_ = isLinearMolecule(myatoms, coords);
        // Symmetrize the atoms
        get_symmetrized_charges(topology_, pd, nullptr, &symmetric_charges_);
    }
    return imm;
}

double ACTMol::bondOrder(int ai, int aj) const
{
    return topology_->findBond(ai, aj).bondOrder();
}

static void printEmap(FILE *fp, const std::map<InteractionType, double> *e)
{
    for(auto ee = e->begin(); ee != e->end(); ee++)
    {
        fprintf(fp, " %s %g", interactionTypeToString(ee->first).c_str(), ee->second);
    }
    fprintf(fp, "\n");
}

void ACTMol::calculateInteractionEnergy(const ForceField                  *pd,
                                        const ForceComputer               *forceComputer,
                                        std::map<InteractionType, double> *einter,
                                        std::vector<gmx::RVec>            *interactionForces,
                                        std::vector<gmx::RVec>            *coords) const
{
    auto &tops = fraghandler_->topologies();
    einter->clear();
    if (tops.size() <= 1)
    {
        return;
    }
    // First, compute the total energy
    gmx::RVec fzero = { 0, 0, 0 };
    interactionForces->resize(coords->size(), fzero);
    if (debug)
    {
        fprintf(debug, "Will compute interaction energy\n");
    }
    (void) forceComputer->compute(pd, topology_, coords, interactionForces, einter);
    double edimer = (*einter)[InteractionType::EPOT];
    if (debug)
    {
        fprintf(debug, "%s: edimer = %g\n", getMolname().c_str(), edimer);
    }
    // Now compute interaction energies if there are fragments
    auto   &astart = fraghandler_->atomStart();
    for(size_t ff = 0; ff < tops.size(); ff++)
    {
        int natom = tops[ff]->atoms().size();
        std::vector<gmx::RVec> forces(natom, fzero);
        std::vector<gmx::RVec> myx(natom);
        int j = 0;
        for (size_t i = astart[ff]; i < astart[ff]+natom; i++)
        {
            copy_rvec((*coords)[i], myx[j]);
            j++;
        }
        std::map<InteractionType, double> energies;
        (void) forceComputer->compute(pd, tops[ff], &myx, &forces, &energies);
        edimer -= energies[InteractionType::EPOT];
        
        if (debug)
        {
            fprintf(debug, "%s initial:", getMolname().c_str());
            printEmap(debug, einter);
        }
        for(const auto &ee : energies)
        {
            auto eptr = einter->find(ee.first);
            if (einter->end() != eptr)
            {
                eptr->second -= ee.second;
            }
            else
            {
                einter->insert({ee.first, -ee.second});
            }
        }
        j = 0;
        for (size_t i = astart[ff]; i < astart[ff]+natom; i++)
        {
            for(int m = 0; m < DIM; m++)
            {
                (*interactionForces)[i][m] -= forces[j][m];
            }
            j++;
        }
        if (debug)
        {
            fprintf(debug, "%s Fragment %zu ", getMolname().c_str(), ff);
            printEmap(debug, &energies);
        }
    }
    if (debug)
    {
        fprintf(debug, "%s result:", getMolname().c_str());
        printEmap(debug, einter);
    }
}

immStatus ACTMol::GenerateAcmCharges(const ForceField       *pd,
                                     const ForceComputer    *forceComp,
                                     std::vector<gmx::RVec> *coords,
                                     std::vector<gmx::RVec> *forces)
{
    std::vector<double> qold;
    fraghandler_->fetchCharges(&qold);
    if (qold.size() != atomsConst().size())
    {
        GMX_THROW(gmx::InternalError(gmx::formatString("Cannot fetch old charges for %s. #atom %lu #qold %zu",
                                                       getMolname().c_str(), atomsConst().size(), qold.size()).c_str()));
    }
    immStatus imm       = immStatus::OK;
    if (fraghandler_->fixedCharges())
    {
        // We do not need to derive charges again, since they were set once already.
        return imm;
    }
    int       iter      = 0;
    bool      converged = true;
    double    EemRms    = 0;
    auto      natom     = atomsConst().size();
    std::map<InteractionType, double> energies;
    do
    {
        EemRms = 0;
        auto eqgen = fraghandler_->generateCharges(debug, getMolname(),
                                                   *coords, pd, atoms());
        if (eQgen::OK == eqgen)
        {
            (void) forceComp->compute(pd, topology_, coords, forces, &energies);
            std::vector<double> qnew;
            fraghandler_->fetchCharges(&qnew);
            apply_symmetrized_charges(&qnew, symmetric_charges_);
            GMX_RELEASE_ASSERT(qold.size()==qnew.size(), "Cannot fetch new charges");
            for (size_t i = 0; i < qnew.size(); i++)
            {
                EemRms  += gmx::square(qnew[i] - qold[i]);
                qold[i]  = qnew[i];
            }
            EemRms   /= natom;
            converged = (EemRms < qTolerance_) || !haveShells();
        }
        else if (eQgen::NOSUPPORT == eqgen)
        {
            imm = immStatus::MissingChargeGenerationParameters;
        }
        else
        {
            imm = immStatus::ChargeGeneration;
        }
        iter++;
    }
    while (imm == immStatus::OK && (!converged) && (iter < maxQiter_));
    if (immStatus::OK == imm)
    {
        if (!converged)
        {
            printf("Alexandria Charge Model did not converge to %g. rms: %g\n",
                   qTolerance_, sqrt(EemRms));
        }
        auto myatoms = atoms();
        for(size_t i = 0; i < natom; i++)
        {
            (*myatoms)[i].setCharge(qold[i]);
        }
        // Loop over qtype properties
        for(auto qp = qProps_.begin(); qp < qProps_.end(); ++qp)
        {
            qp->qPact()->setQ(atomsConst());
            // TODO Is this correct? Should not each qp have it's own coordinates?
            //qp->qPact()->setX(*coords);
        }
    }
    return imm;
}

immStatus ACTMol::GenerateCharges(const ForceField          *pd,
                                  const ForceComputer       *forceComp,
                                  ChargeGenerationAlgorithm  algorithm,
                                  qType                      qtype,
                                  const std::vector<double> &qcustom,
                                  std::vector<gmx::RVec>    *coords,
                                  std::vector<gmx::RVec>    *forces,
                                  bool                       updateQprops)
{
    immStatus imm         = immStatus::OK;
    bool      converged   = false;

    // TODO check whether this needed
    std::map<InteractionType, double> energies;
    auto myatoms = atoms();
    if (algorithm == ChargeGenerationAlgorithm::Custom)
    {
        GMX_RELEASE_ASSERT(atomsConst().size() == qcustom.size(),
                           gmx::formatString("Number of particles (%lu) does not match the number of custom charges (%lu).", atomsConst().size(), qcustom.size()).c_str());
    }
    else if (algorithm == ChargeGenerationAlgorithm::NONE)
    {
        algorithm = pd->chargeGenerationAlgorithm();
        // Check whether there are free charges
        bool allFixed = true;
        for (size_t i = 0; i < myatoms->size(); i++)
        {
            auto atype = (*myatoms)[i].ffType();
            auto ptype = pd->findParticleType(atype);
            auto qff = ptype->parameterConst("charge");
            if (qff.mutability() == Mutability::ACM)
            {
                allFixed = false;
            }
        }    
        if (allFixed)
        {
            algorithm = ChargeGenerationAlgorithm::NONE;
        }
    }
    fraghandler_->setChargeGenerationAlgorithm(algorithm);
    switch (algorithm)
    {
    case ChargeGenerationAlgorithm::NONE:
        {
            if (debug)
            {
                fprintf(debug, "WARNING! Using fixed charges for %s!\n",
                        getMolname().c_str());
            }
            for (size_t i = 0; i < myatoms->size(); i++)
            {
                auto atype = (*myatoms)[i].ffType();
                auto ptype = pd->findParticleType(atype);
                auto qval  = ptype->parameterConst("charge").value();
                (*myatoms)[i].setCharge(qval);
            }
            // If we have shells, we still have to minimize them,
            // but we may want to know the energies anyway.
            (void) forceComp->compute(pd, topology_, coords,
                                      forces, &energies);
            if (haveShells() && updateQprops)
            {
                for(auto qp = qProps_.begin(); qp < qProps_.end(); ++qp)
                {
                    auto qcalc = qp->qPact();
                    qcalc->setQ(*myatoms);
                    auto myx = qcalc->x();
                    (void) forceComp->compute(pd, topology_, &myx, forces, &energies);
                    // TODO, likely we should not change the coordinates here, just the charges
                    qcalc->setX(myx);
                }
            }
            fraghandler_->setCharges(*myatoms);
            
            return immStatus::OK;
        }
        break;
    case ChargeGenerationAlgorithm::Read:
        {
            std::vector<double> qread;
            auto exper = findExperimentConst(JobType::OPT);
            if (nullptr == exper)
            {
                exper  = findExperimentConst(JobType::TOPOLOGY);
            }
            if (nullptr != exper)
            {
                for (auto &ca : exper->calcAtomConst())
                {
                    if (!ca.hasCharge(qtype))
                    {
                        break;
                    }
                    else
                    {
                        qread.push_back(ca.charge(qtype));
                    }
                }
            }
            if (qread.empty())
            {
                return immStatus::NoMolpropCharges;
            }
            size_t j = 0;
            for(size_t i = 0; i < myatoms->size(); i++)
            {
                if ((*myatoms)[i].pType() == eptAtom)
                {
                    (*myatoms)[i].setCharge(qread[j++]);
                }
                else if ((*myatoms)[i].pType() == eptShell)
                {
                    const auto &pId = (*myatoms)[i].ffType();
                    if (!pd->hasParticleType(pId))
                    {
                        return immStatus::AtomTypes;
                    }
                    auto piter = pd->findParticleType(pId);
                    auto q = piter->charge();
                    (*myatoms)[i].setCharge(q);
                    // TODO Do not use -1 here, but particle.core()
                    (*myatoms)[i-1].setCharge(qread[j-1]-q);
                }
            }
            fraghandler_->setCharges(*myatoms);
            return immStatus::OK;
        }
    case ChargeGenerationAlgorithm::Custom:
        {
            for (size_t i = 0; i < myatoms->size(); i++)
            {
                (*myatoms)[i].setCharge(qcustom[i]);
            }
            fraghandler_->setCharges(*myatoms);
            return immStatus::OK;
        }
    case ChargeGenerationAlgorithm::ESP:
        {
            int    maxiter = 5;
            int    iter    = 0;
            // Only fit charges to the first ESP dataset.
            bool   foundESP = false;
            
            // Init Qgresp should be called before this!
            for(auto qp = qProps_.begin(); qp < qProps_.end() && !foundESP; ++qp)
            {
                auto qresp = qp->qgenResp();
                if (qresp->nEsp() == 0)
                {
                    continue;
                }
                foundESP = true;
                qresp->setAtomSymmetry(symmetric_charges_);
                double epsilonr;
                if (!ffOption(*pd, InteractionType::COULOMB, "epsilonr", &epsilonr))
                {
                    epsilonr = 1;
                }
                // Use the coordinate belonging to this RESP data
                auto myx = qresp->coords();
                (void) forceComp->compute(pd, topology_, &myx, forces, &energies);
                qresp->updateAtomCoords(myx);
                qresp->optimizeCharges(epsilonr);
                qresp->calcPot(epsilonr);
                qp->copyRespQ();
                do
                {
                    // Check whether charges have changed
                    double dq2 = 0;
                    for (size_t i = 0; i < myatoms->size(); i++)
                    {
                        // TODO what about shells?
                        double qrq = qresp->getCharge(i);
                        dq2 += gmx::square((*myatoms)[i].charge() - qrq);
                        (*myatoms)[i].setCharge(qrq);
                    }
                    // Copy charges to topology
                    (void) forceComp->compute(pd, topology_, &myx, forces, &energies);
                    qresp->updateAtomCoords(myx);
                    qresp->optimizeCharges(epsilonr);
                    qresp->calcPot(epsilonr);
                    qp->copyRespQ();
                    // If the system is polarizable we need to iterate until convergence since different charges
                    // will give different shell positions which in turn will give different shell positions.
                    converged = (dq2 < qTolerance_) || !pd->polarizable();
                    iter++;
                }
                while ((!converged) && (iter < maxiter));
                if (converged)
                {
                    for (size_t i = 0; i < myatoms->size(); i++)
                    {
                        (*myatoms)[i].setCharge(qresp->getCharge(i));
                    }
                    // TODO not sure whether this is needed but why not.
                    *coords = myx;
                    // Copy charges to fragments
                    fraghandler_->setCharges(*myatoms);
                }
                else
                {
                    imm = immStatus::ChargeGeneration;
                }
            }
        }
        break;
    case ChargeGenerationAlgorithm::EEM:
    case ChargeGenerationAlgorithm::SQE:
        {
            imm = GenerateAcmCharges(pd, forceComp, coords, forces);
        }
        break;
    }
    return imm;
}

void ACTMol::PrintConformation(const char                   *fn,
                               const std::vector<gmx::RVec> &coords,
                               bool                          writeShells,
                               const matrix                  box)
{
    auto title = gmx::formatString("%s processed by ACT - The Alexandria Chemistry Toolkit", getMolname().c_str());

    int        model_nr      = 1;
    char       chain         = ' ';
    gmx_conect conect        = gmx_conect_init();
    auto       itype         = InteractionType::BONDS;
    auto       top           = topology();
    if (top->hasEntry(itype))
    {
        auto &bonds = top->entry(itype);
        for(const auto &b: bonds)
        {
            gmx_conect_add(conect, b->atomIndex(0), b->atomIndex(1));
        }
    }
    auto       epbc          = epbcNONE;
    if (det(box_) > 0)
    {
        epbc = epbcXYZ;
    }
    FILE *fp = gmx_ffopen(fn, "w");
    if (writeShells)
    {
        pdbWriter(fp, title.c_str(), topology_->atoms(), 
                  coords, topology_->residueNames(),
                  epbc, box, chain, model_nr, {},
                  conect);
    }
    else
    {
        pdbWriter(fp, title.c_str(), topology_->atoms(),
                  coords, topology_->residueNames(),
                  epbc, box, chain, model_nr, realAtoms_,
                  conect, true);
    }
    gmx_ffclose(fp);
    gmx_conect_done(conect);
}

static void add_tensor(std::vector<std::string> *commercials,
                       const char               *title,
                       const char               *unit,
                       const std::vector<double> &Q)
{
    double fac = convertFromGromacs(1.0, unit);
    char buf[256];
    snprintf(buf, sizeof(buf), "%s:\n"
             "; ( %6.2f %6.2f %6.2f )\n"
             "; (       %6.2f %6.2f )\n"
             "; (             %6.2f )\n",
             title,
             fac*Q[0], fac*Q[1], fac*Q[2],
             fac*Q[3], fac*Q[4], fac*Q[5]);
    commercials->push_back(buf);
}

void ACTMol::PrintTopology(const char                  *fn,
                          bool                          bVerbose,
                          const ForceField             *pd,
                          const ForceComputer          *forceComp,
                          const CommunicationRecord    *cr,
                          const std::vector<gmx::RVec> &coords,
                          const std::string            &method,
                          const std::string            &basis,
                          bool                          bITP)
{
    char                     buf[256];
    t_mols                   printmol;
    std::vector<std::string> commercials;
    std::vector<double>      vec;
    double                   T = -1;
    std::string              myref;
    auto &qt         = pd->findForcesConst(InteractionType::COULOMB);
    auto iChargeType = name2ChargeType(qt.optionValue("chargetype"));
    std::string              mylot       = makeLot(method, basis);

    FILE *fp = gmx_ffopen(fn, "w");
    
    if (getMolname().size() > 0)
    {
        printmol.name = strdup(getMolname().c_str());
    }
    else if (formula().size() > 0)
    {
        printmol.name = strdup(formula().c_str());
    }
    else
    {
        printmol.name = strdup("Unknown");
    }

    printmol.nr = 1;

    snprintf(buf, sizeof(buf), "Total Mass = %.3f (Da)", totalMass());
    commercials.push_back(buf);
    snprintf(buf, sizeof(buf), "Total Charge = %d (e)", totalCharge());
    commercials.push_back(buf);
    snprintf(buf, sizeof(buf), "Charge Type  = %s\n",
             chargeTypeName(iChargeType).c_str());
    commercials.push_back(buf);

    for(auto qp = qProps_.begin(); qp < qProps_.end(); ++qp)
    {    
        auto qcalc = qp->qPact();
        auto qelec = qp->qPqmConst();
        qcalc->setQ(atomsConst());
        qcalc->setX(coords);
        qcalc->initializeMoments();
        qcalc->calcMoments();
        
        T = -1;
        for(auto &mpo : mpoMultiPoles)
        {
            auto gp = qmProperty(mpo, T, JobType::OPT);
            if (gp)
            {
                if (qelec.hasMultipole(mpo))
                {
                    auto mymu = qelec.getMultipole(mpo);
                    commercials.push_back(gmx::formatString("%s %s (%s)\n",
                                                            mylot.c_str(), mpo_name(mpo), gp->getUnit()));
                    for(auto &fmp : formatMultipole(mpo, mymu))
                    {
                        commercials.push_back(fmp);
                    }
                }
                if (qcalc->hasMultipole(mpo))
                {
                    auto mymu = qcalc->getMultipole(mpo);
                    commercials.push_back(gmx::formatString("Alexandria %s (%s)\n", mpo_name(mpo), gp->getUnit()));
                    for(auto &fmp : formatMultipole(mpo, mymu))
                    {
                        commercials.push_back(fmp);
                    }
                }
            }
        }
    
        if (nullptr != cr)
        {
            qcalc->calcPolarizability(pd, topology(), forceComp);
            auto acalc = qcalc->polarizabilityTensor();
            std::vector<double> ac = { acalc[XX][XX], acalc[XX][YY], acalc[XX][ZZ],
                                       acalc[YY][YY], acalc[YY][ZZ], acalc[ZZ][ZZ] };
            auto unit = mpo_unit2(MolPropObservable::POLARIZABILITY);
            add_tensor(&commercials, "Alexandria Polarizability components (A^3)", unit, ac);
            
            snprintf(buf, sizeof(buf), "Alexandria Isotropic Polarizability: %.2f (A^3)\n",
                     qcalc->isotropicPolarizability());
            commercials.push_back(buf);
        
            snprintf(buf, sizeof(buf), "Alexandria Anisotropic Polarizability: %.2f (A^3)\n",
                     qcalc->anisotropicPolarizability());
            commercials.push_back(buf);
        
            T = -1;
            if (qelec.hasPolarizability())
            {
                auto aelec = qelec.polarizabilityTensor();
                std::vector<double> ae = { aelec[XX][XX], aelec[XX][YY], aelec[XX][ZZ],
                                           aelec[YY][YY], aelec[YY][ZZ], aelec[ZZ][ZZ] };
                snprintf(buf, sizeof(buf), "%s + Polarizability components (A^3)", mylot.c_str());
                add_tensor(&commercials, buf, unit, ae);
                snprintf(buf, sizeof(buf), "%s Isotropic Polarizability: %.2f (A^3)\n",
                         mylot.c_str(), qelec.isotropicPolarizability());
                commercials.push_back(buf);
                snprintf(buf, sizeof(buf), "%s Anisotropic Polarizability: %.2f (A^3)\n",
                         mylot.c_str(), qelec.anisotropicPolarizability());
                commercials.push_back(buf);
            }
        }
    }

    // TODO write a replacement for this function
    print_top_header(fp, pd, bHaveShells_, commercials, bITP);
    write_top(fp, printmol.name, topology_, pd);
    if (!bITP)
    {
        print_top_mols(fp, printmol.name, 1, &printmol);
    }
    if (bVerbose)
    {
        for (auto &entry : *(topology_->entries()))
        {
            auto &fs = pd->findForcesConst(entry.first);
            int ftype = fs.gromacsType();
            if (entry.second.size() > 0)
            {
                printf("There are %4d %s interactions\n", static_cast<int>(entry.second.size()),
                       interaction_function[ftype].name);
            }
        }
        for (auto i = commercials.begin(); (i < commercials.end()); ++i)
        {
            printf("%s\n", i->c_str());
        }
    }

    sfree(printmol.name);
    
    gmx_ffclose(fp);
}

void ACTMol::GenerateCube(const ForceField             *pd,
                          const std::vector<gmx::RVec> &coords,
                          const ForceComputer          *forceComp,
                          real                          spacing,
                          real                          border,
                          const char                   *reffn,
                          const char                   *pcfn,
                          const char                   *pdbdifffn,
                          const char                   *potfn,
                          const char                   *rhofn,
                          const char                   *hisfn,
                          const char                   *difffn,
                          const char                   *diffhistfn,
                          const gmx_output_env_t       *oenv)
{
    auto &qt         = pd->findForcesConst(InteractionType::COULOMB);
    auto iChargeType = name2ChargeType(qt.optionValue("chargetype"));

    if (potfn || hisfn || rhofn || difffn || pdbdifffn)
    {
        char   *act_version = (char *)"act v0.99b";
        double  epsilonr;
        if (!ffOption(*pd, InteractionType::COULOMB, 
                      "epsilonr", &epsilonr))
        {
            epsilonr = 1;
        }
        int index = 0;
        for(auto qc = qProps_.begin(); qc < qProps_.end(); ++qc)
        {
            auto qref  = qc->qPqmConst();
            auto qresp = qc->qgenResp();
            if (0 < qresp->nEsp())
            {
                // We have ESP data!
                auto qcalc = qc->qPact();
                // This code will overwrite files if there are more than one ESP data set
                qcalc->setQ(atomsConst());
                // TODO Check but we should likely not update coords
                // qcalc->setX(coords);
                // Relax shells, position vsites etc. etc.
                auto myx = qresp->coords();
                std::vector<gmx::RVec>            forces(atomsConst().size());
                std::map<InteractionType, double> energies;
                forceComp->compute(pd, topology(), &myx, &forces, &energies);
                qresp->updateAtomCoords(myx);
                qresp->updateAtomCharges(atomsConst());
                qresp->calcPot(epsilonr);
                qresp->potcomp(pcfn, oenv);
                if (pdbdifffn)
                {
                    std::string pdbx = gmx::formatString("%d_%s", index++, pdbdifffn);
                    qresp->writePdbComparison(atomsConst(), pdbx);
                }
                /* This has to be done before the grid is f*cked up by
                   writing a cube file */
                QgenResp qCalc(*qresp);
                QgenResp grref(*qresp);
            
                if (reffn)
                {
                    grref.setAtomInfo(atomsConst(), pd, totalCharge());
                    grref.updateAtomCoords(coords);
                    grref.setAtomSymmetry(symmetric_charges_);
                    grref.readCube(reffn, FALSE);
                }
                else
                {
                    qCalc.makeGrid(spacing, border, coords);
                }
                if (rhofn)
                {
                    std::string buf = gmx::formatString("Electron density generated by %s based on %s charges",
                                                        act_version,
                                                        chargeTypeName(iChargeType).c_str());
                    qCalc.calcRho();
                    qCalc.writeRho(rhofn, buf, oenv);
                }
                if (potfn)
                {
                    std::string buf = gmx::formatString("Potential generated by %s based on %s charges",
                                                        act_version,
                                                        chargeTypeName(iChargeType).c_str());
                    qCalc.calcPot(epsilonr);
                    qCalc.writeCube(potfn, buf, oenv);
                }
                if (hisfn)
                {
                    std::string buf = gmx::formatString("Potential generated by %s based on %s charges",
                                                        act_version,
                                                        chargeTypeName(iChargeType).c_str());
                    qCalc.writeHisto(hisfn, buf, oenv);
                }
                if (difffn || diffhistfn)
                {
                    std::string buf = gmx::formatString("Potential difference generated by %s based on %s charges",
                                                        act_version,
                                                        chargeTypeName(iChargeType).c_str());
                    
                    qCalc.writeDiffCube(&grref, difffn, diffhistfn, buf, oenv, 0);
                }
            }
        }
    }
}
        
void ACTMol::getHarmonics()
{
    for(auto &mpo : { MolPropObservable::FREQUENCY, 
                     MolPropObservable::INTENSITY })
    {
        std::vector<GenericProperty *> harm;
        for (auto &ee : experimentConst())
        {
            if (ee.hasMolPropObservable(mpo))
            {
                harm = ee.propertyConst(mpo);
                break;
            }
        }
        if (!harm.empty())
        {
            for(auto &ff : harm[0]->getVector())
            {
                if (mpo == MolPropObservable::FREQUENCY)
                {
                    ref_frequencies_.push_back(ff);
                }
                else
                {
                    ref_intensities_.push_back(ff);
                }
            }
        }
    }
}

immStatus ACTMol::getExpProps(const ForceField                           *pd,
                              const std::map<MolPropObservable, iqmType> &iqm,
                              real                                        watoms,
                              int                                         maxESP)
{
    std::vector<double> vec;
    std::string         myref;
    std::string         mylot;
    immStatus           imm        = immStatus::OK;
    
    auto &myatoms = atomsConst();
    GMX_RELEASE_ASSERT(myatoms.size() > 0, "No atoms!");
    
    bool foundNothing = true;
    for (const auto &myexp : experimentConst())
    {
        bool     qprop = false;
        ACTQprop actq;
        auto     qelec = actq.qPqm();
        auto     qcalc = actq.qPact();
        auto     props = myexp.propertiesConst();
        auto     xatom = experCoords(myexp.getCoordinates(), topology_);
        std::vector<double> q;
        for (auto prop : props)
        {
            // Check whether this property is in the "Wanted" list
            auto iter = iqm.find(prop.first);
            if (iqm.end() == iter)
            {
                continue;
            }
            // Fetch the property from this experiment
            auto gp = myexp.propertyConst(prop.first);
            switch (prop.first)
            {
            case MolPropObservable::CHARGE:
                {
                    std::string         reference;
                    std::string         lot;
                    if (myexp.getCharges(&q, qType::Calc, &reference, &lot))
                    {
                        qprop = true;
                    }
                    break;
                }
            case MolPropObservable::POTENTIAL:
                {
                    auto qgr         = actq.qgenResp();
                    auto qt          = pd->findForcesConst(InteractionType::COULOMB);
                    auto iChargeType = name2ChargeType(qt.optionValue("chargetype"));
                    qgr->setChargeType(iChargeType);
                    qgr->setAtomInfo(atomsConst(), pd, totalCharge());
                    qgr->updateAtomCoords(xatom);
                    qgr->setAtomSymmetry(symmetric_charges_);
                    qgr->summary(debug);
                    int natoms = nRealAtoms();
                        
                    std::random_device               rd;
                    std::mt19937                     gen(rd());  
                    std::uniform_real_distribution<> uniform(0.0, 1.0);
                    double                           cutoff = 0.01*maxESP;

                    auto eee = prop.second;
                    for (size_t k = 0; k < eee.size(); k++)
                    {
                        auto esp = static_cast<const ElectrostaticPotential *>(eee[k]);
                        auto xyz = esp->xyz();
                        auto V   = esp->V();
                        for(size_t ll = 0; ll < V.size(); ll++)
                        {
                            auto val = uniform(gen);
                            if ((ll >= static_cast<size_t>(natoms) || watoms > 0) && val <= cutoff)
                            {
                                qgr->addEspPoint(xyz[ll][XX], xyz[ll][YY], xyz[ll][ZZ], V[ll]);
                            }
                        }
                    }
                    qprop = true;
                }
                break;
            case MolPropObservable::DIPOLE:
            case MolPropObservable::QUADRUPOLE:
            case MolPropObservable::OCTUPOLE:
            case MolPropObservable::HEXADECAPOLE:
                {
                    for(auto mygp = gp.begin(); mygp < gp.end(); ++mygp)
                    {
                        auto multprop = static_cast<const MolecularMultipole *>(*mygp);
                        qelec->setMultipole(prop.first, multprop->getVector());
                        qprop = true;
                    }
                }
                break;
            case MolPropObservable::POLARIZABILITY:
                {
                    for(auto mygp = gp.begin(); mygp < gp.end(); ++mygp)
                    {
                        auto polprop = static_cast<const MolecularPolarizability *>(*mygp);
                        qelec->setPolarizabilityTensor(polprop->getTensor());
                        qprop = true;
                    }
                }
                break;
            case MolPropObservable::INTERACTIONENERGY:
                {
                    for(auto mygp = gp.begin(); mygp < gp.end(); ++mygp)
                    {
                        auto ieprop = static_cast<const MolecularEnergy *>(*mygp);
                        energy_.insert(std::pair<MolPropObservable, double>(prop.first, ieprop->getValue()));
                        foundNothing = false;
                    }
                }
                break;
            case MolPropObservable::DELTAE0:
            case MolPropObservable::DHFORM:
            case MolPropObservable::DGFORM:
            case MolPropObservable::ZPE:
                {
                    for(auto mygp = gp.begin(); mygp < gp.end(); ++mygp)
                    {
                        auto eprop = static_cast<const MolecularEnergy *>(*mygp);
                        energy_.insert(std::pair<MolPropObservable, double>(prop.first, eprop->getValue()));
                        foundNothing = false;
                    }
                }
                break;
            default:
                break;
            }
        }
        if (qprop)
        {
            if (q.empty())
            {
                q.resize(xatom.size(), 0.0);
                // TODO Check whether this is needed. Likely it is here, since it is the first time.
                qcalc->setQandX(atomsConst(), xatom);
            }
            else
            {
                // TODO Check whether this is needed. Likely it is here, since it is the first time.
                qcalc->setQandX(q, xatom);
            }
            qcalc->initializeMoments();
            qcalc->calcMoments();
            qProps_.push_back(std::move(actq));
            foundNothing = false;
        }
    }
    if (foundNothing)
    {
        ACTQprop actq(topology()->atoms(), xOriginal());
        qProps_.push_back(std::move(actq));
        imm = immStatus::NoData;
    }
    return imm;
}

CommunicationStatus ACTMol::Send(const CommunicationRecord *cr, int dest) const
{
    auto cs = MolProp::Send(cr, dest);
    if (CommunicationStatus::OK == cs)
    {
        cr->send_int(dest, static_cast<int>(dataset_type_));
    }
    return cs;
}

CommunicationStatus ACTMol::BroadCast(const CommunicationRecord *cr, int root, MPI_Comm comm)
{
    auto cs = MolProp::BroadCast(cr, root, comm);
    if (CommunicationStatus::OK == cs)
    {
        int ims = static_cast<int>(dataset_type_);
        cr->bcast(&ims, comm);
        set_datasetType(static_cast<iMolSelect>(ims));
    }
    return cs;
}

CommunicationStatus ACTMol::Receive(const CommunicationRecord *cr, int src)
{
    auto cs = MolProp::Receive(cr, src);
    if (CommunicationStatus::OK == cs)
    {
        set_datasetType(static_cast<iMolSelect>(cr->recv_int(src)));
    }
    return cs;
}

} // namespace alexandria

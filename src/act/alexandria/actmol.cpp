/*
 * This source file is part of the Alexandria Chemistry Toolkit.
 *
 * Copyright (C) 2014-2023
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

#include <map>
#include <random>
#include <string>

#include "act/molprop/molprop_util.h"
#include "act/molprop/multipole_names.h"
#include "act/forcefield/forcefield_parameter.h"
#include "act/forcefield/forcefield_parametername.h"
#include "act/utility/regression.h"
#include "act/utility/units.h"
#include "gromacs/commandline/filenm.h"
#include "gromacs/fileio/pdbio.h"
#include "gromacs/math/vec.h"
#include "gromacs/math/vecdump.h"
#include "gromacs/math/vectypes.h"
#include "gromacs/pbcutil/pbc.h"
#include "gromacs/utility/fatalerror.h"
#include "gromacs/utility/futil.h"
#include "gromacs/utility/strconvert.h"
#include "gromacs/utility/stringcompare.h"
#include "gromacs_top.h"
#include "actmol_low.h"
#include "symmetrize_charges.h"

namespace alexandria
{

ACTMol::ACTMol()
{
    clear_mat(box_);
    for (int m = 0; m < DIM; m++)
    {
        box_[m][m] = 10.0;
    }
}

bool ACTMol::IsSymmetric(real toler) const
{
    real  tm;
    rvec  com, test;
    bool  bSymmAll;
    std::vector<bool> bSymm;
    
    clear_rvec(com);
    tm = 0;
    const auto &myatoms = topology_->atoms();
    std::vector<gmx::RVec> myx = optimizedCoordinates_;
    for (size_t i = 0; i < myatoms.size(); i++)
    {
        real mm  = myatoms[i].mass();
        tm += mm;
        for (int m = 0; (m < DIM); m++)
        {
            com[m] += mm*myx[i][m];
        }
    }
    if (tm > 0)
    {
        for (int m = 0; m < DIM; m++)
        {
            com[m] /= tm;
        }
    }
    for (size_t i = 0; i < myatoms.size(); i++)
    {
        rvec_dec(myx[i], com);
    }

    bSymm.resize(myatoms.size());
    for (size_t i = 0; i < myatoms.size(); i++)
    {
        bSymm[i] = (norm(myx[i]) < toler);
        for (size_t j = i+1; (j < myatoms.size()) && !bSymm[i]; j++)
        {
            rvec_add(myx[i], myx[j], test);
            if (norm(test) < toler)
            {
                bSymm[i] = true;
                bSymm[j] = true;
            }
        }
    }
    bSymmAll = true;
    for (size_t i = 0; i < myatoms.size(); i++)
    {
        bSymmAll = bSymmAll && bSymm[i];
    }

    return bSymmAll;
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
        fprintf(stderr, "atomnumberTotal %d, totalMultiplicity %d, totalCharge %d\n",
                atomnumberTotal, totalMultiplicity(), totalCharge());
        return immStatus::Multiplicity;
    }
    return immStatus::OK;
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
    for (auto &ei : experimentConst())
    {
        gmx::RVec fzero = { 0, 0, 0 };
        // TODO: no need to recompute the energy if we just have
        // done that. Check for OPT being the first calculation.
        const std::vector<gmx::RVec> &xxx = ei.getCoordinates();
        std::vector<gmx::RVec> coords(myatoms.size(), fzero);
        int j = 0;
        for(size_t i = 0; i < myatoms.size(); i++)
        {
            if (myatoms[i].pType() == eptAtom)
            {
                // Read coords from the experimental structure without shells
                copy_rvec(xxx[j], coords[i]);
                j += 1;
            }
        }
        // We need to split the loop to do all the atoms first
        // since we cannot be sure about the order of atoms and shells.
        for(size_t i = 0; i < myatoms.size(); i++)
        {
            if (myatoms[i].pType() == eptAtom)
            {
                // Do nothing
            }
            else if (myatoms[i].pType() == eptShell)
            {
                auto core = myatoms[i].core();
                GMX_RELEASE_ASSERT(core != -1, "Shell without core");
                copy_rvec(coords[core], coords[i]);
            }
            else
            {
                GMX_THROW(gmx::InternalError(gmx::formatString("Don't know how to handle %s particle type", ptype_str[myatoms[i].pType()]).c_str()));
            }
        }
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
        topology_ = new Topology(bondsConst());
    }
    if (immStatus::OK == imm)
    {
        imm = topology_->GenerateAtoms(pd, this,  &optimizedCoordinates_);
    }
    if (immStatus::OK == imm)
    {
        imm = checkAtoms(pd);
    }
    // Create fragments before adding shells!
    if (immStatus::OK == imm)
    {
        fraghandler_ = new FragmentHandler(pd, optimizedCoordinates_, topology_->atoms(),
                                           bondsConst(), fragmentPtr(), missing);
        // Finally, extract frequencies etc.
        getHarmonics();
    }
    
    if (immStatus::OK == imm)
    {
        topology_->build(pd, &optimizedCoordinates_, 175.0, 5.0, missing);
    }
    if (immStatus::OK == imm)
    {
        qProps_.insert({ qType::Calc, QtypeProps(qType::Calc) });
        qProps_.insert({ qType::Elec, QtypeProps(qType::Elec) });
        qProps_.find(qType::Calc)->second.initializeMoments();
    }
    if (immStatus::OK == imm)
    {
        /* Center of charge */
        auto atntot = 0;
        rvec coc    = { 0 };
        auto myatoms = topology_->atoms();
        for (size_t i = 0; i < myatoms.size(); i++)
        {
            auto atn = myatoms[i].atomicNumber();
            atntot  += atn;
            for (auto m = 0; m < DIM; m++)
            {
                coc[m] += optimizedCoordinates_[i][m]*atn;
            }
        }
        /* Center of charge */
        svmul((1.0/atntot), coc, CenterOfCharge_);
        for (auto &qp : qProps_)
        {
            qp.second.setCenterOfCharge(CenterOfCharge_);
        }
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
    return imm;
}

bool ACTMol::linearMolecule() const
{
    const auto myatoms = topology_->atoms();
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
        linear = linear && is_linear(optimizedCoordinates_[core[c-2]],
                                     optimizedCoordinates_[core[c-1]], 
                                     optimizedCoordinates_[core[c]], pbc, th_toler);
        if (!linear)
        {
            break;
        }
    }
    return linear;
}

double ACTMol::bondOrder(int ai, int aj) const
{
    return topology_->findBond(ai, aj)->bondOrder();
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
    double edimer = forceComputer->compute(pd, topology_, coords, interactionForces, einter);
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
        edimer -= forceComputer->compute(pd, tops[ff], &myx, &forces, &energies);
        if (debug)
        {
            fprintf(debug, "%s: edimer = %g\n", getMolname().c_str(), edimer);
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
            fprintf(debug, "%s Fragment %zu Epot %g\n", getMolname().c_str(), ff, 
                    einter->find(InteractionType::EPOT)->second);
        }
    }
    if (debug)
    {
        fprintf(debug, "%s: edimer = %g\n", getMolname().c_str(), edimer);
    }
}

void ACTMol::symmetrizeCharges(const ForceField  *pd,
                               bool               bSymmetricCharges,
                               const char        *symm_string)
{
    if (bSymmetricCharges)
    {
        symmetric_charges_.clear();
        symmetrize_charges(bSymmetricCharges, topology_,
                           pd, symm_string, &symmetric_charges_);
    }
    else
    {
        auto natoms = atomsConst().size();
        for (size_t i = 0; i < natoms; i++)
        {
            symmetric_charges_.push_back(i);
        }
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
    int       iter      = 0;
    bool      converged = false;
    double    EemRms    = 0;
    auto      natom     = atomsConst().size();
    std::map<InteractionType, double> energies;
    do
    {
        if (eQgen::OK == fraghandler_->generateCharges(debug, getMolname(),
                                                       *coords, pd, atoms()))
        {
            (void) forceComp->compute(pd, topology_, coords, forces, &energies);
            EemRms = 0;
            std::vector<double> qnew;
            fraghandler_->fetchCharges(&qnew);
            GMX_RELEASE_ASSERT(qold.size()==qnew.size(), "Cannot fetch new charges");
            for (size_t i = 0; i < qnew.size(); i++)
            {
                EemRms  += gmx::square(qnew[i] - qold[i]);
                qold[i]  = qnew[i];
            }
            EemRms   /= natom;
            converged = (EemRms < qTolerance_) || !haveShells();
            iter++;
        }
        else
        {
            imm = immStatus::ChargeGeneration;
        }
    }
    while (imm == immStatus::OK && (!converged) && (iter < maxQiter_));
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
    auto qcalc = qTypeProps(qType::Calc);
    qcalc->setQ(atomsConst());
    qcalc->setX(*coords);
    return imm;
}

immStatus ACTMol::GenerateCharges(const ForceField          *pd,
                                  const ForceComputer       *forceComp,
                                  ChargeGenerationAlgorithm  algorithm,
                                  qType                      qtype,
                                  const std::vector<double> &qcustom,
                                  std::vector<gmx::RVec>    *coords,
                                  std::vector<gmx::RVec>    *forces)
{
    immStatus imm         = immStatus::OK;
    bool      converged   = false;

    // TODO check whether this needed
    std::map<InteractionType, double> energies;
    auto myatoms = atoms();
    if (algorithm == ChargeGenerationAlgorithm::Custom)
    {
        GMX_RELEASE_ASSERT(atomsConst().size() == qcustom.size(),
                           gmx::formatString("Number of custom charges %lu does not match the number of atoms %lu", qcustom.size(), atomsConst().size()).c_str());
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
            if (haveShells())
            {
                auto qcalc = qTypeProps(qType::Calc);
                qcalc->setQ(*myatoms);
                qcalc->setX(*coords);
            }
            
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

            return immStatus::OK;
        }
    case ChargeGenerationAlgorithm::ESP:
        {
            double chi2[2] = {1e8, 1e8};
            int    cur     = 0;
            int    maxiter = 5;
            int    iter    = 0;
            
            // Init Qgresp should be called before this!
            auto qcalc   = qTypeProps(qType::Calc);
            GMX_RELEASE_ASSERT(qcalc != nullptr, "qType::Calc is not initialized");
            double epsilonr;
            if (!ffOption(*pd, InteractionType::COULOMB, 
                          "epsilonr", &epsilonr))
            {
                epsilonr = 1;
            }

            qcalc->qgenResp()->optimizeCharges(epsilonr);
            qcalc->qgenResp()->calcPot(epsilonr);
            qcalc->copyRespQ();
            
            if (debug)
            {
                fprintf(debug, "RESP: RMS %g\n", chi2[cur]);
            }
            do
            {
                auto qq = qcalc->charge();
                GMX_RELEASE_ASSERT(myatoms->size() == qq.size(),
                                   gmx::formatString("Number of particles (%lu) differs from number of charges (%lu)", 
                                                     myatoms->size(),
                                                     qq.size()).c_str());
                for (size_t i = 0; i < myatoms->size(); i++)
                {
                    (*myatoms)[i].setCharge(qq[i]);
                }
                // Copy charges to topology
                chi2[cur] = forceComp->compute(pd, topology_, coords, forces, &energies);
                qcalc->setX(*coords);
                qcalc->qgenResp()->optimizeCharges(epsilonr);
                qcalc->qgenResp()->calcPot(epsilonr);
                qcalc->copyRespQ();
                if (debug)
                {
                    fprintf(debug, "RESP: RMS %g\n", chi2[cur]);
                }
                converged = (fabs(chi2[cur] - chi2[1-cur]) < qTolerance_) || !pd->polarizable();
                cur       = 1-cur;
                iter++;
            }
            while ((!converged) && (iter < maxiter));
            auto qq = qcalc->charge();
            for (size_t i = 0; i < myatoms->size(); i++)
            {
                (*myatoms)[i].setCharge(qq[i]);
            }
            // Copy charges to fragments
            fraghandler_->setCharges(*myatoms);
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

void ACTMol::CalcPolarizability(const ForceField    *pd,
                                const ForceComputer *forceComp)
{
    auto natoms = atomsConst().size();
    std::vector<gmx::RVec> coordinates(natoms);
    forceComp->calcPolarizability(pd, topology_, &coordinates,
                                  qTypeProps(qType::Calc));
}

void ACTMol::PrintConformation(const char                   *fn,
                               const std::vector<gmx::RVec> &coords,
                               bool                          writeShells,
                               const matrix                  box)
{
    char title[STRLEN];
    
    sprintf(title, "%s processed by ACT - The Alexandria Chemistry Tookit",
            getMolname().c_str());
    int        model_nr      = 1;
    char       chain         = ' ';
    gmx_bool   bTerSepChains = FALSE;
    gmx_conect conect        = gmx_conect_init();
    auto       itype         = InteractionType::BONDS;
    auto       top           = topology();
    if (top->hasEntry(itype))
    {
        auto bonds = top->entry(itype);
        for(const auto &b: bonds)
        {
            auto bb = static_cast<Bond *>(b);
            gmx_conect_add(conect, bb->aI(), bb->aJ());
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
        write_pdbfile(fp, title, nullptr, as_rvec_array(coords.data()),
                      epbc, box, chain, model_nr, conect, bTerSepChains);
    }
    else
    {
        bool usePqrFormat = false;
        write_pdbfile_indexed(fp, title, nullptr, as_rvec_array(coords.data()),
                              epbc, box, chain, model_nr, realAtoms_.size(),
                              realAtoms_.data(), conect, bTerSepChains, usePqrFormat, true);
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

void ACTMol::PrintTopology(const char                   *fn,
                          bool                          bVerbose,
                          const ForceField                *pd,
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
    
    auto qcalc = qTypeProps(qType::Calc);
    auto qelec = qTypeProps(qType::Elec);
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
            auto vec = gp->getVector();
            qelec->setMultipole(mpo, vec);
            auto mymu = qelec->getMultipole(mpo);
            commercials.push_back(gmx::formatString("%s %s (%s)\n",
                                                    mylot.c_str(), mpo_name(mpo), gp->getUnit()));
            for(auto &fmp : formatMultipole(mpo, mymu))
            {
                commercials.push_back(fmp);
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
        CalcPolarizability(pd, forceComp);
        auto qcalc = qTypeProps(qType::Calc);
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
        auto gp = qmProperty(MolPropObservable::POLARIZABILITY,
                             T, JobType::OPT);
        if (gp)
        {
            auto qelec = qTypeProps(qType::Elec);
            auto aelec = qelec->polarizabilityTensor();
            std::vector<double> ae = { aelec[XX][XX], aelec[XX][YY], aelec[XX][ZZ],
                                       aelec[YY][YY], aelec[YY][ZZ], aelec[ZZ][ZZ] };
            snprintf(buf, sizeof(buf), "%s + Polarizability components (A^3)", mylot.c_str());
            add_tensor(&commercials, buf, unit, ae);
            snprintf(buf, sizeof(buf), "%s Isotropic Polarizability: %.2f (A^3)\n",
                     mylot.c_str(), qelec->isotropicPolarizability());
            commercials.push_back(buf);
            snprintf(buf, sizeof(buf), "%s Anisotropic Polarizability: %.2f (A^3)\n",
                     mylot.c_str(), qelec->anisotropicPolarizability());
            commercials.push_back(buf);
        }
    }

    // TODO write a replacement for this function
    print_top_header(fp, pd, bHaveShells_, commercials, bITP);
    write_top(fp, printmol.name, nullptr,
              topology_, nullptr, nullptr, pd);
    if (!bITP)
    {
        print_top_mols(fp, printmol.name, pd->filename().c_str(),
                       nullptr, 0, nullptr, 1, &printmol);
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
        char     *act_version = (char *)"act v0.99b";
        auto qc = qProps_.find(qType::Calc);
        GMX_RELEASE_ASSERT(qc != qProps_.end(), "Cannot find alexandria charge information");
        qc->second.setQ(atomsConst());
        qc->second.setX(coords);
        double epsilonr;
        if (!ffOption(*pd, InteractionType::COULOMB, 
                      "epsilonr", &epsilonr))
        {
            epsilonr = 1;
        }
        qc->second.qgenResp()->calcPot(epsilonr);
        qc->second.qgenResp()->potcomp(pcfn, atomsConst(),
                                       as_rvec_array(coords.data()),
                                       pdbdifffn, oenv);

        /* This has to be done before the grid is f*cked up by
           writing a cube file */
        QgenResp qCalc(*qc->second.qgenResp());
        QgenResp grref(*qc->second.qgenResp());

        if (reffn)
        {
            grref.setAtomInfo(atomsConst(), pd, coords, totalCharge());
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

void ACTMol::calcEspRms(const ForceField             *pd,
                        const std::vector<gmx::RVec> *coords)
{
    int   natoms  = 0;
    double epsilonr;
    if (!ffOption(*pd, InteractionType::COULOMB, 
                  "epsilonr", &epsilonr))
    {
        epsilonr = 1;
    }
    auto &myatoms = atomsConst();
    for (size_t i = 0; i < myatoms.size(); i++)
    {
        if (myatoms[i].pType() == eptAtom)
        {
            natoms++;
        }
    }
    std::vector<gmx::RVec> myx(nRealAtoms());
    std::vector<ActAtom> realAtoms;
    size_t ii = 0;
    for (size_t i = 0; i < myatoms.size(); i++)
    {
        if (myatoms[i].pType() == eptAtom)
        {
            realAtoms.push_back(myatoms[i]);
            copy_rvec((*coords)[i], myx[ii]);
            ii++;
        }
    }
    
    auto qcalc   = qTypeProps(qType::Calc);
    auto qgrcalc = qcalc->qgenResp();
    for(auto &i : qProps_)
    {
        auto qi = i.first;
        if (qType::Calc == qi)
        {
            qgrcalc->updateAtomCharges(atomsConst());
            qgrcalc->updateAtomCoords(*coords);
            qgrcalc->calcPot(epsilonr);
        }
        else if (qType::Elec != qi)
        {
            QgenResp *qgr = i.second.qgenResp();
            qgr->setChargeType(ChargeType::Point);
            qgr->setAtomInfo(realAtoms, pd, myx, totalCharge());
            qgr->updateAtomCharges(i.second.charge());
            for (size_t j = realAtoms.size(); j < qgrcalc->nEsp(); j++)
            {
                auto &ep = qgrcalc->espPoint(j);
                auto &r  = ep.esp();
                qgr->addEspPoint(r[XX], r[YY], r[ZZ], ep.v());
            }
            qgr->calcPot(epsilonr);
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

immStatus ACTMol::getExpProps(const std::map<MolPropObservable, iqmType> &iqm,
                             double                                      T)
{
    int                 natom = 0;
    std::vector<double> vec;
    std::string         myref;
    std::string         mylot;
    immStatus           imm        = immStatus::OK;
    
    auto &myatoms = atomsConst();
    GMX_RELEASE_ASSERT(myatoms.size() > 0, "No atoms!");
    
    // Make a copy of the coordinates without shells
    std::vector<gmx::RVec> xatom(myatoms.size());
    for (size_t i = 0; i < myatoms.size(); i++)
    {
        if (myatoms[i].pType() == eptAtom ||
            myatoms[i].pType() == eptNucleus)
        {
            copy_rvec(xOriginal()[i], xatom[natom]);
            natom++;
        }
    }
    xatom.resize(natom);
    bool foundNothing = true;
    for (const auto &miq : iqm)
    {
        auto mpo = miq.first;
        switch (mpo)
        {
        case MolPropObservable::CHARGE:
        case MolPropObservable::POTENTIAL:
            {
                std::string conf;
                auto ei = findExperimentConst(JobType::OPT);
                if (!ei)
                {
                    ei = findExperimentConst(JobType::TOPOLOGY);
                }
                if (ei)
                {
                    for(auto &i : qTypes())
                    {
                        qType               qi = i.first;
                        std::string         reference;
                        std::string         lot;
                        std::vector<double> q;
                        if (ei->getCharges(&q, qi, &reference, &lot))
                        {
                            auto qp = qProps_.find(qi);
                            if (qp == qProps_.end())
                            {
                                qProps_.insert(std::pair<qType, QtypeProps>(qi, QtypeProps(qi)));
                                qp = qProps_.find(qi);
                                GMX_RELEASE_ASSERT(qp != qProps_.end(), "Could not insert a new QtypeProps in qProps_");
                            }
                            qp->second.setQ(q);
                            qp->second.setX(xatom);
                            qp->second.calcMoments();
                        }
                    }
                    foundNothing = false;
                }
            }
            break;
        case MolPropObservable::DELTAE0:
        case MolPropObservable::DHFORM:
        case MolPropObservable::DGFORM:
        case MolPropObservable::ZPE:
            {
                auto gp = static_cast<const MolecularEnergy *>(qmProperty(mpo, T, JobType::OPT));
                if (gp)
                {
                    energy_.insert(std::pair<MolPropObservable, double>(mpo, gp->getValue()));
                    foundNothing = false;
                }
            }
            break;
        case MolPropObservable::DIPOLE:
        case MolPropObservable::QUADRUPOLE:
        case MolPropObservable::OCTUPOLE:
        case MolPropObservable::HEXADECAPOLE:
            {
                auto gp = static_cast<const MolecularMultipole *>(qmProperty(mpo, T, JobType::OPT));
                if (gp)
                {
                    qProps_.find(qType::Elec)->second.setMultipole(mpo, gp->getVector());
                    foundNothing = false; 
                }
            }
            break;
        case MolPropObservable::POLARIZABILITY:
            {
                auto gp = static_cast<const MolecularPolarizability *>(qmProperty(mpo, T, JobType::OPT));
                if (gp)
                {
                    auto qelec = qTypeProps(qType::Elec);
                    qelec->setPolarizabilityTensor(gp->getTensor());
                    foundNothing = false;
                }
            }
            break;
        default:
            break;
        }
    }
    if (foundNothing)
    {
        imm = immStatus::NoData;
    }
    return imm;
}

void ACTMol::initQgenResp(const ForceField             *pd,
                          const std::vector<gmx::RVec> &coords,
                          real                          watoms,
                          int                           maxESP)
{
    std::string        mylot;
    auto &qt         = pd->findForcesConst(InteractionType::COULOMB);
    auto iChargeType = name2ChargeType(qt.optionValue("chargetype"));
    auto qp          = qTypeProps(qType::Calc);
    QgenResp *qgr = qp->qgenResp();
    qgr->setChargeType(iChargeType);
    qgr->setAtomInfo(atomsConst(), pd, coords, totalCharge());
    qp->setQ(atomsConst());
    qp->setX(coords);
    qgr->setAtomSymmetry(symmetric_charges_);
    qgr->summary(debug);

    int natoms = nRealAtoms();

    std::random_device               rd;
    std::mt19937                     gen(rd());  
    std::uniform_real_distribution<> uniform(0.0, 1.0);
    double                           cutoff = 0.01*maxESP;
 
    auto ci = findExperimentConst(JobType::OPT);
    if (ci)
    {
        int iesp = 0;
        for (auto &epi : ci->electrostaticPotentialConst())
        {
            auto val = uniform(gen);
            if ((iesp >= natoms || watoms > 0) && val <= cutoff)
            {
                auto xu = epi.getXYZunit();
                auto vu = epi.getVunit();
                qgr->addEspPoint(convertToGromacs(epi.getX(), xu),
                                 convertToGromacs(epi.getY(), xu),
                                 convertToGromacs(epi.getZ(), xu),
                                 convertToGromacs(epi.getV(), vu));
            }
            iesp++;
        }
        if (debug)
        {
            fprintf(debug, "%s added %zu out of %zu ESP points to the RESP structure.\n",
                    getMolname().c_str(), qgr->nEsp(),
                    ci->electrostaticPotentialConst().size());
        }
    }
}

QtypeProps *ACTMol::qTypeProps(qType qt)
{
    auto qp = qProps_.find(qt);
    if (qp != qProps_.end())
    {
        return &qp->second;
    }
    return nullptr;
}

const QtypeProps *ACTMol::qTypeProps(qType qt) const
{
    const auto qp = qProps_.find(qt);
    if (qp != qProps_.end())
    {
        return &qp->second;
    }
    return nullptr;
}

void ACTMol::plotEspCorrelation(const ForceField             *pd,
                                const std::vector<gmx::RVec> &coords,
                                const char                   *espcorr,
                                const gmx_output_env_t       *oenv,
                                const ForceComputer          *forceComp)
{
    if (espcorr && oenv)
    {
        auto qgr   = qTypeProps(qType::Calc)->qgenResp();
        qgr->updateAtomCharges(atomsConst());
        qgr->updateAtomCoords(coords);
        std::vector<gmx::RVec> forces(atomsConst().size());
        std::map<InteractionType, double> energies;
        (void) forceComp->compute(pd, topology_, &optimizedCoordinates_,
                                  &forces, &energies);
        qgr->calcPot(1.0);
        qgr->plotLsq(oenv, espcorr);
    }
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

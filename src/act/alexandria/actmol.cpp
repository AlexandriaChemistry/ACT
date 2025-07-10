/*
 * This source file is part of the Alexandria Chemistry Toolkit.
 *
 * Copyright (C) 2014-2025
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
#include <set>
#include <string>

#include "act/alexandria/actmol_low.h"
#include "act/alexandria/gromacs_top.h"
#include "act/alexandria/pdbwriter.h"
#include "act/alexandria/symmetrize_charges.h"
#include "act/basics/msg_handler.h"
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
        bool b = (ActParticle::Atom == aa.pType());
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
                GMX_RELEASE_ASSERT(QgenResp_.natoms() - j > 0, "Atom index out of range copying charges from RESP to Qprops");
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

void ACTMol::checkAtoms(MsgHandler       *msghandler,
                        const ForceField *pd)
{
    int atomnumberTotal = 0;
    auto myatoms = topology_->atoms();
    for (size_t i = 0; i < myatoms.size(); i++)
    {
        auto atype = myatoms[i].ffType();
        if (!pd->hasParticleType(atype))
        {
            msghandler->msg(ACTStatus::Error, ACTMessage::MissingFFParameter,
                            gmx::formatString("Could not find a force field entry for atomtype %s atom %zu in compound '%s'", atype.c_str(), i+1, getMolname().c_str()).c_str());
        }
        else
        {
            atomnumberTotal += pd->findParticleType(atype)->atomnumber();
        }
    }
    if (!msghandler->ok())
    {
        return;
    }
    // Check multiplicity
    int multOk = atomnumberTotal + totalMultiplicity() + totalCharge();
    if (multOk % 2 == 0)
    {
        msghandler->msg(ACTStatus::Warning, ACTMessage::Multiplicity,
                        gmx::formatString("AtomnumberTotal %d, totalMultiplicity %d, totalCharge %d for %s", atomnumberTotal, totalMultiplicity(), totalCharge(),
                                          getMolname().c_str()));
        return;
    }
}

static std::vector<gmx::RVec> experCoords(const std::vector<gmx::RVec> &xxx,
                                          const Topology               *topology)
{
    auto myatoms = topology->atoms();
    if (myatoms.empty())
    {
        // Called this function too early?
        return xxx;
    }
    gmx::RVec fzero = { 0, 0, 0 };
    std::vector<gmx::RVec> coords(myatoms.size(), fzero);
    int j = 0;
    for(size_t i = 0; i < myatoms.size(); i++)
    {
        if (myatoms[i].pType() == ActParticle::Atom)
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
        case ActParticle::Atom:
            // Do nothing
            break;
        case ActParticle::Shell:
        case ActParticle::Vsite:
            {
                auto cores = myatoms[i].cores();
                GMX_RELEASE_ASSERT(!cores.empty(), "Shell or vsite without core");
                copy_rvec(coords[cores[0]], coords[i]);
            }
            break;
        default: // throws
            {
                GMX_THROW(gmx::InternalError(gmx::formatString("Don't know how to handle %s particle type", actParticleToString(myatoms[i].pType()).c_str()).c_str()));
            }
        }
    }
    ForceComputer fcomp;
    fcomp.generateVsites(topology, &coords);
    return coords;
}

bool ACTMol::hasMolPropObservable(MolPropObservable mpo) const
{
    bool hasMPO = false;
    for(const auto &exper : experimentConst())
    {
        hasMPO = hasMPO || exper.hasProperty(mpo);
    }
    return hasMPO;
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
            GMX_THROW(gmx::InvalidInputError(gmx::formatString("No structure at all for '%s'", getMolname().c_str()).c_str()));
        }
    }

    return experCoords(exper->getCoordinates(), topology_);
}

void ACTMol::forceEnergyMaps(MsgHandler                                                        *msghandler,
                             const ForceField                                                  *pd,
                             const ForceComputer                                               *forceComp,
                             std::vector<std::vector<std::pair<double, double> > >             *forceMap,
                             std::vector<ACTEnergy>                                            *energyMap,
                             ACTEnergyMapVector                                                *interactionEnergyMap,
                             std::vector<std::pair<double, std::map<InteractionType, double>>> *energyComponentMap,
                             bool                                                               separateInductionCorrection) const
{
    auto       myatoms = topology_->atoms();
    forceMap->clear();
    energyMap->clear();
    interactionEnergyMap->clear();
    energyComponentMap->clear();
    gmx::RVec fzero = { 0, 0, 0 };
    std::map<MolPropObservable, InteractionType> interE = { 
        { MolPropObservable::INTERACTIONENERGY, InteractionType::EPOT           },
        { MolPropObservable::ELECTROSTATICS,    InteractionType::ELECTROSTATICS },
        { MolPropObservable::DISPERSION,        InteractionType::DISPERSION     },
        { MolPropObservable::EXCHANGE,          InteractionType::EXCHANGE       },
        { MolPropObservable::VDWCORRECTION,     InteractionType::VDWCORRECTION  },
        { MolPropObservable::INDUCTIONCORRECTION, InteractionType::INDUCTIONCORRECTION  },
        { MolPropObservable::INDUCTION,         InteractionType::INDUCTION      },
        { MolPropObservable::CHARGETRANSFER,    InteractionType::CHARGETRANSFER }
    };
    GMX_RELEASE_ASSERT(forceComp, "No force computer supplied");
    for (auto &exper : experimentConst())
    {
        // We compute either interaction energies or normal energies for one experiment
        auto coords  = experCoords(exper.getCoordinates(), topology_);
        bool doInter = false;
        for(auto &ie : interE)
        {
            if (exper.hasProperty(ie.first))
            {
                doInter = true;
            }
        }
        if (doInter)
        {
            // Complicated way of combined multiple terms into one observable
            std::map<InteractionType,  std::set<MolPropObservable>> combined = {
                { InteractionType::ALLELEC, { MolPropObservable::ELECTROSTATICS,
                                              MolPropObservable::INDUCTION,
                                              MolPropObservable::INDUCTIONCORRECTION } },
                { InteractionType::EXCHIND,  { MolPropObservable::EXCHANGE,
                                               MolPropObservable::INDUCTION,
                                               MolPropObservable::INDUCTIONCORRECTION } }
            };
            std::vector<gmx::RVec> interactionForces(myatoms.size(), fzero);
            std::vector<gmx::RVec> mycoords(myatoms.size(), fzero);
            
            std::map<InteractionType, double> einter;
            calculateInteractionEnergy(msghandler,
                                       pd, forceComp, &einter, &interactionForces, &coords,
                                       separateInductionCorrection);
            ACTEnergyMap aemap;
            for (auto &ener : einter)
            {
                std::set<MolPropObservable> mpos;
                // Check whether one is enough or we have to look up multiple
                auto combine = combined.find(ener.first);
                if (combine != combined.end())
                {
                    // Combine multiple QM observable
                    for(const auto &mpo : combine->second)
                    {
                        mpos.insert(mpo);
                    }
                }
                else
                {
                    // Find the correct QM observable
                    for(const auto &ie : interE)
                    {
                        if (ie.second == ener.first)
                        {
                            mpos.insert(ie.first);
                            break;
                        }
                    }
                }
                // Loop over sets of observable to combine QM data
                double my_qm    = 0;
                size_t foundQM  = 0;
                for(const auto &mpo : mpos)
                {
                    if (exper.hasProperty(mpo))
                    {
                        auto propvec = exper.propertyConst(mpo);
                        if (propvec.size() != 1)
                        {
                            msghandler->fatal(gmx::formatString("Expected just one QM value per experiment for %s iso %zu", mpo_name(mpo), propvec.size()).c_str());
                        }
                        my_qm += propvec[0]->getValue();
                        foundQM    += 1;
                    }
                }
                // TODO Store the interaction forces
                ACTEnergy my_all(exper.id());
                my_all.setACT(ener.second);
                if (foundQM == mpos.size() ||
                    (foundQM == mpos.size()-1 && !separateInductionCorrection))
                {
                    my_all.setQM(my_qm);
                }
                aemap.insert({ener.first, my_all});
            }
            if (aemap.size() > 0)
            {
                interactionEnergyMap->push_back(aemap);
            }
        }
        else if (exper.hasProperty(MolPropObservable::DELTAE0))
        {
            std::map<InteractionType, double> energies;
            std::vector<gmx::RVec> forces(myatoms.size(), fzero);
            (void) forceComp->compute(pd, topology_, &coords, &forces, &energies);
            auto eprops = exper.propertyConst(MolPropObservable::DELTAE0);
            if (eprops.size() > 1)
            {
                msghandler->fatal("Multiple energies for this experiment");
            }
            else if (eprops.size() == 1)
            {
                energyMap->push_back(ACTEnergy(exper.id(), eprops[0]->getValue(),
                                               energies[InteractionType::EPOT]));
                energyComponentMap->push_back({ eprops[0]->getValue(), std::move(energies) });
            }
        
            const std::vector<gmx::RVec> &fff = exper.getForces();
            if (!fff.empty())
            {
                size_t ifff = 0;
                std::vector<std::pair<double, double> > thisForce;
                for (size_t i = 0; i < myatoms.size(); i++)
                {
                    if (myatoms[i].pType() == ActParticle::Atom)
                    {
                        if (ifff >= fff.size())
                        {
                            msghandler->fatal(gmx::formatString("Inconsistency: there are %lu atoms and shells, but only %zu forces", myatoms.size(), fff.size()));
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
        if (myatoms[i].pType() == ActParticle::Atom)
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

void ACTMol::GenerateTopology(MsgHandler        *msghandler,
                              ForceField        *pd,
                              missingParameters  missing)
{
    std::string btype1, btype2;

    msghandler->msg(ACTStatus::Debug,
                    gmx::formatString("Generating topology for %s",
                                      getMolname().c_str()));

    generateComposition();
    if (NAtom() <= 0)
    {
        msghandler->msg(ACTStatus::Warning, ACTMessage::AtomTypes, "No atoms");
    }
    /* Store bonds in harmonic potential list first, update type later */
    if (msghandler->ok())
    {
        topology_ = new Topology(*bonds());
    }
    // Get the coordinates.
    std::vector<gmx::RVec> coords;// = xOriginal();

    if (msghandler->ok())
    {
        topology_->GenerateAtoms(msghandler, pd, this,  &coords);
    }
    if (msghandler->ok())
    {
        checkAtoms(msghandler, pd);
    }
    // Create fragments before adding shells!
    if (msghandler->ok())
    {
        auto fptr = fragmentPtr();
        if (fptr->size() > 0)
        {
            fraghandler_ = new FragmentHandler(msghandler,
                                               pd, coords, topology_->atoms(),
                                               bondsConst(), fptr, missing);

            if (fraghandler_->topologies().empty())
            {
                msghandler->msg(ACTStatus::Error, ACTMessage::FragmentHandler,
                                getMolname());
                delete fraghandler_;
            }
            else
            {
                // Finally, extract frequencies etc.
                getHarmonics();
            }
        }
        else
        {
            msghandler->msg(ACTStatus::Error, ACTMessage::FragmentHandler,
                            "No fragments");
        }
    }
    
    if (msghandler->ok())
    {
        topology_->build(msghandler, pd, &coords, 175.0, 5.0, missing);
    }
    if (msghandler->ok())
    {
        auto myatoms = topology()->atoms();
        realAtoms_.clear();
        for(size_t i = 0; i < myatoms.size(); i++)
        {
            if (myatoms[i].pType() == ActParticle::Atom)
            {
                realAtoms_.push_back(i);
            }
        }
    }
    if (msghandler->ok() && missing != missingParameters::Generate)
    {
        std::vector<InteractionType> itUpdate;
        for(auto &entry : *(topology_->entries()))
        {
            itUpdate.push_back(entry.first);
        }
    }
    if (msghandler->ok())
    {
        checkAtoms(msghandler, pd);
    }
    if (msghandler->ok())
    {
        isLinear_ = isLinearMolecule(topology_->atoms(), coords);
        // Symmetrize the atoms
        get_symmetrized_charges(topology_, pd, nullptr, &symmetric_charges_);
    }
    else
    {
        msghandler->msg(ACTStatus::Error, gmx::formatString("Cannot make topology for %s", getMolname().c_str()));
    }
}

double ACTMol::bondOrder(int ai, int aj) const
{
    return topology_->findBond(ai, aj).bondOrder();
}

static void printEmap(MsgHandler *msghandler,
                      const std::map<InteractionType, double> *e)
{
    if (ACTStatus::Debug != msghandler->printLevel())
    {
        return;
    }
    std::string str;
    for(auto ee = e->begin(); ee != e->end(); ee++)
    {
        str += " ";
        str += interactionTypeToString(ee->first);
        str += " ";
        str += gmx_ftoa(ee->second);
    }
    msghandler->writeDebug(str);
}

static void checkEnergies(MsgHandler                              *msghandler,
                          const char                              *info,
                          const std::map<InteractionType, double> &einter)
{
    // Do this kind of check in debug mode only.
    if (ACTStatus::Debug != msghandler->printLevel())
    {
        return;
    }
    std::set<InteractionType> ignore = { InteractionType::ALLELEC, InteractionType::EXCHIND, InteractionType::EPOT };
    double esum = 0;
    for(const auto &e : einter)
    {
        if (ignore.end() == ignore.find(e.first))
        {
            esum += e.second;
        }
    }
    double toler = 1e-3;
    double epot  = einter.find(InteractionType::EPOT)->second;
    double diff  = std::abs(esum - epot);
    double denom = std::abs(esum + epot);
    // Try relative tolerance if the denominator not small,
    // otherwise try absolute tolerance.
    if ((denom > toler && (diff >= toler*denom)) ||
        (denom <= toler && diff >= toler))
    {
        for(const auto &e : einter)
        {
            msghandler->msg(ACTStatus::Info,
                            gmx::formatString("InteractionType %s energy %g\n",
                                              interactionTypeToString(e.first).c_str(), e.second));
        }
        msghandler->fatal(gmx::formatString("%s: found einter %g but sum of terms is %g",
                                            info, epot, esum).c_str());
    }
}

void ACTMol::calculateInteractionEnergy(MsgHandler                        *msghandler,
                                        const ForceField                  *pd,
                                        const ForceComputer               *forceComputer,
                                        std::map<InteractionType, double> *einter,
                                        std::vector<gmx::RVec>            *interactionForces,
                                        std::vector<gmx::RVec>            *coords,
                                        bool                               separateInductionCorrection) const
{
    auto &tops = fraghandler_->topologies();
    einter->clear();
    if (tops.size() != 2)
    {
        return;
    }
    // First, compute the total energy
    gmx::RVec fzero = { 0, 0, 0 };
    interactionForces->resize(coords->size(), fzero);
    msghandler->writeDebug("Will compute interaction energy\n");

    // First compute the relaxed monomer energies
    std::map<InteractionType, double> e_monomer[2];
    auto &astart = fraghandler_->atomStart();
    auto itInduc = InteractionType::INDUCTION;
    auto itICorr = InteractionType::INDUCTIONCORRECTION;
    auto itElec  = InteractionType::ELECTROSTATICS;
    auto itPolar = InteractionType::POLARIZATION;
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
        (void) forceComputer->compute(pd, tops[ff], &myx, &forces, &e_monomer[ff]);
        j = 0;
        for (size_t i = astart[ff]; i < astart[ff]+natom; i++)
        {
            // Copy back coordinates to original array to store relaxed shell positions
            copy_rvec(myx[j], (*coords)[i]);
            // Store interaction forces
            for(int m = 0; m < DIM; m++)
            {
                (*interactionForces)[i][m] -= forces[j][m];
            }
            j++;
        }
        msghandler->writeDebug(gmx::formatString("%s monomer[%zu]:", getMolname().c_str(), ff));
        printEmap(msghandler, &e_monomer[ff]);
        
        // If present, move induction in monomer to electrostatics
        auto ei = e_monomer[ff].find(itInduc);
        if (e_monomer[ff].end() != ei)
        {
            e_monomer[ff][itElec] += ei->second;
            ei->second = 0;
            msghandler->writeDebug(gmx::formatString("%s 'corrected' monomer[%zu]:", getMolname().c_str(), ff));
            printEmap(msghandler, &e_monomer[ff]);
        }
        checkEnergies(msghandler, "Monomer", e_monomer[ff]);
    }
    // TODO forces are overwritten, not trustable!
    std::vector<gmx::RVec> forces(topology_->atoms().size(), fzero);
    std::map<InteractionType, double> e_total;
    auto myatoms = topology_->atoms();
    // Now compute the total relaxed energy
    {
        // Make a copy of the coordinates, since shell positions will change
        std::vector<gmx::RVec> mycoords = *coords;
        // Note that here the shells start in the position of the relaxed
        // monomers such as to only measure the induction due to the dimer
        (void) forceComputer->compute(pd, topology_, &mycoords, &forces, &e_total, fzero, false);
        // Move remaining polarisation energy to the electrostatics.
        // Since we started with the shells in the relaxed monomer position,
        // not the entire polarization will have been moved to induction.
        auto eip = e_total.find(itPolar);
        auto eie = e_total.find(itElec);
        if (e_total.end() != eip &&
            e_total.end() != eie)
        {
            eie->second += eip->second;
            // Erase polarization term!
            e_total.erase(eip);
        }
        checkEnergies(msghandler, "Total", e_total);
    }
    // Now time to compute mutual induction for the dimer.
    // This means, that the shells of one molecule are allowed
    // to relax, but not the other.
    std::map<InteractionType, double> e_dimer[2];
    for(size_t ff = 0; ff < tops.size(); ff++)
    { 
        // Relax molecule ff only
        std::set<int> myshells;
        int natom = tops[ff]->atoms().size();
        for(size_t i = astart[ff]; i < astart[ff]+natom; i++)
        {
            if (myatoms[i].pType() == ActParticle::Shell)
            {
                myshells.insert(i);
            }
        }
        // Make a fresh copy of the coordinates
        std::vector<gmx::RVec> mycoords = *coords;
        // Eqn 6
        (void) forceComputer->compute(pd, topology_, &mycoords, &forces, &e_dimer[ff],
                                      fzero, false, myshells);
        msghandler->writeDebug(gmx::formatString("%s dimer[%lu]:", getMolname().c_str(), ff));
        printEmap(msghandler, &e_dimer[ff]);

        checkEnergies(msghandler, "Dimer", e_dimer[ff]);
    }
    // Now compute interaction energy terms
    *einter = e_total;

    // First subtract everything from the monomers.
    // This gives us correct electrostatics to get Eqn. 8.2 in the manual
    for(size_t ff = 0; ff < tops.size(); ff++)
    {
        for(const auto &ee : e_monomer[ff])
        {
            einter->find(ee.first)->second -= ee.second;
        }
    }
    msghandler->writeDebug(gmx::formatString("%s total:", getMolname().c_str()));
    printEmap(msghandler, &e_total);
    msghandler->writeDebug(gmx::formatString("%s einter:", getMolname().c_str()));
    printEmap(msghandler, einter);

    checkEnergies(msghandler, "Inter 0", *einter);
    // Check, to see whether induction needs to be summed or split
    if (separateInductionCorrection)
    {
        // Compute the induction energy to the second order, Eqn. 8.3 in the manual
        double delta_einduc2 = 0;
        
        for(size_t ff = 0; ff < tops.size(); ff++)
        {
            for(const auto itype : {itInduc})
            {
                auto ee = e_dimer[ff].find(itype);
                if (e_dimer[ff].end() != ee)
                {
                    delta_einduc2 += ee->second;
                    ee->second = 0;
                }
            }
        }
        // Finally, compute the remaining (garbage bin) terms, Eqn. 8.4 in the manual
        // by moving induction energy from the induction term to the
        // induction correction term.
        auto eit = einter->find(itInduc);
        if (einter->end() == eit)
        {
            einter->insert({ itInduc, 0 });
            eit = einter->find(itInduc);
        }

        // Compute higher order induction
        double induc3 = eit->second -  delta_einduc2;
        // Store higher order induction in einter
        einter->insert_or_assign(itICorr, induc3);
        // Correct total induction for this
        eit->second -= induc3;
    }
    else
    {
        // Sum Induction Correction into total Induction
        auto eic = einter->find(itICorr);
        if (einter->end() != eic)
        {
            einter->find(itInduc)->second += eic->second;
            // Erase induction correction term.
            einter->erase(eic);
        }
    }
    checkEnergies(msghandler, "Inter 1", *einter);
    {
        // Gather the terms for ALLELEC output.
        auto eall = einter->find(InteractionType::ALLELEC);
        eall->second = 0;
        if (einter->end() != eall)
        {
            for(auto itype : { itElec, itInduc, itICorr })
            {
                auto ee = einter->find(itype);
                if (einter->end() != ee)
                {
                    eall->second += ee->second;
                }
            }
        }
    }
    checkEnergies(msghandler, "Inter 2", *einter);
    {
        // Gather the terms for EXCHIND output.
        auto eall = einter->find(InteractionType::EXCHIND);
        eall->second = 0;
        if (einter->end() != eall)
        {
            for(auto itype : { InteractionType::EXCHANGE, itInduc, itICorr })
            {
                auto ee = einter->find(itype);
                if (einter->end() != ee)
                {
                    eall->second += ee->second;
                }
            }
        }
    }
    checkEnergies(msghandler, "Inter 3", *einter);
    // Add dimer forces to the interaction forces
    // TODO this is not correct anymore, see above.
    for(size_t i = 0; i < topology_->atoms().size(); i++)
    {
        for(int m = 0; m< DIM; m++)
        {
            (*interactionForces)[i][m] += forces[i][m];
        }
    }
    msghandler->writeDebug(gmx::formatString("%s result:", getMolname().c_str()));
    printEmap(msghandler, einter);
}

ACTMessage ACTMol::GenerateAcmCharges(const ForceField       *pd,
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
    ACTMessage imm       = ACTMessage::OK;
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
                                                   *coords, pd, atoms(),
                                                   symmetric_charges_);
        if (eQgen::OK == eqgen)
        {
            (void) forceComp->compute(pd, topology_, coords, forces, &energies);
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
        }
        else if (eQgen::NOSUPPORT == eqgen)
        {
            imm = ACTMessage::MissingChargeGenerationParameters;
        }
        else
        {
            imm = ACTMessage::ChargeGeneration;
        }
        iter++;
    }
    while (imm == ACTMessage::OK && (!converged) && (iter < maxQiter_));
    if (ACTMessage::OK == imm)
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

void ACTMol::updateQprops(const ForceField          *pd,
                          const ForceComputer       *forceComp,
                          std::vector<gmx::RVec>    *forces)
{
    std::map<InteractionType, double> energies;
    auto myatoms = atoms();
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

void ACTMol::GenerateCharges(MsgHandler                *msghandler,
                             const ForceField          *pd,
                             const ForceComputer       *forceComp,
                             ChargeGenerationAlgorithm  algorithm,
                             qType                      qtype,
                             const std::vector<double> &qcustom,
                             std::vector<gmx::RVec>    *coords,
                             std::vector<gmx::RVec>    *forces,
                             bool                       updateQProps)
{
    bool converged   = false;

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
            msghandler->msg(ACTStatus::Warning,
                            gmx::formatString("WARNING! Using fixed charges for %s!\n", getMolname().c_str()));
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
            fraghandler_->setCharges(*myatoms);
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
                msghandler->msg(ACTStatus::Error, ACTMessage::NoMolpropCharges,
                                getMolname());
                return;
            }
            size_t j = 0;
            for(size_t i = 0; i < myatoms->size(); i++)
            {
                if ((*myatoms)[i].pType() == ActParticle::Atom)
                {
                    (*myatoms)[i].setCharge(qread[j++]);
                }
                else if ((*myatoms)[i].pType() == ActParticle::Shell)
                {
                    const auto &pId = (*myatoms)[i].ffType();
                    if (!pd->hasParticleType(pId))
                    {
                        msghandler->msg(ACTStatus::Error, ACTMessage::AtomTypes,
                                        getMolname());
                        return;
                    }
                    auto piter = pd->findParticleType(pId);
                    auto q = piter->charge();
                    (*myatoms)[i].setCharge(q);
                    // TODO Do not use -1 here, but particle.core()
                    (*myatoms)[i-1].setCharge(qread[j-1]-q);
                }
            }
            fraghandler_->setCharges(*myatoms);
        }
        break;
    case ChargeGenerationAlgorithm::Custom:
        {
            for (size_t i = 0; i < myatoms->size(); i++)
            {
                (*myatoms)[i].setCharge(qcustom[i]);
            }
            fraghandler_->setCharges(*myatoms);
            return;
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
                if (!ffOption(*pd, InteractionType::ELECTROSTATICS, "epsilonr", &epsilonr))
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
                    msghandler->msg(ACTStatus::Error, ACTMessage::ChargeGeneration,
                                    getMolname());
                    return;
                }
            }
        }
        break;
    case ChargeGenerationAlgorithm::EEM:
    case ChargeGenerationAlgorithm::SQE:
        {
            (void) GenerateAcmCharges(pd, forceComp, coords, forces);
        }
        break;
    }
    if (updateQProps)
    {
        updateQprops(pd, forceComp, forces);
    }
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
                       const std::string        &title,
                       const char               *unit,
                       const std::vector<double> &Q)
{
    double fac = convertFromGromacs(1.0, unit);
    commercials->push_back("");
    commercials->push_back(title);
    commercials->push_back(gmx::formatString("  ( %6.2f %6.2f %6.2f )", fac*Q[0], fac*Q[1], fac*Q[2]));
    commercials->push_back(gmx::formatString("  (        %6.2f %6.2f )", fac*Q[3], fac*Q[4]));
    commercials->push_back(gmx::formatString("  (               %6.2f )", fac*Q[5]));
}

std::string ACTMol::levelOfTheory()
{
    std::string s;
    if (!method_.empty() && !basis_.empty())
    {
        s = method_ + "/" + basis_;
    }
    return s;
}

std::vector<std::string> ACTMol::generateCommercials(const ForceField             *pd,
                                                     const ForceComputer          *forceComp,
                                                     const std::vector<gmx::RVec> &coords)
{
    double                   T = -1;
    std::string              myref;
    std::string              mylot       = levelOfTheory();
    std::vector<std::string> commercials;

    commercials.push_back("");
    commercials.push_back(gmx::formatString("Molecule information for %s", getMolname().c_str()));
    commercials.push_back("-----------------------------------------------");
    commercials.push_back(gmx::formatString("Formula      = %s", formula().c_str()));
    commercials.push_back(gmx::formatString("Total Mass   = %.3f (Da)", totalMass()));
    commercials.push_back(gmx::formatString("Total Charge = %d (e)", totalCharge()));
    commercials.push_back(gmx::formatString("Multiplicity = %d", totalMultiplicity()));

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
            if (!gp)
            {
                gp = qmProperty(mpo, T, JobType::TOPOLOGY);
            }
            if (gp)
            {
                if (qelec.hasMultipole(mpo))
                {
                    auto mymu = qelec.getMultipole(mpo);
                    commercials.push_back("");
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
                    commercials.push_back("");
                    commercials.push_back(gmx::formatString("Alexandria %s (%s)\n", mpo_name(mpo), gp->getUnit()));
                    for(auto &fmp : formatMultipole(mpo, mymu))
                    {
                        commercials.push_back(fmp);
                    }
                }
            }
        }
    
        if (pd->polarizable())
        {
            qcalc->calcPolarizability(pd, topology(), forceComp);
            auto acalc = qcalc->polarizabilityTensor();
            std::vector<double> ac = { acalc[XX][XX], acalc[XX][YY], acalc[XX][ZZ],
                                       acalc[YY][YY], acalc[YY][ZZ], acalc[ZZ][ZZ] };
            auto unit = mpo_unit2(MolPropObservable::POLARIZABILITY);
            add_tensor(&commercials, "Alexandria Polarizability components (A^3)", unit, ac);
            double fac = convertFromGromacs(1.0, unit);

            commercials.push_back(gmx::formatString("Alexandria Isotropic Polarizability: %.2f (A^3)\n",
                                                    fac*qcalc->isotropicPolarizability()));
            commercials.push_back(gmx::formatString("Alexandria Anisotropic Polarizability: %.2f (A^3)\n",
                                                    fac*qcalc->anisotropicPolarizability()));
        
            T = -1;
            if (qelec.hasPolarizability())
            {
                auto aelec = qelec.polarizabilityTensor();
                std::vector<double> ae = { aelec[XX][XX], aelec[XX][YY], aelec[XX][ZZ],
                                           aelec[YY][YY], aelec[YY][ZZ], aelec[ZZ][ZZ] };

                add_tensor(&commercials, gmx::formatString("%s + Polarizability components (A^3)", mylot.c_str()),
                           unit, ae);
                commercials.push_back(gmx::formatString("%s Isotropic Polarizability: %.2f (A^3)\n",
                                                        mylot.c_str(), fac*qelec.isotropicPolarizability()));
                commercials.push_back(gmx::formatString("%s Anisotropic Polarizability: %.2f (A^3)\n",
                                                        mylot.c_str(), fac*qelec.anisotropicPolarizability()));
            }
        }
    }
    commercials.push_back("-----------------------------------------------");
    return commercials;
}

void ACTMol::PrintTopology(MsgHandler                   *msg_handler,
                           const char                   *fn,
                           const ForceField             *pd,
                           const ForceComputer          *forceComp,
                           const std::vector<gmx::RVec> &coords,
                           bool                          bITP)
{
    t_mols  printmol;
    FILE   *fp = gmx_ffopen(fn, "w");

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
    auto commercials = generateCommercials(pd, forceComp, coords);

    // TODO write a replacement for this function
    print_top_header(fp, pd, bHaveShells_, commercials, bITP);
    write_top(msg_handler, fp, printmol.name, topology_, pd);
    if (!bITP)
    {
        print_top_mols(fp, printmol.name, 1, &printmol);
    }
    if (msg_handler->verbose())
    {
        for (auto &entry : *(topology_->entries()))
        {
            auto &fs = pd->findForcesConst(entry.first);
            auto ftype = fs.potential();
            if (entry.second.size() > 0)
            {
                printf("There are %4zu %s interactions\n", entry.second.size(),
                       potentialToString(ftype).c_str());
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
    auto &qt          = pd->findForcesConst(InteractionType::ELECTROSTATICS);
    auto  iChargeType = potentialToChargeType(qt.potential());
    
    if (potfn || hisfn || rhofn || difffn || pdbdifffn)
    {
        char   *act_version = (char *)"act v0.99b";
        double  epsilonr;
        if (!ffOption(*pd, InteractionType::ELECTROSTATICS, 
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

void ACTMol::getExpProps(MsgHandler                                 *msghandler,
                         const ForceField                           *pd,
                         const std::map<MolPropObservable, iqmType> &iqm,
                         real                                        watoms,
                         int                                         maxESP)
{
    std::vector<double> vec;
    std::string         myref;
    std::string         mylot;
    
    auto &myatoms = atomsConst();
    if (myatoms.size() == 0)
    {
        msghandler->msg(ACTStatus::Error,
                        gmx::formatString("No atoms to fetch experimental properties for %s.",
                                          getMolname().c_str()));
        return;
    }
    bool foundNothing = true;
    msghandler->msg(ACTStatus::Debug,
                    gmx::formatString("Found %zu experiments for %s",
                                      experimentConst().size(),
                                      getMolname().c_str()));
    for (const auto &myexp : experimentConst())
    {
        bool     qprop = false;
        ACTQprop actq;
        auto     qelec = actq.qPqm();
        auto     qcalc = actq.qPact();
        auto     props = myexp.propertiesConst();
        auto     xatom = experCoords(myexp.getCoordinates(), topology_);
        setBasisset(myexp.getBasisset());
        setMethod(myexp.getMethod());
        std::vector<double> q;
        for (auto prop : props)
        {
            // Check whether this property is in the "Wanted" list
            if (!iqm.empty())
            {
                auto iter = iqm.find(prop.first);
                if (iqm.end() == iter)
                {
                    continue;
                }
            }
            // Fetch the property from this experiment
            auto gp = myexp.propertyConst(prop.first);
            msghandler->msg(ACTStatus::Debug,
                            gmx::formatString("Found %zu %s values for %s",
                                              gp.size(),
                                              mpo_name(prop.first),
                                              getMolname().c_str()));
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
                    msghandler->msg(ACTStatus::Debug,
                                    gmx::formatString("Found %zu charges for %s",
                                                      q.size(), getMolname().c_str()));
                    break;
                }
            case MolPropObservable::POTENTIAL:
                {
                    auto qgr = actq.qgenResp();
                    auto &qt = pd->findForcesConst(InteractionType::ELECTROSTATICS);
                    qgr->setChargeType(potentialToChargeType(qt.potential()));
                    qgr->setAtomInfo(atomsConst(), pd, totalCharge());
                    qgr->updateAtomCoords(xatom);
                    qgr->setAtomSymmetry(symmetric_charges_);
                    qgr->summary(debug);
                    size_t natoms = nRealAtoms();
                        
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
                            if ((ll >= natoms || watoms > 0) && val <= cutoff)
                            {
                                qgr->addEspPoint(xyz[ll][XX], xyz[ll][YY], xyz[ll][ZZ], V[ll]);
                            }
                        }
                        msghandler->msg(ACTStatus::Debug,
                                        gmx::formatString("Found %zu potential points for %s",
                                                          V.size(), getMolname().c_str()));
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
            case MolPropObservable::ELECTROSTATICS:
            case MolPropObservable::INDUCTION:
            case MolPropObservable::EXCHANGE:
            case MolPropObservable::VDWCORRECTION:
            case MolPropObservable::INDUCTIONCORRECTION:
            case MolPropObservable::DISPERSION:
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
            case MolPropObservable::FREQUENCY:
            case MolPropObservable::INTENSITY:
            case MolPropObservable::HF:
            case MolPropObservable::CHARGETRANSFER:
            case MolPropObservable::DSFORM:
            case MolPropObservable::ENTROPY:
            case MolPropObservable::STRANS:
            case MolPropObservable::SROT:
            case MolPropObservable::SVIB:
            case MolPropObservable::CP:
            case MolPropObservable::CV:
                break;
            case MolPropObservable::COORDINATES:
                foundNothing = false;
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
        msghandler->msg(ACTStatus::Warning, ACTMessage::NoData, 
                        gmx::formatString("Did not find any data for %s", getMolname().c_str()));
    }
}

CommunicationStatus ACTMol::Send(const CommunicationRecord *cr, int dest) const
{
    auto cs = MolProp::Send(cr, dest);
    if (CommunicationStatus::OK == cs)
    {
        cr->send(dest, dataset_type_);
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
        iMolSelect ims;
        cr->recv(src, &ims);
        set_datasetType(ims);
    }
    return cs;
}

} // namespace alexandria

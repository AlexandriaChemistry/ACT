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
 * \author Julian Ramon Marrades Furquet <julian@marrad.es>
 * \author Oskar Tegby <oskar.tegby@it.uu.se>
 */

#include "acmfitnesscomputer.h"

#include "act/basics/dataset.h"
#include "act/basics/msg_handler.h"
#include "act/ga/genome.h"
#include "act/utility/communicationrecord.h"

namespace alexandria
{

/* * * * * * * * * * * * * * * * * * * *
* BEGIN: ACMFitnessComputer            *
* * * * * * * * * * * * * * * * * * * */

void ACMFitnessComputer::compute(MsgHandler *msghandler,
                                 ga::Genome *genome,
                                 iMolSelect  trgtFit)
{
    if (nullptr == genome)
    {
        GMX_THROW(gmx::InternalError("Empty genome"));
    }
    // First send around parameters. This code is run on master only.
    distributeTasks(CalcDev::Parameters);
    std::set<int> changed;
    distributeParameters(genome->basesPtr(), changed);
    // Then do the computation
    auto cd = distributeTasks(CalcDev::Compute);
    double fitness = calcDeviation(msghandler, cd, trgtFit);
    genome->setFitness(trgtFit, fitness);
    if (msghandler->tw())
    {
        const auto &targets = sii_->targets();
        if (targets.empty())
        {
            GMX_THROW(gmx::InternalError("No targets whatsoever!"));
        }
        auto ttt = targets.find(trgtFit);
        if (targets.end() != ttt)
        {
            auto fts = ttt->second;
            if (!fts.empty())
            {
                msghandler->tw()->writeLineFormatted("Components of training function for %s set\n", iMolSelectName(trgtFit));
                for (const auto &ft : fts)
                {
                    auto p = ft.second.info();
                    if (!p.empty())
                    {
                        msghandler->tw()->writeLine(p);
                    }
                }
            }
        }
    }
}

CalcDev ACMFitnessComputer::distributeTasks(CalcDev task)
{
    auto cr = sii_->commRec();
    // Send / receive parameters
    if (cr->isHelper())
    {
        // Now who is my middleman?
        int cd;
        cr->recv(cr->superior(), &cd);
        return static_cast<CalcDev>(cd);
    }
    else 
    {
        // Send only to my helpers
        for (auto &dest : cr->helpers())
        {
            cr->send(dest, static_cast<int>(task));
        }
        return task;
    }
}

void ACMFitnessComputer::distributeParameters(const std::vector<double> *params,
                                              const std::set<int>       &changed)
{
    if (debug)
    {
        fprintf(debug, "Starting to distribute parameters\n");
    }
    auto cr = sii_->commRec();
    // Send / receive parameters
    if (cr->isHelper())
    {
        // TODO: Implement broadcast
        // Find out who to talk to
        int src = cr->superior();
        std::vector<double> myparams;
        cr->recv(src, &myparams);
        std::set<int>       mychanged;
        int nchanged;
        cr->recv(src, &nchanged);
        for(int i = 0; i < nchanged; i++)
        {
            int mc;
            cr->recv(src, &mc);
            mychanged.insert(mc);
        }
        sii_->updateForceField(mychanged, myparams);
    }
    else 
    {
        // Send only to my helpers
        for (auto &dest : cr->helpers())
        {
            cr->send(dest, *params);
            cr->send(dest, static_cast<int>(changed.size()));
            for(int iset : changed)
            {
                cr->send(dest, iset);
            }
        }
        sii_->updateForceField(changed, *params);
    }
    if (debug)
    {
        fprintf(debug, "Finished distributing parameters\n");
    }
}

double ACMFitnessComputer::calcDeviation(MsgHandler *msghandler,
                                         CalcDev     task,
                                         iMolSelect  ims)
{
    msghandler->writeDebug("CalcDev starting");

    auto cr = sii_->commRec();
    // Send / receive molselect group.ks
    if (cr->isHelper())
    {
        // Now who is my middleman?
        cr->recv(cr->superior(), &ims);
    }
    else if (CalcDev::ComputeAll != task)
    {
        // Send ims to my helpers
        for (auto &dest : cr->helpers())
        {
            cr->send(dest, ims);
        }
    }
    msghandler->writeDebug(gmx::formatString("CalcDev Going to compute dataset %s\n",
                                             iMolSelectName(ims)));

    // Gather fitting targets
    std::map<eRMS, FittingTarget> *targets = sii_->fittingTargets(ims);
    if (nullptr == targets)
    {
        msghandler->writeDebug(gmx::formatString("Cannot find targets for %s\n", iMolSelectName(ims)));
        return 0;
    }
    // Reset the chi2 in FittingTargets for the given dataset in ims
    sii_->resetChiSquared(ims);

    // If actMaster or actMiddleMan, penalize out of bounds
    if (cr->isMasterOrMiddleMan() && bdc_)
    {
        bdc_->calcDeviation(msghandler, forceComp_, nullptr, nullptr, targets, sii_->forcefield());
    }

    // Loop over molecules
    int ntrain = 0;
    int nlocal = 0;
    auto mymols = molgen_->actmolsPtr();
    for (auto actmol = mymols->begin(); actmol < mymols->end(); ++actmol)
    {
        msghandler->writeDebug(gmx::formatString("CalcDev: mol %s dataset %s\n",
                                                 actmol->getMolname().c_str(),
                                                 iMolSelectName(actmol->datasetType())));
        if (ims != actmol->datasetType())
        {
            continue;
        }
        ntrain++;
        if ((actmol->support() == eSupport::Local) ||
            (task == CalcDev::ComputeAll && actmol->support() == eSupport::Remote))
        {
            std::vector<InteractionType> itUpdate;
            for(auto &io : molgen_->iopt())
            {
                if (io.second)
                {
                    itUpdate.push_back(io.first);
                }
            }
            // Now update the topology
            actmol->topologyPtr()->fillParameters(msghandler, sii_->forcefield(), missingParameters::Error);
            // Fill the fragments too if there are any
            for(auto &ft : actmol->fragmentHandler()->topologiesPtr())
            {
                ft->fillParameters(msghandler, sii_->forcefield(), missingParameters::Error);
            }

            // Run charge generation including shell minimization
            std::vector<gmx::RVec> forces(actmol->atomsConst().size(), { 0, 0, 0 });
            std::vector<gmx::RVec> coords = actmol->xOriginal();
            ACTMessage imm = actmol->GenerateAcmCharges(sii_->forcefield(), forceComp_, &coords, &forces);

            // Check whether we have to disable this compound
            if (ACTMessage::OK != imm && removeMol_)
            {
                actmol->setSupport(eSupport::No);
                continue;
            }

            computeMultipoles(targets, &(*actmol));

            if (devComputers_.size() == 0)
            {
                printf("No devComputers\n");
            }
            for (DevComputer *mydev : devComputers_)
            {
                mydev->calcDeviation(msghandler, forceComp_, &(*actmol), &coords, targets, sii_->forcefield());
            }
            msghandler->writeDebug(gmx::formatString("CalcDev: rank %d mol %s #energies %zu tw %g\n",
                                                     cr->rank(), actmol->getMolname().c_str(), actmol->experimentConst().size(),
                                                     targets->find(eRMS::EPOT)->second.totalWeight()));
            nlocal++;
        }
    }
    // Sum the terms of the chi-squared once we have done calculations
    // for all the molecules.
    sii_->sumChiSquared(task == CalcDev::Compute, ims);
    auto erms = eRMS::TOT;
    auto etot = targets->find(erms);
    if (targets->end() == etot)
    {
        msghandler->msg(ACTStatus::Error, 
                        gmx::formatString("Cannot find %s in targets.",
                                          rmsName(erms)).c_str());
    }
    numberCalcDevCalled_ += 1;
    if (etot->second.chiSquared() == 0 && ntrain > 0)
    {
        std::string msg = gmx::formatString("Zero %s chi squared for %s - this cannot be correct.\n", 
                                            iMolSelectName(ims), rmsName(erms));
        msg += gmx::formatString("There are %d compounds. Task = %s. Nlocal = %d, ntrain = %d.\n",
                                 ntrain, calcDevName(task), nlocal, ntrain);
        msg += "devComputers: ";
        for(auto &d : devComputers_)
        {
            msg += ' ' + d->name();
        }
        msg += "\n";
        for(const auto &ttt: *targets)
        {
            if (ttt.second.weight() > 0 || ttt.second.totalWeight() > 0 )
            {
                msg += gmx::formatString("Weight for %s %g, totalWeight %g chiSquared %g\n", rmsName(ttt.second.erms()),
                                         ttt.second.weight(), ttt.second.totalWeight(),
                                         ttt.second.chiSquared());
            }
        }
        msghandler->msg(ACTStatus::Fatal, msg);
    }
    
    return etot->second.chiSquared();
}

void ACMFitnessComputer::computeMultipoles(std::map<eRMS, FittingTarget> *targets,
                                           ACTMol                        *actmol)
{
    if (targets->find(eRMS::MU)->second.weight() > 0   ||
        targets->find(eRMS::QUAD)->second.weight() > 0 ||
        targets->find(eRMS::OCT)->second.weight() > 0  ||
        targets->find(eRMS::HEXADEC)->second.weight() > 0)
    {
        auto qProps = actmol->qProps();
        for(auto qp = qProps->begin(); qp < qProps->end(); ++qp)
        {
            auto qcalc = qp->qPact();
            qcalc->setQ(actmol->atomsConst());
            qcalc->calcMoments();
        }
    }
}

void ACMFitnessComputer::fillDevComputers(MsgHandler *msghandler,
                                          double      zetaDiff,
                                          bool        haveInductionCorrectionData)
{
    if (sii_->target(iMolSelect::Train, eRMS::BOUNDS)->weight() > 0 ||
        sii_->target(iMolSelect::Train, eRMS::UNPHYSICAL)->weight() > 0)
    {
        bdc_ = new BoundsDevComputer(sii_->optIndexPtr(), zetaDiff);
    }
    if (sii_->target(iMolSelect::Train, eRMS::CHARGE)->weight() > 0 ||
        sii_->target(iMolSelect::Train, eRMS::CM5)->weight() > 0)
    {
        devComputers_.push_back(new ChargeCM5DevComputer());
    }
    if (sii_->target(iMolSelect::Train, eRMS::ESP)->weight() > 0)
    {
        devComputers_.push_back(new EspDevComputer(molgen_->fit("zeta")));
    }
    if (sii_->target(iMolSelect::Train, eRMS::Polar)->weight() > 0)
    {
        if (sii_->forcefield()->polarizable())
        {
            devComputers_.push_back(new PolarDevComputer());
        }
        else
        {
            sii_->target(iMolSelect::Train, eRMS::Polar)->setWeight(0);
        }
    }
    if (sii_->target(iMolSelect::Train, eRMS::MU)->weight() > 0)
    {
        devComputers_.push_back(new MultiPoleDevComputer(MolPropObservable::DIPOLE));
    }
    if (sii_->target(iMolSelect::Train, eRMS::QUAD)->weight() > 0)
    {
        devComputers_.push_back(new MultiPoleDevComputer(MolPropObservable::QUADRUPOLE));
    }
    if (sii_->target(iMolSelect::Train, eRMS::OCT)->weight() > 0)
    {
        devComputers_.push_back(new MultiPoleDevComputer(MolPropObservable::OCTUPOLE));
    }
    if (sii_->target(iMolSelect::Train, eRMS::HEXADEC)->weight() > 0)
    {
        devComputers_.push_back(new MultiPoleDevComputer(MolPropObservable::HEXADECAPOLE));
    }
    if (sii_->target(iMolSelect::Train, eRMS::EPOT)->weight() > 0 ||
        sii_->target(iMolSelect::Train, eRMS::Interaction)->weight() > 0 ||
        sii_->target(iMolSelect::Train, eRMS::Electrostatics)->weight() > 0 ||
        sii_->target(iMolSelect::Train, eRMS::Exchange)->weight() > 0 ||
        sii_->target(iMolSelect::Train, eRMS::Dispersion)->weight() > 0 ||
        sii_->target(iMolSelect::Train, eRMS::Induction)->weight() > 0 ||
        sii_->target(iMolSelect::Train, eRMS::DeltaHF)->weight() > 0 ||
        sii_->target(iMolSelect::Train, eRMS::AllElec)->weight() > 0 ||
        sii_->target(iMolSelect::Train, eRMS::ExchInd)->weight() > 0 ||
        sii_->target(iMolSelect::Train, eRMS::Force2)->weight() > 0)
    {
        auto T = molgen_->enerBoltzTemp();
        std::map<eRMS, double> boltzmann = {
            { eRMS::EPOT, T },
            { eRMS::Interaction, T }
        };
        // Not a nice place to check stuff and throw but we have to pass the information
        // needed to the devComputer.
        double dhfWeight   = sii_->target(iMolSelect::Train, eRMS::DeltaHF)->weight();
        if (dhfWeight > 0 && !haveInductionCorrectionData)
        {
            std::string msg = gmx::formatString("Inconstent input. No induction correction energies present but -fc_deltaHF is %g",
                                                dhfWeight);
            msghandler->msg(ACTStatus::Error, msg);
        }
        auto devcomp     = new ForceEnergyDevComputer(boltzmann);
        if (msghandler->info())
        {
            for(const auto &b : boltzmann)
            {
                msghandler->write(gmx::formatString("Component %s Boltzmann temperature %g",
                                                    rmsName(b.first), b.second));
            }
        }
        devcomp->setSeparateInductionCorrection(dhfWeight >  0 && haveInductionCorrectionData);
        devComputers_.push_back(std::move(devcomp));
    }
    if (sii_->target(iMolSelect::Train, eRMS::FREQUENCY)->weight() > 0)
    {
        devComputers_.push_back(new HarmonicsDevComputer(MolPropObservable::FREQUENCY));
    }
    if (sii_->target(iMolSelect::Train, eRMS::INTENSITY)->weight() > 0)
    {
        devComputers_.push_back(new HarmonicsDevComputer(MolPropObservable::INTENSITY));
    }
}

/* * * * * * * * * * * * * * * * * * * *
* END: ACMFitnessComputer              *
* * * * * * * * * * * * * * * * * * * */

} // namespace alexandria

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
 * \author Julian Ramon Marrades Furquet <julian.marrades@hotmail.es>
 * \author Oskar Tegby <oskar.tegby@it.uu.se>
 */

#include "acmfitnesscomputer.h"

#include "act/basics/dataset.h"
#include "act/ga/genome.h"

namespace alexandria
{

/* * * * * * * * * * * * * * * * * * * *
* BEGIN: ACMFitnessComputer            *
* * * * * * * * * * * * * * * * * * * */

void ACMFitnessComputer::compute(ga::Genome    *genome,
                                 iMolSelect     trgtFit,
                                 bool           verbose)
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
    double fitness = calcDeviation(cd, trgtFit);
    genome->setFitness(trgtFit, fitness);
    if (debug)
    {
        // TODO fix printing
        //tmpInd->printParameters(debug);
    }
    if (verbose && logfile_)
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
                fprintf(logfile_, "Components of training function for %s set\n", iMolSelectName(trgtFit));
                for (const auto &ft : fts)
                {
                    ft.second.print(logfile_);
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
        int src = cr->superior();
        return static_cast<CalcDev>(cr->recv_int(src));
    }
    else 
    {
        // Send only to my helpers
        for (auto &dest : cr->helpers())
        {
            cr->send_int(dest, static_cast<int>(task));
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
        cr->recv_double_vector(src, &myparams);
        std::set<int>       mychanged;
        int nchanged = cr->recv_int(src);
        for(int i = 0; i < nchanged; i++)
        {
            mychanged.insert(cr->recv_int(src));
        }
        sii_->updateForceField(mychanged, myparams);
    }
    else 
    {
        // Send only to my helpers
        for (auto &dest : cr->helpers())
        {
            cr->send_double_vector(dest, params);
            cr->send_int(dest, changed.size());
            for(int iset : changed)
            {
                cr->send_int(dest, iset);
            }
        }
        sii_->updateForceField(changed, *params);
    }
    if (debug)
    {
        fprintf(debug, "Finished distributing parameters\n");
    }
}

double ACMFitnessComputer::calcDeviation(CalcDev    task,
                                         iMolSelect ims)
{
    auto cr = sii_->commRec();
    // Send / receive molselect group.ks
    if (cr->isHelper())
    {
        // Now who is my middleman?
        int src  = cr->superior();
        ims      = cr->recv_iMolSelect(src);
    }
    else if (CalcDev::ComputeAll != task)
    {
        // Send ims to my helpers
        for (auto &dest : cr->helpers())
        {
            cr->send_iMolSelect(dest, ims);
        }
    }
    if (debug)
    {
        fprintf(debug, "Going to do calcDeviation for %s\n",
                iMolSelectName(ims));
    }
    // Gather fitting targets
    std::map<eRMS, FittingTarget> *targets = sii_->fittingTargets(ims);
    if (nullptr == targets)
    {
        if (debug)
        {
            fprintf(debug, "Cannot find targets for %s\n", iMolSelectName(ims));
        }
        return 0;
    }
    // Reset the chi2 in FittingTargets for the given dataset in ims
    sii_->resetChiSquared(ims);

    // If actMaster or actMiddleMan, penalize out of bounds
    if (cr->isMasterOrMiddleMan() && bdc_)
    {
        bdc_->calcDeviation(forceComp_, nullptr, nullptr, targets, sii_->forcefield());
    }

    // Loop over molecules
    int ntrain = 0;
    int nlocal = 0;
    auto mymols = molgen_->actmolsPtr();
    for (auto actmol = mymols->begin(); actmol < mymols->end(); ++actmol)
    {
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
            actmol->topologyPtr()->fillParameters(sii_->forcefield());
            // Fill the fragments too if there are any
            for(auto &ft : actmol->fragmentHandler()->topologiesPtr())
            {
                ft->fillParameters(sii_->forcefield());
            }

            // Run charge generation including shell minimization
            std::vector<gmx::RVec> forces(actmol->atomsConst().size(), { 0, 0, 0 });
            std::vector<gmx::RVec> coords = actmol->xOriginal();
            immStatus imm = actmol->GenerateAcmCharges(sii_->forcefield(), forceComp_, &coords, &forces);

            // Check whether we have to disable this compound
            if (immStatus::OK != imm && removeMol_)
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
                mydev->calcDeviation(forceComp_, &(*actmol), &coords, targets, sii_->forcefield());
            }
            if (debug)
            {
                fprintf(debug, "rank %d mol %s #energies %zu tw %g\n",
                        cr->rank(), actmol->getMolname().c_str(), actmol->experimentConst().size(),
                        targets->find(eRMS::EPOT)->second.totalWeight());
            }
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
        GMX_THROW(gmx::InternalError(gmx::formatString("Cannot find %s in targets.",
                                                       rmsName(erms)).c_str()));
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
            if (ttt.second.weight() > 0)
            {
                msg += gmx::formatString("Weight for %s is %g\n", rmsName(ttt.second.erms()),
                                         ttt.second.weight());
            }
        }
        GMX_THROW(gmx::InvalidInputError(msg.c_str()));
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

void ACMFitnessComputer::fillDevComputers(const bool verbose, double zetaDiff)
{
    if (sii_->target(iMolSelect::Train, eRMS::BOUNDS)->weight() > 0 ||
        sii_->target(iMolSelect::Train, eRMS::UNPHYSICAL)->weight() > 0)
    {
        bdc_ = new BoundsDevComputer(logfile_, verbose, sii_->optIndexPtr(),
                                     zetaDiff);
    }
    if (sii_->target(iMolSelect::Train, eRMS::CHARGE)->weight() > 0 ||
        sii_->target(iMolSelect::Train, eRMS::CM5)->weight() > 0)
    {
        devComputers_.push_back(new ChargeCM5DevComputer(logfile_, verbose));
    }
    if (sii_->target(iMolSelect::Train, eRMS::ESP)->weight() > 0)
    {
        devComputers_.push_back(new EspDevComputer(logfile_, verbose, molgen_->fit("zeta")));
    }
    if (sii_->target(iMolSelect::Train, eRMS::Polar)->weight() > 0)
    {
        if (sii_->forcefield()->polarizable())
        {
            devComputers_.push_back(new PolarDevComputer(logfile_, verbose));
        }
        else
        {
            sii_->target(iMolSelect::Train, eRMS::Polar)->setWeight(0);
        }
    }
    if (sii_->target(iMolSelect::Train, eRMS::MU)->weight() > 0)
    {
        devComputers_.push_back(new MultiPoleDevComputer(logfile_, verbose,
                                                         MolPropObservable::DIPOLE));
    }
    if (sii_->target(iMolSelect::Train, eRMS::QUAD)->weight() > 0)
    {
        devComputers_.push_back(new MultiPoleDevComputer(logfile_, verbose, 
                                                         MolPropObservable::QUADRUPOLE));
    }
    if (sii_->target(iMolSelect::Train, eRMS::OCT)->weight() > 0)
    {
        devComputers_.push_back(new MultiPoleDevComputer(logfile_, verbose, 
                                                         MolPropObservable::OCTUPOLE));
    }
    if (sii_->target(iMolSelect::Train, eRMS::HEXADEC)->weight() > 0)
    {
        devComputers_.push_back(new MultiPoleDevComputer(logfile_, verbose, 
                                                         MolPropObservable::HEXADECAPOLE));
    }
    if (sii_->target(iMolSelect::Train, eRMS::EPOT)->weight() > 0 ||
        sii_->target(iMolSelect::Train, eRMS::Interaction)->weight() > 0 ||
        sii_->target(iMolSelect::Train, eRMS::Electrostatics)->weight() > 0 ||
        sii_->target(iMolSelect::Train, eRMS::Exchange)->weight() > 0 ||
        sii_->target(iMolSelect::Train, eRMS::Dispersion)->weight() > 0 ||
        sii_->target(iMolSelect::Train, eRMS::Induction)->weight() > 0 ||
        sii_->target(iMolSelect::Train, eRMS::AllElec)->weight() > 0 ||
        sii_->target(iMolSelect::Train, eRMS::Force2)->weight() > 0)
    {
        auto T = molgen_->enerBoltzTemp();
        std::map<eRMS, double> boltzmann = {
            { eRMS::EPOT, T },
            { eRMS::Interaction, T }
        };
        devComputers_.push_back(new ForceEnergyDevComputer(logfile_, verbose, boltzmann));
    }
    if (sii_->target(iMolSelect::Train, eRMS::FREQUENCY)->weight() > 0)
    {
        devComputers_.push_back(new HarmonicsDevComputer(logfile_, verbose, MolPropObservable::FREQUENCY));
    }
    if (sii_->target(iMolSelect::Train, eRMS::INTENSITY)->weight() > 0)
    {
        devComputers_.push_back(new HarmonicsDevComputer(logfile_, verbose, MolPropObservable::INTENSITY));
    }
}

/* * * * * * * * * * * * * * * * * * * *
* END: ACMFitnessComputer              *
* * * * * * * * * * * * * * * * * * * */

} // namespace alexandria

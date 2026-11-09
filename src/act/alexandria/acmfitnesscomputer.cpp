/*
 * This source file is part of the Alexandria Chemistry Toolkit.
 *
 * Copyright (C) 2021-2026
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

#include <cstdio>
#include <cstdlib>
#include <functional>
#include <limits>
#include <set>

#include "act/basics/dataset.h"
#include "act/basics/msg_handler.h"
#include "act/ga/genome.h"
#include "act/utility/communicationrecord.h"
#include "act/utility/timings.h"

namespace alexandria
{

/* * * * * * * * * * * * * * * * * * * *
* BEGIN: ACMFitnessComputer            *
* * * * * * * * * * * * * * * * * * * */
bool ACMFitnessComputer::compactCheck()
{
    static const bool check = std::getenv("ACT_GPU_DEVIATION_CHECK") != nullptr;
    return check;
}

ACMFitnessComputer::~ACMFitnessComputer()
{
    if (compactChecks_ > 0 || gpuChecks_ > 0)
    {
        fprintf(stderr, "GPU deviation check: %d evaluations compared, %d differ\n",
                compactChecks_, compactDiffer_);
        fprintf(stderr, "GPU deviation check on the GPU: %d evaluations compared, %d differ\n",
                gpuChecks_, gpuDiffer_);
    }
    for(auto &d : devComputers_)
    {
        delete d;
    }
}

void ACMFitnessComputer::compute(MsgHandler                *msghandler,
                                 ga::Genome                *genome,
                                 iMolSelect                 trgtFit)
{
    if (nullptr == genome)
    {
        GMX_THROW(gmx::InternalError("Empty genome"));
    }
    // First send around parameters. This code is run on master only.
    distributeTasks(CalcDev::Parameters);
    std::set<int> changed;
    distributeParameters(msghandler, genome->basesPtr(), changed);
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
                msghandler->tw()->writeLineFormatted("Components of training function for %s set\n", iMolSelectName(trgtFit).c_str());
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

void ACMFitnessComputer::distributeParameters(MsgHandler                *msghandler,
                                              const std::vector<double> *params,
                                              const std::set<int>       &changed)
{
    static const size_t timer = Timings::instance().index("Parameter distribution");
    ScopedTimer         scoped(timer);
    if (msghandler && msghandler->debug())
    {
        msghandler->writeDebug("Starting to distribute parameters\n");
    }
    auto cr = sii_->commRec();
    // Send / receive parameters
    if (cr->isHelper())
    {
        //! \todo Implement broadcast
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
        sii_->updateForceField(msghandler, mychanged, myparams);
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
        sii_->updateForceField(msghandler, changed, *params);
    }
    if (msghandler && msghandler->debug())
    {
        msghandler->writeDebug("Finished distributing parameters\n");
    }
}

double ACMFitnessComputer::calcDeviation(MsgHandler   *msghandler,
                                         CalcDev       task,
                                         iMolSelect    ims)
{
    DeviationTask dt;
    if (!prepareDeviation(msghandler, task, ims, &dt))
    {
        return 0;
    }
    computeEnergies(msghandler, dt);
    return finishDeviation(msghandler, &dt);
}

std::vector<double> ACMFitnessComputer::calcDeviations(MsgHandler                                     *msghandler,
                                                       const std::vector<const std::vector<double> *> &params,
                                                       iMolSelect                                      ims)
{
    std::vector<double> chi2;
    if (params.empty())
    {
        return chi2;
    }
    const std::set<int> all;
    auto useParameters = [this, msghandler, &all](const std::vector<double> *p)
    {
        sii_->updateForceField(msghandler, all, *p);
    };
    // All sets at once needs a replica for each on the GPU, and no helpers,
    // which would compute their molecules one set at a time
    const size_t nrep = forceComp_->gpuReplicas();
    DeviationTask first;
    bool batch = (params.size() > 1 && params.size() <= nrep &&
                  !sii_->commRec()->isParallelCalc());
    if (batch)
    {
        useParameters(params[0]);
        batch = (prepareDeviation(msghandler, CalcDev::Compute, ims, &first) &&
                 first.gpuCharges && !first.selected.empty());
    }
    if (batch)
    {
        // Each replica its own set, the ones without a set compute the last
        // set again. The force field has the first set already.
        std::vector<DeviationTask> dts(nrep);
        dts[0] = first;
        std::function<void(size_t)> stage = [&](size_t r)
        {
            if (r > 0)
            {
                useParameters(params[std::min(r, params.size() - 1)]);
                (void) prepareDeviation(msghandler, CalcDev::Compute, ims, &dts[r]);
            }
        };
        auto et = energyTask(first);
        setupGpuDeviations(first);
        if (forceComp_->computeRelaxedEnergiesReplicas(msghandler, sii_->forcefield(), first.selected,
                                                       stage, et.separateInductionCorrection))
        {
            for (size_t r = 0; r < params.size(); r++)
            {
                // The deviations from the energies and charges of each set,
                // and the parameters of each, for instance for the bounds
                useParameters(params[r]);
                forceComp_->selectGpuReplica(r);
                sii_->resetChiSquared(ims);
                chi2.push_back(finishDeviation(msghandler, &dts[r]));
            }
            return chi2;
        }
    }
    // One set at a time
    for (auto p : params)
    {
        distributeTasks(CalcDev::Parameters);
        distributeParameters(msghandler, p, all);
        auto cd = distributeTasks(CalcDev::Compute);
        chi2.push_back(calcDeviation(msghandler, cd, ims));
    }
    return chi2;
}

bool ACMFitnessComputer::prepareDeviation(MsgHandler    *msghandler,
                                          CalcDev        task,
                                          iMolSelect     ims,
                                          DeviationTask *dt)
{
    static const size_t timer = Timings::instance().index("Fitness computation overhead");
    ScopedTimer         scoped(timer);
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
                                             iMolSelectName(ims).c_str()));

    // Gather fitting targets
    std::map<eRMS, FittingTarget> *targets = sii_->fittingTargets(ims);
    if (nullptr == targets)
    {
        msghandler->writeDebug(gmx::formatString("Cannot find targets for %s\n", iMolSelectName(ims).c_str()));
        return false;
    }
    // Reset the chi2 in FittingTargets for the given dataset in ims
    sii_->resetChiSquared(ims);


    // Only the parameters that are optimized change, therefore only the
    // topology entries of their interaction types have to be updated.
    // All entries were filled when the topologies were made.
    std::set<InteractionType> itUpdate;
    for(const auto &oi : sii_->optIndex())
    {
        itUpdate.insert(oi.iType());
    }

    // The charge properties are only used for the ESP and the multipoles
    bool needQProps = false;
    for (auto erms : { eRMS::ESP, eRMS::MU, eRMS::QUAD, eRMS::OCT, eRMS::HEXADEC })
    {
        auto t = targets->find(erms);
        needQProps = needQProps || (t != targets->end() && t->second.weight() > 0);
    }
    // With SQE on the GPU, the processor only sets the fixed charges, as
    // long as the energies are computed there and nothing else needs them
    bool gpuCharges = false;
    {
        DeviationTask tmp;
        tmp.targets = targets;
        auto et     = energyTask(tmp);
        gpuCharges  = (!needQProps && et.compute && !et.forces &&
                       forceComp_->gpuChargesOnly(sii_->forcefield()));
    }
    dt->gpuCharges = gpuCharges;
    // Loop over molecules
    int ntrain = 0;
    int nlocal = 0;
    auto mymols = molgen_->actmolsPtr();
    // Nothing has been filled or generated for new molecules
    const size_t never = std::numeric_limits<size_t>::max();
    filledVersion_.resize(mymols->size());
    chargeVersion_.resize(mymols->size(), never);
    fixedChargeVersion_.resize(mymols->size(), never);
    // The molecules in this dataset that this node computes
    std::vector<ACTMol *> selected;
    for (auto actmol = mymols->begin(); actmol < mymols->end(); ++actmol)
    {
        msghandler->writeDebug(gmx::formatString("CalcDev: mol %s dataset %s\n",
                                                 actmol->getMolname().c_str(),
                                                 iMolSelectName(actmol->datasetType()).c_str()));
        if (ims != actmol->datasetType())
        {
            continue;
        }
        ntrain++;
        if ((actmol->support() == eSupport::Local) ||
            (task == CalcDev::ComputeAll && actmol->support() == eSupport::Remote))
        {
            size_t imol = actmol - mymols->begin();
            {
                // Now update the topology, only the interaction types of
                // which parameters changed since the last time.
                static const size_t timer = Timings::instance().index("Topology parameter update");
                ScopedTimer         scoped(timer);
                auto                &filled = filledVersion_[imol];
                std::set<InteractionType> itFill;
                for(auto itype : itUpdate)
                {
                    auto fv = filled.find(itype);
                    if (filled.end() == fv || fv->second != sii_->parameterVersion(itype))
                    {
                        itFill.insert(itype);
                    }
                }
                if (!itFill.empty())
                {
                    actmol->topologyPtr()->fillParameters(msghandler, sii_->forcefield(),
                                                          missingParameters::Error, &itFill);
                    // Fill the fragments too if there are any
                    for(auto &ft : actmol->fragmentHandler()->topologiesPtr())
                    {
                        ft.fillParameters(msghandler, sii_->forcefield(), missingParameters::Error, &itFill);
                    }
                    for(auto itype : itFill)
                    {
                        filled[itype] = sii_->parameterVersion(itype);
                    }
                }
            }
            if (chargeVersion_[imol] != sii_->chargeParameterVersion())
            {
                // Run charge generation including shell minimization,
                // only if parameters changed that the charges depend on.
                static const size_t timer = Timings::instance().index("Charge generation");
                ScopedTimer         scoped(timer);
                if (gpuCharges && actmol->fragmentHandler()->sqeOnGpu(sii_->forcefield()))
                {
                    // The charges of the processor stay outdated, the next
                    // evaluation on the processor generates them
                    if (fixedChargeVersion_[imol] != sii_->chargeParameterVersion())
                    {
                        actmol->updateFixedCharges(sii_->forcefield(), algorithm_);
                        fixedChargeVersion_[imol] = sii_->chargeParameterVersion();
                    }
                }
                else
                {
                    // The coordinates and forces are not used here, and the
                    // charge properties only for the ESP and multipoles
                    std::vector<gmx::RVec> forces(actmol->atomsConst().size(), { 0, 0, 0 });
                    std::vector<gmx::RVec> coords = actmol->xOriginal();
                    actmol->generateCharges(msghandler, sii_->forcefield(), forceComp_,
                                            algorithm_, &coords, &forces, needQProps, false);
                    chargeVersion_[imol] = sii_->chargeParameterVersion();
                }
            }

            // Check whether we have to disable this compound
            if (!msghandler->ok() && removeMol_)
            {
                actmol->setSupport(eSupport::No);
                msghandler->resetStatus();
                continue;
            }

            computeMultipoles(targets, &(*actmol));

            if (devComputers_.size() == 0)
            {
                msghandler->msg(ACTStatus::Warning,
                                "No devComputers: no deviation from reference values will be computed.");
            }
            msghandler->writeDebug(gmx::formatString("CalcDev: rank %d mol %s #energies %zu tw %g\n",
                                                     cr->rank(), actmol->getMolname().c_str(), actmol->experimentConst().size(),
                                                     targets->find(eRMS::EPOT)->second.totalWeight()));
            selected.push_back(&(*actmol));
            nlocal++;
        }
    }
    // If the GPU cannot compute these molecules after all, the processor
    // generates the charges that were left to it, and computes everything
    // from now on, see ForceComputer::prepareGpuResident
    if (gpuCharges && !selected.empty() &&
        !forceComp_->prepareGpuResident(sii_->forcefield(), selected))
    {
        static const size_t timer = Timings::instance().index("Charge generation");
        ScopedTimer         scoped(timer);
        std::vector<ACTMol *> kept;
        for (auto molptr : selected)
        {
            size_t imol = molptr - &(*mymols->begin());
            if (chargeVersion_[imol] != sii_->chargeParameterVersion())
            {
                std::vector<gmx::RVec> forces(molptr->atomsConst().size(), { 0, 0, 0 });
                std::vector<gmx::RVec> coords = molptr->xOriginal();
                molptr->generateCharges(msghandler, sii_->forcefield(), forceComp_,
                                        algorithm_, &coords, &forces, needQProps, false);
                chargeVersion_[imol] = sii_->chargeParameterVersion();
                if (!msghandler->ok() && removeMol_)
                {
                    molptr->setSupport(eSupport::No);
                    msghandler->resetStatus();
                    nlocal--;
                    continue;
                }
            }
            kept.push_back(molptr);
        }
        selected = std::move(kept);
    }
    dt->task     = task;
    dt->ims      = ims;
    dt->targets  = targets;
    dt->selected = std::move(selected);
    dt->ntrain   = ntrain;
    dt->nlocal   = nlocal;
    return true;
}

ForceEnergyDevComputer::EnergyTask ACMFitnessComputer::energyTask(const DeviationTask &dt) const
{
    for (const auto *mydev : devComputers_)
    {
        auto fe = dynamic_cast<const ForceEnergyDevComputer *>(mydev);
        if (fe && dt.targets)
        {
            return fe->energyTask(*dt.targets);
        }
    }
    return ForceEnergyDevComputer::EnergyTask();
}

void ACMFitnessComputer::setupGpuDeviations(const DeviationTask &dt)
{
    if (!compactCheck() || !dt.targets)
    {
        return;
    }
    for (auto *mydev : devComputers_)
    {
        auto fe = dynamic_cast<ForceEnergyDevComputer *>(mydev);
        if (fe)
        {
            (void) fe->setupGpuDeviations(forceComp_, dt.selected, *dt.targets, lossFunction_);
        }
    }
}

void ACMFitnessComputer::computeEnergies(MsgHandler *msghandler, const DeviationTask &dt)
{
    auto et = energyTask(dt);
    setupGpuDeviations(dt);
    if (et.compute)
    {
        forceComp_->computeRelaxedEnergies(msghandler, sii_->forcefield(), dt.selected, true,
                                           et.separateInductionCorrection, et.forces);
    }
}

double ACMFitnessComputer::finishDeviation(MsgHandler *msghandler, DeviationTask *dt)
{
    static const size_t timer = Timings::instance().index("Fitness computation overhead");
    ScopedTimer         scoped(timer);
    auto  task     = dt->task;
    auto  ims      = dt->ims;
    auto *targets  = dt->targets;
    auto &selected = dt->selected;
    int   ntrain   = dt->ntrain;
    int   nlocal   = dt->nlocal;
    // The timer of each DevComputer
    if (devTimers_.size() != devComputers_.size())
    {
        devTimers_.clear();
        for (DevComputer *mydev : devComputers_)
        {
            devTimers_.push_back(Timings::instance().index("Deviation " + mydev->name()));
        }
    }
    for (size_t d = 0; d < devComputers_.size(); d++)
    {
        DevComputer *mydev = devComputers_[d];
        ScopedTimer timer(devTimers_[d]);
        // The energies were computed, see computeEnergies
        auto fe = dynamic_cast<ForceEnergyDevComputer *>(mydev);
        if (fe)
        {
            fe->setComputeEnergies(false);
        }
        // For developers, the same from the energies as the GPU gave them,
        // computed on the processor and on the GPU
        std::map<eRMS, FittingTarget> compact, onGpu;
        bool                          haveCompact = false, haveGpu = false;
        if (fe && compactCheck())
        {
            compact = *targets;
            for (auto &t : compact)
            {
                t.second.reset();
            }
            onGpu       = compact;
            haveCompact = fe->calcDeviationCompact(forceComp_, selected, &compact, lossFunction_);
            haveGpu     = fe->calcDeviationGpu(forceComp_, &onGpu);
        }
        mydev->calcDeviation(msghandler, forceComp_, selected, targets,
                             sii_->forcefield(), lossFunction_);
        // Whether the energy targets are the same as those of calcDeviation
        auto compare = [targets](const std::map<eRMS, FittingTarget> &other, const char *what, int *differ)
        {
            bool same = true;
            for (const auto &t : other)
            {
                if (t.first == eRMS::BOUNDS || t.first == eRMS::TOT)
                {
                    continue;
                }
                const auto &r = targets->find(t.first)->second;
                if (r.chiSquared() != t.second.chiSquared() || r.totalWeight() != t.second.totalWeight())
                {
                    same = false;
                    if (*differ < 10)
                    {
                        fprintf(stderr, "GPU deviation check: %s chi2 %.17g weight %.17g, %s %.17g %.17g\n",
                                rmsName(t.first), r.chiSquared(), r.totalWeight(), what,
                                t.second.chiSquared(), t.second.totalWeight());
                    }
                }
            }
            *differ += same ? 0 : 1;
        };
        if (haveCompact)
        {
            compactChecks_ += 1;
            compare(compact, "compact", &compactDiffer_);
        }
        if (haveGpu)
        {
            gpuChecks_ += 1;
            compare(onGpu, "on the GPU", &gpuDiffer_);
        }
        if (fe)
        {
            fe->setComputeEnergies(true);
        }
    }
    // Sum the terms of the chi-squared once we have done calculations
    // for all the molecules.
    sii_->sumChiSquared(msghandler, task == CalcDev::Compute, ims);
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
                                            iMolSelectName(ims).c_str(), rmsName(erms));
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
                                          bool        haveInductionCorrectionData)
{
    if (sii_->target(iMolSelect::Train, eRMS::BOUNDS)->weight() > 0)
    {
        devComputers_.push_back(new BoundsDevComputer(sii_->optIndexPtr()));
    }
    if (sii_->target(iMolSelect::Train, eRMS::ESP)->weight() > 0)
    {
        devComputers_.push_back(new EspDevComputer());
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

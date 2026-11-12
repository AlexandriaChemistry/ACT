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
    // For developers, how many molecules were computed and skipped, see prepareDeviation
    if (std::getenv("ACT_GPU_SKIP_STATS") && nComputedMol_ + nSkippedMol_ > 0)
    {
        fprintf(stderr, "GPU skip statistics: %zu molecules computed, %zu skipped\n",
                nComputedMol_, nSkippedMol_);
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
    computeEnergies(msghandler, &dt);
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
            if (r >= params.size())
            {
                // A replica without a set of parameters computes nothing
                forceComp_->setGpuActiveMolecules(r, std::vector<char>(first.selected.size(), 0));
            }
            else if (r > 0)
            {
                useParameters(params[std::min(r, params.size() - 1)]);
                dts[r].replica = r;
                (void) prepareDeviation(msghandler, CalcDev::Compute, ims, &dts[r]);
                if (dts[r].skipping)
                {
                    forceComp_->setGpuActiveMolecules(r, dts[r].active);
                }
            }
        };
        auto et = energyTask(first);
        setupGpuDeviations(first);
        // Only the deviations if the GPU can compute them
        bool onlyDeviations = useGpuDeviations(first);
        if (first.skipping)
        {
            forceComp_->setGpuActiveMolecules(0, first.active);
        }
        if (forceComp_->computeRelaxedEnergiesReplicas(msghandler, sii_->forcefield(), first.selected,
                                                       stage, et.separateInductionCorrection, onlyDeviations))
        {
            for (auto &dt : dts)
            {
                dt.gpuDeviations = onlyDeviations;
                if (dt.skipping && !onlyDeviations)
                {
                    msghandler->fatal("Molecules were skipped, but the GPU could not compute the deviations of the others");
                }
            }
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

/*! \brief Whether the processor needs the parameters of an interaction type in the topology
 * When the GPU makes the rows of its parameter tables, the moving shells,
 * the virtual sites and the structures of the charge generation from the
 * force field, see ForceComputer::gpuRecipeRows, the processor only needs
 * the polarization and the virtual sites if it generates the charges of the
 * molecule, which relaxes the shells.
 * \param[in] itype       The interaction type
 * \param[in] hostCharges Whether the processor generates the charges
 */
static bool hostNeedsParameters(InteractionType itype, bool hostCharges)
{
    return hostCharges && (InteractionType::POLARIZATION == itype || isVsite(itype));
}

void ACMFitnessComputer::updateMolecule(MsgHandler                      *msghandler,
                                        ACTMol                          *actmol,
                                        size_t                           imol,
                                        const std::set<InteractionType> &itUpdate,
                                        bool                             gpuCharges,
                                        bool                             needQProps,
                                        bool                             allTypes)
{
    // Whether the processor generates the charges, to be safe also when the
    // charge properties are computed
    bool gpuSqe      = gpuCharges && actmol->fragmentHandler()->sqeOnGpu(sii_->forcefield());
    bool hostCharges = needQProps || !gpuSqe;
    // Without anything to fill, see hostNeedsParameters
    if (allTypes || hostCharges)
    {
        // Now update the topology, only the interaction types of
        // which parameters changed since the last time.
        static const size_t timer = Timings::instance().index("Topology parameter update");
        ScopedTimer         scoped(timer);
        auto                &filled = filledVersion_[imol];
        std::set<InteractionType> itFill;
        for(auto itype : itUpdate)
        {
            if (!allTypes && !hostNeedsParameters(itype, hostCharges))
            {
                continue;
            }
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
    // The GPU generates the charges and makes the fixed ones, see
    // ForceComputer::gpuChargeRecipes. The charges of the processor stay
    // outdated, the next evaluation on the processor makes them.
    if (gpuSqe && !allTypes && forceComp_->gpuChargeRecipes())
    {
        return;
    }
    if (chargeVersion_[imol] != sii_->chargeParameterVersion())
    {
        // Run charge generation including shell minimization,
        // only if parameters changed that the charges depend on.
        static const size_t timer = Timings::instance().index("Charge generation");
        ScopedTimer         scoped(timer);
        if (gpuSqe)
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
}

bool ACMFitnessComputer::needChargeProperties(MsgHandler                          *msghandler,
                                              const std::map<eRMS, FittingTarget> &targets,
                                              iMolSelect                           ims)
{
    static const std::map<eRMS, MolPropObservable> multipoles = {
        { eRMS::MU, MolPropObservable::DIPOLE }, { eRMS::QUAD, MolPropObservable::QUADRUPOLE },
        { eRMS::OCT, MolPropObservable::OCTUPOLE }, { eRMS::HEXADEC, MolPropObservable::HEXADECAPOLE }
    };
    auto &data = qPropsData_[ims];
    bool  need = false;
    for (auto erms : { eRMS::ESP, eRMS::MU, eRMS::QUAD, eRMS::OCT, eRMS::HEXADEC })
    {
        auto t = targets.find(erms);
        if (t == targets.end() || t->second.weight() <= 0)
        {
            continue;
        }
        auto d = data.find(erms);
        if (d == data.end())
        {
            // As EspDevComputer and MultiPoleDevComputer find the data
            bool found = false;
            for (auto &mol : *molgen_->actmolsPtr())
            {
                if (mol.datasetType() != ims)
                {
                    continue;
                }
                for (auto &qp : *mol.qProps())
                {
                    if (eRMS::ESP == erms)
                    {
                        found = found || qp.qgenResp()->nEsp() > 0;
                    }
                    else
                    {
                        found = found || qp.qPqmConst().hasMultipole(multipoles.at(erms));
                    }
                }
                if (found)
                {
                    break;
                }
            }
            d = data.insert({ erms, found }).first;
            if (!found && msghandler)
            {
                msghandler->msg(ACTStatus::Warning,
                                gmx::formatString("The fitting target %s of the %s set has a weight, but none of its molecules has data for it. It does not change the objective function.",
                                                  rmsName(erms), iMolSelectName(ims).c_str()));
            }
        }
        need = need || d->second;
    }
    return need;
}

std::vector<double> ACMFitnessComputer::parameterValues()
{
    std::vector<double> values;
    for (auto &oi : *sii_->optIndexPtr())
    {
        auto p = oi.forceFieldParameter();
        if (!p)
        {
            oi.findForceFieldParameter(sii_->forcefield());
            p = oi.forceFieldParameter();
        }
        values.push_back(p ? p->value() : 0.0);
    }
    return values;
}

void ACMFitnessComputer::buildDependencies(MsgHandler                      *msghandler,
                                           const std::vector<ACTMol *>     &selected,
                                           const std::set<InteractionType> &itUpdate,
                                           bool                             needQProps)
{
    static const size_t timer = Timings::instance().index("Parameter dependencies");
    ScopedTimer         scoped(timer);
    auto pd      = sii_->forcefield();
    auto mymols  = molgen_->actmolsPtr();
    auto refresh = [&]()
    {
        for (auto molptr : selected)
        {
            updateMolecule(msghandler, molptr, molptr - &(*mymols->begin()), itUpdate, true, needQProps);
        }
    };
    // All parameters must be in the force field for updateForceField
    std::vector<double> values;
    {
        std::set<int> all;
        values = parameterValues();
        sii_->updateForceField(msghandler, all, values);
        refresh();
    }
    auto base = forceComp_->gpuMoleculeSignatures(pd);
    molDeps_.assign(selected.size(), {});
    if (base.size() != selected.size())
    {
        return;
    }
    size_t ndep    = 0;
    auto  &optIndex = *sii_->optIndexPtr();
    for (size_t k = 0; k < values.size(); k++)
    {
        // A small change within the bounds of the parameter
        auto   p  = optIndex[k].forceFieldParameter();
        double v  = values[k];
        double d  = 1e-3*(std::abs(v) + 1e-3);
        double nv = v + d;
        if (p && nv > p->maximum())
        {
            nv = v - d;
        }
        if (p && nv < p->minimum())
        {
            nv = 0.5*(p->minimum() + p->maximum());
        }
        if (!p || nv == v)
        {
            // It cannot be changed here, to be safe all molecules depend on it
            for (auto &deps : molDeps_)
            {
                deps.push_back(k);
            }
            ndep += selected.size();
            continue;
        }
        auto changed = values;
        changed[k]   = nv;
        sii_->updateForceField(msghandler, { static_cast<int>(k) }, changed);
        refresh();
        auto sig = forceComp_->gpuMoleculeSignatures(pd);
        for (size_t m = 0; m < selected.size(); m++)
        {
            if (sig[m] != base[m])
            {
                molDeps_[m].push_back(k);
                ndep += 1;
            }
        }
        sii_->updateForceField(msghandler, { static_cast<int>(k) }, values);
    }
    refresh();
    depsMols_.assign(selected.begin(), selected.end());
    depsVersion_ = forceComp_->gpuVersion();
    molCache_.clear();
    if (msghandler->info() && !values.empty())
    {
        msghandler->write(gmx::formatString("A change of one parameter changes %.1f of the %zu molecules on average, only those are computed again",
                                            static_cast<double>(ndep)/values.size(), selected.size()));
    }
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
    bool needQProps = needChargeProperties(msghandler, *targets, ims);
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
    // Molecules whose parameters did not change can be skipped, if the GPU
    // computes only the deviations, see useGpuDeviations, and the
    // dependencies of the molecules on the parameters are known. For
    // developers, ACT_NO_GPU_SKIP turns this off, read each time such that
    // a test can compare with and without, see TrainFFMpiTest.
    const bool        noSkip = std::getenv("ACT_NO_GPU_SKIP") != nullptr;
    static const bool noDev  = std::getenv("ACT_NO_GPU_DEVIATIONS") != nullptr;
    bool skipCandidate = (!noSkip && !noDev && !compactCheck() && task == CalcDev::Compute && gpuCharges);
    bool haveDeps      = (skipCandidate && depsVersion_ != 0 && depsVersion_ == forceComp_->gpuVersion());
    bool canSkip       = haveDeps && forceComp_->gpuDeviationsReady();
    if (haveDeps && cacheVersion_ != forceComp_->gpuDeviationsVersion())
    {
        // The deviations of the molecules are those of other targets
        molCache_.clear();
        cacheVersion_ = forceComp_->gpuDeviationsVersion();
    }
    if (haveDeps && molCache_.size() <= dt->replica)
    {
        molCache_.resize(dt->replica + 1, std::vector<std::array<MolCacheEntry, 2>>(molDeps_.size()));
    }
    std::vector<double>   values = haveDeps ? parameterValues() : std::vector<double>();
    std::vector<ACTMol *> skipped;
    // When the GPU makes the rows of its parameter tables from the force
    // field, the topologies only need the parameters that the processor uses,
    // see updateMolecule. If the GPU does not after all, for instance since
    // its molecules are made again, the others are filled below.
    const bool allTypes = !(gpuCharges && forceComp_->gpuRecipeRows());
    dt->active.clear();
    dt->keys.clear();
    dt->hit.clear();
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
            // Molecules whose parameters did not change since they were
            // computed are not computed again, see buildDependencies
            size_t m      = selected.size();
            bool   active = true;
            int    hit    = -1;
            std::vector<double> key;
            if (haveDeps && m < molDeps_.size())
            {
                for (auto k : molDeps_[m])
                {
                    key.push_back(values[k]);
                }
                for (int e = 0; e < 2; e++)
                {
                    const auto &entry = molCache_[dt->replica][m][e];
                    if (entry.valid && entry.key == key)
                    {
                        hit = e;
                    }
                }
                active = !(canSkip && hit >= 0 &&
                           actmol->fragmentHandler()->sqeOnGpu(sii_->forcefield()));
            }
            if (active)
            {
                updateMolecule(msghandler, &(*actmol), imol, itUpdate, gpuCharges, needQProps, allTypes);
            }
            else
            {
                skipped.push_back(&(*actmol));
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
            dt->active.push_back(active ? 1 : 0);
            dt->keys.push_back(std::move(key));
            dt->hit.push_back(hit);
            selected.push_back(&(*actmol));
            nlocal++;
        }
    }
    // The molecules can only be skipped if they are those of the
    // dependencies, otherwise they are updated after all
    bool consistent = (haveDeps && selected.size() == depsMols_.size() &&
                       std::equal(selected.begin(), selected.end(), depsMols_.begin()));
    auto updateSkipped = [&]()
    {
        for (auto molptr : skipped)
        {
            updateMolecule(msghandler, molptr, molptr - &(*mymols->begin()), itUpdate, gpuCharges, needQProps, allTypes);
        }
        skipped.clear();
        std::fill(dt->active.begin(), dt->active.end(), 1);
    };
    if (!consistent)
    {
        updateSkipped();
    }
    // If the GPU cannot compute these molecules after all, the processor
    // generates the charges that were left to it, and computes everything
    // from now on, see ForceComputer::prepareGpuResident
    bool gpuOk = true;
    if (gpuCharges && !selected.empty() &&
        !forceComp_->prepareGpuResident(sii_->forcefield(), selected))
    {
        gpuOk = false;
        updateSkipped();
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
    if (!allTypes && !(gpuOk && forceComp_->gpuRecipeRows()))
    {
        for (auto molptr : selected)
        {
            updateMolecule(msghandler, molptr, molptr - &(*mymols->begin()), itUpdate, gpuCharges, needQProps);
        }
    }
    dt->skipping = canSkip && consistent && gpuOk;
    if (!dt->skipping)
    {
        std::fill(dt->active.begin(), dt->active.end(), 1);
    }
    if (!(haveDeps && consistent && gpuOk))
    {
        // Without keys nothing is kept, see finishDeviation
        for (auto &key : dt->keys)
        {
            key.clear();
        }
        std::fill(dt->hit.begin(), dt->hit.end(), -1);
    }
    dt->haveKeys = haveDeps && consistent && gpuOk;
    // Which molecules depend on which parameters, once for the molecules on
    // the GPU, but not while the parameters of the replicas are staged, see calcDeviations
    if (skipCandidate && gpuOk && dt->replica == 0 && !selected.empty() &&
        depsVersion_ != forceComp_->gpuVersion())
    {
        buildDependencies(msghandler, selected, itUpdate, needQProps);
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

bool ACMFitnessComputer::useGpuDeviations(const DeviationTask &dt)
{
    // For developers, to compare with computing the energies of the molecules
    static const bool off = std::getenv("ACT_NO_GPU_DEVIATIONS") != nullptr;
    if (off || compactCheck() || dt.task != CalcDev::Compute || !dt.gpuCharges || !dt.targets ||
        dt.selected.empty())
    {
        return false;
    }
    for (auto *mydev : devComputers_)
    {
        auto fe = dynamic_cast<ForceEnergyDevComputer *>(mydev);
        if (fe)
        {
            return fe->setupGpuDeviations(forceComp_, dt.selected, *dt.targets, lossFunction_);
        }
    }
    return false;
}

void ACMFitnessComputer::computeEnergies(MsgHandler *msghandler, DeviationTask *dt)
{
    auto et = energyTask(*dt);
    setupGpuDeviations(*dt);
    dt->gpuDeviations = false;
    if (et.compute)
    {
        // Only the deviations if the GPU can compute them, of the molecules
        // whose parameters changed, see prepareDeviation
        if (dt->skipping)
        {
            forceComp_->setGpuActiveMolecules(-1, dt->active);
        }
        if (useGpuDeviations(*dt) &&
            forceComp_->computeGpuDeviations(msghandler, sii_->forcefield(), dt->selected))
        {
            dt->gpuDeviations = true;
            return;
        }
        if (dt->skipping)
        {
            msghandler->fatal("Molecules were skipped, but the GPU could not compute the deviations of the others");
        }
        forceComp_->computeRelaxedEnergies(msghandler, sii_->forcefield(), dt->selected, true,
                                           et.separateInductionCorrection, et.forces);
    }
}

double ACMFitnessComputer::finishDeviation(MsgHandler *msghandler, DeviationTask *dt)
{
    static const size_t timer = Timings::instance().index("Fitness computation overhead");
    ScopedTimer         scoped(timer);
    if (!gpuWhyNotReported_ && !forceComp_->gpuResidentWhyNot().empty())
    {
        // Once, such that it is clear why the processor computes everything
        msghandler->msg(ACTStatus::Warning,
                        gmx::formatString("The GPU does not compute the molecules of the training, since %s. The processor computes them.",
                                          forceComp_->gpuResidentWhyNot().c_str()));
        gpuWhyNotReported_ = true;
    }
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
    // The charges of the molecules that were skipped are not those of
    // these parameters, but those of the replica staged last or of the
    // other set kept, therefore their bounds are those kept with their
    // deviations, see MolCacheEntry.
    std::vector<double> molBound;
    if (dt->gpuDeviations && targets->find(eRMS::BOUNDS) != targets->end())
    {
        molBound.resize(selected.size(), 0);
        for (size_t m = 0; m < selected.size(); m++)
        {
            if (dt->skipping && !dt->active[m])
            {
                molBound[m] = molCache_[dt->replica][m][dt->hit[m]].bound;
            }
            else
            {
                molBound[m] = BoundsDevComputer::moleculeBound(selected[m], sii_->forcefield(), lossFunction_);
            }
        }
    }
    for (size_t d = 0; d < devComputers_.size(); d++)
    {
        DevComputer *mydev = devComputers_[d];
        ScopedTimer timer(devTimers_[d]);
        if (!molBound.empty() && dynamic_cast<BoundsDevComputer *>(mydev))
        {
            auto bounds = targets->find(eRMS::BOUNDS);
            for (auto b : molBound)
            {
                bounds->second.increase(1, b);
            }
            continue;
        }
        // The energies were computed, see computeEnergies
        auto fe = dynamic_cast<ForceEnergyDevComputer *>(mydev);
        if (fe && dt->gpuDeviations)
        {
            // Only the deviations, the energies of the molecules were not
            // computed. Those of the molecules that were skipped are those
            // kept for the same parameters, see prepareDeviation.
            size_t nmol = selected.size();
            std::vector<GpuMoleculeDeviations>         fresh(nmol);
            std::vector<const GpuMoleculeDeviations *> devs(nmol, nullptr);
            auto  &cache = molCache_.size() > dt->replica ? molCache_[dt->replica] : molCache_.emplace_back();
            for (size_t m = 0; m < nmol; m++)
            {
                if (dt->skipping && !dt->active[m])
                {
                    devs[m] = &cache[m][dt->hit[m]].dev;
                    nSkippedMol_ += 1;
                }
                else if (forceComp_->gpuMoleculeDeviations(m, &fresh[m]))
                {
                    devs[m] = &fresh[m];
                    nComputedMol_ += 1;
                }
                else
                {
                    msghandler->fatal("The GPU did not compute the deviations of the energies");
                }
            }
            fe->addGpuDeviations(devs, targets);
            // Keep the new ones, and the last other one of each molecule
            for (size_t m = 0; m < nmol && dt->haveKeys && m < cache.size(); m++)
            {
                auto &entries = cache[m];
                if (dt->skipping && !dt->active[m])
                {
                    if (dt->hit[m] == 1)
                    {
                        std::swap(entries[0], entries[1]);
                    }
                }
                else
                {
                    entries[1] = std::move(entries[0]);
                    entries[0] = { true, dt->keys[m], std::move(fresh[m]), m < molBound.size() ? molBound[m] : 0.0 };
                }
            }
            continue;
        }
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

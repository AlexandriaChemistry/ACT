/*
 * This source file is part of the Alexandria Chemistry Toolkit.
 *
 * Copyright (C) 2021-2026
 *
 * this program is free software; you can redistribute it and/or
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
#include "actmiddleman.h" 

#include "act/alexandria/acminitializer.h"
#include "act/alexandria/mcmcmutator.h"
#include "act/alexandria/percentmutator.h"
#include "act/alexandria/staticindividualinfo.h"
#include "act/forces/forcecomputerstatistics.h"
#include "act/utility/portable_random.h"
#include "act/utility/communicationrecord.h"
#include "gromacs/utility/gmxassert.h"
#include "gromacs/utility/textwriter.h"
#include "act/utility/timings.h"

namespace alexandria
{

ACTMiddleMan::~ACTMiddleMan()
{
    delete mutator_;
    for (auto m : lockstepMutators_)
    {
        delete m;
    }
}

void sendQueueWork(const CommunicationRecord             *cr,
                   int                                    dest,
                   const std::vector<const ga::Genome *> &genomes,
                   TrainFFMiddlemanMode                   mode,
                   const std::vector<int>                &seeds,
                   int                                    generation)
{
    GMX_RELEASE_ASSERT(genomes.size() == seeds.size(), "A seed for each genome");
    // II. Signify the middleman to continue, as in the loop of ACTMiddleMan::run
    cr->send_data(dest);
    cr->send(dest, static_cast<int>(genomes.size()));
    for (size_t i = 0; i < genomes.size(); i++)
    {
        genomes[i]->Send(cr, dest);
        // III. What to do, with the random numbers
        cr->send(dest, mode);
        cr->send(dest, seeds[i]);
    }
    cr->send(dest, generation);
}

ACTMiddleMan::ACTMiddleMan(MsgHandler                *msghandler,
                           MolGen                    *mg,
                           StaticIndividualInfo      *sii,
                           GAConfigHandler           *gach,
                           BayesConfigHandler        *bch,
                           gmx_output_env_t          *oenv,
                           ChargeGenerationAlgorithm  algorithm,
                           bool                       openConvFiles)
: gach_(gach), id_(sii->commRec()->middleManOrdinal())
{
    // This ica
    int seed;
    sii->commRec()->recv(sii->commRec()->superior(), &seed);
    // Standard mersenne_twister_engine seeded with what we received
    std::mt19937 gen(seed);
    
    ACMInitializer initializer(sii, gach->randomInit(), portableSeed(&gen));
    
    // Create and initialize the individual
    ind_.reset(std::move(initializer.initialize()));

    // Create force computer
    double epsilonr;
    epsilonr = sii->forcefield()->epsilonr();
    forceComp_.init(bch->shellToler(), bch->shellMaxIter(), epsilonr);
    forceComp_.setKeepShells(bch->shellKeep());
    forceComp_.setGpuResident(bch->gpuResident());

    //! \todo Fitness computer. What about those false flags?
    fitComp_.init(msghandler, sii, mg, false, &forceComp_, algorithm, gach->lossFunction());
    
    if (gach->optimizer() == OptimizerAlg::GA)
    {
        mutator_ = new alexandria::PercentMutator(sii, portableSeed(&gen), algorithm, gach->percent());
    }
    else
    {
        //! \todo we need to make some logfiles for the middlemen, because they apparently cannot write to the global logfile
        auto mut = new alexandria::MCMCMutator(portableSeed(&gen),
                                               bch, &fitComp_, sii,
                                               bch->evaluateTestset(),
                                               algorithm, 
                                               gach->maxGenerations());
        mut->setNstFitness(gach->nstFitness());
        if (openConvFiles)
        {
            mut->openParamConvFiles(oenv);
            mut->openChi2ConvFile(oenv);
        }
        mutator_ = std::move(mut);
        // The other individuals that are computed at once, which get their
        // seeds with their work, without convergence files. MCMCMutator::mutateLockstep
        // does not evaluate the test set.
        if (sii->commRec()->queue() && !bch->evaluateTestset())
        {
            for (int i = 1; i < gach->lockstep(); i++)
            {
                auto other = new alexandria::MCMCMutator(portableSeed(&gen), bch, &fitComp_, sii,
                                                         bch->evaluateTestset(), algorithm,
                                                         gach->maxGenerations());
                other->setNstFitness(gach->nstFitness());
                lockstepMutators_.push_back(other);
            }
        }
    }
    // A copy of the molecules on the GPU for each individual computed at
    // once, see ACMFitnessComputer::calcDeviations. Helpers would compute
    // their molecules one individual at a time.
    if (sii->commRec()->queue() && gach->lockstep() > 1 && !sii->commRec()->isParallelCalc())
    {
        forceComp_.setGpuReplicas(gach->lockstep());
    }
}

void ACTMiddleMan::computeQueueWork(MsgHandler *msghandler)
{
    auto cr     = ind_->sii()->commRec();
    int  master = cr->superior();
    // As sendQueueWork sends it. The individuals may be any, their random
    // numbers and the generation come with them, such that the results do
    // not depend on which middleman computes them.
    int n;
    cr->recv(master, &n);
    std::vector<ga::Genome> genomes(n);
    std::vector<int>        seeds(n);
    TrainFFMiddlemanMode    mode = TrainFFMiddlemanMode::FITNESS;
    for (int i = 0; i < n; i++)
    {
        genomes[i].Receive(cr, master);
        TrainFFMiddlemanMode m;
        cr->recv(master, &m);
        GMX_RELEASE_ASSERT(i == 0 || m == mode, "All individuals sent at once have the same task");
        mode = m;
        cr->recv(master, &seeds[i]);
    }
    int generation;
    cr->recv(master, &generation);
    std::vector<const std::vector<double> *> params;
    for (auto &g : genomes)
    {
        params.push_back(g.basesPtr());
    }
    auto train = iMolSelect::Train;
    auto test  = iMolSelect::Test;
    switch (mode)
    {
    case TrainFFMiddlemanMode::STATISTICS:
        {
            alexandria::ForceComputerStatistics fcStats;
            (void) fcStats.statistics(cr, fitComp_.forceComputer(), 0, false);
        }
        break;
    case TrainFFMiddlemanMode::FITNESS:
        {
            auto chi2 = fitComp_.calcDeviations(msghandler, params, train);
            for (int i = 0; i < n; i++)
            {
                cr->send(master, chi2[i]);
            }
        }
        break;
    case TrainFFMiddlemanMode::INITIAL:
        {
            auto chi2 = fitComp_.calcDeviations(msghandler, params, train);
            for (int i = 0; i < n; i++)
            {
                genomes[i].setFitness(train, chi2[i]);
                if (gach_->evaluateTestset())
                {
                    fitComp_.compute(msghandler, &genomes[i], test);
                }
                genomes[i].Send(cr, master);
            }
        }
        break;
    case TrainFFMiddlemanMode::MUTATION:
        {
            // Each starts a new round of mutations from its genome
            std::vector<ga::Genome> best = genomes;
            std::vector<bool>       minimum(n, false);
            auto mcmc = dynamic_cast<MCMCMutator *>(mutator_);
            // In lockstep only if the GPU computes the individuals together.
            // Otherwise ACMFitnessComputer::calcDeviations computes them one
            // after another, each time with all parameters changed, such that
            // all molecules are computed again rather than those that the
            // changed parameter affects. One after another then gives the same
            // results, much faster.
            if (mcmc && n > 1 && lockstepMutators_.size() + 1 >= static_cast<size_t>(n) &&
                forceComp_.gpuCharges(ind_->sii()->forcefield()))
            {
                std::vector<MCMCMutator *> mutators;
                std::vector<ga::Genome *>  gptr, bptr;
                for (int i = 0; i < n; i++)
                {
                    auto m = i == 0 ? mcmc : lockstepMutators_[i-1];
                    m->reseed(seeds[i]);
                    m->setGeneration(generation);
                    mutators.push_back(m);
                    gptr.push_back(&genomes[i]);
                    bptr.push_back(&best[i]);
                }
                MCMCMutator::mutateLockstep(msghandler, mutators, gptr, bptr);
                for (int i = 0; i < n; i++)
                {
                    minimum[i] = mutators[i]->foundMinimum();
                }
            }
            else
            {
                for (int i = 0; i < n; i++)
                {
                    mutator_->reseed(seeds[i]);
                    mutator_->setGeneration(generation);
                    mutator_->mutate(msghandler, &genomes[i], &best[i], gach_->prMut());
                    if (mcmc)
                    {
                        minimum[i] = mcmc->foundMinimum();
                    }
                    else
                    {
                        fitComp_.compute(msghandler, &genomes[i], train);
                    }
                }
            }
            for (int i = 0; i < n; i++)
            {
                // IV. Send the mutated genome, or with Hybrid the best one
                // found if it is better
                auto *result = minimum[i] ? &best[i] : &genomes[i];
                if (gach_->evaluateTestset())
                {
                    fitComp_.compute(msghandler, result, test);
                }
                result->Send(cr, master);
            }
        }
        break;
    }
}
    
void ACTMiddleMan::run(MsgHandler *msghandler)
{
    static const size_t timer = Timings::instance().index("Middleman loop");
    const ScopedTimer   scoped(timer);
    auto cr = ind_->sii()->commRec();
    GMX_RELEASE_ASSERT(cr->isMiddleMan(), "I thought I was the middle man...");
    // Check whether we need to accept a genome from the master
    int read;
    int master = cr->superior();
    cr->recv(master, &read);
    if (read == 1)
    {
        cr->recv(master, ind_->genomePtr()->basesPtr());
    }
    // ForceComputer Statistics
    alexandria::ForceComputerStatistics fcStats;
    
    // Start by computing my own fitness
    fitComp_.compute(msghandler, ind_->genomePtr(), iMolSelect::Train);
    if (gach_->evaluateTestset())
    {
        fitComp_.compute(msghandler, ind_->genomePtr(), iMolSelect::Test);
    }
    // Copy initial genome to bestgenome
    ind_->setBestGenome(ind_->genome());
    // I.
    // Send my initial genome and fitness to the master
    ind_->genome().Send(cr, master);
    if (msghandler->debug())
    {
        msghandler->writeDebug(ind_->bestGenomePtr()->print("Initial genome"));
    }
    auto cont = CommunicationStatus::OK;

    // II. Do we need to continue?
    cont = cr->recv_data(master);
    // Two iterations in this loop corresponds to one in the master.
    while (CommunicationStatus::RECV_DATA == cont)
    {
        if (cr->queue())
        {
            computeQueueWork(msghandler);
            cont = cr->recv_data(master);
            continue;
        }
        // Now get the new genome
        ind_->genomePtr()->Receive(cr, master);
        // Copy current genome to bestgenome for new round of mutations
        ind_->setBestGenome(ind_->genome());
        // III. Receive assignment from master.  
        TrainFFMiddlemanMode mode;
        cr->recv(master, &mode);
        switch (mode)
        {
        case TrainFFMiddlemanMode::FITNESS:
            fitComp_.compute(msghandler, ind_->genomePtr(), iMolSelect::Train);
            cr->send(master, ind_->genome().fitness(iMolSelect::Train));
            break;
        case TrainFFMiddlemanMode::INITIAL:
            // As for my own first individual above
            fitComp_.compute(msghandler, ind_->genomePtr(), iMolSelect::Train);
            if (gach_->evaluateTestset())
            {
                fitComp_.compute(msghandler, ind_->genomePtr(), iMolSelect::Test);
            }
            ind_->genome().Send(cr, master);
            break;
        case TrainFFMiddlemanMode::STATISTICS:
            (void) fcStats.statistics(cr, fitComp_.forceComputer(), 0, false);
            break;
        case TrainFFMiddlemanMode::MUTATION:
            mutator_->mutate(msghandler, ind_->genomePtr(), ind_->bestGenomePtr(), 
                             gach_->prMut());

            if (gach_->optimizer() == OptimizerAlg::GA)
            {
                fitComp_.compute(msghandler, ind_->genomePtr(), iMolSelect::Train);
                if (gach_->evaluateTestset())
                {
                    fitComp_.compute(msghandler, ind_->genomePtr(), iMolSelect::Test);
                }
                // IV. Send the mutated genome
                ind_->genomePtr()->Send(cr, master);
            }
            else
            {
                // IV. If we are working with Hybrid or MCMC, send 
                // the best genome found to the MASTER.
                if (mutator_->foundMinimum())
                {
                    if (gach_->evaluateTestset())
                    {
                        fitComp_.compute(msghandler, ind_->bestGenomePtr(), iMolSelect::Test);
                    }
                    ind_->bestGenomePtr()->Send(cr, master);
                }
                else
                {
                    if (gach_->evaluateTestset())
                    {
                        fitComp_.compute(msghandler, ind_->genomePtr(), iMolSelect::Test);
                    }
                    ind_->genomePtr()->Send(cr, master);
                }
            }
            // With a queue, the master asks for the statistics of all
            // individuals at once, see TrainFFMiddlemanMode::STATISTICS
            if (msghandler->info() && !cr->queue())
            {
                (void) fcStats.statistics(cr, fitComp_.forceComputer(), 0, false);
            }
            break;
        }
        // II.
        // Do we need to continue?
        cont = cr->recv_data(master);
    }
    // Stop my helpers too.
    stopHelpers();
}

void ACTMiddleMan::stopHelpers()
{
    (void) fitComp_.distributeTasks(CalcDev::Stop);
}

} // namespace alexandria

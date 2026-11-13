/*
 * This source file is part of the Alexandria Chemistry Toolkit.
 *
 * Copyright (C) 2014-2026
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
#include "mcmcmutator.h"

#include <filesystem>
#include <sys/types.h>
#include <sys/stat.h>

#include "gromacs/utility/basedefinitions.h"

#include "act/alexandria/sensitivity.h"
#include "act/utility/memory_check.h"
#include "act/utility/jsontree.h"
#include "gromacs/fileio/xvgr.h"
#include "gromacs/utility/textwriter.h"
#include "act/utility/timings.h"

namespace alexandria
{

MCMCMutator::MCMCMutator(int                        seed,
                         BayesConfigHandler        *bch,
                         ACMFitnessComputer        *fitComp,
                         StaticIndividualInfo      *sii,
                         bool                       evaluateTestSet,
                         ChargeGenerationAlgorithm  algorithm,
                         int                        maxGenerations)
    : Mutator(seed, algorithm), evaluateTestSet_(evaluateTestSet)
{
    bch_     = bch;
    fitComp_ = fitComp;
    sii_     = sii;
    size_t nParam = sii_->nParam();
    pMean_.resize(nParam, 0.0);
    pSigma_.resize(nParam, 0.0);
    attemptedMoves_.resize(nParam, 0);
    acceptedMoves_.resize(nParam, 0);
    maxGenerations_ = maxGenerations;
}

bool MCMCMutator::checkMutate(MsgHandler       *msghandler,
                              const ga::Genome *genome,
                              MutateState      *state)
{
    size_t nParam = genome->nBase();
    state->sum.assign(nParam, 0);
    state->sumOfSquares.assign(nParam, 0);
    state->nsum = 0;
    if (sii_->xvgConv().empty() || sii_->xvgEpot().empty())
    {
        if (msghandler)
        {
            msghandler->msg(ACTStatus::Warning,
                            "You forgot to call setOutputFiles. Back to the drawing board.");
        }
        return false;
    }
    if (genome->nBase() == 0)
    {
        fprintf(stderr, "No parameters to train. Did you run alexandria geometry_ff?\n");
        return false;
    }
    return true;
}

void MCMCMutator::startMutate(MsgHandler  *msghandler,
                              MutateState *state,
                              bool         messages)
{
    // Save initial evaluation
    state->initEval = state->prevEval;
    if (msghandler)
    {
        msghandler->writeDebug(memory_usage());
    }
    state->beta0 = bch_->computeBeta(0, maxGenerations_, 0);
    if (messages && bch_->globalAnnealBegin() < 1 && myGeneration_ == 0)
    {
        msghandler->msg(ACTStatus::Info,
                        gmx::formatString("Global annealing starting from T = %g, cooling over %d generations",
                                          bch_->temperature(), maxGenerations_));
    }
    if (messages && bch_->annealing() && std::abs(state->beta0*bch_->temperature()-1) > 1e-6)
    {
        msghandler->msg(ACTStatus::Info,
                        gmx::formatString("Local annealing starting from T = %g", 1.0/state->beta0));
    }
    state->iterOffset = myGeneration_*bch_->maxIter();
}

void MCMCMutator::collectStep(const ga::Genome *genome, MutateState *state, int iter)
{
    // For the second half of the optimization, collect data
    // to find the mean and standard deviation of each
    // parameter
    //! \todo make optional
    if (iter >= bch_->maxIter()/2)
    {
        size_t k = 0;
        for (auto &p : genome->bases())
        {
            state->sum[k]          += p;
            state->sumOfSquares[k] += gmx::square(p);
            k                      += 1;
        }
        state->nsum++;
    }
}

void MCMCMutator::endMutate(const ga::Genome *bestGenome, MutateState *state)
{
    myGeneration_ += 1;
    // OPTIMIZATION IS COMPLETE
    //! \todo Make optional
    computeMeanSigma(state->sum, state->nsum, &state->sumOfSquares);

    // Assume no better minimum was found
    bMinimum_ = false;
    auto ims  = iMolSelect::Train;
    if (bestGenome->hasFitness(ims) &&
        bestGenome->fitness(ims) < state->initEval[ims])
    {
        bMinimum_ = true;
    }
}

void MCMCMutator::mutate(MsgHandler                *msghandler,
                         ga::Genome                *genome,
                         ga::Genome                *bestGenome,
                         gmx_unused double          prMut)
{
    static const size_t timer = Timings::instance().index("MCMC steps");
    ScopedTimer         scoped(timer);
    MutateState state;
    if (!checkMutate(msghandler, genome, &state))
    {
        return;
    }
    // Gather initial chi2 evaluations for training and test sets
    // Update all parameters in forcefield
    std::set<int> changed;
    auto cd = CalcDev::Parameters;
    fitComp_->distributeTasks(cd);
    fitComp_->distributeParameters(msghandler, genome->basesPtr(), changed);

    cd        = CalcDev::Compute;
    auto ims  = iMolSelect::Train;
    fitComp_->distributeTasks(cd);
    auto chi2 = fitComp_->calcDeviation(msghandler, cd, ims);
    state.prevEval[ims] = chi2;
    genome->setFitness(ims, chi2);
    if (evaluateTestSet_)
    {
        auto ims_test  = iMolSelect::Test;
        auto chi2_test = fitComp_->calcDeviation(msghandler, cd, ims_test);
        state.prevEval[ims_test] = chi2_test;
        genome->setFitness(ims_test, chi2_test);
    }
    startMutate(msghandler, &state, true);
    if (!evaluateTestSet_ && fitComp_->chainReady())
    {
        mutateOnGpu(msghandler, genome, bestGenome, &state);
        endMutate(bestGenome, &state);
        return;
    }
    // Optimization loop
    const auto imsTrain = iMolSelect::Train;
    const auto imsTest  = iMolSelect::Test;
    size_t     nParam   = genome->nBase();
    for (int iter = 0; iter < bch_->maxIter(); iter++)
    {
        for (size_t pp = 0; pp < nParam; pp++)
        {
            auto proposal = proposeStep(genome);
            std::set<int> changedOne = { static_cast<int>(proposal.paramIndex) };
            // Update FF parameter data structure with
            // the new value of parameter j
            (void) fitComp_->distributeTasks(CalcDev::Parameters);
            fitComp_->distributeParameters(msghandler, genome->basesPtr(), changedOne);

            std::map<iMolSelect, double> currEval;
            // Evaluate the energy on training set
            (void) fitComp_->distributeTasks(CalcDev::Compute);
            currEval[imsTrain] = fitComp_->calcDeviation(msghandler, CalcDev::Compute, imsTrain);
            // Evaluate the energy on the test set only on whole steps!
            if (evaluateTestSet_)
            {
                if (pp == 0)
                {
                    // Recompute for test set
                    currEval[imsTest] = fitComp_->calcDeviation(msghandler, CalcDev::Compute, imsTest);
                }
                else
                {
                    // Keep previous value
                    currEval[imsTest] = state.prevEval.find(imsTest)->second;
                }
            }
            if (!finishStep(msghandler, genome, bestGenome, &state, proposal, currEval, pp, iter, false))
            {
                // forcefield on helpers needs to know if we undo the change
                fitComp_->distributeTasks(CalcDev::Parameters);
                fitComp_->distributeParameters(msghandler, genome->basesPtr(), changedOne);
            }
            collectStep(genome, &state, iter);
        }
    }
    endMutate(bestGenome, &state);
}

std::vector<gpu::ChainStep> MCMCMutator::chainSteps(size_t nParam)
{
    std::vector<gpu::ChainStep> steps;
    for (int iter = 0; iter < bch_->maxIter(); iter++)
    {
        double beta = bch_->computeBeta(myGeneration_, maxGenerations_, iter);
        for (size_t pp = 0; pp < nParam; pp++)
        {
            gpu::ChainStep st;
            st.param  = randIndex();
            st.change = randNum();
            st.accept = randNum();
            st.beta   = beta;
            steps.push_back(st);
        }
    }
    return steps;
}

void MCMCMutator::replayStep(MsgHandler              *msghandler,
                             ga::Genome              *genome,
                             ga::Genome              *bestGenome,
                             MutateState             *state,
                             const gpu::ChainStep    &step,
                             const gpu::ChainResults &results,
                             size_t                   s,
                             size_t                   pp,
                             int                      iter,
                             bool                     lockstep)
{
    const auto imsTrain = iMolSelect::Train;
    Proposal proposal;
    proposal.paramIndex   = step.param;
    proposal.storeParam   = genome->bases()[proposal.paramIndex];
    proposal.acceptRandom = step.accept;
    attemptedMoves_[proposal.paramIndex] += 1;
    genome->setBase(proposal.paramIndex, results.trial[s]);
    std::map<iMolSelect, double> currEval = { { imsTrain, results.chi2[s] } };
    (void) finishStep(msghandler, genome, bestGenome, state, proposal, currEval, pp, iter, lockstep,
                      results.accepted[s] ? 1 : 0);
    collectStep(genome, state, iter);
}

void MCMCMutator::mutateOnGpu(MsgHandler  *msghandler,
                              ga::Genome  *genome,
                              ga::Genome  *bestGenome,
                              MutateState *state)
{
    const auto imsTrain = iMolSelect::Train;
    size_t     nParam   = genome->nBase();
    auto       steps    = chainSteps(nParam);
    gpu::ChainResults results;
    bool haveBest = bestGenome->hasFitness(imsTrain);
    fitComp_->runChain(msghandler, genome->bases(), steps, state->prevEval[imsTrain], haveBest,
                       haveBest ? bestGenome->fitness(imsTrain) : 0.0, &results);
    // The steps again on the processor, as the GPU did them, for the
    // statistics, the output and the best genome
    size_t s = 0;
    for (int iter = 0; iter < bch_->maxIter(); iter++)
    {
        for (size_t pp = 0; pp < nParam; pp++, s++)
        {
            replayStep(msghandler, genome, bestGenome, state, steps[s], results, s, pp, iter, true);
        }
    }
    if (genome->bases() != results.params)
    {
        msghandler->fatal("The parameters of the steps on the processor differ from those on the GPU");
    }
    // Writing a checkpoint above changed the force field
    sii_->updateForceField(msghandler, std::set<int>(), genome->bases());
}

void MCMCMutator::mutateLockstep(MsgHandler                       *msghandler,
                                 const std::vector<MCMCMutator *> &mutators,
                                 const std::vector<ga::Genome *>  &genomes,
                                 const std::vector<ga::Genome *>  &bestGenomes)
{
    static const size_t timer = Timings::instance().index("MCMC steps");
    ScopedTimer         scoped(timer);
    size_t n = mutators.size();
    if (n == 0)
    {
        return;
    }
    GMX_RELEASE_ASSERT(genomes.size() == n && bestGenomes.size() == n,
                       "A genome and a best genome for each mutator");
    std::vector<MutateState>                 states(n);
    std::vector<const std::vector<double> *> params;
    for (size_t i = 0; i < n; i++)
    {
        GMX_RELEASE_ASSERT(!mutators[i]->evaluateTestSet_, "Cannot evaluate the test set in lockstep");
        if (!mutators[i]->checkMutate(msghandler, genomes[i], &states[i]))
        {
            return;
        }
        params.push_back(genomes[i]->basesPtr());
    }
    auto       fitComp  = mutators[0]->fitComp_;
    const auto imsTrain = iMolSelect::Train;
    auto       chi2     = fitComp->calcDeviations(msghandler, params, imsTrain);
    for (size_t i = 0; i < n; i++)
    {
        states[i].prevEval[imsTrain] = chi2[i];
        genomes[i]->setFitness(imsTrain, chi2[i]);
        mutators[i]->startMutate(msghandler, &states[i], i == 0);
    }
    // All genomes have the same number of parameters and iterations
    const auto *bch    = mutators[0]->bch_;
    size_t      nParam = genomes[0]->nBase();
    if (fitComp->chainsReady(n))
    {
        // Each genome a chain of its own on the GPU, see mutateOnGpu
        std::vector<std::vector<gpu::ChainStep>> steps(n);
        std::vector<char>                        haveBest(n);
        std::vector<double>                      prevChi2(n), bestChi2(n);
        for (size_t i = 0; i < n; i++)
        {
            steps[i]    = mutators[i]->chainSteps(nParam);
            haveBest[i] = bestGenomes[i]->hasFitness(imsTrain) ? 1 : 0;
            bestChi2[i] = haveBest[i] ? bestGenomes[i]->fitness(imsTrain) : 0.0;
            prevChi2[i] = states[i].prevEval[imsTrain];
        }
        std::vector<gpu::ChainResults> results;
        fitComp->runChains(msghandler, params, steps, prevChi2, haveBest, bestChi2, &results);
        // The steps again on the processor in the order of the loop below
        size_t s = 0;
        for (int iter = 0; iter < bch->maxIter(); iter++)
        {
            for (size_t pp = 0; pp < nParam; pp++, s++)
            {
                for (size_t i = 0; i < n; i++)
                {
                    mutators[i]->replayStep(msghandler, genomes[i], bestGenomes[i], &states[i], steps[i][s],
                                            results[i], s, pp, iter, true);
                }
            }
        }
        for (size_t i = 0; i < n; i++)
        {
            if (genomes[i]->bases() != results[i].params)
            {
                msghandler->fatal("The parameters of the steps on the processor differ from those on the GPU");
            }
            mutators[i]->endMutate(bestGenomes[i], &states[i]);
        }
        return;
    }
    std::vector<Proposal> proposals(n);
    for (int iter = 0; iter < bch->maxIter(); iter++)
    {
        for (size_t pp = 0; pp < nParam; pp++)
        {
            for (size_t i = 0; i < n; i++)
            {
                proposals[i] = mutators[i]->proposeStep(genomes[i]);
            }
            chi2 = fitComp->calcDeviations(msghandler, params, imsTrain);
            for (size_t i = 0; i < n; i++)
            {
                std::map<iMolSelect, double> currEval = { { imsTrain, chi2[i] } };
                (void) mutators[i]->finishStep(msghandler, genomes[i], bestGenomes[i], &states[i],
                                               proposals[i], currEval, pp, iter, true);
                mutators[i]->collectStep(genomes[i], &states[i], iter);
            }
        }
    }
    for (size_t i = 0; i < n; i++)
    {
        mutators[i]->endMutate(bestGenomes[i], &states[i]);
    }
}

MCMCMutator::Proposal MCMCMutator::proposeStep(ga::Genome *genome)
{
    Proposal proposal;
    // Pick a random parameter index
    proposal.paramIndex = randIndex();
    // Store the original value of the parameter
    proposal.storeParam = genome->bases()[proposal.paramIndex];
    // Change the parameter
    changeParam(genome, proposal.paramIndex);
    attemptedMoves_[proposal.paramIndex] += 1;
    // Also when the step goes downhill, see Proposal
    proposal.acceptRandom = randNum();
    return proposal;
}

bool MCMCMutator::finishStep(MsgHandler                         *msghandler,
                             ga::Genome                         *genome,
                             ga::Genome                         *bestGenome,
                             MutateState                        *state,
                             const Proposal                     &proposal,
                             const std::map<iMolSelect, double> &currEval,
                             size_t                              pp,
                             int                                 iter,
                             bool                                lockstep,
                             int                                 decided)
{
    const auto imsTrain   = iMolSelect::Train;
    const auto imsTest    = iMolSelect::Test;
    auto       prevEval   = &state->prevEval;
    size_t     paramIndex = proposal.paramIndex;
    double     deltaEval  = currEval.find(imsTrain)->second - prevEval->find(imsTrain)->second;

    // Accept any downhill move
    bool accept = (deltaEval < 0);

    // For an uphill move apply the Metropolis Criteria
    // to decide whether to accept or reject the new parameter,
    // unless the GPU decided already, see mutateOnGpu
    if (decided >= 0)
    {
        accept = decided != 0;
    }
    else if (!accept)
    {
        // Get beta for the MCMC algorithm.
        state->beta0 = bch_->computeBeta(myGeneration_, maxGenerations_, iter);
        const double mcProbability   = std::exp( - ( (state->beta0) / (sii_->weightedTemperature())[paramIndex] ) * deltaEval );
        accept = (mcProbability > proposal.acceptRandom);
    }

    // Fractional iteration taking into account the inner loop with <pp> over <nParam>
    const double xiter = state->iterOffset + iter + (1.0*pp)/genome->nBase();
    if (accept)
    {
        // If the parameter change is accepted
        prevEval->find(imsTrain)->second = currEval.find(imsTrain)->second;
        genome->setFitness(imsTrain, currEval.find(imsTrain)->second);
        if (evaluateTestSet_)
        {
            prevEval->find(imsTest)->second = currEval.find(imsTest)->second;
            genome->setFitness(imsTest, currEval.find(imsTest)->second);
        }
        acceptedMoves_[paramIndex] += 1;
        // and if a new minimum was found
        if (!bestGenome->hasFitness(imsTrain) ||
            currEval.find(imsTrain)->second < bestGenome->fitness(imsTrain))
        {
            // If pointer to log file exists, write information about new minimum
            printNewMinimum(msghandler, currEval, xiter);

            *bestGenome = *genome;
            if (bch_->checkPoint())
            {
                if (msghandler)
                {
                    msghandler->write(gmx::formatString("Will write FF to checkpoint file %s",
                                                        sii_->outputFile().c_str()));
                }
                if (lockstep)
                {
                    // The force field may have the parameters of another genome
                    sii_->updateForceField(msghandler, std::set<int>(), genome->bases());
                }
                sii_->saveState(false);
            }
        }
    }
    else
    {
        // If the parameter change is not accepted
        // Set the old value of the parameter back
        genome->setBase(paramIndex, proposal.storeParam);
    }

    if (nStep_ % nstFitness_ == 0)
    {
        printParameterStep(genome, xiter);
        printChi2Step(*prevEval, xiter);
    }
    nStep_ += 1;
    return accept;
}

void MCMCMutator::computeMeanSigma(const std::vector<double>    &sum,
                                   const int                     nsum,
                                         std::vector<double>    *sum_of_sq)
{
    if (nsum <= 0) return;
    double ps2 = 0.0;
    for (size_t k = 0; k < sii_->nParam(); k++)
    {
        pMean_[k]        = (sum[k]/nsum);
        (*sum_of_sq)[k] /= nsum;
        ps2              = std::max(0.0, (*sum_of_sq)[k]-gmx::square(pMean_[k]));
        pSigma_[k]       = std::sqrt(ps2);
    }
}                                         

void MCMCMutator::changeParam(ga::Genome *genome,
                              size_t      j)
{
    std::vector<double> *param = genome->basesPtr();
    
    GMX_RELEASE_ASSERT(j < genome->nBase(), "Parameter out of range");
    // One random number for each step, see Proposal
    double rnd = randNum();
    if (rnd == 0.5)
    {
        // Make sure the parameter changes!
        rnd = 1;
    }
    (*param)[j] += (2*rnd-1) * sii_->stepSize()[j];
    if (sii_->mutability()[j] == Mutability::Bounded)
    {
        if ((*param)[j] < sii_->lowerBound()[j])
        {
            (*param)[j] = sii_->lowerBound()[j];
        }
        else if ((*param)[j] > sii_->upperBound()[j])
        {
            (*param)[j] = sii_->upperBound()[j];
        }
    }
}

void MCMCMutator::printMonteCarloStatistics(gmx::TextWriter  *tw,
                                            const ga::Genome &initialGenome,
                                            const ga::Genome &bestGenome)
{
    if (!tw)
    {
        return;
    }

    const auto paramNames = sii_->paramNames();
    const auto ntrain = sii_->nTrain();
    const auto weightedTemperature = sii_->weightedTemperature();

    tw->writeStringFormatted("Monte Carlo statistics of parameters after optimization\n");
    tw->writeStringFormatted("#best %zu #mean %zu #sigma %zu #param %zu\n",
                             bestGenome.nBase(), pMean_.size(), 
                             pSigma_.size(), paramNames.size());
    tw->writeStringFormatted("Parameter                     Ncopies Initial   Best    Mean    Sigma Attempt  Acceptance  T-Weight\n");
    for (size_t k = 0; k < bestGenome.nBase(); k++)
    {
        double acceptance_ratio = 0;
        if (attemptedMoves_[k] > 0)
        {
            acceptance_ratio = 100*(double(acceptedMoves_[k])/attemptedMoves_[k]);
        }
        tw->writeStringFormatted("%-30s  %5d  %6.3f  %6.3f  %6.3f  %6.3f    %4d %5.1f%%  %10.5f\n",
                                 paramNames[k].c_str(), ntrain[k],
                                 initialGenome.base(k), bestGenome.base(k), 
                                 pMean_[k], pSigma_[k],
                                 attemptedMoves_[k], acceptance_ratio, weightedTemperature[k]);
    }
}

void MCMCMutator::printNewMinimum(MsgHandler                         *msghandler,
                                  const std::map<iMolSelect, double> &chi2,
                                  double                              xiter)
{
    if (!msghandler->verbose())
    {
        return;
    }
    auto tw = msghandler->tw();
    if (!tw)
    {
        tw = msghandler->twDebug();
    }
    if (tw)
    {
        tw->writeStringFormatted("Middleman %i ", sii_->id());
        if (evaluateTestSet_)
        {
            tw->writeStringFormatted("iter %10g. Found new minimum at %10g. Corresponding energy on the test set: %g\n",
                                     xiter, chi2.find(iMolSelect::Train)->second,
                                     chi2.find(iMolSelect::Test)->second);
        }
        else
        {
            tw->writeStringFormatted("iter %10g. Found new minimum at %10g\n",
                                     xiter, chi2.find(iMolSelect::Train)->second);
        }
    }
}             

void MCMCMutator::printParameterStep(ga::Genome *genome,
                                     double      xiter)
{
    const auto &paramClassIndex = sii_->paramClassIndex();

    auto bases = genome->bases();
    if (fpc_.size() == sii_->paramClass().size())
    {
        // Write iteration number to each parameter convergence 
        // surveillance file
        for (auto &fp: fpc_) 
        {
            if (fp.get())
            {
                fprintf(fp.get(), "%8f", xiter);
            }
        }
        // Write value of each parameter to its respective surveillance file
        for (size_t k = 0; k < genome->nBase(); k++)  
        {
            if (fpc_[paramClassIndex[k]].get())
            {
                fprintf(fpc_[paramClassIndex[k]].get(), "  %10g", bases[k]);
            }
        }
        for (auto &fp: fpc_)
        {
            if (fp.get())
            {
                fprintf(fp.get(), "\n");
            }
        }
    }
}                             

void MCMCMutator::printChi2Step(const std::map<iMolSelect, double> &chi2,
                                double                              xiter)
{
    if (fpe_.get() == nullptr)
    {
        // If fpe is a null pointer, return
        return;
    }

    if (evaluateTestSet_)
    {
        fprintf(fpe_.get(), "%8f  %10g  %10g\n",
                xiter, chi2.find(iMolSelect::Train)->second,
                chi2.find(iMolSelect::Test)->second);
    }
    else
    {
        fprintf(fpe_.get(), "%8f  %10g\n",
                xiter, chi2.find(iMolSelect::Train)->second);
    }
}                    

void MCMCMutator::openParamConvFiles(const gmx_output_env_t *oenv)
{
    const std::vector<std::string> &pClass = sii_->paramClass();
    for (size_t i = 0; i < pClass.size(); i++)
    {
        // The files of each individual are in its directory
        auto fileName = gmx::formatString("%sind%d-%s-%s",
                                          sii_->prefix().c_str(), sii_->id(),
                                          pClass[i].c_str(),
                                          std::filesystem::path(sii_->xvgConv()).filename().c_str());
        gmx::FilePtr fp(xvgropen(fileName.c_str(), "Parameter convergence",
                                 "iteration", "", oenv));
        fpc_.push_back(std::move(fp));

        std::vector<const char*> tmpParamNames;
        for (size_t j = 0; j < sii_->paramNames().size(); j++)
        {
            if ( (sii_->paramClassIndex())[j] == i )
            {
                tmpParamNames.push_back( (sii_->paramNames())[j].c_str() );
            }
        }
        xvgr_legend(fpc_[i].get(), tmpParamNames.size(),
                    tmpParamNames.data(), oenv);
    }
}

void MCMCMutator::openChi2ConvFile(const gmx_output_env_t *oenv)
{
    auto fileName = gmx::formatString("%sind%d-%s", sii_->prefix().c_str(), sii_->id(),
                                      std::filesystem::path(sii_->xvgEpot()).filename().c_str());
    fpe_.reset(xvgropen(fileName.c_str(),
                        "Chi squared",
                        "Iteration",
                        "Unknown units",
                        oenv));
    if (evaluateTestSet_)
    {
        std::vector<std::string> legend;
        legend.push_back(iMolSelectName(iMolSelect::Train));
        legend.push_back(iMolSelectName(iMolSelect::Test));
        xvgrLegend(fpe_.get(), legend, oenv);
    }
}

} //namespace alexandria

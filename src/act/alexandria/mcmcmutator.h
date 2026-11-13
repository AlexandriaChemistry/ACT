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

#ifndef ALEXANDRIA_MCMCMUTATOR_H
#define ALEXANDRIA_MCMCMUTATOR_H

#include <algorithm>
#include <cmath>

#include "act/forces/forcecomputer.h"
#include "act/ga/genome.h"
#include "act/ga/mutator.h"
#include "gromacs/utility/fileptr.h"
#include "confighandler.h"
#include "acmfitnesscomputer.h"
#include "acmindividual.h"

namespace alexandria
{

class JsonTree;

/*!
 * \brief Mutator which applies MCMC optimization to an ACMIndividual.
 */
class MCMCMutator : public ga::Mutator
{

private:
    //! Did we find a minimum yet?
    bool                       bMinimum_  = false;
    //! Pointer to BayesConfigHandler
    BayesConfigHandler        *bch_       = nullptr;
    //! Pointer to ACMFitnessComputer
    ACMFitnessComputer        *fitComp_   = nullptr;
    //! Pointer to Individual
    StaticIndividualInfo      *sii_       = nullptr;
    //! Attempted changes for each parameter
    std::vector<int>           attemptedMoves_;
    //! Accepted changes for each parameter
    std::vector<int>           acceptedMoves_;
    //! Mean of each parameter
    std::vector<double>        pMean_;
    //! Standard deviation of each parameter
    std::vector<double>        pSigma_;
    //! Internal counter for generation number, see setGeneration
    int                        myGeneration_   = 0;
    //! Number of steps between writing to the convergence files
    int                        nstFitness_     = 1;
    //! Number of steps done
    long                       nStep_          = 0;
    //! Copy of max generations from GA config handler
    int                        maxGenerations_ = 1;
    //! Convergence file for each parameter type
    std::vector<gmx::FilePtr>  fpc_;
    //! Convergence file for Chi2
    gmx::FilePtr               fpe_;
    //! Evaluate the test set at each iteration?
    bool                       evaluateTestSet_ = false;
    //! All used randIndices
    std::vector<int>           randIndices_;
    /*!
     * \brief Change force field parameter at a given index for an individual
     * \param[in] genome   pointer to the genome
     * \param[in] j     index of the parameter to change
     */
    void changeParam(ga::Genome *genome,
                     size_t      j);

    /*!
     * \brief Print new minimum to log file and, if necessary, print params to debug file
     * \param[in] msghandler        The message handler, what else?
     * \param[in] chi2              the new minimum
     * \param[in] xiter             fractional iteration. E.g., if we are halfway through iteration 3, it is 3.5
     */
    void printNewMinimum(MsgHandler                         *msghandler,
                         const std::map<iMolSelect, double> &chi2,
                         double                              xiter);

    /*!
     * \brief Print parameter values to their respective surveillance files
     * \param[in] genome Pointer to the genome
     * \param[in] xiter  Fractional iteration. E.g., if we are halfway through iteration 3, it is 3.5
     */
    void printParameterStep(ga::Genome *genome,
                            double      xiter);                                          

    /*!
     * \brief Write \f$ \chi^2 \f$ value of a genome to its convergence file, if it exists
     * \param[in] chi2                  chi2 for train and test sets
     * \param[in] xiter                 fractional iteration. E.g., if we are halfway through iteration 3, it is 3.5
     */
    void printChi2Step(const std::map<iMolSelect, double> &chi2,
                       double                              xiter);

    //! \return a random index of the force field parameter vector, see indexOfRandom
    int randIndex()
    {
        int ri = static_cast<int>(indexOfRandom(randNum(), sii_->nParam()));
        randIndices_.push_back(ri);
        return ri;
    }

    /*!
     * \brief Compute mean and standard deviation for each force field parameter
     * \param[in] sum           over "nsum" iterations, the sum of each parameter
     * \param[in] nsum          number of iterations to compute statistics over
     * \param[in] sum_of_sq     over "nsum" iterations, the sum of each parameter squared
     */
    void computeMeanSigma(const std::vector<double>    &sum,
                          const int                     nsum,
                                std::vector<double>    *sum_of_sq);

    //! The state of mutate between its steps
    struct MutateState
    {
        //! The chi squared of the current parameters, and that of the first ones
        std::map<iMolSelect, double> prevEval, initEval;
        //! Sums of the parameters and of their squares, for computeMeanSigma
        std::vector<double>          sum, sumOfSquares;
        //! The number of terms in the sums
        int                          nsum       = 0;
        //! Beta for annealing
        double                       beta0      = 0;
        //! Total iteration number in all generations prior to this one
        int                          iterOffset = 0;
    };

    //! A change of one parameter, see proposeStep
    /*! \brief A step of the Monte Carlo chain
     * Each step draws the same three random numbers, the parameter, the
     * change and the one of the Metropolis criterion, also when the latter
     * is not needed. The random numbers then only depend on the seed and the
     * number of steps, not on the values of the objective function, such
     * that they can be drawn before the steps are computed.
     */
    struct Proposal
    {
        //! The parameter
        size_t paramIndex = 0;
        //! Its value before the change
        double storeParam = 0;
        //! The random number of the Metropolis criterion, in [0, 1)
        double acceptRandom = 0;
    };

    /*!
     * \brief Check whether the genome can be mutated, and prepare the sums
     * \param[in]  msghandler The message and status handler
     * \param[in]  genome     The genome
     * \param[out] state      The state
     * \return whether there is something to do
     */
    bool checkMutate(MsgHandler *msghandler, const ga::Genome *genome, MutateState *state);

    /*!
     * \brief Start the steps, after the chi squared of the genome was computed
     * \param[in]    msghandler The message and status handler
     * \param[inout] state      The state, with the chi squared of the genome
     * \param[in]    messages   Whether to print the annealing settings
     */
    void startMutate(MsgHandler *msghandler, MutateState *state, bool messages);

    /*!
     * \brief Change a random parameter, the first half of a step of MCMC
     * \param[inout] genome The genome
     * \return the change
     */
    Proposal proposeStep(ga::Genome *genome);

    /*!
     * \brief Accept or reject the change, the second half of a step of MCMC
     * \param[in]    msghandler The message and status handler
     * \param[inout] genome     The genome, its parameter is restored if rejected
     * \param[inout] bestGenome The best genome, to be filled
     * \param[inout] state      The state
     * \param[in]    proposal   The change of proposeStep
     * \param[in]    currEval   The chi squared with the change
     * \param[in]    pp         Index of inner loop over number of parameters
     * \param[in]    iter       Current iteration number
     * \param[in]    lockstep   Whether the force field may have the parameters
     *                          of another genome, see mutateLockstep
     * \param[in]    decided    Whether the step is accepted, -1 to decide here
     * \return whether the change was accepted
     */
    bool finishStep(MsgHandler                         *msghandler,
                    ga::Genome                         *genome,
                    ga::Genome                         *bestGenome,
                    MutateState                        *state,
                    const Proposal                     &proposal,
                    const std::map<iMolSelect, double> &currEval,
                    size_t                              pp,
                    int                                 iter,
                    bool                                lockstep,
                    int                                 decided = -1);

    /*! \brief The steps of mutate on the GPU, see ACMFitnessComputer::runChain
     * The steps are replayed on the processor with the decisions of the
     * GPU, for the statistics, the output and the best genome.
     * \param[in]    msghandler For messages
     * \param[inout] genome     The genome
     * \param[inout] bestGenome The best genome
     * \param[inout] state      The state of mutate
     */
    void mutateOnGpu(MsgHandler  *msghandler,
                     ga::Genome  *genome,
                     ga::Genome  *bestGenome,
                     MutateState *state);

    /*! \brief The steps of the chain on the GPU with their random numbers
     * In the order of proposeStep and changeParam, which do not depend on
     * whether the steps are accepted.
     * \param[in] nParam The number of parameters
     * \return The steps of all iterations
     */
    std::vector<gpu::ChainStep> chainSteps(size_t nParam);

    /*! \brief A step of the chain on the GPU again on the processor, see mutateOnGpu
     * \param[in]    msghandler For messages
     * \param[inout] genome     The genome
     * \param[inout] bestGenome The best genome
     * \param[inout] state      The state of mutate
     * \param[in]    step       The step
     * \param[in]    results    The results of the steps on the GPU
     * \param[in]    s          The index of the step
     * \param[in]    pp         The parameter step of the iteration
     * \param[in]    iter       The iteration
     * \param[in]    lockstep   Whether the force field may have the parameters of another genome
     */
    void replayStep(MsgHandler              *msghandler,
                    ga::Genome              *genome,
                    ga::Genome              *bestGenome,
                    MutateState             *state,
                    const gpu::ChainStep    &step,
                    const gpu::ChainResults &results,
                    size_t                   s,
                    size_t                   pp,
                    int                      iter,
                    bool                     lockstep);

    /*!
     * \brief Add the parameters to the sums in the second half of the iterations
     * \param[in]    genome The genome
     * \param[inout] state  The state
     * \param[in]    iter   Current iteration number
     */
    void collectStep(const ga::Genome *genome, MutateState *state, int iter);

    /*!
     * \brief Finish the mutation
     * \param[in]    bestGenome The best genome
     * \param[inout] state      The state
     */
    void endMutate(const ga::Genome *bestGenome, MutateState *state);

public:
    //! Set the generation, see ga::Mutator::setGeneration
    void setGeneration(int generation) { myGeneration_ = generation; }


    /*!
     * \brief Constructor
     * \param[in] seed            seed coming from the middle man creating this mutator
     * \param[in] bch             pointer to BayesConfigHandler object
     * \param[in] fitComp         pointer to ACMFitnessComputer object
     * \param[in] sii             pointer to StaticIndividualInfo object
     * \param[in] evaluateTestSet Whether or not to evaluate the test set
     *                            every once in a while
     * \param[in] algorithm       The charge generation algorithm
     * \param[in] maxGenerations  If this is part of the HYBRID training this may be needed for annealing
     */
    MCMCMutator(int                        seed,
                BayesConfigHandler        *bch,
                ACMFitnessComputer        *fitComp,
                StaticIndividualInfo      *sii,
                bool                       evaluateTestSet,
                ChargeGenerationAlgorithm  algorithm,
                int                        maxGenerations);
 
    /*!
     * \brief Run the Markov chain Monte carlo (MCMC) simulation
     * \param[in]  msghandler The message and status handler
     * \param[in]  genome     pointer to the genome
     * \param[out] bestGenome pointer to the best genome
     * \param[in]  prMut      Probability for mutation. 
     *                   Abused as a boolean for evaluating the test set here, if > 0.
     */
    virtual void mutate(MsgHandler                *msghandler,
                        ga::Genome                *genome,
                        ga::Genome                *bestGenome,
                        double                     prMut);

    //! \return whether the test set is evaluated during the steps, which mutateLockstep does not do
    bool evaluateTestSet() const { return evaluateTestSet_; }

    /*!
     * \brief Run MCMC for several genomes at once, each with its own mutator
     * The steps of all genomes are computed at once, see
     * ACMFitnessComputer::calcDeviations. Each mutator draws the same random
     * numbers as mutate, such that the result of each genome is that of
     * mutate, except for where the relaxation of the shells starts from,
     * see ForceComputer::setKeepShells. The mutators must not evaluate the
     * test set.
     * \param[in]    msghandler  The message and status handler
     * \param[in]    mutators    The mutator of each genome, they share the fitness computer
     * \param[inout] genomes     The genomes
     * \param[out]   bestGenomes The best genome of each
     */
    static void mutateLockstep(MsgHandler                       *msghandler,
                               const std::vector<MCMCMutator *> &mutators,
                               const std::vector<ga::Genome *>  &genomes,
                               const std::vector<ga::Genome *>  &bestGenomes);

    //! \return the number of calls to the objective function MCMCMutator::MCMC() routine
    size_t numberObjectiveFunctionCalls() const
    {
        return 1 + bch_->maxIter() * sii_->nParam();
    }

    //! \return whether a minimum was found
    bool foundMinimum() { return bMinimum_; }

    //! \return the number of steps in one mutation
    size_t numberOfSteps() const
    {
        return bch_->maxIter() * sii_->nParam();
    }

    /*! \brief Set how often to write to the convergence files
     * \param[in] nstFitness Number of steps between writing
     */
    void setNstFitness(int nstFitness) { nstFitness_ = nstFitness; }

    /*!
     * \brief Print the MC statistics to a file.
     * \param[in] tw            TextWriter
     * \param[in] initialGenome The initial genome
     * \param[in] bestGenome    The genome with highest fitness
     */
    void printMonteCarloStatistics(gmx::TextWriter  *tw,
                                   const ga::Genome &initialGenome,
                                   const ga::Genome &bestGenome);

    /*!
     * \brief Open parameter convergence files
     * \param[in] oenv the GROMACS output environment
     */
    void openParamConvFiles(const gmx_output_env_t *oenv);

    /*!
     * \brief Open a \f$ \chi^2 \f$ convergence file
     * \param[in] oenv              the GROMACS output environment
     */
    void openChi2ConvFile(const gmx_output_env_t *oenv);
    //! \return all random indices produced
    const std::vector<int> &randIndices() const { return randIndices_; }

    /*! \brief The index of a parameter for a random number
     * Each of the n indices has the same probability, the interval
     * [k/n, (k+1)/n) of the random number gives index k.
     * \param[in] r A random number in [0, 1)
     * \param[in] n The number of parameters, larger than zero
     * \return The index, from 0 to n-1
     */
    static size_t indexOfRandom(double r, size_t n)
    {
        return std::min(n - 1, static_cast<size_t>(std::floor(r*n)));
    }
};


} //namespace alexandria


#endif //ALEXANDRIA_MCMCMUTATOR_H

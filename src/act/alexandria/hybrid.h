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
 */
#ifndef ALEXANDRIA_HYBRID_GAMC_H
#define ALEXANDRIA_HYBRID_GAMC_H

#include <cstdio>

#include <vector>

#include "act/ga/crossover.h"
#include "act/ga/fitness_computer.h"
#include "act/ga/genetic_algorithm.h"
#include "act/ga/initializer.h"
#include "act/ga/probability_computer.h"
#include "act/ga/selector.h"
#include "act/ga/terminator.h"

#include "act/utility/communicationrecord.h"
#include "confighandler.h"
#include "staticindividualinfo.h"

struct gmx_output_env;

namespace ga
{

class MsgHandler;

/*! \brief Genetic algorithm with selection, crossover and mutation,
 * used for both the GA and the HYBRID optimizer. In the HYBRID case
 * the mutation consists of MCMC steps on each individual.
 */
class HybridGAMC : public GeneticAlgorithm
{
private:
    //! Who am I?
    alexandria::StaticIndividualInfo *sii_;
    //! GAConfigHandler pointer
    alexandria::GAConfigHandler      *gach_;
    //! Whether to do adaptive steps
    bool                              adaptiveSteps_;
    //! Output filename for fitness files
    const char                       *fitnessFile_;
    //! Gene pool input (may be null pointer)
    const char                       *gpin_;
    //! Gene pool output
    const char                       *gpout_;
    //! seed for random numbers
    int                               seed_;

    //! The work of a middleman for an individual, see runQueue
    struct QueueTask
    {
        //! The individual in the gene pool
        size_t                           individual;
        //! What to do with it
        alexandria::TrainFFMiddlemanMode mode;
    };

    /*! \brief Hand out the tasks to the middlemen, each gets the next as soon as it is done
     *
     * See CommunicationRecord::queue. Each task gets the seed of its
     * individual, drawn from gen for all individuals in their order, and the
     * generation, such that the results do not depend on which middleman
     * does which task, nor on the number of middlemen.
     * \param[inout] pool       The gene pool, the individuals are updated
     * \param[in]    tasks      The tasks
     * \param[in]    generation The generation, from zero
     * \param[inout] gen        The random number generator
     */
    void runQueue(GenePool                     *pool,
                  const std::vector<QueueTask> &tasks,
                  int                           generation,
                  std::mt19937                 &gen);

    /*! \brief Routine to manage crossovers
     * \param[in]  msghandler For debugging
     * \param[in]  pold       The old population
     * \param[out] pnew       The new population
     * \param[in]  dis        Random number distributor
     * \param[in]  gen        Random number generator
     */
    void doCrossOver(alexandria::MsgHandler                 *msghandler,
                     const GenePool                         *pold,
                     GenePool                               *pnew,
                     std::uniform_real_distribution<double> &dis,
                     std::mt19937                           &gen);

public:
    /*!
     * \brief Constructor for self-building
     */
    HybridGAMC(Initializer                         *initializer,
               FitnessComputer                     *fitnessComputer,
               ProbabilityComputer                 *probComputer,
               Selector                            *selector,
               Crossover                           *crossover,
               Mutator                             *mutator,
               std::vector<Terminator*>            *terminators,
               std::vector<Penalizer*>             *penalizers,
               alexandria::StaticIndividualInfo    *sii,
               alexandria::GAConfigHandler         *gach,
               const char                          *fitnessFileName,
               const char                          *genePoolIn,
               const char                          *genePoolOut,
               int                                  seed,
               bool                                 adaptiveSteps)
    : GeneticAlgorithm(initializer, fitnessComputer, probComputer, selector, crossover,
                       mutator, terminators, penalizers),
      sii_(sii), gach_(gach), adaptiveSteps_(adaptiveSteps), fitnessFile_(fitnessFileName),
      gpin_(genePoolIn), gpout_(genePoolOut), seed_(seed)
    {}

    //! \copydoc ga::GeneticAlgorithm::evolve
    bool evolve(alexandria::MsgHandler       *msghandler,
                std::map<iMolSelect, Genome> *bestGenome) override;

};

} // namespace ga

#endif

﻿/*! \internal \brief
 * Implements part of the alexandria program.
 * \author Mohammad Mehdi Ghahremanpour <mohammad.ghahremanpour@icm.uu.se>
 * \author David van der Spoel <david.vanderspoel@icm.uu.se>
 * \author Julian Ramon Marrades Furquet <julian.marrades@hotmail.es>
 */
#ifndef ALEXANDRIA_ACM_GA_H
#define ALEXANDRIA_ACM_GA_H

#include <cstdio>

#include <vector>

#include "act/ga//Crossover.h"
#include "act/ga//FitnessComputer.h"
#include "act/ga//GeneticAlgorithm.h"
#include "act/ga//Initializer.h"
#include "act/ga//ProbabilityComputer.h"
#include "act/ga//Selector.h"
#include "act/ga//Terminator.h"

#include "confighandler.h"
#include "staticindividualinfo.h"

struct gmx_output_env;

namespace ga
{

class HybridGAMC : public GeneticAlgorithm
{
private:
    //! Who am I?
    alexandria::StaticIndividualInfo *sii_;
    //! GAConfigHandler pointer
    alexandria::GAConfigHandler      *gach_;
    //! logFile
    FILE                             *logFile_;
    //! seed for random numbers
    int                               seed_;
public:
    /*!
     * \brief Constructor for self-building
     */
    HybridGAMC(FILE                                *logFile,
               Initializer                         *initializer,
               FitnessComputer                     *fitnessComputer,
               ProbabilityComputer                 *probComputer,
               Selector                            *selector,
               Crossover                           *crossover,
               Mutator                             *mutator,
               std::vector<Terminator*>             terminators,
               alexandria::StaticIndividualInfo    *sii,
               alexandria::GAConfigHandler         *gach,
               int                                  seed)
    : GeneticAlgorithm(initializer, fitnessComputer, probComputer, selector, crossover,
                       mutator, terminators, gach->popSize()),
      sii_(sii), gach_(gach), logFile_(logFile), seed_(seed) {}

    //! \copydocs ga::GeneticAlgorithm::evolve
    virtual bool evolve(ga::Genome *bestGenome);

};

class MCMC : public GeneticAlgorithm
{
private:
    //! Who am I?
    alexandria::StaticIndividualInfo *sii_;
    //! GAConfigHandler pointer
    alexandria::GAConfigHandler      *gach_;
    //! logFile
    FILE                             *logFile_;
public:
    /*!
     * \brief Constructor for self-building
     */
    MCMC(FILE                                *logFile,
         Initializer                         *initializer,
         FitnessComputer                     *fitnessComputer,
         Mutator                             *mutator,
         alexandria::StaticIndividualInfo    *sii,
         alexandria::GAConfigHandler         *gach)
    : GeneticAlgorithm(initializer, fitnessComputer, nullptr, nullptr, nullptr,
                       mutator, std::vector<Terminator*>(), gach->popSize()),
      sii_(sii), gach_(gach), logFile_(logFile) {}

    //! \copydocs ga::GeneticAlgorithm::evolve
    virtual bool evolve(ga::Genome *bestGenome);

};

} // namespace ga

#endif

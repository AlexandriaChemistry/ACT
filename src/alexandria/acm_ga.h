/*! \internal \brief
 * Implements part of the alexandria program.
 * \author Mohammad Mehdi Ghahremanpour <mohammad.ghahremanpour@icm.uu.se>
 * \author David van der Spoel <david.vanderspoel@icm.uu.se>
 * \author Julian Ramon Marrades Furquet <julian.marrades@hotmail.es>
 */
#ifndef ALEXANDRIA_ACM_GA_H
#define ALEXANDRIA_ACM_GA_H

#include <cstdio>

#include <vector>

#include "ga/Crossover.h"
#include "ga/FitnessComputer.h"
#include "ga/GeneticAlgorithm.h"
#include "ga/Initializer.h"
#include "ga/Mutator.h"
#include "ga/ProbabilityComputer.h"
#include "ga/Selector.h"
#include "ga/Sorter.h"
#include "ga/Terminator.h"

#include "confighandler.h"
#include "staticindividualinfo.h"

struct gmx_output_env;
struct t_commrec;
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
public:
    /*!
     * \brief Constructor for self-building
     */
    HybridGAMC(FILE                                *logFile,
               Initializer                         *initializer,
               FitnessComputer                     *fitnessComputer,
               Sorter                              *sorter,
               ProbabilityComputer                 *probComputer,
               Selector                            *selector,
               Crossover                           *crossover,
               Mutator                             *mutator,
               Terminator                          *terminator,
               alexandria::StaticIndividualInfo    *sii,
               alexandria::GAConfigHandler         *gach)
    : GeneticAlgorithm(initializer, fitnessComputer,
                       sorter, probComputer, selector, crossover, mutator, terminator,
                       gach->popSize()), sii_(sii), gach_(gach), logFile_(logFile) {}
 
    //! \brief Evolve the initial population
    virtual void evolve(ga::Genome *bestGenome);
    
};

class MCMC : public GeneticAlgorithm
{
private:
    //! Who am I?
    alexandria::StaticIndividualInfo *sii_;
    //! logFile
    FILE                             *logFile_;
    //! Should we regularly evaluate the test set?
    bool                               evaluateTestSet_;
public:
    /*!
     * \brief Constructor for self-building
     */
    MCMC(FILE                                *logFile,
         Initializer                         *initializer,
         FitnessComputer                     *fitnessComputer,
         Sorter                              *sorter,
         ProbabilityComputer                 *probComputer,
         Selector                            *selector,
         Crossover                           *crossover,
         Mutator                             *mutator,
         Terminator                          *terminator,
         alexandria::StaticIndividualInfo    *sii,
         alexandria::GAConfigHandler         *gach,
         bool                                 evaluateTestSet)
    : GeneticAlgorithm(initializer, fitnessComputer,
                       sorter, probComputer, selector, crossover, mutator, terminator,
                       gach->popSize()),
      sii_(sii), logFile_(logFile), evaluateTestSet_(evaluateTestSet) {}
    
    //! \brief Evolve the initial population
    virtual void evolve(ga::Genome *bestGenome);
    
};

} // namespace ga

#endif

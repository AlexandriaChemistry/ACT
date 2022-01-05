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

#include "confighandler.h"
#include "ga/Crossover.h"
#include "ga/FitnessComputer.h"
#include "ga/GeneticAlgorithm.h"
#include "ga/Initializer.h"
#include "ga/Mutator.h"
#include "ga/ProbabilityComputer.h"
#include "ga/Selector.h"
#include "ga/Sorter.h"
#include "ga/Terminator.h"

struct gmx_output_env;
struct t_commrec;
namespace ga
{

class PureGA : public GeneticAlgorithm
{
private:
    //! GAConfigHandler pointer
    alexandria::GAConfigHandler *gach_;
public:
    /*!
     * \brief Constructor for self-building
     */
    PureGA(FILE                                *logFile,
           struct gmx_output_env_t             *oenv,
           Initializer                         *initializer,
           FitnessComputer                     *fitnessComputer,
           Sorter                              *sorter,
           ProbabilityComputer                 *probComputer,
           Selector                            *selector,
           Crossover                           *crossover,
           Mutator                             *mutator,
           Terminator                          *terminator,
           alexandria::GAConfigHandler         *gach,
           bool                                 evaluateTestSet,
           const  std::string                   &outputFile)
    : GeneticAlgorithm(logFile, oenv, initializer, fitnessComputer,
                       sorter, probComputer, selector, crossover, mutator, terminator,
                       gach->popSize(), evaluateTestSet, outputFile), gach_(gach) {}
 
    //! \brief Evolve the initial population
    virtual void evolve();
    
};

class HybridGA : public GeneticAlgorithm
{
public:
    /*!
     * \brief Constructor for self-building
     */
    HybridGA(FILE                                *logFile,
             struct gmx_output_env_t             *oenv,
             Initializer                         *initializer,
             FitnessComputer                     *fitnessComputer,
             Sorter                              *sorter,
             ProbabilityComputer                 *probComputer,
             Selector                            *selector,
             Crossover                           *crossover,
             Mutator                             *mutator,
             Terminator                          *terminator,
             alexandria::GAConfigHandler         *gach,
             bool                                 evaluateTestSet,
             const  std::string                  &outputFile)
    : GeneticAlgorithm(logFile, oenv, initializer, fitnessComputer,
                       sorter, probComputer, selector, crossover, mutator, terminator,
                       gach->popSize(), evaluateTestSet, outputFile) {}
    
    //! \brief Evolve the initial population
    virtual void evolve();
    
};

} // namespace ga

#endif

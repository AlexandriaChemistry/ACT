#ifndef GA_GENETICALGORITHM_H
#define GA_GENETICALGORITHM_H


// #include "aliases.h"

#include "Initializer.h"
#include "FitnessComputer.h"
#include "Sorter.h"
#include "ProbabilityComputer.h"
#include "Selector.h"
#include "Crossover.h"
#include "Mutator.h"
#include "Terminator.h"

#include "alexandria/confighandler.h"
#include "alexandria/sharedindividualinfo.h"


namespace ga
{


/*!
 * Class which encapsulates a genetic algorithm
 */
class GeneticAlgorithm
{

private:

    //! BayesConfigHandler pointer
    alexandria::BayesConfigHandler *bch_;
    //! GAConfigHandler pointer
    alexandria::GAConfigHandler *gach_;

    // FIXME: swap to a vector of individual pointers, and decide whether probability
    // lies within the individual or outside of it (in a vector)
    //! Old population
    std::vector<Individual*> oldPop_;
    //! New population, which emerges from the old population
    std::vector<Individual*> newPop_;
    //! Temporal storage to swap "oldPop" and "newPop" after each generation
    std::vector<Individual*> tmpPop_;
    //! Probability of selection for each individual
    std::vector<double> probability_;

    //! Initializes each individual in the population
    Initializer            *initializer_;
    //! Computes fitness for each individual in the population
    FitnessComputer        *fitComputer_;
    //! Sorts the individuals based on their fitness
    Sorter                 *sorter_;
    //! Computes the probability of selection of each individual
    ProbabilityComputer    *probComputer_;
    //! Selects an individual from the population based on its probability
    Selector               *selector_;
    //! Grabs 2 individuals and crosses their genes to generate 2 new individuals
    Crossover              *crossover_;
    //! Mutates the genes of the individuals
    Mutator                *mutator_;
    //! Checks if the evolution should continue or be terminated
    Terminator             *terminator_;

public:

    //! Default constructor
    GeneticAlgorithm() {}

    /*!
     * Constructor
     */
    GeneticAlgorithm(const int mindata_,
                     alexandria::SharedIndividualInfo *sii,
                     const bool randInit_,
                     const std::string) {}

    /*!
     * Evolve the initial population
     * @return          a tuple containing the final population, the final fitness, the best individual, the fitness
     *                  of the best individual, and the number of generations
     */
    const ga_result_t evolve();

};


} //namespace ga


#endif //GA_GENETICALGORITHM_H

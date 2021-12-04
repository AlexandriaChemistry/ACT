#ifndef ACT_GENETICALGORITHM_H
#define ACT_GENETICALGORITHM_H


#include "aliases.h"

#include "Initializer.h"
#include "FitnessComputer.h"
#include "Sorter.h"
#include "ProbabilityComputer.h"
#include "Selector.h"
#include "Crossover.h"
#include "Mutator.h"
#include "Terminator.h"


namespace ga
{


/*!
 * Class which encapsulates a genetic algorithm
 */
class GeneticAlgorithm
{

    //! Amount of individuals in the population
    int popSize;
    //! Amount of genes in each individual
    int chromosomeLength;
    //! Amount of top individuals to be propagated, unchanged, to the next generation
    int nElites;

    //! Old population
    matrix oldPop;
    //! New population, which emerges from the old population
    matrix newPop;
    //! Temporal storage to swap "oldPop" and "newPop" after each generation
    matrix tmpPop;
    //! Fitness score for each individual
    vector fitness;
    //! Probability of selection for each individual
    vector probability;

    //! Initializes each individual in the population
    Initializer            *initializer;
    //! Computes fitness for each individual in the population
    FitnessComputer        *fitComputer;
    //! Sorts the individuals based on their fitness
    Sorter                 *sorter;
    //! Computes the probability of selection of each individual
    ProbabilityComputer    *probComputer;
    //! Selects an individual from the population based on its probability
    Selector               *selector;
    //! Grabs 2 individuals and crosses their genes to generate 2 new individuals
    Crossover              *crossover;
    //! Mutates the genes of the individuals
    Mutator                *mutator;
    //! Checks if the evolution should continue or be terminated
    Terminator             *terminator;

public:
    /*!
     * Create a new GeneticAlgorithm object
     * @param popSize               size of the population
     * @param chromosomeLength      length of each individual
     * @param nElites               the amount of best individuals to move unchanged to the next generation
     * @param initializer           Initializer object
     * @param fitComputer           FitnessComputer object
     * @param sorter                Sorter object
     * @param probComputer          ProbabilityComputer object
     * @param selector              Selector object
     * @param crossover             Crossover object
     * @param mutator               Mutator object
     * @param terminator            Terminator object
     */
    GeneticAlgorithm(const int                  popSize,
                     const int                  chromosomeLength,
                     const int                  nElites,
                           Initializer         *initializer,
                           FitnessComputer     *fitComputer,
                           Sorter              *sorter,
                           ProbabilityComputer *probComputer,
                           Selector            *selector,
                           Crossover           *crossover,
                           Mutator             *mutator,
                           Terminator          *terminator);

    /*!
     * Evolve the initial population
     * @param prCross   the probability of crossover
     * @param prMut     the probability of mutation
     * @param verbose   The higher the value, the more debug prints
     * @return          a tuple containing the final population, the final fitness, the best individual, the fitness
     *                  of the best individual, and the number of generations
     */
    const ga_result_t evolve(const double   prCross,
                             const double   prMut,
                             const int      verbose);

    /*!
     * Set a new value for \p nElites
     * @param nElites   the new value
     */
    void setnElites(const int nElites) { this->nElites = nElites; }

    /*!
     * Get the value of \p nElites
     * @return  the value of \p nElites
     */
    int getnElites() { return nElites; }

};


} //namespace ga


#endif //ACT_GENETICALGORITHM_H

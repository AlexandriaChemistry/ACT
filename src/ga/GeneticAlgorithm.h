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


class GeneticAlgorithm {

    // Scalars
    int popSize;
    int chromosomeLength;

    // Vectors and matrices
    matrix oldPop;
    matrix newPop;
    matrix tmpPop;
    vector fitness;
    vector probability;

    // Object pointers
    Initializer initializer;
    FitnessComputer fitComputer;
    Sorter sorter;
    ProbabilityComputer probComputer;
    Selector selector;
    Crossover crossover;
    Mutator mutator;
    Terminator terminator;

public:
    /*!
     * Create a new GeneticAlgorithm object
     * @param popSize               size of the population
     * @param chromosomeLength      length of each individual
     * @param initializer           Initializer object
     * @param fitComputer           FitnessComputer object
     * @param sorter                Sorter object
     * @param probComputer          ProbabilityComputer object
     * @param selector              Selector object
     * @param crossover             Crossover object
     * @param mutator               Mutator object
     * @param terminator            Terminator object
     */
    GeneticAlgorithm(const int              popSize,
                     const int              chromosomeLength,
                     Initializer            initializer,
                     FitnessComputer        fitComputer,
                     Sorter                 sorter,
                     ProbabilityComputer    probComputer,
                     Selector               selector,
                     Crossover              crossover,
                     Mutator                mutator,
                     Terminator             terminator);

    /*!
     * Evolve the initial population
     * @param prCross   the probability of crossover
     * @param prMut     the probability of mutation
     * @return          a tuple containing the final population, the final fitness, the best individual, the fitness
     *                  of the best individual, and the number of generations
     */
    const ga_result_t evolve(const double   prCross,
                             const double   prMut);

};


#endif //ACT_GENETICALGORITHM_H

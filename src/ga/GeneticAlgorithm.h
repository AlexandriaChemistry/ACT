#ifndef ACT_GENETICALGORITHM_H
#define ACT_GENETICALGORITHM_H

#include "aliases.h"

#include "Initializer.h"
#include "FitnessComputer.h"
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
    double** oldPop;
    double** newPop;
    double** tmpPop;
    double* fitness;
    double* probability;

    // Function pointers
    void (*initialize) (double* const, const int);
    void (*computeFitness) (double* const, double* const, const int);
    void (*computeProbability) (double* const, double* const, const int);
    double* const (*select) (double** const, double* const, const int);
    void (*crossover) (double* const, double* const, double* const, double* const);
    void (*mutate) (double* const);
    bool (*terminate) (double** const, double* const, const int, const int);

public:
    /*!
     * Default constructor for the GeneticAlgorithm class
     * @param popSize               the number of individuals in the population. Assumed to be even
     * @param chromosomeLength      the number of genes in each individuals
     * @param initialize            pointer to a function which takes an individual (1) (and its length (2))
     *                              and initializes it
     * @param computeFitness        pointer to a function which takes an individual (1) (and its length(3)) and
     *                              computes its fitness, writing the result using the given pointer (2)
     * @param computeProbability    pointer to a function which takes the fitness array (1) (and its length (3)) and
     *                              normalizes it writing to the probability array (2)
     * @param select                pointer to a function which takes the population (1) (and its size (3)) and
     *                              the probabilities (2); and returns a pointer to the selected individual
     * @param crossover             pointer to a function which takes two parents (1, 2) and writes the resulting
     *                              crossover to two children (3, 4)
     * @param mutate                pointer to a function which takes a pointer to a gene and mutates it
     * @param terminate             pointer to a function which takes the population (1) (and its size (4)), the
     *                              fitness array (2), and the generation number (3); and decides whether the evolution
     *                              is over
     */
    GeneticAlgorithm(const int popSize,
                     const int chromosomeLength,
                     Initializer initializer,
                     FitnessComputer fitComputer,
                     ProbabilityComputer probComputer,
                     Selector selector,
                     Crossover crossover,
                     Mutator mutator,
                     Terminator terminator);

    /*!
     * Evolve the initial population
     * @param prCross   the probability of crossover
     * @param prMut     the probability of mutation
     * @return          a tuple containing the final population, the final fitness, the best individual, the fitness
     *                  of the best individual, and the number of generations
     */
    const ga_result_t evolve(const double prCross, const double prMut);

};


#endif //ACT_GENETICALGORITHM_H

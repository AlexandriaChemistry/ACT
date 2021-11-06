#ifndef ACT_GENETICALGORITHM_H
#define ACT_GENETICALGORITHM_H

#include <tuple>

#include "Initializer.h"
#include "FitnessComputer.h"
#include "ProbabilityComputer.h"
#include "Selector.h"
#include "Crossover.h"


using ga_result = std::tuple<double** const, double* const, double* const, const double, const int>


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

    // Object pointers
    Initializer initializer;
    FitnessComputer ftComputer;
    ProbabilityComputer probComputer;
    Selector selector;
    Crossover crossover;
    void (*mutate) (double* const);
    bool (*terminate) (double** const, double* const, const int, const int);

public:


    GeneticAlgorithm(const int popSize,
                     const int chromosomeLength,
                     Initializer initializer,
                     FitnessComputer ftComputer,
                     ProbabilityComputer probComputer,
                     Selector selector,
                     Crossover crossover,
                     void (*const mutate) (double* const),
                     bool (*const terminate) (double** const, double* const, const int, const int));

    /*!
     * Evolve the initial population
     * @param prCross   the probability of crossover
     * @param prMut     the probability of mutation
     * @return          a tuple containing the final population, the final fitness, the best individual, the fitness
     *                  of the best individual, and the number of generations
     */
    const ga_result evolve(const double prCross, const double prMut);

};


#endif //ACT_GENETICALGORITHM_H

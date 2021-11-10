#include "GeneticAlgorithm.h"

#include <random>
#include <stdio.h>

#include "aliases.h"

#include "Initializer.h"
#include "FitnessComputer.h"
#include "ProbabilityComputer.h"
#include "Crossover.h"
#include "Mutator.h"
#include "Terminator.h"

#include "helpers.h"


GeneticAlgorithm::GeneticAlgorithm(const int                popSize,
                                   const int                chromosomeLength,
                                   Initializer              initializer,
                                   FitnessComputer          fitComputer,
                                   Sorter                   sorter,
                                   ProbabilityComputer      probComputer,
                                   Selector                 selector,
                                   Crossover                crossover,
                                   Mutator                  mutator,
                                   Terminator               terminator) {

    this->popSize = popSize;
    this->chromosomeLength = chromosomeLength;
    this->initializer = initializer;
    this->fitComputer = fitComputer;
    this->sorter = sorter;
    this->probComputer = probComputer;
    this->selector = selector;
    this->crossover = crossover;
    this->mutator = mutator;
    this->terminator = terminator;

    // Initialize the data structures
    this->oldPop = allocateMatrix(popSize, chromosomeLength);
    this->newPop = allocateMatrix(popSize, chromosomeLength);
    this->fitness = vector(popSize);
    this->probability = vector(popSize);

}


const ga_result_t GeneticAlgorithm::evolve(const double     prCross,
                                           const double     prMut
                                           const bool       verbose) {

    if (verbose) printf("Starting evolution...\n");

    // Random number generation
    std::random_device rd;  // Will be used to obtain a seed for the random number engine
    std::mt19937 gen(rd()); // Standard mersenne_twister_engine seeded with rd()
    std::uniform_real_distribution<> dis(0.0, 1.0);

    // Iteration variables
    int i, j, k;

    // Chromosomes
    vector parent1;
    vector parent2;
//    vector child1;
//    vector child2;

    // Generations
    int generation = 0;
    if (verbose) printf("Generation: %i\n", generation);

    // Initialize the population
    if (verbose) printf("Initializing individuals...\n");
    for (vector ind : oldPop) initializer.initialize(ind, chromosomeLength);

    // Compute fitness
    if (verbose) printf("Computing initial fitness...\n");
    for (i = 0; i < popSize; i++) {
        fitComputer.compute(oldPop[i], &fitness[i], chromosomeLength);
    }

    // If verbose, print best individual
    if (verbose) {
        const int index = findMaximumIndex(fitness, popSize);
        printf("Best individual: ");
        printVector(oldPop[index]);
        printf("Max fitness: %f", fitness[index]);
    }

    // Iterate and create new generations
    do {

        // Increase generation counter
        generation++;
        if (verbose) printf("Generation: %i\n", generation);

        // Sort individuals based on fitness
        if (verbose) printf("Sorting... (if needed)\n");
        sorter.sort(oldPop, fitness, popSize);

        // Normalize the fitness into a probability
        if (verbose) printf("Computing probabilities...\n");
        probComputer.compute(fitness, probability, popSze);

        // Generate new population
        if (verbose) printf("Generating new population...\n");
        for (i = 0; i < popSize; i+=2) {

            // Select parents
            parent1 = selector.select(oldPop, probability, popSize);
            parent2 = selector.select(oldPop, probability, popSize);

            // Do crossover
            if (dis(gen) <= prCross) {
                crossover.offspring(parent1, parent2, newPop[i], newPop[i+1],
                                    chromosomeLength);
            } else {
                copyVectorValues(parent1, newPop[i], 0, chromosomeLength);
                copyVectorValues(parent2, newPop[i+1], 0, chromosomeLength);
            }

            // Do mutation in each child, and compute fitness to avoid another traversal
            for (k = 0; k < 2; k++) {
                for (j = 0; j < chromosomeLength; j++) {
                    if (dis(gen) <= prMut) mutator.mutate(&newPop[i+k][j]);
                }
                // Compute fitness
                fitComputer.compute(newPop[i], &fitness[i], chromosomeLength);
            }

            // Swap oldPop and newPop
            tmpPop = oldPop;
            oldPop = newPop;
            newPop = tmpPop;

            // If verbose, print best individual
            if (verbose) {
                const int index = findMaximumIndex(fitness, popSize);
                printf("Best individual: ");
                printVector(oldPop[index]);
                printf("Max fitness: %f", fitness[index]);
            }

        }

        if (verbose) printf("Checking termination conditions...\n");

    } while(!terminator.terminate(oldPop, fitness, generation, popSize));

    if (verbose) printf("Evolution is done!\n");

    int bestFitIndex = findMaximumIndex(fitness, popSize);

    return {oldPop, fitness, oldPop[bestFitIndex], fitness[bestFitIndex], generation};

}

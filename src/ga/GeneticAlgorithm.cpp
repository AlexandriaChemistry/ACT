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


GeneticAlgorithm::GeneticAlgorithm(const int                 popSize,
                                   const int                 chromosomeLength,
                                   Initializer*              initializer,
                                   FitnessComputer*          fitComputer,
                                   Sorter*                   sorter,
                                   ProbabilityComputer*      probComputer,
                                   Selector*                 selector,
                                   Crossover*                crossover,
                                   Mutator*                  mutator,
                                   Terminator*               terminator) {

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
    oldPop = allocateMatrix(popSize, chromosomeLength);
    newPop = allocateMatrix(popSize, chromosomeLength);
    fitness = vector(popSize);
    probability = vector(popSize);

}


const ga_result_t GeneticAlgorithm::evolve(const double     prCross,
                                           const double     prMut,
                                           const bool       verbose) {

    if (verbose) printf("\nStarting evolution...\n");

    // Random number generation
    std::random_device rd;  // Will be used to obtain a seed for the random number engine
    std::mt19937 gen(rd()); // Standard mersenne_twister_engine seeded with rd()
    std::uniform_real_distribution<> dis(0.0, 1.0);

    // Iteration variables
    int i, j, k;

    // Indices for parents
    int parent1;
    int parent2;

    // Generations
    int generation = 0;
    if (verbose) printf("Generation: %i\n", generation);

    // Initialize the population
    if (verbose) printf("Initializing individuals and compute initial fitness...\n");
    for (i = 0; i < popSize; i++) {
        (*initializer).initialize(oldPop[i], chromosomeLength);
        (*fitComputer).compute(oldPop[i], fitness, i, chromosomeLength);
    }

    // If verbose, print best individual
    if (verbose) {
        const int index = findMaximumIndex(fitness, popSize);
        printf("Best individual: ");
        printVector(oldPop[index]);
        printf("Max fitness: %f\n", fitness[index]);
    }

    // Iterate and create new generations
    do {

        // Increase generation counter
        generation++;
        if (verbose) printf("Generation: %i\n", generation);

        // Sort individuals based on fitness
        if (verbose) printf("Sorting... (if needed)\n");
        (*sorter).sort(oldPop, fitness, popSize);

        // Normalize the fitness into a probability
        if (verbose) printf("Computing probabilities...\n");
        (*probComputer).compute(fitness, probability, popSize);

        // Generate new population
        if (verbose) printf("Generating new population...\n");
        for (i = 0; i < popSize; i+=2) {
            if (verbose) printf("i = %i, %i\n", i, i+1);

            // Select parents
            parent1 = (*selector).select(probability, popSize);
            parent2 = (*selector).select(probability, popSize);

            // Do crossover
            if (dis(gen) <= prCross) {
                printf("Doing crossover...\n");
                (*crossover).offspring(oldPop[parent1], oldPop[parent2], newPop[i],
                                       newPop[i+1], chromosomeLength);
            } else {
                printf("Omitting crossover...\n");
                copyVectorValues(oldPop[parent1], newPop[i], 0, chromosomeLength);
                copyVectorValues(oldPop[parent2], newPop[i+1], 0, chromosomeLength);
            }

            // Do mutation in each child, and compute fitness to avoid another traversal
            if (verbose) printf("Doing mutation...\n");
            for (k = 0; k < 2; k++) {
                for (j = 0; j < chromosomeLength; j++) {
                    if (dis(gen) <= prMut) (*mutator).mutate(newPop[i+k], j);
                }
                // Compute fitness
                (*fitComputer).compute(newPop[i], fitness, i, chromosomeLength);
            }

        }

        // Swap oldPop and newPop
        if (verbose) printf("Swapping oldPop and newPop...\n");
        tmpPop = oldPop;
        oldPop = newPop;
        newPop = tmpPop;

        // If verbose, print best individual
        if (verbose) {
            const int index = findMaximumIndex(fitness, popSize);
            printf("Best individual: ");
            printVector(oldPop[index]);
            printf("Max fitness: %f\n", fitness[index]);
        }

        if (verbose) printf("Checking termination conditions...\n");

    } while(!(*terminator).terminate(oldPop, fitness, generation, popSize, chromosomeLength));

    if (verbose) printf("Evolution is done!\n\n");

    int bestFitIndex = findMaximumIndex(fitness, popSize);

    return {oldPop, fitness, oldPop[bestFitIndex], fitness[bestFitIndex], generation};

}

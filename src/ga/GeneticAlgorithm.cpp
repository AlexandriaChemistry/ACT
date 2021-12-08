#define GEN_PRINTS 500

#include "GeneticAlgorithm.h"

#include <stdio.h>

#include "aliases.h"

#include "Initializer.h"
#include "FitnessComputer.h"
#include "ProbabilityComputer.h"
#include "Crossover.h"
#include "Mutator.h"
#include "Terminator.h"

#include "ga_helpers.h"


namespace ga
{


GeneticAlgorithm::GeneticAlgorithm(const int                    popSize,
                                   const int                    chromosomeLength,
                                   const int                    nElites,
                                         Initializer           *initializer,
                                         FitnessComputer       *fitComputer,
                                         Sorter                *sorter,
                                         ProbabilityComputer   *probComputer,
                                         Selector              *selector,
                                         Crossover             *crossover,
                                         Mutator               *mutator,
                                         Terminator            *terminator)
{

    // TODO: Make sure that there is an even number of individuals in the population
    // TODO: Make sure nElites is even
    // TODO: Make sure an actual sorter is given when nElites > 0
    // TODO: Make sure an actual sorter is given when using RankProbabilityComputer

    this->popSize           = popSize;
    this->chromosomeLength  = chromosomeLength;
    this->nElites           = nElites;
    this->initializer       = initializer;
    this->fitComputer       = fitComputer;
    this->sorter            = sorter;
    this->probComputer      = probComputer;
    this->selector          = selector;
    this->crossover         = crossover;
    this->mutator           = mutator;
    this->terminator        = terminator;

    // Initialize the data structures
    oldPop      = allocateMatrix(popSize, chromosomeLength);
    newPop      = allocateMatrix(popSize, chromosomeLength);
    fitness     = vector(popSize);
    probability = vector(popSize);

}


const ga_result_t GeneticAlgorithm::evolve(const double     prCross,
                                           const double     prMut,
                                           const int        verbose)
{

    if (verbose >= 1) printf("\nStarting evolution...\n");

    // Random number generation
    std::random_device rd;  // Will be used to obtain a seed for the random number engine
    std::mt19937 gen(rd()); // Standard mersenne_twister_engine seeded with rd()
    std::uniform_real_distribution<double> dis(0.0, 1.0);

    // Iteration variables
    int i, k;

    // Indices for parents
    int parent1;
    int parent2;

    // Generations
    int generation = 0;
    if (verbose >= 2 or (verbose >= 1 and generation % GEN_PRINTS == 0))
    {
        printf("\nGeneration: %i\n", generation);
    }

    // Initialize the population and compute fitness
    if (verbose >= 2) printf("Initializing individuals and computing initial fitness...\n");
    for (i = 0; i < popSize; i++)
    {
        (*initializer).initialize(&(oldPop[i]), chromosomeLength);
        (*fitComputer).compute(oldPop[i], &fitness, i, chromosomeLength);
    }

    // If verbose, print best individual
    if (verbose >= 1)
    {
        const int index = findMaximumIndex(fitness, popSize);
        printf("Best individual: ");
        printVector(oldPop[index]);
        printf("Max fitness: %f\n", fitness[index]);
    }
    if (verbose >= 2)
    {
        printf("Population:\n");
        printMatrix(oldPop);

        printf("Fitness vector: ");
        printVector(fitness);
    }

    // Iterate and create new generations
    do
    {

        // Increase generation counter
        generation++;
        if (verbose >= 2 or (verbose >= 1 and generation % GEN_PRINTS == 0))
        {
            printf("\nGeneration: %i\n", generation);
        }

        // Sort individuals based on fitness
        if (verbose >= 2) printf("Sorting... (if needed)\n");
        (*sorter).sort(&oldPop, &fitness, popSize);
        if (verbose >= 2)
        {
            printf("Population after sorting:\n");
            printMatrix(oldPop);
            printf("Fitness vector after sorting: ");
            printVector(fitness);
        }

        // Normalize the fitness into a probability
        if (verbose >= 2) printf("Computing probabilities...\n");
        (*probComputer).compute(fitness, &probability, popSize);
        if (verbose >= 2)
        {
            printf("Probabilities: ");
            printVector(probability);
        }

        // Move the "nElites" best individuals (unchanged) into the new population (assuming population is sorted)
        if (verbose >= 2) printf("Moving the %i best individual(s) into the new population...\n", nElites);
        for (i = 0; i < nElites; i++) newPop[i] = oldPop[i];

        // Generate new population after the elitism
        if (verbose >= 2) printf("Generating the rest of the new population...\n");
        for (i = nElites; i < popSize; i += 2)
        {
            if (verbose >= 3) printf("i = %i, %i\n", i, i + 1);

            // Select parents
            parent1 = (*selector).select(probability, popSize);
            parent2 = (*selector).select(probability, popSize);
            if (verbose >= 3) printf("parent1: %i; parent2: %i\n", parent1, parent2);

            // Do crossover
            if (dis(gen) <= prCross)
            {
                if (verbose >= 3) printf("Doing crossover...\n");
                (*crossover).offspring(oldPop[parent1], oldPop[parent2], &(newPop[i]),
                                        &(newPop[i+1]), chromosomeLength);
            }
            else
            {
                if (verbose >= 3) printf("Omitting crossover...\n");
                newPop[i] = oldPop[parent1];
                newPop[i+1] = oldPop[parent2];
            }

            // Do mutation in each child
            if (verbose >= 3) printf("Doing mutation...\n");
            for (k = 0; k < 2; k++)
            {
                (*mutator).mutate(&(newPop[i + k]), chromosomeLength, prMut);
            }

        }

        // Swap oldPop and newPop
        if (verbose >= 2) printf("Swapping oldPop and newPop...\n");
        tmpPop = oldPop;
        oldPop = newPop;
        newPop = tmpPop;

        // Compute fitness
        for (i = 0; i < popSize; i++) {
            // Compute fitness
            (*fitComputer).compute(oldPop[i], &fitness, i, chromosomeLength);
        }

        // If verbose, print best individual
        if (verbose >= 2 or (verbose >= 1 and generation % GEN_PRINTS == 0))
        {
            const int index = findMaximumIndex(fitness, popSize);
            printf("Best individual: ");
            printVector(oldPop[index]);
            printf("Max fitness: %f\n", fitness[index]);
        }
        if (verbose >= 2)
        {
            printf("Population:\n");
            printMatrix(oldPop);

            printf("Fitness vector: ");
            printVector(fitness);
        }

        if (verbose >= 2) printf("Checking termination conditions...\n");

    } while (!(*terminator).terminate(oldPop, fitness, generation, popSize, chromosomeLength));

    if (verbose >= 1) printf("\nEvolution is done!\n");

    int bestFitIndex = findMaximumIndex(fitness, popSize);

    return {oldPop, fitness, oldPop[bestFitIndex], fitness[bestFitIndex], generation};

}


}  //namespace ga

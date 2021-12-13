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

    popSize_           = popSize;
    chromosomeLength_  = chromosomeLength;
    nElites_           = nElites;
    initializer_       = initializer;
    fitComputer_       = fitComputer;
    sorter_            = sorter;
    probComputer_      = probComputer;
    selector_          = selector;
    crossover_         = crossover;
    mutator_           = mutator;
    terminator_        = terminator;

    // Initialize the data structures
    oldPop_      = allocateMatrix(popSize, chromosomeLength);
    newPop_      = allocateMatrix(popSize, chromosomeLength);
    fitness_     = vector(popSize);
    probability_ = vector(popSize);

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
    for (i = 0; i < popSize_; i++)
    {
        initializer_->initialize(&(oldPop_[i]), chromosomeLength_);
        fitComputer_->compute(oldPop_[i], &fitness_, i, chromosomeLength_);
    }

    // If verbose, print best individual
    if (verbose >= 1)
    {
        const int index = findMaximumIndex(fitness_, popSize_);
        printf("Best individual: ");
        printVector(oldPop_[index]);
        printf("Max fitness: %f\n", fitness_[index]);
    }
    if (verbose >= 2)
    {
        printf("Population:\n");
        printMatrix(oldPop_);

        printf("Fitness vector: ");
        printVector(fitness_);
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
        sorter_->sort(&oldPop_, &fitness_, popSize_);
        if (verbose >= 2)
        {
            printf("Population after sorting:\n");
            printMatrix(oldPop_);
            printf("Fitness vector after sorting: ");
            printVector(fitness_);
        }

        // Normalize the fitness into a probability
        if (verbose >= 2) printf("Computing probabilities...\n");
        probComputer_->compute(fitness_, &probability_, popSize_);
        if (verbose >= 2)
        {
            printf("Probabilities: ");
            printVector(probability_);
        }

        // Move the "nElites" best individuals (unchanged) into the new population (assuming population is sorted)
        if (verbose >= 2) printf("Moving the %i best individual(s) into the new population...\n", nElites_);
        for (i = 0; i < nElites_; i++) newPop_[i] = oldPop_[i];

        // Generate new population after the elitism
        if (verbose >= 2) printf("Generating the rest of the new population...\n");
        for (i = nElites_; i < popSize_; i += 2)
        {
            if (verbose >= 3) printf("i = %i, %i\n", i, i + 1);

            // Select parents
            parent1 = selector_->select(probability_, popSize_);
            parent2 = selector_->select(probability_, popSize_);
            if (verbose >= 3) printf("parent1: %i; parent2: %i\n", parent1, parent2);

            // Do crossover
            if (dis(gen) <= prCross)
            {
                if (verbose >= 3) printf("Doing crossover...\n");
                crossover_->offspring(oldPop_[parent1], oldPop_[parent2], &(newPop_[i]),
                                      &(newPop_[i+1]), chromosomeLength_);
            }
            else
            {
                if (verbose >= 3) printf("Omitting crossover...\n");
                newPop_[i] = oldPop_[parent1];
                newPop_[i+1] = oldPop_[parent2];
            }

            // Do mutation in each child
            if (verbose >= 3) printf("Doing mutation...\n");
            for (k = 0; k < 2; k++)
            {
                mutator_->mutate(&(newPop_[i + k]), prMut);
            }

        }

        // Swap oldPop and newPop
        if (verbose >= 2) printf("Swapping oldPop and newPop...\n");
        tmpPop_ = oldPop_;
        oldPop_ = newPop_;
        newPop_ = tmpPop_;

        // Compute fitness
        for (i = 0; i < popSize_; i++) {
            // Compute fitness
            fitComputer_->compute(oldPop_[i]);
        }

        // If verbose, print best individual
        if (verbose >= 2 or (verbose >= 1 and generation % GEN_PRINTS == 0))
        {
            const int index = findMaximumIndex(fitness_, popSize_);
            printf("Best individual: ");
            printVector(oldPop_[index]);
            printf("Max fitness: %f\n", fitness_[index]);
        }
        if (verbose >= 2)
        {
            printf("Population:\n");
            printMatrix(oldPop_);

            printf("Fitness vector: ");
            printVector(fitness_);
        }

        if (verbose >= 2) printf("Checking termination conditions...\n");

    } while (!terminator_->terminate(oldPop_, fitness_, generation, popSize_, chromosomeLength_));

    if (verbose >= 1) printf("\nEvolution is done!\n");

    int bestFitIndex = findMaximumIndex(fitness_, popSize_);

    return {oldPop_, fitness_, oldPop_[bestFitIndex], fitness_[bestFitIndex], generation};

}


}  //namespace ga

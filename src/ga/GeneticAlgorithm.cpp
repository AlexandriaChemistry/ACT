#include "GeneticAlgorithm.h"

#include "Initializer.h"
#include "FitnessComputer.h"
#include "ProbabilityComputer.h"
#include "Crossover.h"
#include "Mutator.h"
#include "Terminator.h"

#include "helpers.h"

GeneticAlgorithm::GeneticAlgorithm(const int popSize,
                                   const int chromosomeLength,
                                   Initializer initializer,
                                   FitnessComputer fitComputer,
                                   ProbabilityComputer probComputer,
                                   Selector selector,
                                   Crossover crossover,
                                   Mutator mutator,
                                   Terminator terminator) {

    this->popSize = popSize;
    this->chromosomeLength = chromosomeLength;
    this->initializer = initializer;
    this->fitComputer = fitComputer;
    this->probComputer = probComputer;
    this->selector = selector;
    this->crossover = crossover;
    this->mutator = mutator;
    this->terminator = terminator;

    // Initialize the data structures
    this->oldPop = allocateMatrix(popSize, chromosomeLength);
    this->newPop = allocateMatrix(popSize, chromosomeLength);
    this->fitness = allocateArray(popSize);
    this->probability = allocateArray(popSize);

}


const int GeneticAlgorithm::evolve(const double prCross, const double prMut) {

    // Iteration variables
    int i, j, k;

    // Chromosomes
    double* parent1;
    double* parent2;
    double* child1;
    double* child2;

    // Generations
    int generation = 0;

    // Initialize the population
    initializer.initialize(oldPop, popSize);

    // Compute fitness
    for (i = 0; i < popSize; i++) {
        fitComputer.compute(oldPop[i], &fitness[i], chromosomeLength);
    }

    // Iterate and create new generations
    do {

        // Increase generation counter
        generation++;

        // Normalize the fitness into a probability
        probComputer.compute(fitness, probability, popSize);

        // Generate new population
        for (i = 0; i < popSize; i+=2) {

            // Select parents
            parent1 = selector.select(oldPop, probability, popSize);
            parent2 = selector.select(oldPop, probability, popSize);

            // Do crossover
            if (rand01() <= prCross) {
                crossover.offspring(parent1, parent2, newPop[i], newPop[i+1]);
            } else {
                copyArrayValues(parent1, newPop[i], chromosomeLength);
                copyArrayValues(parent2, newPop[i+1], chromosomeLength);
            }

            // Do mutation in each child, and compute fitness to avoid another traversal
            for (k = 0; k < 2; k++) {
                for (j = 0; j < chromosomeLength; j++) {
                    if (rand01() <= prMut) mutator.mutate(&newPop[i+k][j]);
                }
                // Compute fitness
                fitComputer.compute(newPop[i], &fitness[i], chromosomeLength);
            }

            // Swap oldPop and newPop
            tmpPop = oldPop;
            oldPop = newPop;
            newPop = tmpPop;

        }

    } while(!terminator.terminate(oldPop, fitness, generation, popSize));

}

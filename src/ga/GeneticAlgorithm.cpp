#include "GeneticAlgorithm.h"

#include "Initializer.h"
#include "FitnessComputer.h"
#include "ProbabilityComputer.h"
#include "Crossover.h"

#include "helpers.h"

GeneticAlgorithm::GeneticAlgorithm(const int popSize,
                                   const int chromosomeLength,
                                   Initializer initializer,
                                   FitnessComputer ftComputer,
                                   ProbabilityComputer probComputer,
                                   Selector selector,
                                   Crossover crossover,
                                   void (*const mutate)(double *const),
                                   bool (*const terminate)(double **const, double *const, const int, const int)) {

    this->popSize = popSize;
    this->chromosomeLength = chromosomeLength;
    this->initializer = initializer;
    this->ftComputer = ftComputer;
    this->probComputer = probComputer;
    this->selector = selector;
    this->crossover = crossover;
    this->mutate = mutate;
    this->terminate = terminate;

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
        ftComputer.compute(oldPop[i], &fitness[i], chromosomeLength);
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
                    if (rand01() <= prMut) mutate(&newPop[i+k][j]);
                }
                // Compute fitness
                computeFitness(newPop[i], &fitness[i], chromosomeLength);
            }

            // Swap oldPop and newPop
            tmpPop = oldPop;
            oldPop = newPop;
            newPop = tmpPop;

        }

    } while(!terminate(oldPop, fitness, generation, popSize));

}

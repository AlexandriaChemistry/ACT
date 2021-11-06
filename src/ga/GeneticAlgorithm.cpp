#include "GeneticAlgorithm.h"
#include "helpers.h"

GeneticAlgorithm::GeneticAlgorithm(const int popSize, const int chromosomeLength,
                                   void (*const initialize)(double *const, const int),
                                   void (*const computeFitness)(double *const, double *const, const int),
                                   void (*const computeProbability)(double *const, double *const, const int),
                                   double* const (*const select)(double **const, double *const, const int),
                                   void (*const crossover)(double *const, double *const, double *const, double *const),
                                   void (*const mutate)(double *const),
                                   bool (*const terminate)(double **const, double *const, const int, const int)) {

    this->popSize = popSize;
    this->chromosomeLength = chromosomeLength;
    this->initialize = initialize;
    this->computeFitness = computeFitness;
    this->computeProbability = computeProbability;
    this->select = select;
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
    initialize(oldPop, popSize);

    // Compute fitness
    for (i = 0; i < popSize; i++) {
        computeFitness(oldPop[i], &fitness[i], chromosomeLength);
    }

    // Iterate and create new generations
    do {

        // Increase generation counter
        generation++;

        // Normalize the fitness into a probability
        computeProbability(fitness, probability, popSize);

        // Generate new population
        for (i = 0; i < popSize; i+=2) {

            // Select parents
            parent1 = select(oldPop, probability, popSize);
            parent2 = select(oldPop, probability, popSize);

            // Do crossover
            if (rand01() <= prCross) {
                crossover(parent1, parent2, newPop[i], newPop[i+1]);
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

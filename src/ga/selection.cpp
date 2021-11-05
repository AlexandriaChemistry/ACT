#include <time.h>

#include "selection.h"


/*!
 * Select parents from the population with roulette wheel selection.
 *
 * @param parents parent indexes for each new individual. Shape: (2*popSize, gene size)
 * @param population the population
 * @param fitness fitness for the population.
 * @param popSize the number of individuals in the population. Assumed to be even.
 */
void rouletteWheel(double parents[][], double population[][], double fitness[], const int popSize) {
    // Loop variables
    int i;
    // Normalize the fitness (change to LAPACK)
    double normalizedFitness[popSize] = {};
    double sum = 0;
    for (i = 0; i < popSize; i++) sum += fitness[i];
    for (i = 0; i < popSize; i++) normalizedFitness[i] = fitness[i] / sum;
    for (i = 0; i < popSize; i++) {
        parents[2*i] = selectIndFromPropFitness(population, normalizedFitness, popSize)
        parents[2*i+1] = selectIndFromPropFitness(population, normalizedFitness, popSize)
    }
}


/*!
 * Select parents from the population with Boltzmann selection.
 *
 * @param parents parent indexes for each new individual. Shape: (2*popSize, gene size)
 * @param population the population
 * @param fitness fitness for the population.
 * @param popSize the number of individuals in the population. Assumed to be even.
 * @param T temperature
 */
void boltzmann(double parents[][], double population[][], double fitness[], const int popSize, const double T) {
    // Loop variables
    int i;
    // Normalize the fitness (change to LAPACK)
    double normalizedFitness[popSize] = {};
    double sum = 0;
    for (i = 0; i < popSize; i++) sum += exp(fitness[i]/T);
    for (i = 0; i < popSize; i++) normalizedFitness[i] = exp(fitness[i]/T) / sum;
    for (i = 0; i < popSize; i++) {
        parents[2*i] = selectIndFromPropFitness(population, normalizedFitness, popSize)
        parents[2*i+1] = selectIndFromPropFitness(population, normalizedFitness, popSize)
    }
}


/*!
 * Select an individual based on its normalized fitness
 *
 * @param population the population
 * @param normalizedFitness fitness which adds up to 1
 * @param popSize the size of the population
 * @return a pointer to the first element of the selected individual
 */
double* selectIndFromPropFitness(double population[][], double normalizedFitness[], const int popSize) {
    //    srand((unsigned) time(NULL));
    double randNumb = (double) rand()/RAND_MAX;
    int i = 0;
    while (randNumb > 0) {
        randNumb -= normalizedFitness[i];
        i += 1;
    }
    return population[i-1];
}

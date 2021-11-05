#include <time.h>
#include <tuple>

#include "selection.h"

/*!
 * Select an individual from the population with roulette wheel selection.
 *
 * @param population matrix in which each row is an individual.
 * @param normalizedFitness normalized fitness for the population. Adds uo to 1.
 * @param popSize the number of individuals in the population.
 * @return a pointer to the selected individual.
 */
double* rouletteWheel(double population[][], double normalizedFitness[], const int popSize) {

//    srand((unsigned) time(NULL));

    double randNumb = (double) rand()/RAND_MAX;
    int i = 0;

    while (randNumb > 0) {
        randNumb -= normalizedFitness[i];
        i += 1;
    }

    return population[i-1];

}

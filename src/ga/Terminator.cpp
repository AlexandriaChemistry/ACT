#include "Terminator.h"

#include "helpers.h"
#include "aliases.h"


bool SimpleTerminator::terminate(const matrix   population,
                                 const vector&  fitness,
                                 const int      generationNumber,
                                 const int      popSize,
                                 const int      chromosomeLength) {
    double maximumFitness = fitness[findMaximumIndex(fitness, popSize)];
    return maximumFitness >= 1 / (0.01 * chromosomeLength);
}

#ifndef ACT_PROBABILITYCOMPUTER_H
#define ACT_PROBABILITYCOMPUTER_H

#include "aliases.h"

/*!
 * Abstract class for computing the selection probability of each individual in the population
 */
class ProbabilityComputer {

public:
    /*!
     * Compute the probability of each individual in the population
     * @param fitness   the fitness of each individual
     * @param prob      structure to store the probabilities
     * @param popSize   size of the population
     */
    virtual void compute(const vector fitness, const vector prob, const int popSize);

};

#endif //ACT_PROBABILITYCOMPUTER_H

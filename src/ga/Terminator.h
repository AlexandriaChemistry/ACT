#ifndef ACT_TERMINATOR_H
#define ACT_TERMINATOR_H

#include "aliases.h"

/*!
 * Abstract class for termination conditions
 */
class Terminator {

public:
    /*!
     * Check whether the evolution should be terminated
     * @param population            each row is an individual
     * @param fitness               fitness of each individual
     * @param generationNumber      generation number
     * @param popSize               size of the population
     * @return                      true if the evolution is complete, false otherwise
     */
    virtual bool terminate(const matrix population, const vector fitness, const int generationNumber,
                           const int popSize);

};

#endif //ACT_TERMINATOR_H

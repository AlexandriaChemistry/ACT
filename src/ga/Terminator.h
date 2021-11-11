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
     * @param chromosomeLength      length of each individual
     * @return                      true if the evolution is complete, false otherwise
     */
    virtual bool terminate(const matrix     population,
                           const vector&    fitness,
                           const int        generationNumber,
                           const int        popSize,
                           const int        chromosomeLength) { return true; };

};


/*!
 * Toy terminator class which returns true when the sum of squared values in the vector is no larger than
 * tolerance * <chromosomeLength>
 */
class SimpleTerminator : public Terminator {

private:
    double tolerance;

public:
    SimpleTerminator(const double     tolerance);
    bool terminate(const matrix     population,
                   const vector&    fitness,
                   const int        generationNumber,
                   const int        popSize,
                   const int        chromosomeLength);

};

/*!
 * Toy terminator class which returns true when the sum of squared values in the vector is no larger than
 * tolerance * <chromosomeLength>
 */
class GenerationTerminator : public Terminator {

private:
    int maxGenerations;

public:
    GenerationTerminator(const int  maxGenerations);
    bool terminate(const matrix     population,
                   const vector&    fitness,
                   const int        generationNumber,
                   const int        popSize,
                   const int        chromosomeLength);

};

#endif //ACT_TERMINATOR_H

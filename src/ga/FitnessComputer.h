#ifndef ACT_FITNESSCOMPUTER_H
#define ACT_FITNESSCOMPUTER_H

#include "aliases.h"

/*!
 * Abstract class for computing the fitness of an individual
 */
class FitnessComputer {

public:
    /*!
     * Compute the fitness of an individual
     * @param individual    the individual
     * @param fitness       the fitness vector
     * @param indIndex      index of the individual in the population
     * @param length        length of the chromosome
     */
    virtual void compute(const vector&       individual,
                               vector&       fitness,
                         const int           indIndex,
                         const int           length) {};

};

class SimpleFitnessComputer : public FitnessComputer {

public:
    void compute(const vector&       individual,
                       vector&       fitness,
                 const int           indIndex,
                 const int           length);
};

#endif //ACT_FITNESSCOMPUTER_H

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
     * @param ftPtr         pointer to where the fitness should be written
     * @param length        length of the chromosome
     */
    virtual void compute(const vector       individual,
                         double* const      ftPtr,
                         const int          length);

};

class SimpleFitnessComputer : public FitnessComputer {

public:
    void compute(const vector       individual,
                   double* const      ftPtr,
                   const int          length);
};

#endif //ACT_FITNESSCOMPUTER_H

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
    double compute(const vector       individual,
                   double* const      ftPtr,
                   const int          length) {
        double sum = 0;

        for (int i = 0; i < length; i++)
            sum += individual[i]*individual[i];

        return 1 / sum;
    }
};

#endif //ACT_FITNESSCOMPUTER_H

#ifndef ACT_FITNESSCOMPUTER_H
#define ACT_FITNESSCOMPUTER_H


#include "aliases.h"


namespace ga
{


/*!
 * Abstract class for computing the fitness of an individual
 */
class FitnessComputer
{

public:

    /*!
     * Compute the fitness of an individual
     * @param individual    the individual
     * @param fitness       pointer to the fitness vector
     * @param indIndex      index of the individual in the population
     * @param length        length of the chromosome
     */
    virtual void compute(const vector  &individual,
                                vector  *fitness,
                            const int      indIndex,
                            const int      length) = 0;

};


/*!
 * Fitness computer that gives higher fitness to individuals close to the 0 vector.
 */
class SimpleFitnessComputer : public FitnessComputer
{

public:

    /*!
     * Compute the fitness of the \f$i\f$-th individual according to the formula
     * \f[
     *      f_i = \frac{ 1 }{ \sum \limits_{ j=1 }^{ m } { x_j }^2 },
     * \f]
     * where \f$x_j\f$ is the \f$j\f$-th gene of the \f$i\f$-th individual.
     *
     * @param individual    the individual
     * @param fitness       pointer to the fitness vector
     * @param indIndex      index of the individual in the population
     * @param length        length of the chromosome
     */
    virtual void compute(const vector  &individual,
                               vector  *fitness,
                         const int      indIndex,
                         const int      length);
};


} //namespace ga


#endif //ACT_FITNESSCOMPUTER_H

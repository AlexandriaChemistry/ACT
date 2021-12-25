/*! \internal \brief
 * Implements part of the alexandria program.
 * \author Julian Ramon Marrades Furquet <julianramon.marradesfurquet.8049@student.uu.se>
 * \author Oskar Tegby <oskar.tegby@it.uu.se>
 */


#ifndef GA_MUTATOR_H
#define GA_MUTATOR_H


#include <random>
#include <time.h>

#include "Individual.h"


namespace ga
{


/*!
 * Abstract class for gene mutation
 */
class Mutator
{

private:

    std::random_device                      rd_base;
    std::mt19937                            gen_base;
    std::uniform_real_distribution<double>  dis_base;

protected:

    /*!
     * Constructor which initializes the random number generator
     */
    Mutator()
    : gen_base(rd_base()), dis_base(std::uniform_real_distribution<double>(0.0, 1.0))
    {
        gen_base.seed(::time(NULL));
    };

    /*!
     * Obtain a random number in \f$ [0, 1] \f$
     * @returns     a double-precision floating point number in \f$ [0, 1] \f$
     */
    double randNum() const { return dis_base(gen_base); }

public:

    /*!
     * Mutate an individual's genes (in place)
     * @param individual        pointer to the individual to mutate
     * @param prMut             probability of mutating a gene
     */
    virtual void mutate(      Individual   *individual,
                        const double        prMut) = 0;

};


} //namespace ga


#endif //GA_MUTATOR_H

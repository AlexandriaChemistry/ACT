#ifndef GA_CROSSOVER_H
#define GA_CROSSOVER_H


#include <random>
#include <time.h>
#include <vector>

#include "Individual.h"


namespace ga
{


/*!
 * Abstract class to perform crossover
 */
class Crossover
{

private:

    std::random_device                     rd_base;
    std::mt19937                           gen_base;
    std::uniform_int_distribution<size_t>  dis_base;

protected:

    /*!
     * Create a new crossover object.
     * @param chromosomeLength  length of the chromosome
     */
    Crossover(const size_t chromosomeLength)
    : gen_base(rd_base()), dis_base(std::uniform_int_distribution<size_t>(1, chromosomeLength - 1))
    {
        gen_base.seed(::time(NULL));
    }

    /*!
     * Pick a random gene index
     * @returns     an integer representing an index
     */
    size_t randIndex() { return dis_base(gen_base); }

public:

    /*!
     * Perform crossover operation
     * @param parent1   the first parent
     * @param parent2   the second parent
     * @param child1    the first child to write to
     * @param child2    the second child to write to
     * @param length    length of each individual
     */
    virtual void offspring(Individual  *parent1,
                           Individual  *parent2,
                           Individual  *child1,
                           Individual  *child2) = 0;

};


} //namespace ga


#endif //GA_CROSSOVER_H

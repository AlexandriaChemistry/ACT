/*! \internal \brief
 * Implements part of the alexandria program.
 * \author Julian Ramon Marrades Furquet <julian.marrades@hotmail.es>
 */


#ifndef GA_CROSSOVER_H
#define GA_CROSSOVER_H


#include <random>
#include <time.h>
#include <vector>

#include "Genome.h"


namespace ga
{


/*!
 * \brief Abstract class to perform crossover.
 * Given two \ref Individual (parents), it will mix their genomes to generate two new Individual (children)
 */
class Crossover
{

private:

    // Random number stuff
    std::random_device                     rd_base;
    std::mt19937                           gen_base;
    std::uniform_int_distribution<size_t>  dis_base;

protected:

    /*!
     * \brief Constructor
     * \param[in] chromosomeLength  length of the chromosome
     */
    Crossover(const size_t chromosomeLength)
    : gen_base(rd_base()), dis_base(std::uniform_int_distribution<size_t>(1, chromosomeLength - 1))
    {
        gen_base.seed(::time(NULL));
    }

    /*!
     * \brief Pick a random gene index
     * \return the selected index
     */
    size_t randIndex() { return dis_base(gen_base); }

public:

    /*!
     * \brief Perform crossover operation
     * \param[in] parent1   the first parent
     * \param[in] parent2   the second parent
     * \param[in] child1    the first child to write to
     * \param[in] child2    the second child to write to
     */
    virtual void offspring(Genome *parent1,
                           Genome *parent2,
                           Genome *child1,
                           Genome *child2) = 0;

};


} //namespace ga


#endif //GA_CROSSOVER_H

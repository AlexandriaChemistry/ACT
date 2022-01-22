/*! \internal \brief
 * Implements part of the alexandria program.
 * \author Julian Ramon Marrades Furquet <julian.marrades@hotmail.es>
 */


#ifndef GA_SELECTOR_H
#define GA_SELECTOR_H

#include <random>
#include <time.h>
#include <vector>

#include "Genome.h"

namespace ga
{


/*!
 * \brief Abstract class to select an Individual from the population based on its selection probability
 */
class Selector
{

public:

    /*!
     * Select an individual (by index) from the population
     * \param[in] pop   Pointer to the genes
     * \return          the index of the selected individual
     */
    virtual int select(const std::vector<Genome> *pop) = 0;

};


/*!
 * \brief Class for roulette-based selection. Uses cumulative probability to perform selection.
 */
class RouletteSelector : public Selector
{

private:

    // Random number stuff
    std::random_device                      rd;
    std::mt19937                            gen;
    std::uniform_real_distribution<double>  dis;

public:

    //! \brief Constructor
    RouletteSelector()
    : gen(rd()), dis(std::uniform_real_distribution<>(0.0, 1.0))
    {
        gen.seed(::time(NULL));
    }

    virtual int select(const std::vector<Genome> *pop);

};


} //namespace ga


#endif //GA_SELECTOR_H

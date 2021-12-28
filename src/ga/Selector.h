/*! \internal \brief
 * Implements part of the alexandria program.
 * \author Julian Ramon Marrades Furquet <julian.marrades@hotmail.es>
 */


#ifndef GA_SELECTOR_H
#define GA_SELECTOR_H


#include <random>
#include <time.h>
#include <vector>

#include "Individual.h"


namespace ga
{


/*!
 * Abstract class to select an individual from the population based on its selection probability
 */
class Selector
{

public:

    /*!
     * Select an individual (by index) from the population
     * @param pop               the population
     * @return                  the index of the selected individual
     */
    virtual int select(const std::vector<Individual*> &pop) = 0;

};


/*!
 * Class for roulette-based selection. Uses cumulative probability to perform selection.
 */
class RouletteSelector : public Selector
{

private:

    std::random_device                      rd;
    std::mt19937                            gen;
    std::uniform_real_distribution<double>  dis;

public:

    /*!
     * Create a new instance of RouletteSelector.
     */
    RouletteSelector()
    : gen(rd()), dis(std::uniform_real_distribution<>(0.0, 1.0))
    {
        gen.seed(::time(NULL));
    }

    virtual int select(const std::vector<Individual*> &pop);

};


} //namespace ga


#endif //GA_SELECTOR_H

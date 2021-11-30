#ifndef ACT_SELECTOR_H
#define ACT_SELECTOR_H

#include "aliases.h"

#include <random>
#include <time.h>


namespace ga
{

    /*!
     * Abstract class to select an individual from the population based on its selection probability
     */
    class Selector
    {

    public:
        /*!
         * Select an individual from the population
         * @param probability       selection probability of each individual
         * @param popSize           size of the population
         * @return                  the selected individual
         */
        virtual int select(const vector  &probability,
                           const int      popSize) = 0;

    };


/*!
 * Class for roulette-based selection. Uses cumulative probability to perform selection.
 */
    class RouletteSelector : public Selector
    {

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

        int select(const vector  &probability,
                   const int      popSize);
    };

}


#endif //ACT_SELECTOR_H

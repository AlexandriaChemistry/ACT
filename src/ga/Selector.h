#ifndef ACT_SELECTOR_H
#define ACT_SELECTOR_H

#include "aliases.h"

#include <random>
#include <time.h>

/*!
 * Abstract class to select an individual from the population
 */
class Selector {

public:
    /*!
     * Select an individual from the population
     * @param probability       probability of each individual
     * @param popSize           size of the population
     * @return                  the selected individual
     */
    virtual const int select(const vector&    probability,
                             const int        popSize) { return 0; };

};


/*!
 * Class for roulette-based selection
 */
 class RouletteSelector: public Selector {

     std::random_device                     rd;
     std::mt19937                           gen;
     std::uniform_real_distribution<double> dis;

 public:
     /*!
      * Create a new instance of RouletteSelector
      */
     RouletteSelector()
     : gen(rd()), dis(std::uniform_real_distribution<>(0.0, 1.0)) {
         gen.seed(::time(NULL));
     }

     const int select(const vector&   probability,
                      const int       popSize);
 };


#endif //ACT_SELECTOR_H

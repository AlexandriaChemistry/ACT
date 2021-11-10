#ifndef ACT_SELECTOR_H
#define ACT_SELECTOR_H

#include "aliases.h"

#include <random>

/*!
 * Abstract class to select an individual from the population
 */
class Selector {

public:
    /*!
     * Select an individual from the population
     * @param population        each row is an individual
     * @param probability       probability of each individual
     * @param popSize           size of the population
     * @return                  the selected individual
     */
    virtual const vector select(const matrix    population,
                                const vector    probability,
                                const int       popSize) { return {}; };

};


/*!
 * Class for roulette-based selection
 */
 class RouletteSelector: public Selector {

     std::random_device rd;  // Will be used to obtain a seed for the random number engine
     std::mt19937 gen;
     std::uniform_real_distribution<> dis;

 public:
     /*!
      * Create a new instance of RouletteSelector
      */
     RouletteSelector();

     const vector select(const matrix   population,
                         const vector   probability,
                         const int      popSize);
 };


#endif //ACT_SELECTOR_H

#ifndef ACT_CROSSOVER_H
#define ACT_CROSSOVER_H

#include <random>

#include "aliases.h"


/*!
 * Abstract class to perform crossover
 */
class Crossover {

protected:
    std::random_device rd;  // Will be used to obtain a seed for the random number engine
    std::mt19937 gen; // Standard mersenne_twister_engine seeded with rd()
    std::uniform_int_distribution<> dis;

public:
    /*!
     * Create a new crossover object.
     * @param chromosomeLength  length of the chromosome
     */
    Crossover(const int chromosomeLength);

    /*!
     * Perform crossover operation
     * @param parent1   the first parent
     * @param parent2   the second parent
     * @param child1    the first child to write to
     * @param child2    the second child to write to
     * @param length    length of each individual
     */
    virtual void offspring(const vector     parent1,
                           const vector     parent2,
                           const vector     child1,
                           const vector     child2,
                           const int        length);

};


/*!
 * Class for single-point crossover operation.
 */
class SinglePointCrossover : public Crossover {

public:
    SinglePointCrossover(const int chromosomeLength): Crossover(chromosomeLength);
    void offspring(const vector     parent1,
                   const vector     parent2,
                   const vector     child1,
                   const vector     child2,
                   const int        length);

};


/*!
 * Class for double-point crossover operation.
 */
class DoublePointCrossover : public Crossover {

public:
    DoublePointCrossover(const int chromosomeLength): Crossover(chromosomeLength);
    void offspring(const vector     parent1,
                   const vector     parent2,
                   const vector     child1,
                   const vector     child2,
                   const int        length);

};

#endif //ACT_CROSSOVER_H

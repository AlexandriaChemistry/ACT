#ifndef ACT_CROSSOVER_H
#define ACT_CROSSOVER_H


#include <random>
#include <time.h>

#include "aliases.h"


namespace ga
{


/*!
 * Abstract class to perform crossover
 */
class Crossover
{

private:

    std::random_device                  rd;
    std::mt19937                        gen;
    std::uniform_int_distribution<int>  dis;

protected:

    /*!
     * Create a new crossover object.
     * @param chromosomeLength  length of the chromosome
     */
    Crossover(const int chromosomeLength)
    : gen(rd()), dis(std::uniform_int_distribution<>(1, chromosomeLength - 1))
    {
        gen.seed(::time(NULL));
    }

    /*!
     * Pick a random gene index
     * @returns     an integer representing an index
     */
    int randIndex();

public:

    /*!
     * Perform crossover operation
     * @param parent1   the first parent
     * @param parent2   the second parent
     * @param child1    the first child to write to
     * @param child2    the second child to write to
     * @param length    length of each individual
     */
    virtual void offspring(const vector &parent1,
                           const vector &parent2,
                                 vector *child1,
                                 vector *child2,
                           const int     length) = 0;

};


/*!
 * Class for the single-point crossover operation.
 */
class SinglePointCrossover : public Crossover
{

public:

    /*!
     * Create a new SinglePointCrossover object
     * @param chromosomeLength  amount of genes in each individual
     */
    SinglePointCrossover(const int chromosomeLength)
    : Crossover(chromosomeLength) {}

    virtual void offspring(const vector &parent1,
                           const vector &parent2,
                                 vector *child1,
                                 vector *child2,
                           const int     length);

};


/*!
 * Class for the double-point crossover operation.
 */
class DoublePointCrossover : public Crossover
{

public:

    /*!
     * Create a new DoublePointCrossover object
     * @param chromosomeLength  amount of genes in each individual
     */
    DoublePointCrossover(const int chromosomeLength)
    : Crossover(chromosomeLength) {};

    virtual void offspring(const vector &parent1,
                           const vector &parent2,
                                 vector *child1,
                                 vector *child2,
                           const int     length);

};


/*!
 * Class for the n-point crossover operation.
 */
class NPointCrossover : public Crossover
{

private:

    int numberOfCrossovers;
    vector crossoverIndices;

public:

    /*!
     * Create a new NPointCrossover object
     * @param chromosomeLength       amount of genes in each individual
     * @param numberOfCrossovers     the number of places the genes swap
     */
    NPointCrossover(const int chromosomeLength,
                    const int numberOfCrossovers)
    : Crossover(chromosomeLength)
    {
        this->numberOfCrossovers = numberOfCrossovers;
    };

    virtual void offspring(const vector &parent1,
                           const vector &parent2,
                                 vector *child1,
                                 vector *child2,
                           const int     length);

};


} //namespace ga


#endif //ACT_CROSSOVER_H

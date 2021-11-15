#ifndef ACT_CROSSOVER_H
#define ACT_CROSSOVER_H

#include <random>
#include <time.h>

#include "aliases.h"

namespace ga {

    /*!
     * Abstract class to perform crossover
     */
    class Crossover {

        std::random_device                  rd;
        std::mt19937                        gen;
        std::uniform_int_distribution<int>  dis;

    public:
        /*!
         * Create a new crossover object.
         * @param chromosomeLength  length of the chromosome
         */
        Crossover(const int chromosomeLength)
        : gen(rd()), dis(std::uniform_int_distribution<>(1, chromosomeLength - 2)) {
            gen.seed(::time(NULL));
        }

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
                               vector &child1,
                               vector &child2,
                               const int length) {};

        /*!
         * Return random index
         */
        int randIndex();

    };


    /*!
     * Class for single-point crossover operation.
     */
    class SinglePointCrossover : public Crossover {

    public:
        SinglePointCrossover(const int chromosomeLength) : Crossover(chromosomeLength) {}

        void offspring(const vector &parent1,
                       const vector &parent2,
                       vector &child1,
                       vector &child2,
                       const int length);

    };


    /*!
     * Class for double-point crossover operation.
     */
    class DoublePointCrossover : public Crossover {

    public:
        DoublePointCrossover(const int chromosomeLength) : Crossover(chromosomeLength) {};

        void offspring(const vector &parent1,
                       const vector &parent2,
                       vector &child1,
                       vector &child2,
                       const int length);

    };

}

#endif //ACT_CROSSOVER_H

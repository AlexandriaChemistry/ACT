#ifndef ACT_TERMINATOR_H
#define ACT_TERMINATOR_H

#include "aliases.h"


namespace ga
{

    /*!
     * Abstract class to check for evolution termination conditions
     */
    class Terminator
    {

    public:
        /*!
         * Check whether the evolution should be terminated
         * @param population            each row is an individual
         * @param fitness               fitness of each individual
         * @param generationNumber      generation number
         * @param popSize               size of the population
         * @param chromosomeLength      length of each individual
         * @return                      true if the evolution is complete, false otherwise
         */
        virtual bool terminate(const matrix    &population,
                               const vector    &fitness,
                               const int        generationNumber,
                               const int        popSize,
                               const int        chromosomeLength) { return true; }

    };


    /*!
     * Toy terminator class which returns true when the genes in the best individual are sufficiently close to 0
     */
    class SimpleTerminator : public Terminator
    {

        double tolerance;

    public:
        /*!
         * Create a new SimpleTerminator object
         * @param tolerance     the tolerance value
         */
        SimpleTerminator(const double tolerance);

        /*!
         * Check whether the evolution should be terminated.<br>
         * Assume the \f$i\f$-th individual has the highest fitness. We will terminate when
         * \f[
         *      f_i \geq \frac{ 1 }{ tolerance \cdot m },
         * \f]
         * where \f$m\f$ is the amount of genes in each individual.
         *
         * @param population            each row is an individual
         * @param fitness               fitness of each individual
         * @param generationNumber      generation number
         * @param popSize               size of the population
         * @param chromosomeLength      length of each individual
         * @return                      true if the evolution is complete, false otherwise
         */
        bool terminate(const matrix    &population,
                       const vector    &fitness,
                       const int        generationNumber,
                       const int        popSize,
                       const int        chromosomeLength);

    };

    /*!
     * Terminator which stops evolution after a given amount of generations.
     */
    class GenerationTerminator : public Terminator
    {

        int maxGenerations;

    public:
        /*!
         * Create a new GenerationTerminator object
         * @param maxGenerations    the maximum amount of generations
         */
        GenerationTerminator(const int maxGenerations);

        /*!
         * Will return true when \p generationNumber \f$\geq\f$ \p maxGenerations, and false otherwise.
         * @param population            each row is an individual
         * @param fitness               fitness of each individual
         * @param generationNumber      generation number
         * @param popSize               size of the population
         * @param chromosomeLength      length of each individual
         * @return                      true if the evolution is complete, false otherwise
         */
        bool terminate(const matrix    &population,
                       const vector    &fitness,
                       const int        generationNumber,
                       const int        popSize,
                       const int        chromosomeLength);

    };

}

#endif //ACT_TERMINATOR_H

#ifndef ACT_MUTATOR_H
#define ACT_MUTATOR_H

#include "aliases.h"

#include <random>
#include <time.h>


namespace ga
{

    /*!
     * Abstract class for gene mutation
     */
    class Mutator
    {

        std::random_device                      rd_base;
        std::mt19937                            gen_base;
        std::uniform_real_distribution<double>  dis_base;

    protected:
        /*!
         * Constructor which initializes the random number generator
         */
        Mutator()
        : gen_base(rd_base()), dis_base(std::uniform_real_distribution<>(0.0, 1.0))
        {
            gen_base.seed(::time(NULL));
        };

        /*!
         * Obtain a random number in \f$ [0, 1] \f$
         * @returns     a double-precision floating point number in \f$ [0, 1] \f$
         */
        double randNum() { return dis_base(gen_base); }

    public:

        /*!
         * Mutate an individual's genes (in place)
         * @param individual        pointer to the individual to mutate
         * @param chromosomeLength  number of genes in an individual
         * @param prMut             probability of mutating a gene
         */
        virtual void mutate(      vector   *individual,
                            const int       chromosomeLength,
                            const double    prMut) {};

    };

    /*!
     * Class for percentual gene mutation
     * Modifies a gene by a maximum of (\p frac*100)% of its value
     */
    class PercentMutator : public Mutator
    {

        std::random_device                      rd;
        std::mt19937                            gen;
        std::uniform_real_distribution<double>  dis;

    public:
        /*!
         * Create a new PercentMutation object
         * @param frac  the fraction of change in [0, 1]
         */
        PercentMutator(const double frac)
        : Mutator(), gen(rd()), dis(std::uniform_real_distribution<>(1 - frac, 1 + frac))
        {
            gen.seed(::time(NULL));
        }

        /*!
         * Mutate genes (in place) by multiplying them by a random number in \f$[1-frac, 1+frac]\f$
         * @param individual        pointer to the individual to mutate
         * @param chromosomeLength  number of genes in an individual
         * @param prMut             probability of mutating a gene
         */
        void mutate(      vector   *individual,
                    const int       chromosomeLength,
                    const double    prMut);

    };


    /*!
     * Class for range gene mutation
     * Modifies a gene by a maximum of 2 * \p range.
     */
    class RangeMutator : public Mutator
    {

        std::random_device                      rd;
        std::mt19937                            gen;
        std::uniform_real_distribution<double>  dis;

    public:
        /*!
         * Create a new RangeMutator object
         * @param range  the range for maximum change
         */
        RangeMutator(const double range)
        : Mutator(), gen(rd()), dis(std::uniform_real_distribution<>(-range, range))
        {
            gen.seed(::time(NULL));
        }

        /*!
         * Mutate genes (in place) by generating a random number in \f$[-range, range]\f$ and adding it to the genes.
         * @param individual        pointer to the individual to mutate
         * @param chromosomeLength  number of genes in an individual
         * @param prMut             probability of mutating a gene
         */
        void mutate(      vector   *individual,
                    const int       chromosomeLength,
                    const double    prMut);

    };

}

#endif //ACT_MUTATOR_H

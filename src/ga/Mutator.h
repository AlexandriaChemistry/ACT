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

    public:
        /*!
         * Mutate a gene (in place)
         * @param individual    pointer to the individual to mutate
         * @param indGene       index of the gene to alter
         */
        virtual void mutate(      vector   *individual,
                            const int       indGene) {};

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
        : gen(rd()), dis(std::uniform_real_distribution<>(1 - frac, 1 + frac))
        {
            gen.seed(::time(NULL));
        }

        /*!
         * Mutate a gene (in place) by multiplying it by a random number in \f$[1-frac, 1+frac]\f$
         * @param individual    pointer to individual to mutate
         * @param indGene       index of the gene to alter
         */
        void mutate(      vector   *individual,
                    const int       indGene);

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
        : gen(rd()), dis(std::uniform_real_distribution<>(-range, range))
        {
            gen.seed(::time(NULL));
        }

        /*!
         * Mutate a gene (in place) by generating a random number in \f$[-range, range]\f$ and adding it to the gene.
         * @param individual    pointer to individual to mutate
         * @param indGene       index of the gene to alter
         */
        void mutate(      vector   *individual,
                    const int       indGene);

    };

}

#endif //ACT_MUTATOR_H

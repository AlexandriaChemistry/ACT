#ifndef ACT_MUTATOR_H
#define ACT_MUTATOR_H

#include "aliases.h"

#include <random>
#include <time.h>

/*!
 * Abstract class for gene mutation
 */
class Mutator {

public:
    /*!
     * Mutate a gene
     * @param individual    individual to mutate
     * @param indGene       index of the gene to alter
     */
    virtual void mutate(      vector&   individual,
                        const int       indGene) {};

};

/*!
 * Class for percentual gene mutation
 * Modifies a gene by a maximum of <frac*100>% of its value
 */
class PercentMutator:public Mutator {

    std::random_device rd;  // Will be used to obtain a seed for the random number engine
    std::mt19937 gen; // Standard mersenne_twister_engine seeded with rd()
    std::uniform_real_distribution<double> dis;

public:
    /*!
     * Create a new PercentMutation object
     * @param frac  the fraction of change in [0, 1]
     */
    PercentMutator(const double frac)
    : gen(rd()), dis(std::uniform_real_distribution<>(1-frac, 1+frac)) {
        gen.seed(::time(NULL));
    }

    void mutate(      vector&   individual,
                const int       indGene);

};

#endif //ACT_MUTATOR_H

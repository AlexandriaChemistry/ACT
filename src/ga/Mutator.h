#ifndef ACT_MUTATOR_H
#define ACT_MUTATOR_H

#include <random>

/*!
 * Abstract class for gene mutation
 */
class Mutator {

public:
    /*!
     * Mutate a gene
     * @param gene  pointer to the gene
     */
    virtual void mutate(double* const gene);

};

/*!
 * Class for percentual gene mutation
 * Modifies a gene by a maximum of <frac*100>% of its value
 */
class PercentMutator:public Mutator {

    std::random_device rd;  // Will be used to obtain a seed for the random number engine
    std::mt19937 gen; // Standard mersenne_twister_engine seeded with rd()
    std::uniform_real_distribution<> dis;

public:
    /*!
     * Create a new PercentMutation object
     * @param frac  the fraction of change in [0, 1]
     */
    PercentMutator(const double frac);
    void mutate(double* const gene);

};

#endif //ACT_MUTATOR_H

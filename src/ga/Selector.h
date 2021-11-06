#ifndef ACT_SELECTOR_H
#define ACT_SELECTOR_H

/*!
 * Abstract class to select an individual from the population
 */
class Selector {

public:
    /*!
     * Select an individual from the population
     * @param population        each row is an individual
     * @param probabilities     probability of each individual
     * @param popSize           size of the population
     * @return                  a pointer to the selected individual
     */
    virtual double* const select(double** const population, double* const probabilities, const int popSize);

};


#endif //ACT_SELECTOR_H

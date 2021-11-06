#ifndef ACT_CROSSOVER_H
#define ACT_CROSSOVER_H

/*!
 * Abstract class to perform crossover
 */
class Crossover {

public:
    /*!
     * Perform crossover operation
     * @param parent1   the first parent
     * @param parent2   the second parent
     * @param child1    the first child to write to
     * @param child2    the second child to write to
     */
    virtual void offspring(double* const parent1, double* const parent2, double* const child1, double* const child2);

};

#endif //ACT_CROSSOVER_H

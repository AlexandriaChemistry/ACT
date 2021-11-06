#ifndef ACT_MUTATOR_H
#define ACT_MUTATOR_H

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

#endif //ACT_MUTATOR_H

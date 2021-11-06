#ifndef ACT_INITIALIZER_H
#define ACT_INITIALIZER_H

/*!
 * Abstract class for initializing individuals
 */
class Initializer {

public:
    /*!
     * Initialize an individual
     * @param individual    the individual to initialize
     * @param length        length of the chromosome
     */
    virtual void initialize(double* const individual, const int length);

};

#endif //ACT_INITIALIZER_H

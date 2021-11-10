#ifndef ACT_INITIALIZER_H
#define ACT_INITIALIZER_H

#include "aliases.h"

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
    virtual void initialize(const vector    individual,
                            const int       length);

};

#endif //ACT_INITIALIZER_H

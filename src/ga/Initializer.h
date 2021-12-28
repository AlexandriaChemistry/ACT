/*! \internal \brief
 * Implements part of the alexandria program.
 * \author Julian Ramon Marrades Furquet <julian.marrades@hotmail.es>
 */


#ifndef GA_INITIALIZER_H
#define GA_INITIALIZER_H


#include "Individual.h"

#include <time.h>
#include <random>


namespace ga
{


/*!
* Abstract class for initializing individuals
*/
class Initializer
{

public:

    /*!
     * Initialize an individual
     * @param individual    pointer to pointer to the individual to initialize
     */
    virtual void initialize(Individual **individual) = 0;
};


} //namespace ga


#endif //GA_INITIALIZER_H

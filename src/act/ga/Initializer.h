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
* \brief Abstract class for initializing Individual instances
*/
class Initializer
{

public:

    /*!
     * \brief Initialize an Individual
     * \return   pointer to a new individual
     */
    virtual Individual *initialize() = 0;
};


} //namespace ga


#endif //GA_INITIALIZER_H

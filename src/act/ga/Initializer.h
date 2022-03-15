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
* \brief Abstract class for initializing and randomizing Individual/Genome instances
*/
class Initializer
{

public:

    /*!
     * \brief Initialize an Individual
     * \return   pointer to a new individual
     */
    virtual Individual *initialize() = 0;

    /*!
     * \brief Randomize a Genome object
     * \param[in] genome the Genome to randomize
     */
    virtual void randomizeGenome(Genome *genome) = 0;
};


} //namespace ga


#endif //GA_INITIALIZER_H

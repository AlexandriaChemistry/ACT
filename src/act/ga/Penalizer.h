#ifndef GA_PENALIZER_H
#define GA_PENALIZER_H

/*! \internal \brief
 * Implements part of the alexandria program.
 * \author Julian Ramon Marrades Furquet <julian.marrades@hotmail.es>
 */

#include "GenePool.h"

namespace ga
{

/*!
 * \brief Penalizes a population
 */
class Penalizer
{

protected:

    //! \brief Default constructor
    Penalizer() {}

public:

    /*!
     * \brief Penalize a population
     * \param[inout] pool       the collection of genomes
     * \param[in]    generation the current generation number
     */
    virtual void penalize(      GenePool *pool,
                          const int       generation) = 0;

}; // class Penalizer

} // namespace ga


#endif // GA_PENALIZER_H
/*! \internal \brief
 * Implements part of the alexandria program.
 * \author Julian Ramon Marrades Furquet <julian.marrades@hotmail.es>
 */
#ifndef GA_INDIVIDUAL_H
#define GA_INDIVIDUAL_H

#include <cstdio>

#include <map>

#include "Dataset.h"
#include "Genome.h"

namespace ga
{

//! \brief Abstract individual for genetic algorithms
class Individual
{

protected:

    //! \brief Default constructor
    Individual() {}

public:

    /* * * * * * * * * * * * * * * * * * * * * *
    * BEGIN: Cloning                           *
    * * * * * * * * * * * * * * * * * * * * * */

    //! \return a copy of this Individual
    virtual Individual *clone() = 0;

    /*!
     * \brief Copy the genome from another Individual
     * \param[in] other pointer to another Individual
     */
    virtual void copyGenome(Genome *other) = 0;

    /* * * * * * * * * * * * * * * * * * * * * *
    * END: Cloning                             *
    * * * * * * * * * * * * * * * * * * * * * */
};


} //namespace ga


#endif //GA_INDIVIDUAL_H

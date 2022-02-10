/*! \internal \brief
 * Implements part of the alexandria program.
 * \author Julian Ramon Marrades Furquet <julian.marrades@hotmail.es>
 */


#include "Terminator.h"

#include "gromacs/utility/basedefinitions.h"

#include "GenePool.h"

namespace ga
{


/* * * * * * * * * * * * * * * * * * * * * *
* BEGIN: GenerationTerminator              *
* * * * * * * * * * * * * * * * * * * * * */

bool GenerationTerminator::terminate(gmx_unused const GenePool *pool,
                                                const int       generationNumber)
{
    return generationNumber >= maxGenerations_;
}

/* * * * * * * * * * * * * * * * * * * * * *
* END: GenerationTerminator                *
* * * * * * * * * * * * * * * * * * * * * */


} //namespace ga

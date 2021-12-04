#include "Terminator.h"

#include "ga_helpers.h"
#include "aliases.h"


namespace ga
{


/* * * * * * * * * * * * * * * * * * * * * *
* BEGIN: SimpleTerminator                  *
* * * * * * * * * * * * * * * * * * * * * */

bool SimpleTerminator::terminate(const matrix  &population,
                                 const vector  &fitness,
                                 const int      generationNumber,
                                 const int      popSize,
                                 const int      chromosomeLength)
{
    double maximumFitness = fitness[findMaximumIndex(fitness, popSize)];
    return maximumFitness >= 1 / (tolerance * chromosomeLength);
}

/* * * * * * * * * * * * * * * * * * * * * *
* END: SimpleTerminator                    *
* * * * * * * * * * * * * * * * * * * * * */

/* * * * * * * * * * * * * * * * * * * * * *
* BEGIN: GenerationTerminator              *
* * * * * * * * * * * * * * * * * * * * * */

bool GenerationTerminator::terminate(const matrix  &population,
                                     const vector  &fitness,
                                     const int      generationNumber,
                                     const int      popSize,
                                     const int      chromosomeLength)
{
    return generationNumber >= maxGenerations;
}

/* * * * * * * * * * * * * * * * * * * * * *
* END: GenerationTerminator                *
* * * * * * * * * * * * * * * * * * * * * */


} //namespace ga

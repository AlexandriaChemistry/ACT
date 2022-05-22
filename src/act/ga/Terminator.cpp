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
    if (generationNumber >= maxGenerations_)
    {
        if (outfile_)
        {
            fprintf(
                outfile_,
                "GenerationTerminator: evolution will be terminated as the maximum number of generations (%d) has been reached.\n",
                maxGenerations_
            );
        }
        return true;
    }
    else
    {
        return false;
    }
}

/* * * * * * * * * * * * * * * * * * * * * *
* END: GenerationTerminator                *
* * * * * * * * * * * * * * * * * * * * * */

/* * * * * * * * * * * * * * * * * * * * * *
* BEGIN: TestGenTerminator                 *
* * * * * * * * * * * * * * * * * * * * * */

bool TestGenTerminator::terminate(           const GenePool *pool,
                                  gmx_unused const int       generationNumber)
{
    remaining_ -= 1;
    
    // Get best genome in population
    const auto bestGenome = pool->getBest(iMolSelect::Train);
    // Check if its test fitness is better than the best one so far
    const double tmpTestFit = bestGenome.fitness(iMolSelect::Test);
    if (tmpTestFit < bestFitness_)
    {
        bestFitness_ = tmpTestFit;
        remaining_ = generations_;
    }

    if (remaining_ <= 0)
    {
        if (outfile_)
        {
            fprintf(
                outfile_,
                "TestGenTerminator: evolution will be terminated as the best test fitness (%lf) has not improved in the last (%d) generation(s).\n",
                bestFitness_,
                generations_
            );
        }
        return true;
    }
    else
    {
        return false;
    }
}

/* * * * * * * * * * * * * * * * * * * * * *
* END: TestGenTerminator                   *
* * * * * * * * * * * * * * * * * * * * * */


} //namespace ga

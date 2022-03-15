/*! \internal \brief
 * Implements part of the alexandria program.
 * \author Julian Ramon Marrades Furquet <julian.marrades@hotmail.es>
 */

#include "Penalizer.h"

namespace ga
{

/* * * * * * * * * * * * * * * * * * * * * *
* BEGIN: VolumeFractionPenalizer           *
* * * * * * * * * * * * * * * * * * * * * */

VolumeFractionPenalizer::VolumeFractionPenalizer(      FILE        *logfile,
                                                 const double       totalVolume,
                                                 const double       volFracLimit,
                                                 const double       popFrac,
                                                       Initializer *initializer)
: Penalizer(logfile), totalVolume_(totalVolume), volFracLimit_(volFracLimit),
  popFrac_(popFrac), initializer_(initializer)
{}

void VolumeFractionPenalizer::penalize(                 GenePool *pool,
                                       gmx_unused const int       generation)
{
    // Get the volume of the population
    const double poolVolume = getPoolVolume(*pool);
    const double volFrac_ = poolVolume/totalVolume_;
    if (volFrac_ < volFracLimit_)
    {
        if (logfile_)
        {
            fprintf(
                logfile_,
                "VolumeFractionPenalizer: population volume is %lf, total volume is %lf, fraction %lf is below the limit %lf, thus randomizing the worst %lf%% of genomes...\n",
                poolVolume,
                totalVolume_,
                volFrac_,
                volFracLimit_,
                popFrac_ * 100
            );
        }
        // Randomize
        size_t i = static_cast<size_t>(
            (1-popFrac_) * static_cast<double>(pool->popSize())
        );
        for (; i < pool->popSize(); i++)
        {
            initializer_->randomizeGenome(pool->genomePtr(i));
        }
    }
}

/* * * * * * * * * * * * * * * * * * * * * *
* END: VolumeFractionPenalizer             *
* * * * * * * * * * * * * * * * * * * * * */

} // namespace ga

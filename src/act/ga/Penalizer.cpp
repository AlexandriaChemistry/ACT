/*! \internal \brief
 * Implements part of the alexandria program.
 * \author Julian Ramon Marrades Furquet <julian.marrades@hotmail.es>
 */

#include <limits>

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

double VolumeFractionPenalizer::getPoolVolume(const GenePool &pool) const
{
    double volume = 1;
    for (size_t i = 0; i < pool.genomeSize(); i++)  // For each parameter
    {
        double maximum = std::numeric_limits<double>::min();
        double minimum = std::numeric_limits<double>::max();
        for (auto genome : pool.genePool())  // For each genome
        {
            if (genome.base(i) > maximum)
            {
                maximum = genome.base(i);
            }
            else if (genome.base(i) < minimum)
            {
                minimum = genome.base(i);
            }
        }
        volume *= maximum - minimum;
    }
    GMX_RELEASE_ASSERT(
        volume >= 0,
        "VolumeFractionPenalizer: the volume of the population is negative. Something went really wrong"
    );
    return volume;
}

/* * * * * * * * * * * * * * * * * * * * * *
* END: VolumeFractionPenalizer             *
* * * * * * * * * * * * * * * * * * * * * */

} // namespace ga

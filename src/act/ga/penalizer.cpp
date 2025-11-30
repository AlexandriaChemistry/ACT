/*
 * This source file is part of the Alexandria Chemistry Toolkit.
 *
 * Copyright (C) 2021-2025
 *
 * Developers:
 *             Mohammad Mehdi Ghahremanpour,
 *             Julian Marrades,
 *             Marie-Madeleine Walz,
 *             Paul J. van Maaren,
 *             David van der Spoel (Project leader)
 *
 * This program is free software; you can redistribute it and/or
 * modify it under the terms of the GNU General Public License
 * as published by the Free Software Foundation; either version 2
 * of the License, or (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this program; if not, write to the Free Software
 * Foundation, Inc., 51 Franklin Street, Fifth Floor,
 * Boston, MA  02110-1301, USA.
 */
/*! \internal \brief
 * Implements part of the alexandria program.
 * \author Julian Ramon Marrades Furquet <julian@marrad.es>
 */

#include <limits>
#include <cmath>

#include "penalizer.h"

#include "gromacs/utility/textwriter.h"

namespace ga
{

/* * * * * * * * * * * * * * * * * * * * * *
* BEGIN: VolumeFractionPenalizer           *
* * * * * * * * * * * * * * * * * * * * * */

VolumeFractionPenalizer::VolumeFractionPenalizer(      gmx_output_env_t  *oenv,
                                                 const bool               logVolume,
                                                 const double             totalVolume,
                                                 const double             volFracLimit,
                                                 const double             popFrac,
                                                       Initializer       *initializer)
: logVolume_(logVolume), totalVolume_(totalVolume),
  volFracLimit_(volFracLimit), popFrac_(popFrac), initializer_(initializer)
{
    outfile_ = xvgropen(
       logVolume ? "vfp_log_volume.xvg" : "vfp_volume.xvg",
       "",
       "generation",
       logVolume ? "log_volume" : "volume",
       oenv
    );
    std::vector<const char*> tmpParamNames = {
        logVolume ? "pop_log_volume" : "pop_volume",
        "volume_fraction"
    };
    xvgr_legend(outfile_, tmpParamNames.size(), tmpParamNames.data(), oenv);
}

bool VolumeFractionPenalizer::penalize(gmx::TextWriter *tw,
                                       GenePool        *pool,
                                       const int        generation)
{
    // Get the volume of the population
    const double poolVolume = getPoolVolume(*pool);
    const double volFrac_ = poolVolume/totalVolume_;
    // Print to output file
    fprintf(outfile_, "%d %lf %lf\n", generation, poolVolume, volFrac_);
    if (volFrac_ < volFracLimit_)
    {
        if (tw)
        {
            tw->writeStringFormatted("Generation %d -> VolumeFractionPenalizer %s: population volume is %g, total volume is %g, fraction %g is below the limit %g, thus randomizing the worst %g%% of genomes...\n",
                                     generation,
                                     logVolume_ ? "(log scale)" : "",
                                     poolVolume,
                                     totalVolume_,
                                     volFrac_,
                                     volFracLimit_,
                                     popFrac_ * 100
                                     );
        }
        // Randomize
        size_t i = static_cast<size_t>(
            lround(
                (1-popFrac_) * static_cast<double>(pool->popSize())
            )
        );
        for (; i < pool->popSize(); i++)
        {
            initializer_->randomizeGenome(pool->genomePtr(i));
        }
        return true;
    }
    else
    {
        return false;
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
    return logVolume_ ? log(volume) : volume;
}

VolumeFractionPenalizer::~VolumeFractionPenalizer()
{
    xvgrclose(outfile_);
}

/* * * * * * * * * * * * * * * * * * * * * *
* END: VolumeFractionPenalizer             *
* * * * * * * * * * * * * * * * * * * * * */

/* * * * * * * * * * * * * * * * * * * * * *
* BEGIN: CatastrophePenalizer              *
* * * * * * * * * * * * * * * * * * * * * */

CatastrophePenalizer::CatastrophePenalizer(const int          seed,
                                           const int          genInterval,
                                           const double       popFrac,
                                                 Initializer *initializer,
                                           const size_t       popSize)
: gen(rd()), genInterval_(genInterval), popFrac_(popFrac),
  initializer_(initializer), availableIndices_(popSize)
{
    gen.seed(seed);
    for (size_t i = 0; i < popSize; i++)
    {
        availableIndices_[i] = i;
    }
}

/*! \brief Print cause of rapid decline in population
 * \param[in] tw          Text writer
 * \param[in] generation  Generation number
 * \param[in] percentage  Percentage of population that will die
 * \param[in] gen         Part of the random number generation
 */
static void printCause(gmx::TextWriter *tw,
                       int              generation,
                       double           percentage,
                       std::mt19937    &gen)
{
    std::vector<const char *> ccc = {
        "Earth was hit by a strong burst of gamma radiation wiping out %g percent\nof the population. New mutants will take over.\n",
        "Forest fires wiped out %g percent of the population.\nSeedlings will replace them.\n",
        "Anthrax bacteria were released from thawing permafrost, killing %g percent of\nthe population. Breeding of new individuals commenced.\n",
        "Collapse of the Thwaites glacier in Antarctica leads to flooding globally\nkilling %g percent of the population. New life forms will repopulate.\n",
        "Melting of the Greenland ice cap released so much fresh water that the\nNorth Atlantic Current stopped and Europe became uninhabitable.\n%g percent of the population will be replaced.",
        "Deforestation leading to catastrophic loss of biodiversity reduced fitness.\n%g percent new species will be created.\n",
        "Alien invaders exterminate %g percent of the population and replaced them\nwith teenage Mutant Ninja turtles.\n",
        "An asteroid hits earth and kills %g%% of all living organism.\nEvolution starts anew.\n",
        "A tsunami hits the western pacific rim wiping out %g%% of the population.\nRebuilding will commence.\n",
        "Strict reproduction laws in patriarchic countries lead to an inbreeding depression.\nTo address the reduced fitness of the population %g%% will be replaced.\n",
        "Tornadoes accross North America wipe out most major cities.\n%g percent of the population is replaced.\n",
        "Hurricanes wipe out the population in south-east Asia.\n%g percent of the population will be rejuvenated.\n",
        "Collapsed hydrological dams devastate land in South America.\n%g percent of the population is reborn.\n",
        "Drought, exacerbated by climate change, wipes out the population of\nsub-Saharan Africa.\n%g%% of the world population is replaced.\n",
        "A solar flare wipes out all communication on earth, killing %g percent.\nNewborns on the way.\n",
        "Major rivers on the Indian subcontent dry up due to dwindling glaciers in the Himalayas.\n%g percent of the world population is replaced.\n",
        "The empire fires off the death star, scraping the surface of the earth.\n%g percent of individuals are exterminated.\n",
        "Desertification makes agriculture in  southern Europe impossible.\n%g percent of the population will be regrown.\n",
        "Mutant fungi wipe out industrial agriculture in North America.\n%g percent of the population has to start from scratch.\n",
        "Tricholoma fungi adapt to oil palm trees in south east Asia, turning them into matsutake mushrooms.\n%g percent of the population are impoverished and replaced.\n"
    };
    std::uniform_int_distribution<int>  distr(0, ccc.size());
    int csel = distr(gen) % ccc.size();
    tw->writeStringFormatted("\nGeneration %d: ", generation);
    tw->writeStringFormatted(ccc[csel], percentage);
}

bool CatastrophePenalizer::penalize(gmx::TextWriter *tw,
                                    GenePool        *pool,
                                    const int        generation)
{
    if (generation > 0 && generation % genInterval_ == 0)
    {
        if (tw)
        {
            printCause(tw, generation, popFrac_ * 100, gen);
        }
        // Shuffle the indices vector and randomize the last (100 - popFrac_*100)%
        std::shuffle(availableIndices_.begin(), availableIndices_.end(), gen);
        size_t i = static_cast<size_t>(
            lround(
                (1-popFrac_) * static_cast<double>(availableIndices_.size())
            )
        );
        for (; i < availableIndices_.size(); i++)
        {
            initializer_->randomizeGenome(pool->genomePtr(availableIndices_[i]));
        }
        return true;
    }
    else
    {
        return false;
    }
}

/* * * * * * * * * * * * * * * * * * * * * *
* END: CatastrophePenalizer                *
* * * * * * * * * * * * * * * * * * * * * */

} // namespace ga

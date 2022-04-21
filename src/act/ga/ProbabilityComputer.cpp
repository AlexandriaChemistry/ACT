/*! \internal \brief
 * Implements part of the alexandria program.
 * \author Julian Ramon Marrades Furquet <julian.marrades@hotmail.es>
 */
#include "ProbabilityComputer.h"

#include <math.h>

namespace ga
{

/* * * * * * * * * * * * * * * * * * * * * *
* BEGIN: FitnessProbabilityComputer        *
* * * * * * * * * * * * * * * * * * * * * */

void FitnessProbabilityComputer::compute(                 std::vector<Genome> *pop,
                                         gmx_unused const int                  generation)
{
    double       total = 0;
    const double epsilon = 1e-4;
    for (size_t i = 0; i < pop->size(); i++)
    {
        inverses_[i] = 1 / ( epsilon + (*pop)[i].fitness(iMolSelect::Train) );
        total += inverses_[i];
    }
    for (size_t i = 0; i < pop->size(); i++)
    {
        (*pop)[i].setProbability( inverses_[i] / total );
    }
}

/* * * * * * * * * * * * * * * * * * * * * *
* END: FitnessProbabilityComputer          *
* * * * * * * * * * * * * * * * * * * * * */

/* * * * * * * * * * * * * * * * * * * * * *
* BEGIN: BoltzmannProbabilityComputer      *
* * * * * * * * * * * * * * * * * * * * * */

double BoltzmannProbabilityComputer::computeTemperature(const int generation)
{
    double temp = temperature_;
    if (generation >= maxGenerations_)
    {
        return 1e-6;
    }
    else if (generation >= boltzAnneal_ * maxGenerations_)
    {
        temp = ( temperature_ / ( boltzAnneal_ * maxGenerations_ - maxGenerations_ ) ) * generation + ( ( temperature_ / ( maxGenerations_ - boltzAnneal_ * maxGenerations_ ) ) * maxGenerations_ );
    }
    return temp;
}

void BoltzmannProbabilityComputer::compute(      std::vector<Genome> *pop,
                                           const int                  generation)
{
    double       total   = 0;
    const double epsilon = 1e-4;
    const double temp    = boltzAnneal_ < 1 ? computeTemperature(generation) : temperature_;
    for (size_t i = 0; i < pop->size(); i++)
    {
        exponentials_[i] = exp(( 1 / ( epsilon + (*pop)[i].fitness(iMolSelect::Train) ) ) / temp);
        total += exponentials_[i];
    }
    for (size_t i = 0; i < pop->size(); i++)
    {
        (*pop)[i].setProbability(exponentials_[i] / total);
    }
}

/* * * * * * * * * * * * * * * * * * * * * *
* END: BoltzmannProbabilityComputer        *
* * * * * * * * * * * * * * * * * * * * * */

/* * * * * * * * * * * * * * * * * * * * * *
* BEGIN: RankProbabilityComputer           *
* * * * * * * * * * * * * * * * * * * * * */

void RankProbabilityComputer::compute(                 std::vector<Genome> *pop,
                                      gmx_unused const int                  generation)
{
    for (size_t i = 0; i < pop->size(); i++)
    {
        (*pop)[i].setProbability((pop->size() - i) / sumOfRanks_);
    }
}

/* * * * * * * * * * * * * * * * * * * * * *
* END: RankProbabilityComputer             *
* * * * * * * * * * * * * * * * * * * * * */


} //namespace ga


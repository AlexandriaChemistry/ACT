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

void FitnessProbabilityComputer::compute(std::vector<Genome> *pop)
{
    double       total = 0;
    const double epsilon = 1e-4;
    for (size_t i = 0; i < pop->size(); i++)
    {
        total += 1 / ( epsilon + (*pop)[i].fitness(iMolSelect::Train) );
    }
    for (size_t i = 0; i < pop->size(); i++)
    {
        (*pop)[i].setProbability( ( 1 / ( epsilon + (*pop)[i].fitness(iMolSelect::Train) ) ) / total );
    }
}

/* * * * * * * * * * * * * * * * * * * * * *
* END: FitnessProbabilityComputer          *
* * * * * * * * * * * * * * * * * * * * * */

/* * * * * * * * * * * * * * * * * * * * * *
* BEGIN: BoltzmannProbabilityComputer      *
* * * * * * * * * * * * * * * * * * * * * */

void BoltzmannProbabilityComputer::compute(std::vector<Genome> *pop)
{
    double       total = 0;
    const double epsilon = 1e-4;
    for (size_t i = 0; i < pop->size(); i++)
    {
        exponentials_[i] = exp(( 1 / ( epsilon + (*pop)[i].fitness(iMolSelect::Train) ) ) / temperature_);
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

void RankProbabilityComputer::compute(std::vector<Genome> *pop)
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


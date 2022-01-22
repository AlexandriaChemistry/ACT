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
    double total = 0;
    for (size_t i = 0; i < pop->size(); i++)
    {
        total += (*pop)[i].fitness(iMolSelect::Train);
    }
    for (size_t i = 0; i < pop->size(); i++)
    {
        (*pop)[i].setProbability((*pop)[i].fitness(iMolSelect::Train) / total);
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
    //! Stores \f$ e^{f_i} \f$ for each Individual \f$ i \f$
    std::vector<double> exponentials;
    double              total = 0;
    for (size_t i = 0; i < pop->size(); i++)
    {
        exponentials.push_back(exp((*pop)[i].fitness(iMolSelect::Train) / temperature_));
        total += exponentials[i];
    }
    for (size_t i = 0; i < pop->size(); i++)
    {
        (*pop)[i].setProbability(exponentials[i] / total);
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


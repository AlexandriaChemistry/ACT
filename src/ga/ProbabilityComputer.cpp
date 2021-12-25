#include "ProbabilityComputer.h"

#include <math.h>


namespace ga
{


/* * * * * * * * * * * * * * * * * * * * * *
* BEGIN: FitnessProbabilityComputer        *
* * * * * * * * * * * * * * * * * * * * * */

void FitnessProbabilityComputer::compute(std::vector<Individual*> *pop)
{
    double total = 0;
    int i;
    for (i = 0; i < pop->size(); i++)
    {
        total += (*pop)[i]->fitnessTrain();
    }
    for (i = 0; i < pop->size(); i++)
    {
        (*pop)[i]->setProbability((*pop)[i]->fitnessTrain() / total);
    }
}

/* * * * * * * * * * * * * * * * * * * * * *
* END: FitnessProbabilityComputer          *
* * * * * * * * * * * * * * * * * * * * * */

/* * * * * * * * * * * * * * * * * * * * * *
* BEGIN: BoltzmannProbabilityComputer      *
* * * * * * * * * * * * * * * * * * * * * */

void BoltzmannProbabilityComputer::compute(std::vector<Individual*> *pop)
{
    double total = 0;
    int i;
    for (i = 0; i < pop->size(); i++)
    {
        exponentials_[i] = exp((*pop)[i]->fitnessTrain() / temperature_);
        total += exponentials_[i];
    }
    for (i = 0; i < pop->size(); i++)
    {
        (*pop)[i]->setProbability(exponentials_[i] / total);
    }
}

/* * * * * * * * * * * * * * * * * * * * * *
* END: BoltzmannProbabilityComputer        *
* * * * * * * * * * * * * * * * * * * * * */

/* * * * * * * * * * * * * * * * * * * * * *
* BEGIN: RankProbabilityComputer           *
* * * * * * * * * * * * * * * * * * * * * */

void RankProbabilityComputer::compute(std::vector<Individual*> *pop)
{
    for (int i = 0; i < pop->size(); i++)
    {
        (*pop)[i]->setProbability((pop->size() - i) / sumOfRanks_);
    }
}

/* * * * * * * * * * * * * * * * * * * * * *
* END: RankProbabilityComputer             *
* * * * * * * * * * * * * * * * * * * * * */


} //namespace ga


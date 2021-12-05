#include "ProbabilityComputer.h"

#include "aliases.h"

#include <math.h>


namespace ga
{


/* * * * * * * * * * * * * * * * * * * * * *
* BEGIN: FitnessProbabilityComputer        *
* * * * * * * * * * * * * * * * * * * * * */

void FitnessProbabilityComputer::compute(const vector  &fitness,
                                               vector  *prob,
                                         const int      popSize)
{
    double total = 0;
    int i;
    for (i = 0; i < popSize; i++)
    {
        total += fitness[i];
    }
    for (i = 0; i < popSize; i++)
    {
        (*prob)[i] = fitness[i] / total;
    }
}

/* * * * * * * * * * * * * * * * * * * * * *
* END: FitnessProbabilityComputer          *
* * * * * * * * * * * * * * * * * * * * * */

/* * * * * * * * * * * * * * * * * * * * * *
* BEGIN: BoltzmannProbabilityComputer      *
* * * * * * * * * * * * * * * * * * * * * */

void BoltzmannProbabilityComputer::compute(const vector    &fitness,
                                                 vector    *prob,
                                           const int        popSize)
{
    double total = 0;
    int i;
    for (i = 0; i < popSize; i++)
    {
        exponentials[i] = exp(fitness[i] / temperature);
        total += exponentials[i];
    }
    for (i = 0; i < popSize; i++)
    {
        (*prob)[i] = exponentials[i] / total;
    }
}

/* * * * * * * * * * * * * * * * * * * * * *
* END: BoltzmannProbabilityComputer        *
* * * * * * * * * * * * * * * * * * * * * */

/* * * * * * * * * * * * * * * * * * * * * *
* BEGIN: RankProbabilityComputer           *
* * * * * * * * * * * * * * * * * * * * * */

void RankProbabilityComputer::compute(const vector &fitness,
                                            vector *prob,
                                      const int     popSize)
{
    for (int i = 0; i < popSize; i++)
    {
        (*prob)[i] = (popSize - i) / sumOfRanks;
    }
}

/* * * * * * * * * * * * * * * * * * * * * *
* END: RankProbabilityComputer             *
* * * * * * * * * * * * * * * * * * * * * */


} //namespace ga


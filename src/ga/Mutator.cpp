#include "Mutator.h"


namespace ga
{


/* * * * * * * * * * * * * * * * * * * * * *
* BEGIN: PercentMutator                    *
* * * * * * * * * * * * * * * * * * * * * */

void PercentMutator::mutate(      vector   *individual,
                            const int       chromosomeLength,
                            const double    prMut)
{
    for (int i = 0; i < chromosomeLength; i++)
    {
        if (randNum() <= prMut) (*individual)[i] *= dis(gen);
    }
}

/* * * * * * * * * * * * * * * * * * * * * *
* END: PercentMutator                      *
* * * * * * * * * * * * * * * * * * * * * */

/* * * * * * * * * * * * * * * * * * * * * *
* BEGIN: RangeMutator                      *
* * * * * * * * * * * * * * * * * * * * * */

void RangeMutator::mutate(      vector   *individual,
                          const int       chromosomeLength,
                          const double    prMut)
{
    for (int i = 0; i < chromosomeLength; i++)
    {
        if (randNum() <= prMut) (*individual)[i] += dis(gen);
    }
}

/* * * * * * * * * * * * * * * * * * * * * *
* END: RangeMutator                        *
* * * * * * * * * * * * * * * * * * * * * */


} //namespace ga

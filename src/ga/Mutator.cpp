#include "Mutator.h"


namespace ga
{

    void PercentMutator::mutate(      vector   *individual,
                                const int       chromosomeLength,
                                const double    prMut)
    {
        for (int i = 0; i < chromosomeLength; i++)
        {
            if (randNum()) (*individual)[i] *= dis(gen);
        }
    }


    void RangeMutator::mutate(      vector   *individual,
                              const int       chromosomeLength,
                              const double    prMut)
    {
        for (int i = 0; i < chromosomeLength; i++)
        {
            if (randNum()) (*individual)[i] += dis(gen);
        }
    }

}

#include "ProbabilityComputer.h"

#include "aliases.h"

#include <math.h>


namespace ga
{

    void FitnessProbabilityComputer::compute(const vector&  fitness,
                                                   vector&  prob,
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
            prob[i] = fitness[i] / total;
        }
    }


    BoltzmannProbabilityComputer::BoltzmannProbabilityComputer(const int    popSize,
                                                               const double temperature)
    {
        exponentials = vector(popSize);
        this->temperature = temperature;
    }


    void BoltzmannProbabilityComputer::compute(const vector&    fitness,
                                                     vector&    prob,
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
            prob[i] = exponentials[i] / total;
        }
    }


    RankProbabilityComputer::RankProbabilityComputer(const int popSize)
    {
        sumOfRanks = popSize * (popSize + 1) / 2;
    }


    void RankProbabilityComputer::compute(const vector& fitness,
                                                vector& prob,
                                          const int     popSize)
    {
        for (int i = 0; i < popSize; i++)
        {
            prob[i] = (popSize - i) / sumOfRanks;
        }
    }

}


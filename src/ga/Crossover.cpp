#include "Crossover.h"

#include "helpers.h"
#include <random>
#include <cmath>
#include <algorithm>

namespace ga
{

    int Crossover::randIndex()
    {
        return dis(gen);
    }

    void SinglePointCrossover::offspring(const vector &parent1,
                                         const vector &parent2,
                                               vector *child1,
                                               vector *child2,
                                         const int     length)
    {

        const int index = randIndex();

        int i;
        for (i = 0; i < index; i++)
        {
            (*child1)[i] = parent1[i];
            (*child2)[i] = parent2[i];
        }
        for (i = index; i < length; i++)
        {
            (*child1)[i] = parent2[i];
            (*child2)[i] = parent1[i];
        }

    }


    void DoublePointCrossover::offspring(const vector &parent1,
                                         const vector &parent2,
                                               vector *child1,
                                               vector *child2,
                                         const int     length)
    {

        const int tmp1 = randIndex();
        int tmp2 = randIndex();
        while (tmp1 == tmp2)
        {
            tmp2 = randIndex();
        }
        const int diff = std::abs(tmp1 - tmp2);

        const int index1 = std::min(tmp1, tmp2);

        int i;
        for (i = 0; i < index1; i++)
        {
            (*child1)[i] = parent1[i];
            (*child2)[i] = parent2[i];
        }
        for (i = index1; i < index1 + diff; i++)
        {
            (*child1)[i] = parent2[i];
            (*child2)[i] = parent1[i];
        }
        for (i = index1 + diff; i < length; i++)
        {
            (*child1)[i] = parent1[i];
            (*child2)[i] = parent2[i];
        }

    }

    void NPointCrossover::offspring(const vector &parent1,
                                    const vector &parent2,
                                          vector *child1,
                                          vector *child2,
                                    const int     length) {

        assert(numberOfCrossovers < length);

        int index;
        int crossover;
        int remainingCrossovers = numberOfCrossovers;
        vector indexCandidates(length - 2);
        vector crossoverIndices(numberOfCrossovers + 2);

        if (numberOfCrossovers == length - 1) {
            for (index = 0; index < length; index++)
                crossoverIndices[index] = index;
            crossover = 0;
        } else {
            for (index = 1; index < length - 1; index++)
                indexCandidates[index - 1] = index;

            crossoverIndices[0] = 0;
            crossoverIndices[numberOfCrossovers + 1] = length - 1;

            crossover = 1;

            while (crossover < numberOfCrossovers + 1) {
                index = randIndex() - 1;
                while (indexCandidates[index] == -1 || index > length - 3)
                    index = randIndex() - 1;
                crossoverIndices[crossover] = indexCandidates[index];
                indexCandidates[index] = -1;
                crossover++;
            }

            crossover = 0;

            sort(crossoverIndices.begin(), crossoverIndices.end());
        }

        while (remainingCrossovers > 0) {
              for (index = crossoverIndices[crossover]; index < crossoverIndices[crossover + 1]; index++) {
                  (*child1)[index] = parent1[index];
                  (*child2)[index] = parent2[index];
              } for (index = crossoverIndices[crossover + 1]; index < crossoverIndices[crossover + 2]; index++) {
                  (*child1)[index] = parent2[index];
                  (*child2)[index] = parent1[index];
              }
              crossover += 2;
              remainingCrossovers -= 2;
        }

        if (numberOfCrossovers % 2 == 0) {
            for (index = crossoverIndices[crossover]; index < length; index++) {
                (*child1)[index] = parent1[index];
                (*child2)[index] = parent2[index];
            }
        } else {
            (*child1)[length - 1] = parent2[length - 1];
            (*child2)[length - 1] = parent1[length - 1];
        }
    }
}

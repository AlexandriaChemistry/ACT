#include "FitnessComputer.h"

#include "aliases.h"

void SimpleFitnessComputer::compute(const vector&     individual,
                                          vector&     fitness,
                                    const int         indIndex,
     			              		const int         length) {
    double sum = 0;

    for (int i = 0; i < length; i++)
        sum += individual[i] * individual[i];

    fitness[indIndex] = 1 / sum;
}
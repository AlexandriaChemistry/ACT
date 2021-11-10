#include "FitnessComputer.h"

#include "aliases.h"

void SimpleFitnessComputer::compute(const vector       individual,
			                		double* const      ftPtr,
     			              		const int          length) {
    double sum = 0;

    for (int i = 0; i < length; i++)
        sum += individual[i]*individual[i];

    return 1 / sum;
}
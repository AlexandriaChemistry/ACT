#include "ProbabilityComputer.h"


void FitnessProbabilityComputer::compute(const vector fitness, const vector prob, const int popSize) {
    double total = 0;
    int i;
    for (i = 0; i < popSize; i++) {
        total += fitness[i];
    }
    for (i = 0; i < popSize; i++) {
        prob[i] = fitness[i] / total;
    }
}

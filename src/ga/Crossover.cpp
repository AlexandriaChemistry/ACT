#include "Crossover.h"

#include <random>
#include <cmath>
#include <algorithm>


Crossover::Crossover(const int chromosomeLength) {
    gen = std::mt19937(rd());
    // Exclude first and last indices. Otherwise, there is no crossover.
    dis = std::uniform_int_distribution<>(1, chromosomeLength - 2);
}


int Crossover::randIndex() {
    return dis(gen);
}


void SinglePointCrossover::offspring(const vector   parent1,
                                     const vector   parent2,
                                           vector   child1,
                                           vector   child2,
                                     const int      length) {

    const int index = randIndex();

    int i;
    for (i = 0; i < index; i++) {
        child1[i] = parent1[i];
        child2[i] = parent2[i];
    }
    for (i = index; i < length; i++) {
        child1[i] = parent2[i];
        child2[i] = parent1[i];
    }

}


void DoublePointCrossover::offspring(const vector   parent1,
                                     const vector   parent2,
                                           vector   child1,
                                           vector   child2,
                                     const int      length) {

    const int tmp1 = randIndex();
    int tmp2 = randIndex();
    while (tmp1 == tmp2) {
        tmp2 = randIndex();
    }
    const int diff = std::abs(tmp1 - tmp2);

    const int index1 = std::min(tmp1, tmp2);

    int i;
    for (i = 0; i < index1; i++) {
        child1[i] = parent1[i];
        child2[i] = parent2[i];
    }
    for (i = index1; i < index1 + diff; i++) {
        child1[i] = parent2[i];
        child2[i] = parent1[i];
    }
    for (i = index1 + diff; i < length; i++) {
        child1[i] = parent1[i];
        child2[i] = parent2[i];
    }

}
#include <tuple>

#include "crossover.h"

void singlePoint(double parent1[], double parent2[], double child1[], double child2[], const int length) {
    const int index = rand() % length;
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
#include <tuple>

#include "crossover.h"

std::tuple<double*, double*> singlePoint(double parent1[], double parent2[], const int length) {
    const int index = rand() % length;
    double child1[length] = {};
    double child2[length] = {};
    int i;
    for (i = 0; i < index; i++) {
        child1[i] = parent1[i];
        child2[i] = parent2[i];
    }
    for (i = index; i < length; i++) {
        child1[i] = parent2[i];
        child2[i] = parent1[i];
    }
    std::make_tuple(child1, child2);
}
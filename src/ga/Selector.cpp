#include "Selector.h"

#include "aliases.h"

#include <random>


RouletteSelector::RouletteSelector() {

    gen = std::mt19937(rd);
    dis = std::uniform_real_distribution<>(0.0, 1.0);

}


const vector RouletteSelector::select(const matrix population, const vector probability, const int popSize) {

    double num = dis(gen);
    int i = 0;
    while (num > 0) {
        num += probability[i];
        i++;
    }
    return population[i-1];

}

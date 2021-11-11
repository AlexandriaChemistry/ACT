#include "Selector.h"

#include "aliases.h"

#include <random>


const int RouletteSelector::select(const vector&  probability,
                                   const int      popSize) {

    double num = dis(gen);
    int i = 0;
    while (num > 0 and i < popSize) {
        num += probability[i];
        i++;
    }
    return i-1;

}

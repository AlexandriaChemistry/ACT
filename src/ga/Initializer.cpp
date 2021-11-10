#include "Initializer.h"

SimpleInitializer::SimpleInitializer(const double min, const double max) {
    gen = std::mt19937(rd());
    dis = std::uniform_real_distribution<>(min, max);
}

void SimpleInitializer::initialize(      vector     individual,
                                   const int        length) {
    for (int i = 0; i < length; i++) {
        individual[i] = dis(gen);
    }
}

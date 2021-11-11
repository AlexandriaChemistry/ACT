#include "Initializer.h"

#include <stdio.h>


void SimpleInitializer::initialize(      vector&    individual,
                                   const int        length) {
    for (int i = 0; i < length; i++) {
        individual[i] = dis(gen);
    }
}

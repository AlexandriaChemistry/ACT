#include "Mutator.h"


void PercentMutator::mutate(      vector&   individual,
                            const int       indGen) {
    individual[indGen] *= dis(gen);
}

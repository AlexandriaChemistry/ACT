#include "Mutator.h"


void PercentMutator::mutate(double *const gene) {
    (*gene) *= dis(gen);
}

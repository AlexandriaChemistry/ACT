#include "Mutator.h"


PercentMutator::PercentMutator(const double frac) {
    gen = std::mt19937(rd());
    dis = std::uniform_real_distribution<>(1-frac, 1+frac);
}


void PercentMutator::mutate(double *const gene) {
    (*gene) = dis(gen) * (*gene);
}

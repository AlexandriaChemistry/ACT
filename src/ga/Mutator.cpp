#include "Mutator.h"


namespace ga {

    void PercentMutator::mutate(      vector&   individual,
                                const int       indGene) {
        individual[indGene] *= dis(gen);
    }


    void RangeMutator::mutate(      vector& individual,
                              const int     indGene) {
        individual[indGene] += dis(gen);
    }

}

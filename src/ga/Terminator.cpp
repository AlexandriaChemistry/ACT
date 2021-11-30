#include "Terminator.h"

#include "ga_helpers.h"
#include "aliases.h"


namespace ga
{

    SimpleTerminator::SimpleTerminator(const double tolerance)
    {
        this->tolerance = tolerance;
    }


    GenerationTerminator::GenerationTerminator(const int maxGenerations)
    {
        this->maxGenerations = maxGenerations;
    }


    bool SimpleTerminator::terminate(const matrix  &population,
                                     const vector  &fitness,
                                     const int      generationNumber,
                                     const int      popSize,
                                     const int      chromosomeLength)
    {
        double maximumFitness = fitness[findMaximumIndex(fitness, popSize)];
        return maximumFitness >= 1 / (tolerance * chromosomeLength);
    }


    bool GenerationTerminator::terminate(const matrix  &population,
                                         const vector  &fitness,
                                         const int      generationNumber,
                                         const int      popSize,
                                         const int      chromosomeLength)
    {
        return generationNumber >= maxGenerations;
    }

}

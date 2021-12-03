#include "aliases.h"
#include "ga_includes.h"

#include "ga_helpers.h"

#include <stdio.h>
#include <stdlib.h>

using namespace ga;

int main(int argc, char const *argv[])
{

    if (argc != 8)
    {
        printf("\nUsage: ./test <nElites> <popSize> <chromLen> <tolerance> <verbose> <nrep> <ncrossovers>\n\n");
        return -1;
    }

	const int       popSize         = atoi(argv[2]); // This
	const int       chromLen        = atoi(argv[3]);
//	const int       maxGenerations  = 5;
//    const double    mutFrac         = 0.1;
    const int       nElites         = atoi(argv[1]); // This
    const double    tolerance       = atof(argv[4]);
    const int       verbose         = atoi(argv[5]);
    const int       nrep            = atoi(argv[6]);
    const int       ncrossovers     = atoi(argv[7]); // This

    // Print given configuration
    printf("\nRun configuration:\n");
    printf("  - popSize: %i\n", popSize);
    printf("  - chromLen: %i\n", chromLen);
    printf("  - nElites: %i\n", nElites);
    printf("  - tolerance: %lf\n", tolerance);
    printf("  - verbose: %i\n", verbose);
    printf("  - nrep: %i\n", nrep);
    printf("  - ncrossovers: %i\n", ncrossovers);


	SimpleInitializer           init(-10.0, 10.0);
	SimpleFitnessComputer       fit;
//    EmptySorter                 sort;
//    MergeSorter                 sort(popSize, chromLen, true);
	QuickSorter                 sort(popSize, true);
//    FitnessProbabilityComputer  procomp;
	RankProbabilityComputer     procomp(popSize);
	RouletteSelector            select;
	NPointCrossover             npoint(chromLen, ncrossovers);
	// SinglePointCrossover        singlepoint(chromLen);
//  PercentMutator              mutate(mutFrac);
    RangeMutator                mutate(0.5);
    SimpleTerminator            terminate(tolerance);
//	GenerationTerminator        terminate(maxGenerations);

	GeneticAlgorithm ga = GeneticAlgorithm(popSize,
										   chromLen,
                                           nElites,
										   &init,
										   &fit,
										   &sort,
										   &procomp,
										   &select,
										   &npoint,
										   &mutate,
										   &terminate);

    vector gens(nrep);

    for (int i = 0; i < nrep; i++) {
        const ga_result_t result = ga.evolve(0.35, 0.01, verbose);
        gens[i] = result.generations;
    }

    const double avg = vectorMEAN(gens, nrep);
    const double std = vectorSTD(gens, avg, nrep);

    printf("\nGenerations required: %f +- %f\n\n", avg, std);

//    printf("\nEvolution took %i generations\n", result.generations);
//    printf("Best individual: ");
//    printVector(result.bestIndividual);
//    printf("Best fitness: %f\n\n", result.bestFitness);

	return 0;
}

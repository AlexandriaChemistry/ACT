#include "aliases.h"
#include "includes.h"

#include "helpers.h"

#include <stdio.h>

using namespace ga;

int main(int argc, char const *argv[]) {

	const int       popSize         = 200;
	const int       chromLen        = 10;
//	const int       maxGenerations  = 5;
//    const double    mutFrac         = 0.1;
    const double    tolerance       = 0.000001;

	SimpleInitializer           init(-10.0, 10.0);
	SimpleFitnessComputer       fit;
//    EmptySorter                 sort;
//    MergeSorter                 sort(popSize, chromLen, false);
	QuickSorter                 sort(popSize, true);
//    FitnessProbabilityComputer  procomp;
	RankProbabilityComputer     procomp(popSize);
	RouletteSelector            select;
	SinglePointCrossover        singlepoint(chromLen);
//  PercentMutator              mutate(mutFrac);
    RangeMutator                mutate(0.5);
    SimpleTerminator            terminate(tolerance);
//	GenerationTerminator        terminate(maxGenerations);

	GeneticAlgorithm ga = GeneticAlgorithm(popSize,
										   chromLen,
										   &init,
										   &fit,
										   &sort,
										   &procomp,
										   &select,
										   &singlepoint,
										   &mutate,
										   &terminate);

	const ga_result_t result = ga.evolve(0.35, 0.01, 1);

    printf("\nEvolution took %i generations\n", result.generations);
    printf("Best individual: ");
    printVector(result.bestIndividual);
    printf("Best fitness: %f\n\n", result.bestFitness);

	return 0;
}

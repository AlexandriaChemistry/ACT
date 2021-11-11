#include "aliases.h"
#include "includes.h"

#include "helpers.h"

#include <stdio.h>

int main(int argc, char const *argv[]) {
	
	const int popSize = 10;
	const int chromLen = 5;

	SimpleInitializer init(-10.0, 10.0);
	SimpleFitnessComputer fit;
	MergeSorter sort(popSize, chromLen);
	RankProbabilityComputer procomp(popSize);
	RouletteSelector select;
	SinglePointCrossover singlepoint(chromLen);
    PercentMutator mutate(0.1);
	SimpleTerminator terminate;

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

	ga_result_t result = ga.evolve(0.35, 0.01, true);

    printf("\nEvolution took %i generations\n", result.generations);
    printf("Best individual: ");
    printVector(result.bestIndividual);
    printf("Best fitness: %f\n\n", result.bestFitness);

//    vector vec = vector(chromLen);
//    init.initialize(vec, chromLen);
//    printVector(vec);
//
//    vector fitness = vector(5);
//    fit.compute(vec, fitness, 1, chromLen);
//    printVector(fitness);

	return 0;
}

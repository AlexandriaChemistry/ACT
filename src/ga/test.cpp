#include "aliases.h"
#include "includes.h"

#include "helpers.h"

#include <stdio.h>

int main(int argc, char const *argv[]) {
	
	const int popSize = 100;
	const int chromLen = 10;
	const int maxGenerations = 3;
	// const double tolerance = 0.0001;

	SimpleInitializer init(-10.0, 10.0);
	SimpleFitnessComputer fit;
	EmptySorter sort;
	FitnessProbabilityComputer procomp;
	RouletteSelector select;
	SinglePointCrossover singlepoint(chromLen);
    PercentMutator mutate(0.2);
	// SimpleTerminator terminate(tolerance);
	GenerationTerminator terminate(maxGenerations);

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

	ga_result_t result = ga.evolve(0.35, 0.01, false);

    printf("\nEvolution took %i generations", result.generations);
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

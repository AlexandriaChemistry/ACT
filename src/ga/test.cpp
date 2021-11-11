#include "aliases.h"
#include "includes.h"

#include "helpers.h"

int main(int argc, char const *argv[]) {
	
	const int popSize = 10;
	const int chromLen = 5;

	SimpleInitializer init(-10.0, 10.0);
	SimpleFitnessComputer fit;
	EmptySorter sort;
	FitnessProbabilityComputer procomp;
	RouletteSelector select;
	SinglePointCrossover singlepoint(chromLen);
    PercentMutator mutate(0.2);
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

	ga.evolve(0.5, 0.01, true);

//    vector vec = vector(chromLen);
//    init.initialize(vec, chromLen);
//    printVector(vec);
//
//    vector fitness = vector(5);
//    fit.compute(vec, fitness, 1, chromLen);
//    printVector(fitness);

	return 0;
}

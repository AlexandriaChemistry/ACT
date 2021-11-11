#include "aliases.h"
#include "includes.h"

#include "helpers.h"

int main(int argc, char const *argv[]) {
	
//	const int popSize = 10;
	const int chromLen = 5;

	SimpleInitializer init(-10.0, 10.0);
//	FitnessComputer fit = SimpleFitnessComputer();
//	Sorter sort = EmptySorter();
//	ProbabilityComputer procomp = FitnessProbabilityComputer();
//	Selector select = RouletteSelector();
//	SinglePointCrossover singlepoint(chromLen);
//	Mutator mutate = PercentMutator(0.2);
//	Terminator terminate = SimpleTerminator();
//
//	GeneticAlgorithm ga = GeneticAlgorithm(popSize,
//										   chromLen,
//										   &init,
//										   &fit,
//										   &sort,
//										   &procomp,
//										   &select,
//										   &singlepoint,
//										   &mutate,
//										   &terminate);
//
//	ga.evolve(0.5, 0.01, true);

    vector vec = vector(chromLen);
    init.initialize(vec, chromLen);
    printVector(vec);

	return 0;
}

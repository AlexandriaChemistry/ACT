#include "aliases.h"
#include "includes.h"

int main(int argc, char const *argv[]) {
	
	const int popSize = 100;
	const int chromLen = 10;

	Initializer init = SimpleInitializer(-10, 10);
	FitnessComputer fit = SimpleFitnessComputer();
	Sorter sort = EmptySorter();
	ProbabilityComputer procomp = FitnessProbabilityComputer();
	Selector select = RouletteSelector();
	SinglePointCrossover singlepoint(chromLen);
	Mutator mutate = PercentMutator(0.2);
	Terminator terminate = SimpleTerminator();

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

	return 0;
}

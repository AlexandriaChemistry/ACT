#include "aliases.h"
#include "includes.h"

#include "helpers.h"

#include <stdio.h>

using namespace ga;

int main(int argc, char const *argv[]) {

	const int       popSize         = 12;
	const int       chromLen        = 4;
	const int       maxGenerations  = 5;
	const int				crossovers 			= 1;
  const double    mutFrac         = 0.1;
  const double    tolerance       = 0.000001;

		// vector parent1(chromLen);
		// vector parent2(chromLen);
		// vector child1(chromLen);
		// vector child2(chromLen);

		// for (int index = 0; index < chromLen; index++) {
		// 		parent1[index] = 1;
		// 		parent2[index] = 2;
		// }

	SimpleInitializer           init(-10.0, 10.0);
	SimpleFitnessComputer       fit;
  // EmptySorter                 sort;
	// MergeSorter                 sort(popSize, chromLen);
	QuickSorter                 sort(popSize);
  // FitnessProbabilityComputer  procomp;
	RankProbabilityComputer     procomp(popSize);
	RouletteSelector            select;
	NPointCrossover        			npoint(chromLen, crossovers);
  PercentMutator              mutate(mutFrac);
  // RangeMutator                mutate(0.5);
  // SimpleTerminator            terminate(tolerance);
	GenerationTerminator        terminate(maxGenerations);

	GeneticAlgorithm ga = GeneticAlgorithm(popSize,
										   chromLen,
										   &init,
										   &fit,
										   &sort,
										   &procomp,
										   &select,
										   &npoint,
										   &mutate,
										   &terminate);

	const ga_result_t result = ga.evolve(0.35, 0.01, 5);

  printf("\nEvolution took %i generations\n", result.generations);
  printf("Best individual: ");
  printVector(result.bestIndividual);
  printf("Best fitness: %f\n\n", result.bestFitness);

	return 0;
}

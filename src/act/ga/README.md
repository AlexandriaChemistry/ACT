# ga

This directory provides a generic, template-based genetic algorithm (GA) framework used for evolutionary optimisation of force-field parameters and other continuous search problems.

## Contents

### Population management
- **gene_pool.h / gene_pool.cpp** — `GenePool`: maintains the current population of individuals across generations; handles selection, replacement, and elitism.
- **genome.h / genome.cpp** — `Genome`: the floating-point parameter vector carried by each individual.
- **individual.h** — Abstract `Individual` base class; concrete subclasses add problem-specific fitness data.
- **genetic_algorithm.h / genetic_algorithm.cpp** — `GeneticAlgorithm`: main loop that iterates selection → crossover → mutation → evaluation until a termination criterion is met.

### Operators
- **crossover.h** — Abstract `Crossover` interface.
- **npointcrossover.h / npointcrossover.cpp** — N-point crossover operator that recombines two parent genomes at *n* randomly chosen cut points.
- **mutator.h** — Abstract `Mutator` interface; concrete implementations live in `src/act/alexandria/` (e.g. `PercentMutator`, `McmcMutator`).
- **selector.h / selector.cpp** — Selection strategies (tournament, roulette-wheel, etc.) that choose parents for the next generation.
- **initializer.h** — Abstract `Initializer` interface for creating the initial population.

### Fitness and probability
- **fitness_computer.h** — Abstract `FitnessComputer` interface; problem-specific implementations compute scalar fitness scores from a genome.
- **probability_computer.h / probability_computer.cpp** — Converts raw fitness scores to selection probabilities.
- **penalizer.h / penalizer.cpp** — Applies constraint penalties to out-of-bounds parameter values.
- **sorter.h** — Sorts individuals by fitness to support elitist strategies.
- **terminator.h / terminator.cpp** — Termination criteria (maximum generations, fitness threshold, stagnation detection).

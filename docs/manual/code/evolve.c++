void evolve(int popSize, int nElites, int prCross, int prMut)
{
    Genome oldPop[popSize] = initialize(popSize); // Initialize
    Genome newPop[];
    computeDeviation(oldPop); // Deviation from data
    int generation = 0;
    do
    {
        sort(oldPop); // Sorting
        if (penalize(oldPop, generation)) // Penalty?
        {
            // Recompute deviation and sort again
            computeDeviation(oldPop);
            sort(oldPop);
        }
        generation++;
        computeProbabilities(oldPop); // Selection probabilities
        // Elitism
        for (int i = 0; i < nElites; i++)
            newPop.add(oldPop[i]);
        // Rest of population
        for (int i = nElites; i < popSize; i += 2)
        {
            Genome parent1, parent2 = select(oldPop); // Selection
            Genome child1, child2;
            if (random() <= prCross) // Crossover?
                child1, child2 = crossover(parent1, parent2);
            else
                child1, child2 = parent1, parent2;
            // Mutation
            for (Genome child : {child1, child2})
            {
                mutate(child, prMut);
                newPop.add(child); // Add to new population
            }
        }
        computeDeviation(newPop); // Deviation from data
        oldPop = newPop; // Swap populations
        newPop.clear();  // Erase all genomes in the population
    }
    // Termination
    while (!terminate(oldPop, generation));
}
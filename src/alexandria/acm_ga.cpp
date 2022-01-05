#include "acm_ga.h"
#include "mcmcmutator.h"

namespace ga
{

void HybridGA::evolve()
{
    if (populationSize() < 1)
    {
        fprintf(stderr, "No individuals at all!\n");
        return;
    }
    // Simplify syntax
    using alexandria::ACMIndividual;
    using alexandria::MCMCMutator;
    
    auto oldPop = oldPopPtr();
    oldPop->clear();
    // Initialize population/s
    for (int i = 0; i < populationSize(); i++)
    {
        oldPop->push_back(initializer()->initialize());
    }
    
    // Cast each individual to ACMIndividual for easier evolution
    std::vector<ACMIndividual*> acmPop;
    for (Individual *ind : *oldPop)
    {
        acmPop.push_back(static_cast<ACMIndividual*>(ind));
    }
    if (acmPop[0]->sii()->nParam() < 1)
    {
        fprintf(stderr, "Cannot evolve a chromosome without genes.\n");
        return;
    }
    
    // Open files of each individual
    for (ACMIndividual *ind : acmPop)
    {
        ind->openParamConvFiles(oenv());
        ind->openChi2ConvFile(oenv(), evaluateTestSet());
    }
    
    // Evolve each individual
    for (ACMIndividual *ind : acmPop)
    {
        mutator()->mutate(ind, evaluateTestSet());
    }

    // Collect results into the best individual
    int bestIndex = 0;
    double bestFitness = (*oldPop)[0]->fitnessTrain();
    for (int i = 0; i < populationSize(); i++)
    {
        if (acmPop[i]->fitnessTrain() < bestFitness)
        {
            bestFitness = acmPop[i]->fitnessTrain();
            bestIndex = i;
        }
    }
    setBestIndividual(acmPop[bestIndex]);
    
    // Close files of each individual
    for (ACMIndividual *ind : acmPop)
    {
        ind->closeConvFiles();
    }
}

void PureGA::evolve()
{
    if (gach_->popSize() < 2)
    {
        fprintf(stderr, "Need at least two individuals in the population.\n");
        return;
    }
    if (logFile())
    {
        fprintf(logFile(), "\nStarting GA/HYBRID evolution\n");
    }
    // Open surveillance files for fitness
    openFitnessFiles();
    
    // Random number generation
    std::random_device rd;  // Will be used to obtain a seed for the random number engine
    std::mt19937 gen(rd()); // Standard mersenne_twister_engine seeded with rd()
    std::uniform_real_distribution<double> dis(0.0, 1.0);
    
    // Iteration variables
    size_t i, k;
    
    // Indices for parents
    int parent1;
    int parent2;
    
    // Generations
    int generation = 0;
    if (logFile())
    {
        fprintf(logFile(), "\nGeneration %i\n", generation);
    
        // Initialize the population and compute fitness
        fprintf(logFile(), "Initializing individuals and computing initial fitness...\n");
    }
    auto oldPop = oldPopPtr();
    // Initialize population/s
    for (int i = 0; i < populationSize(); i++)
    {
        oldPop->push_back(initializer()->initialize());
        fitnessComputer()->compute((*oldPop)[i], Target::Train);
    }
    auto firstInd = static_cast<alexandria::ACMIndividual*>((*oldPop)[0]);
    if (firstInd->sii()->nParam() < 1)
    {
        fprintf(stderr, "Cannot evolve a chromosome without genes.\n");
        return;
    }
    fprintFitness();

    // FIXME: THIS IS NOT GENERAL. Open files of each individual
    for (Individual *ind : *oldPop)
    {
        alexandria::ACMIndividual *tmpInd = static_cast<alexandria::ACMIndividual*>(ind);
        tmpInd->openParamConvFiles(oenv());
        tmpInd->openChi2ConvFile(oenv(), evaluateTestSet());
    }
    
    // Copy individuals into newPop_
    copyOldToNewPopulations();
    
    // Initialize best individual
    auto bestIndex = findBestIndex();
    setBestIndividual((*oldPop)[bestIndex]);
    
    fprintPop();
    fprintBestIndInPop();
    fprintBestInd();
    
    // Iterate and create new generation
    do
    {
        // Increase generation counter
        generation++;
        if (logFile())
        {
            fprintf(logFile(), "\nGeneration %i\n", generation);
        
            // Sort individuals based on fitness
            fprintf(logFile(), "Sorting... (if needed)\n");
        }
        sorter()->sort(oldPop);
        fprintPop();
        
        // Normalize the fitness into a probability
        if (logFile())
        {
            fprintf(logFile(), "Computing probabilities...\n");
        }
        probabilityComputer()->compute(oldPop);
        if (logFile())
        {
            fprintProbability();
        
            // Move the "nElites" best individuals (unchanged) into the new population (assuming population is sorted)
            fprintf(logFile(), "Moving the %i best individual(s) into the new population...\n", gach_->nElites());
        }
        auto newPop = newPopPtr();
        for (i = 0; i < (size_t) gach_->nElites(); i++)
        {
            (*newPop)[i]->copyGenome(((*oldPop)[i]));
        }
        
        // Generate new population after the elitism
        if (logFile())
        {
            fprintf(logFile(), "Generating the rest of the new population...\n");
        }
        for (i = gach_->nElites(); i < oldPop->size(); i += 2)
        {
            if (logFile())
            {
                fprintf(logFile(), "i = %zu, %zu\n", i, i + 1);
            }
            // Select parents
            parent1 = selector()->select(*oldPop);
            parent2 = selector()->select(*oldPop);
            if (logFile())
            {
                fprintf(logFile(), "parent1: %i; parent2: %i\n", parent1, parent2);

                // Do crossover
                fprintf(logFile(), "Before crossover\n");
                fprintf(logFile(), "Parent 1:\n");
            }
            (*oldPop)[parent1]->fprintSelf(logFile());
            if (logFile())
            {
                fprintf(logFile(), "Parent 2:\n");
            }
            (*oldPop)[parent2]->fprintSelf(logFile());
            if (dis(gen) <= gach_->prCross())  // If crossover is to be performed
            {
                if (logFile())
                {
                    fprintf(logFile(), "Doing crossover...\n");
                }
                crossover()->offspring((*oldPop)[parent1], (*oldPop)[parent2], (*newPop)[i], (*newPop)[i+1]);
            }
            else
            {
                if (logFile())
                {
                    fprintf(logFile(), "Omitting crossover...\n");
                }
                (*newPop)[i]->copyGenome((*oldPop)[parent1]);
                (*newPop)[i+1]->copyGenome((*oldPop)[parent2]);
            }
            if (logFile())
            {
                fprintf(logFile(), "Child 1:\n");
            }
            (*newPop)[i]->fprintSelf(logFile());
            if (logFile())
            {
                fprintf(logFile(), "Child 2:\n");
            }
            (*newPop)[i+1]->fprintSelf(logFile());

            // Do mutation in each child
            if (logFile())
            {
                fprintf(logFile(), "Doing mutation...\n");
            }
            for (k = 0; k < 2; k++)
            {
                mutator()->mutate((*newPop)[i + k], gach_->prMut());
            }
            if (logFile())
            {
                fprintf(logFile(), "Child 1:\n");
            }
            (*newPop)[i]->fprintSelf(logFile());
            if (logFile())
            {
                fprintf(logFile(), "Child 2:\n");
            }
            (*newPop)[i+1]->fprintSelf(logFile());
        }

        // Swap oldPop and newPop
        if (logFile())
        {
            fprintf(logFile(), "Swapping oldPop and newPop...\n");
        }
        swapOldNewPopulations();
        
        // Compute fitness
        if (logFile())
        {
            fprintf(logFile(), "Computing fitness of new generation...\n");
        }
        for (i = 0; i < oldPop->size(); i++)
        {
            fitnessComputer()->compute((*oldPop)[i], Target::Train);
        }
        fprintFitness();

        fprintPop();
        fprintBestIndInPop();

        // Check if a better individual was found, and update if so
        Individual *tmpBest = (*oldPop)[findBestIndex()];
        if (tmpBest->fitnessTrain() < bestInd()->fitnessTrain())  // If we have a new best
        {
            if (logFile())
            {
                fprintf(logFile(), "A new best individual has been found!\n");
                fprintf(logFile(), "Previous best:\n");
            }
            bestInd()->fprintSelf(logFile());
            if (logFile())
            {
                fprintf(logFile(), "New best:\n");
            }
            setBestIndividual(tmpBest->clone());
            if (logFile())
            {
                bestInd()->fprintSelf(logFile());
            }
        }
        else
        {
            // New best not found...
            if (logFile())
            {
                fprintf(logFile(), "HAVE NOT FOUND a new best individual...\n");
            }
        }

        if (logFile())
        {
            fprintf(logFile(), "Checking termination conditions...\n");
        }
    }
    while (!terminator()->terminate(*oldPop, generation));
    
    // FIXME: THIS IS NOT GENERAL. Close files of each individual
    for (Individual *ind : *oldPop)
    {
        static_cast<alexandria::ACMIndividual*>(ind)->closeConvFiles();
    }
    
    // Close surveillance files for fitness
    closeFitnessFiles();
    
    if (logFile())
    {
        fprintf(logFile(), "\nGA/HYBRID Evolution is done!\n");
    }
    fprintBestInd();
    
}

} // namespace ga

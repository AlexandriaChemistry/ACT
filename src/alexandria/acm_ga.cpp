#include "acm_ga.h"

#include <algorithm>

#include "ga/GenePool.h"

#include "gmx_simple_comm.h"
#include "mcmcmutator.h"
#include "tune_ff.h"

namespace ga
{

bool MCMC::evolve(ga::Genome *bestGenome)
{
    if (populationSize() != 1)
    {
        fprintf(stderr, "MCMC needs to have exactly one individual!\n");
        return false;
    }
    if (sii_->nParam() < 1)
    {
        fprintf(stderr, "Cannot evolve a chromosome without genes.\n");
        return false;
    }
    // Simplify syntax, create individual
    auto *ind = static_cast<alexandria::ACMIndividual *>(initializer()->initialize());
    
    mutator()->mutate(ind->genomePtr(), bestGenome, evaluateTestSet_);
    mutator()->stopHelpers();
    mutator()->finalize();
    
    delete ind;
    
    return mutator()->foundMinimum();
}

bool HybridGAMC::evolve(ga::Genome *bestGenome)
{
    if (gach_->popSize() < 2)
    {
        fprintf(stderr, "Need at least two individuals in the population.\n");
        return false;
    }
    auto cr = sii_->commrec();
    if (cr->nmiddlemen < 1)
    {
        fprintf(stderr, "Need at least two cores/processes to run the genetic algorithm.\n");
        return false; 
    }
    if (logFile_)
    {
        fprintf(logFile_, "\nStarting GA/HYBRID evolution\n");
    }
    // Open surveillance files for fitness
    openFitnessFiles();
    
    // Random number generation
    std::random_device rd;  // Will be used to obtain a seed for the random number engine
    std::mt19937 gen(rd()); // Standard mersenne_twister_engine seeded with rd()
    std::uniform_real_distribution<double> dis(0.0, 1.0);
    
    // Indices for parents
    size_t parent1, parent2;
    
    // Generations
    int generation = 0;
    if (logFile_)
    {
        fprintf(logFile_, "\nGeneration %i\n", generation);
    
        // Initialize the population and compute fitness
        fprintf(logFile_, "Initializing individuals and computing initial fitness...\n");
    }
#ifdef OLD
    auto oldPop = oldPopPtr();
    // Initialize population only for this node
    for (int i = 0; i < populationSize(); i++)
    {
        oldPop->push_back(initializer()->initialize());
        fitnessComputer()->compute((*oldPop)[i], Target::Train);
    }
    auto firstInd = static_cast<alexandria::ACMIndividual*>((*oldPop)[0]);
#endif

    // Create the gene pools
    GenePool *pool[2];
    int       pold = 0;
#define pnew (1-pold)
    pool[pold] = new GenePool(sii_->nParam());
    pool[pnew] = new GenePool(sii_->nParam()); 
    
    // Load the initial genomes from the middlemen. 
    // This is needed since they have read their own parameters
    // from the Poldata structures.
    for(int i = 0; i < cr->nmiddlemen; i++)
    {
        int src = middleManGlobalIndex(cr, i);
        ga::Genome genome;
        genome.Receive(cr, src);
        pool[pold]->addGenome(genome);
        pool[pnew]->addGenome(genome);
    }
    // Now we have filled the gene pool and initial fitness values
    pool[pold]->print(logFile_);
    fprintFitness(*(pool[pold]));

    auto bestIndex = pool[pold]->findBestIndex();
    // TODO Check whether we need to update this at all here
    *bestGenome    = pool[pold]->genome(bestIndex);
    
    // When random initialization, assume a better minimum has been found no matter what
    bool bMinimum = gach_->randomInit() ? true : false;
    // Iterate and create new generation
    do
    {
        // Increase generation counter
        generation++;
        fprintf(logFile_, "\nGeneration %i\n", generation);
        
        // Sort individuals in increasing order of fitness
        auto gp = pool[pold]->genePoolPtr();
        if (gach_->sort())
        {
            fprintf(logFile_, "Sorting old population...\n");
            std::sort(gp->begin(), gp->end(), 
                      [](const Genome &a, const Genome &b) -> bool
                      { return a.fitness(iMolSelect::Train) < b.fitness(iMolSelect::Train); });
        }
        
        // Normalize the fitness into a probability
        fprintf(logFile_, "Computing probabilities...\n");
        probabilityComputer()->compute(gp);
        
        pool[pold]->print(logFile_);
        
        if (gach_->nElites() > 0)
        {
            // Move the "nElites" best individuals (unchanged) into the new population 
            // (assuming population is sorted)
            fprintf(logFile_, "Moving the %i best individual(s) into the new population...\n", gach_->nElites());
            for (int i = 0; i < gach_->nElites(); i++)
            {
                pool[pnew]->replaceGenome(i, pool[pold]->genome(i));
            }
            
            // Generate new population after the elitism
            fprintf(logFile_, "Generating the rest of the new population...\n");
        }
        for (size_t i = gach_->nElites(); i < pool[pold]->popSize(); i += 2)
        {
            // Select parents
            parent1 = selector()->select(gp);
            parent2 = selector()->select(gp);
            fprintf(logFile_, "parent1: %zu; parent2: %zu\n", parent1, parent2);
            
            // If crossover is to be performed
            if (dis(gen) <= gach_->prCross())  
            {
                // Do crossover
                fprintf(logFile_, "Before crossover\n");
                pool[pold]->genome(parent1).print("Parent 1:", logFile_);
                pool[pold]->genome(parent2).print("Parent 2:", logFile_);
                
                fprintf(logFile_, "Doing crossover...\n");
                crossover()->offspring(pool[pold]->genomePtr(parent1),
                                       pool[pold]->genomePtr(parent2),
                                       pool[pnew]->genomePtr(i),
                                       pool[pnew]->genomePtr(i+1));
            }
            else
            {
                fprintf(logFile_, "Omitting crossover...\n");
                pool[pnew]->replaceGenome(i,   pool[pold]->genome(parent1));
                pool[pnew]->replaceGenome(i+1, pool[pold]->genome(parent2));
            }
            pool[pnew]->genomePtr(i)->unsetFitness(iMolSelect::Train);
            pool[pnew]->genomePtr(i+1)->unsetFitness(iMolSelect::Train);
            pool[pnew]->genome(i).print("Child 1:", logFile_);
            pool[pnew]->genome(i+1).print("Child 2:", logFile_);

            // Do mutation in each child
            fprintf(logFile_, "Doing mutation...\n");
            for (size_t k = 0; k < 2; k++)
            {
                mutator()->mutate(pool[pnew]->genomePtr(i + k), bestGenome, gach_->prMut());
            }
            pool[pnew]->genome(i).print("Child 1:", logFile_);
            pool[pnew]->genome(i+1).print("Child 2:", logFile_);
        }

        // Swap oldPop and newPop
        fprintf(logFile_, "Swapping oldPop and newPop...\n");
        pold = pnew;
        
        // Now time to send out the new genomes to the individuals
        for(int i = 0; i < cr->nmiddlemen; i++)
        {
            int dest = middleManGlobalIndex(cr, i);
            // Send a 1 to tell the middlemen to continue
            gmx_send_int(cr, dest, 1);
            // Send the data set, 0 for iMolSelect::Train
            gmx_send_int(cr, dest, 0);
            // Now send the new bases
            gmx_send_double_vector(cr, dest, pool[pold]->genomePtr(i)->basesPtr());
        }
        fprintf(logFile_, "Computing fitness of new generation...\n");
        // And receive back the updated fitnesses
        for(int i = 0; i < cr->nmiddlemen; i++)
        {
            int src = middleManGlobalIndex(cr, i);
            auto fitness = gmx_recv_double(cr, src);
            pool[pold]->genomePtr(i)->setFitness(iMolSelect::Train, fitness);
        }
        fprintFitness(*(pool[pold]));
        
        // Check if a better genome was found, and update if so
        size_t newBest = pool[pold]->findBestIndex();
        if (pool[pold]->genome(newBest).fitness(iMolSelect::Train) < 
            bestGenome->fitness(iMolSelect::Train))  // If we have a new best
        {
            bestGenome->print("A new best individual has been found!\nPrevious best:\n",
                              logFile_);
            *bestGenome = pool[pold]->genome(newBest); 
            bestGenome->print("New best:\n", logFile_);
            bMinimum = true;
        }
    }
    while (!terminator()->terminate(pool[pold], generation));
    
    // Close surveillance files for fitness
    closeFitnessFiles();
    
    if (logFile_)
    {
        fprintf(logFile_, "\nGA/HYBRID Evolution is done!\n");
    }
    bestGenome->print("Best: ", logFile_);

    return bMinimum;    
}

} // namespace ga

#include "acm_ga.h"

#include <algorithm>

#include "act/ga//GenePool.h"

#include "act/utility/communicationrecord.h"
#include "mcmcmutator.h"
#include "tune_ff.h"

namespace ga
{

bool MCMC::evolve(ga::Genome *bestGenome)
{
    if (sii_->nParam() < 1)
    {
        fprintf(stderr, "Cannot evolve a chromosome without genes.\n");
        return false;
    }
    
    auto cr = sii_->commRec();

    // Create a gene pool
    GenePool pool(sii_->nParam());
    // Create and add our own individual (Will be the first one in the pool)
    auto *ind = static_cast<alexandria::ACMIndividual *>(initializer()->initialize());
    // Compute its fitness
    fitnessComputer()->compute(ind->genomePtr(), iMolSelect::Train);
    pool.addGenome(ind->genome());
    // Receive initial genomes from middlemen
    for (auto &src : cr->middlemen())
    {
        if (src != cr->rank())
        {
            ga::Genome genome;
            genome.Receive(cr, src);
            pool.addGenome(genome);
        }
    }
    GMX_RELEASE_ASSERT(static_cast<int>(pool.popSize()) == gach_->popSize(),
                       "The initial population does not match the specified population size...");
    // Print the genomes to the logfile
    pool.print(logFile_);

    // Update best genome
    auto bestIndex = pool.findBestIndex();
    *bestGenome    = pool.genome(bestIndex);

    // When random initialization, assume a better minimum has been found no matter what
    bool bMinimum = gach_->randomInit() ? true : false;

    // Resend the genomes back to the middlemen (they expect them anyway...)
    for (size_t i = 1; i < pool.popSize(); i++)
    {
        int dest = cr->middlemen()[i-1];
        // Tell the middle man to continue
        cr->send_data(dest);
        // Send the data set
        cr->send_iMolSelect(dest, iMolSelect::Train);
        // Now resend the bases
        cr->send_double_vector(dest, pool.genomePtr(i)->basesPtr());
    }

    // Mutate my own genome
    mutator()->mutate(ind->genomePtr(), ind->bestGenomePtr(), gach_->prMut());
    // Bring it into the population
    // FIXME: what if -bForceOutput? Make it sensitive to the flag
    pool.replaceGenome(0, ind->bestGenome());

    // Fetch the mutated genomes and their fitness. FIXME: use them when -bForceOutput
    for (size_t i = 1; i < pool.popSize(); i++)
    {
        int src      = cr->middlemen()[i-1];
        cr->recv_double_vector(src, pool.genomePtr(i)->basesPtr());  // Receiving the mutated parameters
        auto fitness = cr->recv_double(src);  // Receiving the new training fitness
        pool.genomePtr(i)->setFitness(iMolSelect::Train, fitness);
    }

    // Fetch the best genomes FIXME: use them when -nobForceOutput
    for (size_t i = 1; i < pool.popSize(); i++)
    {
        int src = cr->middlemen()[i-1];
        pool.genomePtr(i)->Receive(cr, src);
    }

    // Print the genomes to the logfile
    pool.print(logFile_);

    // Check if a better genome was found, and update if so
    auto imstr = iMolSelect::Train;
    size_t newBest = pool.findBestIndex();
    if (newBest < pool.popSize() &&
        pool.genome(newBest).hasFitness(imstr))
    {
        if (pool.genome(newBest).fitness(imstr) < 
            bestGenome->fitness(imstr))  // If we have a new best
        {
            bestGenome->print("A new best individual has been found!\nPrevious best:\n",
                              logFile_);
            *bestGenome = pool.genome(newBest); 
            bestGenome->print("New best:\n", logFile_);
            fprintf(logFile_, "\nMCMC Statistics for the master node only\n");
            auto mymut = reinterpret_cast<alexandria::MCMCMutator *>(mutator());
            mymut->printMonteCarloStatistics(logFile_, ind->initialGenome(),
                                             *bestGenome);

            bMinimum = true;
        }
    }
    else
    {
        fprintf(stderr, "No best genome in pool. WTF?\n");
    }

    // delete ind;  // FIXME: compilation warning about non-virtual destructor
    
    return bMinimum;
}

bool HybridGAMC::evolve(ga::Genome *bestGenome)
{
    if (gach_->popSize() < 2)  // FIXME: This is already checked in GAConfigHandler::check_pargs()
    {
        fprintf(stderr, "Need at least two individuals in the population.\n");
        return false;
    }
    auto cr = sii_->commRec();
    if (cr->nmiddlemen() < 2)  // FIXME: have we already checked that the number of processors is the correct one?
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
    gen.seed(seed_);
    
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

    // Create the gene pools
    GenePool *pool[2];
    int       pold = 0;
#define pnew (1-pold)
    pool[pold] = new GenePool(sii_->nParam());
    pool[pnew] = new GenePool(sii_->nParam()); 
    
    // Create and add our own individual (Will be the first one in the pool)
    auto *ind = static_cast<alexandria::ACMIndividual *>(initializer()->initialize());
    // Compute its fitness
    fitnessComputer()->compute(ind->genomePtr(), iMolSelect::Train);
    pool[pold]->addGenome(ind->genome());
    pool[pnew]->addGenome(ind->genome());
    // Load the initial genomes from the middlemen. 
    // This is needed since they have read their own parameters
    // from the Poldata structures.
    for(auto &src : cr->middlemen())
    {
        ga::Genome genome;
        genome.Receive(cr, src);
        pool[pold]->addGenome(genome);
        pool[pnew]->addGenome(genome);
    }
    // Now we have filled the gene pool and initial fitness values
    pool[pold]->print(logFile_);
    // Print fitness to surveillance files
    fprintFitness(*(pool[pold]));

    auto bestIndex = pool[pold]->findBestIndex();
    // TODO: Check whether we need to update this at all here
    *bestGenome    = pool[pold]->genome(bestIndex);
    
    // When random initialization, assume a better minimum has been found no matter what
    bool bMinimum = gach_->randomInit() ? true : false;
    // Iterate and create new generation
    do
    {
        // Increase generation counter
        generation++;
        if (debug)
        {
            fprintf(debug, "\nGeneration %i\n", generation);
        }
        // Sort individuals in increasing order of fitness
        auto gp = pool[pold]->genePoolPtr();
        if (gach_->sort())
        {
            if (debug)
            {
                fprintf(debug, "Sorting old population...\n");
            }
            pool[pold]->sort(iMolSelect::Train);
        }
        
        // Normalize the fitness into a probability
        if (debug)
        {
            fprintf(debug, "Computing probabilities...\n");
        }
        probabilityComputer()->compute(gp);
        
        if (debug)
        {
            pool[pold]->print(debug);
        }
        if (gach_->nElites() > 0)  // FIXME: can we remove the > 0?
        {
            // Move the "nElites" best individuals (unchanged) into the new population 
            // (assuming population is sorted)
            if (debug)
            {
                fprintf(debug, "Moving the %i best individual(s) into the new population...\n", gach_->nElites());
            }
            for (int i = 0; i < gach_->nElites(); i++)
            {
                pool[pnew]->replaceGenome(i, pool[pold]->genome(i));
            }
            
            // Generate new population after the elitism
            if (debug)
            {
                fprintf(logFile_, "Generating the rest of the new population...\n");
            }
        }
        for (size_t i = gach_->nElites(); i < pool[pold]->popSize(); i += 2)
        {
            // Select parents
            parent1 = selector()->select(gp);
            parent2 = selector()->select(gp);
            if (debug)
            {
                fprintf(debug, "parent1: %zu; parent2: %zu\n", parent1, parent2);
            }
            
            // If crossover is to be performed
            if (dis(gen) <= gach_->prCross())  
            {
                // Do crossover
                if (debug)
                {
                    fprintf(debug, "Before crossover\n");
                    pool[pold]->genome(parent1).print("Parent 1:", debug);
                    pool[pold]->genome(parent2).print("Parent 2:", debug);
                }
                if (debug)
                {
                    fprintf(debug, "Doing crossover...\n");
                }
                crossover()->offspring(pool[pold]->genomePtr(parent1),
                                       pool[pold]->genomePtr(parent2),
                                       pool[pnew]->genomePtr(i),
                                       pool[pnew]->genomePtr(i+1));
            }
            else
            {
                if (debug)
                {
                    fprintf(debug, "Omitting crossover...\n");
                }
                pool[pnew]->replaceGenome(i,   pool[pold]->genome(parent1));
                pool[pnew]->replaceGenome(i+1, pool[pold]->genome(parent2));
            }
            pool[pnew]->genomePtr(i)->unsetFitness(iMolSelect::Train);
            pool[pnew]->genomePtr(i+1)->unsetFitness(iMolSelect::Train);
            if (debug)
            {
                pool[pnew]->genome(i).print("Child 1:", debug);
                pool[pnew]->genome(i+1).print("Child 2:", debug);
            }
        }
        if (debug)
        {
            fprintf(debug, "Sending for mutation...\n");
        }
        for (size_t i = std::max(1, gach_->nElites()); i < pool[pnew]->popSize(); i++)
        {
            int dest = cr->middlemen()[i-1];
            // Signify the middlemen to continue
            cr->send_data(dest);
            // Send the data set 
            cr->send_iMolSelect(dest, iMolSelect::Train);
            // Now send the new bases
            cr->send_double_vector(dest, pool[pnew]->genomePtr(i)->basesPtr());
        }
        // Mutate the MASTER's genome if no elitism
        if (gach_->nElites() == 0)  // FIXME: can we just negate instead of comparing?
        {
            if (debug)
            {
                fprintf(debug, "Mutating the MASTER's genome...\n");
            }
            mutator()->mutate(pool[pnew]->genomePtr(0), ind->bestGenomePtr(), gach_->prMut());
            if (gach_->optimizer() == alexandria::OptimizerAlg::GA)
            {
                fitnessComputer()->compute(pool[pnew]->genomePtr(0), iMolSelect::Train);
            }
        }
        if (debug)
        {
            fprintf(debug, "Fetching mutated children and fitness from new generation...\n");
        }
        // Receive the new children (parameters + fitness) from the middle men for the non elitist
        // FIXME: if we end up sending more stuff, it might be worth it to just send the entire genome
        for (size_t i = std::max(1, gach_->nElites()); i < pool[pnew]->popSize(); i++)
        {
            int src      = cr->middlemen()[i-1];
            cr->recv_double_vector(src, pool[pnew]->genomePtr(i)->basesPtr());  // Receiving the mutated parameters
            auto fitness = cr->recv_double(src);  // Receiving the new training fitness
            pool[pnew]->genomePtr(i)->setFitness(iMolSelect::Train, fitness);
        }

        // Swap oldPop and newPop
        if (debug)
        {
            fprintf(debug, "Swapping oldPop and newPop...\n");
        }
        pold = pnew;

        // Print population again!
        if (debug)
        {
            pool[pold]->print(debug);
        }

        // Print fitness to surveillance files
        fprintFitness(*(pool[pold]));
        
        // Check if a better genome was found, and update if so
        auto   imstr   = iMolSelect::Train;
        size_t newBest = pool[pold]->findBestIndex();
        if (newBest < pool[pold]->popSize() &&
            pool[pold]->genome(newBest).hasFitness(imstr))
        {
            if (pool[pold]->genome(newBest).fitness(imstr) < 
                bestGenome->fitness(imstr))  // If we have a new best
            {
                bestGenome->print("A new best individual has been found!\nPrevious best:\n",
                                  logFile_);
                *bestGenome = pool[pold]->genome(newBest); 
                bestGenome->print("New best:\n", logFile_);
                bMinimum = true;
            }
        }
        else
        {
            fprintf(stderr, "No best genome in pool. Que?\n");
        }
    }
    while (!terminate(pool[pold], generation));
    
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

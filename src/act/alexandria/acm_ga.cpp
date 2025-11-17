/*
 * This source file is part of the Alexandria Chemistry Toolkit.
 *
 * Copyright (C) 2021-2025
 *
 * Developers:
 *             Mohammad Mehdi Ghahremanpour, 
 *             Julian Marrades,
 *             Marie-Madeleine Walz,
 *             Paul J. van Maaren, 
 *             David van der Spoel (Project leader)
 *
 * This program is free software; you can redistribute it and/or
 * modify it under the terms of the GNU General Public License
 * as published by the Free Software Foundation; either version 2
 * of the License, or (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this program; if not, write to the Free Software
 * Foundation, Inc., 51 Franklin Street, Fifth Floor, 
 * Boston, MA  02110-1301, USA.
 */
#include "acm_ga.h"

#include <algorithm>
#include <ctime>
#include <regex>

#include "act/basics/msg_handler.h"
#include "act/ga/gene_pool.h"
#include "act/utility/communicationrecord.h"
#include "mcmcmutator.h"
#include "train_ff.h"

namespace ga
{

bool MCMC::evolve(alexandria::MsgHandler       *msghandler,
                  std::map<iMolSelect, Genome> *bestGenome)
{
    if (sii_->nParam() < 1)
    {
        msghandler->msg(alexandria::ACTStatus::Error,
                        "Cannot evolve a chromosome without genes.\n");
        return false;
    }
    
    auto cr = sii_->commRec();

    // Tell the middleman no genepool was read.
    // TODO: Implement genepool reading in MCMC
    int read = 0;
    for(auto &ii : cr->middlemen())
    {
        cr->send(ii, read);
    }
    // Dataset(s)
    const auto imstr = iMolSelect::Train;
    const auto imste = iMolSelect::Test;

    // Create a gene pool
    GenePool pool(sii_->nParam());
    // Create and add our own individual (Will be the first one in the pool)
    auto *ind = std::move(static_cast<alexandria::ACMIndividual *>(initializer()->initialize()));

    // Compute its fitness
    auto tw = msghandler->tw();
    msghandler->msg(alexandria::ACTStatus::Info,
                    "MASTER's initial parameter vector chi2 components:");

    fitnessComputer()->compute(msghandler, ind->genomePtr(), imstr);
    // Not really needed but just to print the components
    fitnessComputer()->compute(msghandler, ind->genomePtr(), imste);
    
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
    // Print the genomes to the debug stream if requested.
    if (msghandler->debug())
    {
        for(const auto &p: pool.print())
        {
            msghandler->tw()->writeString(p);
        }
    }

    // Update best genome
    (*bestGenome)[imstr] = pool.getBest(imstr);

    // When random initialization, assume a better minimum has been found no matter what
    bool bMinimum = gach_->randomInit() ? true : false;

    // Resend the genomes back to the middlemen (they expect them anyway...)
    int i = 1;
    for (auto &dest : cr->middlemen())
    {
        if (dest != cr->rank())
        {
            // Tell the middle man to continue
            cr->send_data(dest);
            // Send the data set
            cr->send(dest, imstr);
            // Now resend the bases
            cr->send(dest, pool.genomePtr(i)->bases());
            // Tell the middleman to carry the MUTATION mode
            cr->send(dest, alexandria::TrainFFMiddlemanMode::MUTATION);
            i += 1;
        }
    }

    // Mutate my own genome
    mutator()->mutate(msghandler, ind->genomePtr(), ind->bestGenomePtr(), gach_->prMut());
    // Bring it into the population
    pool.replaceGenome(0, ind->bestGenome());

    // Fetch the mutated genomes and their fitness.
    for (size_t i = 1; i < pool.popSize(); i++)
    {
        int src      = cr->middlemen()[i-1];
        // Receiving the mutated parameters
        cr->recv(src, pool.genomePtr(i)->basesPtr());
        // Receiving the new training fitness
        double fitness;
        cr->recv(src, &fitness);
        pool.genomePtr(i)->setFitness(imstr, fitness);
    }
    // Print the genomes to the logfile
    if (msghandler->verbose())
    {
        for(const auto &p : pool.print())
        {
            msghandler->write(p);
        }
    }

    // Check if a better genome was found, and update if so
    const auto tmpGenome = pool.getBest(imstr);
    const auto tmpBest   = bestGenome->find(imstr)->second;
    if (tmpGenome.fitness(imstr) < tmpBest.fitness(imstr))  // If we have a new best
    {
        (*bestGenome)[imstr] = tmpGenome;
        msghandler->write(tmpGenome.print("New best for train"));
        msghandler->write("MCMC Statistics for the master node only");
        auto mymut = reinterpret_cast<alexandria::MCMCMutator *>(mutator());
        mymut->printMonteCarloStatistics(tw, ind->initialGenome(), tmpGenome);

        bMinimum = true;
    }
    else
    {
        msghandler->msg(alexandria::ACTStatus::Error,
                        "No best genome in pool. WTF?\n");
    }

    // Save last population
    lastPop_ = pool;
    // Clean
    delete ind;
    return bMinimum;
}

bool HybridGAMC::evolve(alexandria::MsgHandler       *msghandler,
                        std::map<iMolSelect, Genome> *bestGenome)
{
    auto cr = sii_->commRec();
    // FIXME: have we already checked that the number of processors is the correct one?
    if (cr->nmiddlemen() < 2)
    {
        fprintf(stderr, "Need at least two cores/processes to run the genetic algorithm.\n");
        return false; 
    }
    msghandler->write("\nStarting GA/HYBRID evolution\n");
    // Random number generation
    std::random_device rd;  // Will be used to obtain a seed for the random number engine
    std::mt19937 gen(rd()); // Standard mersenne_twister_engine seeded with rd()
    std::uniform_real_distribution<double> dis(0.0, 1.0);
    gen.seed(seed_);
    
    // Datasets
    const auto imstr = iMolSelect::Train;
    const auto imste = iMolSelect::Test;

    // Generations
    int generation = 0;
    // Initialize the population and compute fitness
    msghandler->write("Initializing individuals and computing initial fitness...");

    // Create the gene pools
    GenePool *pool[2];
    int       pold = 0;
#define pnew (1-pold)
    pool[pold] = new GenePool(sii_->nParam());
    pool[pnew] = new GenePool(sii_->nParam()); 
    
    // Create and add our own individual (Will be the first one in the pool)
    auto *ind = static_cast<alexandria::ACMIndividual *>(initializer()->initialize());
    // Check whether we need to read a gene pool
    int read = 0;
    if (gpin_)
    {
        pool[pold]->read(gpin_);
        if (pool[pold]->popSize() - gach_->popSize() != 0)
        {
            GMX_THROW(gmx::InvalidInputError(gmx::formatString("Read a gene pool with %zu individuals but that does not match the population size of %d", pool[pold]->popSize(), gach_->popSize()).c_str()));
        }
        for(int ii = 0; ii < gach_->popSize(); ++ii)
        {
            if (pool[pold]->genome(ii).nBase() != sii_->nParam())
            {
                GMX_THROW(gmx::InvalidInputError(gmx::formatString("Genome %d in genepool read from '%s' has %zu elements, but I expected %zu", ii, gpin_, pool[pold]->genome(ii).nBase(), sii_->nParam()).c_str()));
            }
        }
        *pool[pnew] = *pool[pold];
        msghandler->write(gmx::formatString("\nRead gene pool with %zu individuals and %zu bases from file %s\n\n", 
                                            pool[pold]->popSize(), pool[pold]->genome(0).nBase(), gpin_));

        ind->copyGenome(pool[pold]->genome(0));
        read = 1;
        msghandler->msg(alexandria::ACTStatus::Debug,
                        gmx::formatString("Will send genomes to %zu middlemen\n", cr->middlemen().size()));
    }
    for(auto &ii : cr->middlemen())
    {
        cr->send(ii, read);
        if (read == 1)
        {
            // Genome index is not the same thing as middleman index...
            // https://github.com/dspoel/ACT/issues/560
            int genome_index = ii/(1+cr->nhelper_per_middleman());
            cr->send(ii, pool[pold]->genomePtr(genome_index)->bases());
            msghandler->msg(alexandria::ACTStatus::Debug,
                            gmx::formatString("Sent genome to middleman %d\n", ii));
        }
    }

    // Compute its fitness
    msghandler->write("MASTER's initial parameter vector chi2 components:");

    fitnessComputer()->compute(msghandler, ind->genomePtr(), imstr);
    // Maybe not really needed but just to print the components
    fitnessComputer()->compute(msghandler, ind->genomePtr(), imste);

    if (read == 0)
    {
        pool[pold]->addGenome(ind->genome());
        pool[pnew]->addGenome(ind->genome());
    
        // I. Match numbers to corresponding ones in actmiddleman.cpp
        // Load the initial genomes from the middlemen. 
        // This is needed since they have read their own parameters
        // from the ForceField structures.
        for(auto &src : cr->middlemen())
        {
            ga::Genome genome;
            genome.Receive(cr, src);
            pool[pold]->addGenome(genome);
            pool[pnew]->addGenome(genome);
        }
    }
    if (msghandler->info())
    {
        // Now we have filled the gene pool and initial fitness values
        for(const auto &p : pool[pold]->print())
        {
            msghandler->write(p);
        }
    }
    // Initialize bestGenome
    (*bestGenome)[imstr] = pool[pold]->getBest(imstr);
    if (gach_->evaluateTestset())
    {
        (*bestGenome)[imste] = pool[pold]->getBest(imstr);
    }
    // Open surveillance files for fitness, after initialization of bestGenome
    for(const auto &bg : *bestGenome)
    {
        openFitnessFiles(fitnessFile_, bg.first);
    }
    // Print original fitness to surveillance files
    fprintFitness(*(pool[pold]));

    // When random initialization, assume a better minimum has been found no matter what
    bool bMinimum = gach_->randomInit() ? true : false;
    msghandler->write(gmx::formatString("Starting %d generations of force field training.\n",
                                        gach_->maxGenerations()));

    // Store starting time
    time_t start_time = std::time(nullptr);
    // Iterate and create new generation
    do
    {
        // Sort individuals in increasing order of fitness
        auto gp = pool[pold]->genePoolPtr();
        if (gach_->sort())
        {
            pool[pold]->sort(imstr);
            if (msghandler->debug())
            {
                msghandler->msg(alexandria::ACTStatus::Debug,
                                "Sorting old population...\n");
                for(const auto &p : pool[pold]->print())
                {
                    msghandler->tw()->writeString(p);
                }
            }
        }

        // Penalize
        if (penalize(msghandler->tw(), pool[pold], generation))
        {
            // Send to middlemen for fitness recomputation
            for (size_t i = 1; i < pool[pold]->popSize(); i++)
            {
                int dest = cr->middlemen()[i-1];
                // II.
                // Signify the middlemen to continue
                cr->send_data(dest);
                // Send the data set FIXME: is this necessary? It will always be train?
                cr->send(dest, imstr);
                // Now send the new bases
                cr->send(dest, pool[pold]->genomePtr(i)->bases());
                // III.
                // Tell the middleman to carry the MUTATION mode
                cr->send(dest, alexandria::TrainFFMiddlemanMode::FITNESS);
            }
            // Recompute my fitness
            {
                auto genome0 = pool[pold]->genomePtr(0);
                genome0->unsetFitness(imstr);
                fitnessComputer()->compute(msghandler, genome0, imstr);
            }
            // Receive fitness from middlemen
            for (size_t i = 1; i < pool[pold]->popSize(); i++)
            {
                int src = cr->middlemen()[i-1];
                double fitness;
                cr->recv(src, &fitness);  // Receiving the new training fitness
                pool[pold]->genomePtr(i)->setFitness(imstr, fitness);
            }
            // Sort again if needed
            if (gach_->sort())
            {
                pool[pold]->sort(imstr);
            }
            // Print population to debug if we have penalized the population
            if (msghandler->debug())
            {
                msghandler->msg(alexandria::ACTStatus::Debug,
                                "Population has been penalized!\n");
                for(const auto &p : pool[pold]->print())
                {
                    msghandler->write(p);
                }
            }
        }

        // Increase generation counter
        generation++;

        // Normalize the fitness into a probability
        if (debug)
        {
            fprintf(debug, "Computing probabilities...\n");
        }
        probabilityComputer()->compute(gp, generation);
        
        if (msghandler->debug())
        {
            for(const auto &p : pool[pold]->print())
            {
                msghandler->write(p);
            }
        }
        if (gach_->nElites() > 0)
        {
            // Move the "nElites" best individuals (unchanged) into the new population 
            // (assuming population is sorted)
            if (msghandler->debug())
            {
                msghandler->msg(alexandria::ACTStatus::Debug,
                                gmx::formatString("Moving the %i best individual(s) into the new population...\n", gach_->nElites()));
            }
            for (int i = 0; i < gach_->nElites(); i++)
            {
                pool[pnew]->replaceGenome(i, pool[pold]->genome(i));
            }
            
            // Generate new population after the elitism
            msghandler->msg(alexandria::ACTStatus::Debug,
                            "Generating the rest of the new population...\n");
        }
        for (size_t i = gach_->nElites(); i < pool[pold]->popSize(); i += 2)
        {
            // Select parents
            auto parent1 = selector()->select(gp);
            // We do not want the same two parents for a child.
            // https://github.com/dspoel/ACT/issues/549
            auto parent2 = parent1;
            int  mytry   = 0;
            int  maxtry  = 5000;
            while (parent1 == parent2 && mytry < maxtry)
            {
                parent2 = selector()->select(gp);
                mytry++;
            }
            if (mytry == maxtry)
            {
                parent2 = (parent1 + 1) % pool[pold]->popSize();
            }
            auto child1  = i;
            auto child2  = i+1;
            if (debug)
            {
                fprintf(debug, "parent1: %d parent2: %d child1: %zu child2: %zu\n", parent1, parent2, child1, child2);
            }
            
            // If crossover is to be performed
            if (dis(gen) <= gach_->prCross())  
            {
                // Do crossover
                if (msghandler->debug())
                {
                    msghandler->write("Before crossover\n");
                    msghandler->write(pool[pold]->genome(parent1).print("Parent 1:"));
                    msghandler->write(pool[pold]->genome(parent2).print("Parent 2:"));
                    msghandler->write("Doing crossover...\n");
                }
                crossover()->offspring(pool[pold]->genomePtr(parent1),
                                       pool[pold]->genomePtr(parent2),
                                       pool[pnew]->genomePtr(child1),
                                       pool[pnew]->genomePtr(child2));
            }
            else
            {
                if (msghandler->debug())
                {
                    msghandler->write("Omitting crossover...\n");
                }
                pool[pnew]->replaceGenome(child1, pool[pold]->genome(parent1));
                pool[pnew]->replaceGenome(child2, pool[pold]->genome(parent2));
            }
            pool[pnew]->genomePtr(child1)->unsetFitness(imstr);
            pool[pnew]->genomePtr(child2)->unsetFitness(imstr);
            if (gach_->evaluateTestset())
            {
                pool[pnew]->genomePtr(child1)->unsetFitness(imste);
                pool[pnew]->genomePtr(child2)->unsetFitness(imste);
            }
            if (msghandler->debug())
            {
                msghandler->write(pool[pnew]->genome(child1).print("Child 1:"));
                msghandler->write(pool[pnew]->genome(child2).print("Child 2:"));
            }
        }
        if (msghandler->debug())
        {
            msghandler->msg(alexandria::ACTStatus::Debug, "Sending for mutation...");
        }

        for (size_t i = std::max(1, gach_->nElites()); i < pool[pnew]->popSize(); i++)
        {
            int dest = cr->middlemen()[i-1];
            // II.
            // Signify the middlemen to continue
            cr->send_data(dest);
            // Send the data set FIXME: is this necessary? It will always be train?
            cr->send(dest, imstr);
            // Now send the new bases
            cr->send(dest, pool[pnew]->genomePtr(i)->bases());
            // III.
            // Tell the middleman to carry the MUTATION mode
            cr->send(dest, alexandria::TrainFFMiddlemanMode::MUTATION);
        }
        // Mutate the MASTER's genome if no elitism
        if (gach_->nElites() == 0)
        {
            if (debug)
            {
                fprintf(debug, "Mutating the MASTER's genome...\n");
            }
            auto g0ptr = pool[pnew]->genomePtr(0);
            mutator()->mutate(msghandler, g0ptr, ind->bestGenomePtr(), gach_->prMut());
            if (mutator()->foundMinimum())
            {
                // Store master's best genome in the pool at position 0
                pool[pnew]->replaceGenome(0, ind->bestGenome());
            }
            if (gach_->optimizer() == alexandria::OptimizerAlg::GA)
            {
                // For HYBRID the fitness is already computed by the mutator
                fitnessComputer()->compute(msghandler, g0ptr, imstr);
            }
            if (gach_->evaluateTestset())
            {
                fitnessComputer()->compute(msghandler, g0ptr, imste);
            }
        }
        if (msghandler->debug())
        {
            msghandler->write("Fetching mutated children and fitness from new generation...\n");
        }
        // Receive the new children (parameters + fitness) from the middle men for the
        // non elitist.
        // FIXME: if we end up sending more stuff, it might be worth it to just send the entire genome
        for (size_t i = std::max(1, gach_->nElites()); i < pool[pnew]->popSize(); i++)
        {
            int src      = cr->middlemen()[i-1];
            // IV.
            // Receive the mutated parameters
            cr->recv(src, pool[pnew]->genomePtr(i)->basesPtr());
            // V.
            // and the fitness.
            double fitness;
            cr->recv(src, &fitness);
            pool[pnew]->genomePtr(i)->setFitness(imstr, fitness);
            if (gach_->evaluateTestset())
            {
                double fitnessTest;
                cr->recv(src, &fitnessTest);  // Receiving the new test fitness
                pool[pnew]->genomePtr(i)->setFitness(imste, fitnessTest);
            }
        }

        // Swap oldPop and newPop
        if (msghandler->debug())
        {
            msghandler->write("Swapping oldPop and newPop...\n");
        }
        pold = pnew;

        // Print population again!
        if (msghandler->debug())
        {
            for(const auto &p : pool[pold]->print())
            {
                msghandler->write(p);
            }
        }

        // Print fitness to surveillance files
        fprintFitness(*(pool[pold]));

        // Check if a better genome (for train) was found, and update if so
        const auto tmpGenome   = pool[pold]->getBest(imstr);
        const auto tmpBest     = (*bestGenome)[imstr];
        time_t my_time     = std::time(nullptr);
        time_t diff_time   = my_time - start_time;
        time_t finish_time = my_time + (diff_time / generation) * (gach_->maxGenerations() - generation);
        std::string t_now(std::ctime(&my_time));
        std::string t_finish(std::ctime(&finish_time));
        // regex code from https://www.systutorials.com/how-to-remove-newline-characters-from-a-string-in-c/
        std::regex newlines_re("\n+");
        auto mess_start = gmx::formatString("Generation %d/%d. Time: %s, expect to finish at %s",
                                            generation, gach_->maxGenerations(),
                                            std::regex_replace(t_now, newlines_re, "").c_str(),
                                            std::regex_replace(t_finish, newlines_re, "").c_str());

        if (tmpGenome.fitness(imstr) < tmpBest.fitness(imstr))  // If we have a new best
        {
            auto mess = gmx::formatString("\n%s. New best individual for train", mess_start.c_str());
            msghandler->write(tmpGenome.print(mess.c_str()));
            (*bestGenome)[imstr] = tmpGenome;
            bMinimum = true;

            if (gpout_)
            {
                pool[pold]->write(gpout_);
            }
        }
        else
        {
            msghandler->write(mess_start);
        }
        if (gach_->evaluateTestset())
        {
            const auto tmpBestTest = (*bestGenome)[imste];
            if (tmpGenome.fitness(imste) < tmpBestTest.fitness(imste))
            {
                auto mess = gmx::formatString("%s. New best individual for test", mess_start.c_str());
                (*bestGenome)[imste] = tmpGenome;
                msghandler->write(tmpGenome.print(mess.c_str()));
            }
        }
        msghandler->flush();
    }
    while (!terminate(msghandler->tw(), pool[pold], generation));

    // Close surveillance files for fitness
    closeFitnessFiles();

    msghandler->write("GA/HYBRID Evolution is done!");

    msghandler->write((*bestGenome)[imstr].print("Best (Train): "));
    if (gach_->evaluateTestset())
    {
        msghandler->write((*bestGenome)[imste].print("Best (Test): "));
    }

    // Save last population
    lastPop_ = *(pool[pold]);
    // Clean
    for (int i = 0; i < 2; i++)
    {
        delete pool[i];
    }

    return bMinimum;
}

} // namespace ga

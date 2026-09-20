/*
 * This source file is part of the Alexandria Chemistry Toolkit.
 *
 * Copyright (C) 2021-2026
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
#include "mcmc.h"

#include <algorithm>
#include <ctime>
#pragma GCC diagnostic ignored "-Wpragmas"
#pragma GCC diagnostic ignored "-Wunknown-warning-option"
#pragma GCC diagnostic ignored "-Wmaybe-uninitialized"
#include <regex>
#pragma GCC diagnostic pop

#include "act/basics/msg_handler.h"
#include "act/forces/forcecomputerstatistics.h"
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
    //! \todo: Implement genepool reading in MCMC
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
    auto ind = static_cast<alexandria::ACMIndividual *>(initializer()->initialize());

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
        // Tell the middle man to continue
        cr->send_data(dest);
        // Now resend the genome
        pool.genomePtr(i)->Send(cr, dest);
        // Tell the middleman to carry the MUTATION mode
        cr->send(dest, alexandria::TrainFFMiddlemanMode::MUTATION);
        i += 1;
    }

    // Mutate my own genome
    ind->setBestGenome(ind->genome());
    mutator()->mutate(msghandler, ind->genomePtr(), ind->bestGenomePtr(), gach_->prMut());
    // Bring it into the population
    pool.replaceGenome(0, ind->bestGenome());

    // Fetch the mutated genomes and their fitness.
    for (size_t i = 1; i < pool.popSize(); i++)
    {
        int src      = cr->middlemen()[i-1];
        // Receiving the mutated genome
        pool.genomePtr(i)->Receive(cr, src);
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
        msghandler->msg(alexandria::ACTStatus::Warning,
                        "No better genome found. Please check your input and output.\n");
        bMinimum = false;
    }
    // ForceComputer Statistics
    alexandria::ForceComputerStatistics fcStats;
    if (msghandler->info())
    {
        std::string stats = fcStats.statistics(cr,
                                               static_cast<const alexandria::ACMFitnessComputer *>(fitnessComputer())->forceComputer(),
                                               0,
                                               true);
        msghandler->msg(alexandria::ACTStatus::Info, stats);
    }

    // Save last population
    lastPop_ = pool;
    // Clean
    delete ind;
    return bMinimum;
}

} // namespace ga

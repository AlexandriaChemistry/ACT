/*
 * This source file is part of the Alexandria Chemistry Toolkit.
 *
 * Copyright (C) 2014-2025
 *
 * this program is free software; you can redistribute it and/or
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
#include "actmiddleman.h" 

#include "act/alexandria/acminitializer.h"
#include "act/alexandria/mcmcmutator.h"
#include "act/alexandria/percentmutator.h"
#include "act/alexandria/staticindividualinfo.h"
#include "act/utility/communicationrecord.h"
#include "gromacs/utility/gmxassert.h"
#include "gromacs/utility/textwriter.h"

namespace alexandria
{
 
ACTMiddleMan::ACTMiddleMan(MsgHandler                *msghandler,
                           MolGen                    *mg,
                           StaticIndividualInfo      *sii,
                           GAConfigHandler           *gach,
                           BayesConfigHandler        *bch,
                           gmx_output_env_t          *oenv,
                           ChargeGenerationAlgorithm  algorithm,
                           bool                       openConvFiles)
: gach_(gach), id_(sii->commRec()->middleManOrdinal())
{
    // This ica
    int seed;
    sii->commRec()->recv(sii->commRec()->superior(), &seed);
    // Standard mersenne_twister_engine seeded with what we received
    std::mt19937 gen(seed);
    // Default constructor to cover all available (positive) range
    std::uniform_int_distribution<int> dis(0);
    
    auto initializer = new ACMInitializer(sii, gach->randomInit(), dis(gen));
    
    // Create and initialize the individual
    ind_ = static_cast<ACMIndividual *>(initializer->initialize());

    // Create force computer
    forceComp_ = new ForceComputer(bch->shellToler(), bch->shellMaxIter());

    // Fitness computer FIXME: what about those false flags?
    fitComp_ = new ACMFitnessComputer(msghandler, sii, mg, false, forceComp_, algorithm);
    
    if (gach->optimizer() == OptimizerAlg::GA)
    {
        mutator_ = new alexandria::PercentMutator(sii, dis(gen), algorithm, gach->percent());
    }
    else
    {
        // FIXME: we need to make some logfiles for the middlemen, because they apparently cannot write to the global logfile
        auto mut = new alexandria::MCMCMutator(dis(gen),
                                               bch, fitComp_, sii,
                                               bch->evaluateTestset(),
                                               algorithm, 
                                               gach->maxGenerations());
        if (openConvFiles)
        {
            mut->openParamConvFiles(oenv);
            mut->openChi2ConvFile(oenv);
        }
        mutator_ = mut;
    }
}
    
void ACTMiddleMan::run(MsgHandler *msghandler)
{
    auto cr = ind_->sii()->commRec();
    GMX_RELEASE_ASSERT(cr->isMiddleMan(), "I thought I was the middle man...");
    // Check whether we need to accept a genome from the master
    int read;
    int master = cr->superior();
    cr->recv(master, &read);
    if (read == 1)
    {
        cr->recv(master, ind_->genomePtr()->basesPtr());
    }
    // Start by computing my own fitness
    fitComp_->compute(msghandler, ind_->genomePtr(), iMolSelect::Train);
    if (gach_->evaluateTestset())
    {
        fitComp_->compute(msghandler, ind_->genomePtr(), iMolSelect::Test);
    }
    // I.
    // Send my initial genome and fitness to the master if needed
    if (read == 0)
    {
        ind_->genome().Send(cr, master);
    }
    auto cont = CommunicationStatus::OK;

    // II.
    // Do we need to continue?
    cont = cr->recv_data(master);
    while (CommunicationStatus::RECV_DATA == cont)
    {
        // Get the dataset
        // FIXME: is this really necessary?
        iMolSelect ims;
        cr->recv(master, &ims);
        
        // Now get the new bases.
        cr->recv(master, ind_->genomePtr()->basesPtr());
        
        // III.
        // Receive assignment from master.  
        TrainFFMiddlemanMode mode;
        cr->recv(master, &mode);
        if (mode == TrainFFMiddlemanMode::MUTATION)
        {
            mutator_->mutate(msghandler, ind_->genomePtr(), ind_->bestGenomePtr(), 
                             gach_->prMut());
            
            if (gach_->optimizer() == OptimizerAlg::GA)
            {
                fitComp_->compute(msghandler, ind_->genomePtr(), ims);
                // IV.
                // Send the mutated vector
                cr->send(master, *ind_->genomePtr()->basesPtr());

                // V.
                // Send the new train fitness
                cr->send(master, ind_->genome().fitness(ims));
                if (gach_->evaluateTestset())
                {
                    fitComp_->compute(msghandler, ind_->genomePtr(), iMolSelect::Test);
                    cr->send(master, ind_->genome().fitness(iMolSelect::Test));
                }
            }
            else
            {
                // IV.
                // If we are working with Hybrid or MCMC, send 
                // the best genome found to the MASTER.
                cr->send(master, *ind_->bestGenomePtr()->basesPtr());
                //ind_->bestGenome().Send(cr, master);
                // V.
                // Send the new train fitness
                cr->send(master, ind_->bestGenome().fitness(ims));
                if (gach_->evaluateTestset())
                {
                    fitComp_->compute(msghandler, ind_->bestGenomePtr(), iMolSelect::Test);
                    cr->send(master, ind_->bestGenome().fitness(iMolSelect::Test));
                }
            }
        }
        else if (mode == TrainFFMiddlemanMode::FITNESS)
        {
            fitComp_->compute(msghandler, ind_->genomePtr(), ims);
            cr->send(master, ind_->genome().fitness(ims));
        }
        cont = cr->recv_data(master);
    }
    
    // Stop my helpers too.
    stopHelpers();
}

void ACTMiddleMan::printStatistics(gmx::TextWriter *tw)
{
   if (gach_->optimizer() == OptimizerAlg::MCMC && tw)
   {       
       tw->writeStringFormatted("Middle man %i\n", id_);
       static_cast<MCMCMutator *>(mutator_)->printMonteCarloStatistics(tw, ind_->initialGenome(),
                                                                       ind_->bestGenome());
   }
}

void ACTMiddleMan::stopHelpers()
{
    (void) fitComp_->distributeTasks(CalcDev::Stop);
}

} // namespace alexandria

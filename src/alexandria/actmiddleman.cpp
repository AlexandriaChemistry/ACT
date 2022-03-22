#include "actmiddleman.h" 

#include "gromacs/utility/gmxassert.h"

#include "acminitializer.h"
#include "act/utility/communicationrecord.h"
#include "mcmcmutator.h"
#include "percentmutator.h"
#include "staticindividualinfo.h"

namespace alexandria
{
 
ACTMiddleMan::ACTMiddleMan(MolGen               *mg,
                           StaticIndividualInfo *sii,
                           GAConfigHandler      *gach,
                           BayesConfigHandler   *bch,
                           bool                  flush,
                           gmx_output_env_t     *oenv)
: gach_(gach), sii_(sii), id_(sii->commRec()->middleManOrdinal())
{
    // This ica
    int seed = bch->seed();
    // Create random number generator and feed it the global seed
    std::random_device rd;  // Will be used to obtain a seed for the random number engine
    std::mt19937 gen(rd()); // Standard mersenne_twister_engine seeded with rd()
    std::uniform_int_distribution<int> dis(0); // Default constructor to cover all available (positive) range
    gen.seed(seed);
    // Use the random number generator to get a seed for this processor based on the global seed
    // Skip the first "nodeid" numbers and grab the next one
    for (int i = 0; i < sii->commRec()->middleManOrdinal(); i++)
    {
        dis(gen);
    }
    seed = dis(gen);
    
    auto initializer = new ACMInitializer(sii, gach->randomInit(), seed);
    
    // Create and initialize the individual
    ind_ = static_cast<ACMIndividual *>(initializer->initialize());

    // Fitness computer FIXME: what about those false flags?
    fitComp_ = new ACMFitnessComputer(nullptr, false, sii, mg, false);
    // Create and initialize the mutator
    sii->makeIndividualDir();  // We need to call this before opening working files!
    if (gach->optimizer() == OptimizerAlg::GA)
    {
        mutator_ = new alexandria::PercentMutator(sii, seed, gach->percent());
    }
    else
    {
        // FIXME: we need to make some logfiles for the middlemen, because they apparently cannot write to the global logfile
        auto mut = new alexandria::MCMCMutator(nullptr, false, flush, seed, bch, fitComp_, sii, bch->evaluateTestset());
        mut->openParamConvFiles(oenv);
        mut->openChi2ConvFile(oenv);
        mutator_ = mut;
    }
}
    
void ACTMiddleMan::run()
{
    auto cr = ind_->sii()->commRec();
    GMX_RELEASE_ASSERT(cr->isMiddleMan(), "I thought I was the middle man...");
    // Start by computing my own fitness
    fitComp_->compute(ind_->genomePtr(), iMolSelect::Train);
    if (gach_->evaluateTestset())
    {
        fitComp_->compute(ind_->genomePtr(), iMolSelect::Test);
    }
    // The send my initial genome and fitness to the master
    ind_->genome().Send(cr, 0);
    auto cont = CommunicationStatus::OK;
    
    cont = cr->recv_data(0);
    while (CommunicationStatus::RECV_DATA == cont)
    {
        // Get the dataset
        // FIXME: is this really necessary?
        iMolSelect ims = cr->recv_iMolSelect(0);
        
        // Now get the parameters
        cr->recv_double_vector(0, ind_->genomePtr()->basesPtr());
        
        TuneFFMiddlemanMode mode = cr->recv_ff_middleman_mode(0);
        if (mode == TuneFFMiddlemanMode::MUTATION)
        {
            mutator_->mutate(ind_->genomePtr(), ind_->bestGenomePtr(), gach_->prMut());
            
            if (gach_->optimizer() == OptimizerAlg::GA)
            {
                fitComp_->compute(ind_->genomePtr(), ims);
            }
            
            // Send the mutated vector
            cr->send_double_vector(0, ind_->genomePtr()->basesPtr());
            
            // Send the new train fitness
            cr->send_double(0, ind_->genome().fitness(ims));
            if (gach_->evaluateTestset() && gach_->optimizer() != OptimizerAlg::MCMC)
            {
                fitComp_->compute(ind_->genomePtr(), iMolSelect::Test);
                cr->send_double(0, ind_->genome().fitness(iMolSelect::Test));
            }
            
            // If we are working with MCMC, send the best found to the MASTER
            if (gach_->optimizer() == OptimizerAlg::MCMC)
            {
                ind_->bestGenome().Send(cr, 0);
            }
        }
        else if (mode == TuneFFMiddlemanMode::FITNESS)
        {
            fitComp_->compute(ind_->genomePtr(), ims);
            cr->send_double(0, ind_->genome().fitness(ims));
        }
        cont = cr->recv_data(0);
    }
    
    // Stop my helpers too.
    stopHelpers();
    // Save the last genome
    sii_->saveState(true, sii_->outputFileLast());
}

void ACTMiddleMan::printStatistics(FILE *logFile)
{
   if (gach_->optimizer() == OptimizerAlg::MCMC && logFile)
   {       
       fprintf(logFile, "Middle man %i\n", id_);
       static_cast<MCMCMutator *>(mutator_)->printMonteCarloStatistics(logFile, ind_->initialGenome(),
                                                                       ind_->bestGenome());
   }
}

void ACTMiddleMan::stopHelpers()
{
    std::vector<double> dummy;
    fitComp_->calcDeviation(&dummy,
                            alexandria::CalcDev::Final, iMolSelect::Train);
}

} // namespace alexandria

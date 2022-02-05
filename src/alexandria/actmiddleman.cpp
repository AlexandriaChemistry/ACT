#include "actmiddleman.h" 

#include "gromacs/utility/gmxassert.h"

#include "acminitializer.h"
#include "communicationrecord.h"
#include "mcmcmutator.h"
#include "percentmutator.h"
#include "staticindividualinfo.h"

namespace alexandria
{
 
ACTMiddleMan::ACTMiddleMan(FILE                 *logFile,
                           MolGen               *mg,
                           StaticIndividualInfo *sii,
                           GAConfigHandler      *gach,
                           BayesConfigHandler   *bch,
                           bool                  verbose,
                           gmx_output_env_t     *oenv)
: gach_(gach), logFile_(logFile), id_(sii->commRec()->middleManOrdinal())
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
    fitComp_ = new ACMFitnessComputer(nullptr, sii, mg, 
                                      false, false, false);
    // Create and initialize the mutator
    if (gach->optimizer() == OptimizerAlg::GA)
    {
        mutator_ = new alexandria::PercentMutator(sii, seed, gach->percent());
    }
    else
    {
        auto mut = new alexandria::MCMCMutator(nullptr, verbose, seed, bch, fitComp_, sii, bch->evaluateTestset());
        mut->makeWorkDir();  // We need to call this before opening working files!
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
    // The send my initial genome and fitness to the master
    ind_->genome().Send(cr, 0);
    auto cont = CommunicationStatus::OK;
    
    do
    {
        cont = cr->recv_data(0);
        if (cont == CommunicationStatus::RECV_DATA)  // FIXME: reverse condition and save an identation
        {
            // Get the dataset
            iMolSelect ims = cr->recv_iMolSelect(0);

            // Now get the parameters
            cr->recv_double_vector(0, ind_->genomePtr()->basesPtr());
            
            mutator_->mutate(ind_->genomePtr(), ind_->bestGenomePtr(), gach_->prMut());

            if (gach_->optimizer() == OptimizerAlg::GA)
            {
                fitComp_->compute(ind_->genomePtr(), ims);
            }

            // Send the mutated vector
            cr->send_double_vector(0, ind_->genomePtr()->basesPtr());

            // Send the new fitness
            cr->send_double(0, ind_->genome().fitness(ims));

            // If we are working with MCMC, send the best found to the MASTER
            if (gach_->optimizer() == OptimizerAlg::MCMC)
            {
                ind_->bestGenome().Send(cr, 0);
            }
        }
    }
    while (CommunicationStatus::RECV_DATA == cont);
    // TODO: Print Monte Carlo statistics if necessary
    // FIXME: in HYBRID this won't be correct
    // FIXME: do this in MASTER too
    // if (gach_->optimizer() != OptimizerAlg::GA)
    // {
    //     fprintf(logFile_, "Middle man %i\n", id_);
    //     static_cast<MCMCMutator *>(mutator_)->printMonteCarloStatistics(logFile_, ind_->initialGenome(),
    //                                                                     ind_->bestGenome());
    // }
    // Stop my helpers too.
    stopHelpers();
}

void ACTMiddleMan::stopHelpers()
{
    std::vector<double> dummy;
    fitComp_->calcDeviation(&dummy,
                            alexandria::CalcDev::Final, iMolSelect::Train);
}

} // namespace alexandria

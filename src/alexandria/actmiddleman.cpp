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
                           const std::string    &outputFile,
                           bool                  verbose,
                           gmx_output_env_t     *oenv)
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
    for (int i = 0; i < sii->commRec()->rank(); i++)
    {
        dis(gen);
    }
    seed = dis(gen);
    
    auto initializer = new ACMInitializer(sii, gach->randomInit(), outputFile, seed);
    
    // Create and initialize the individual
    ind_ = static_cast<ACMIndividual *>(initializer->initialize());

    // Fitness computer
    fitComp_ = new ACMFitnessComputer(nullptr, sii, mg, 
                                      false, false, false);
    // Create and initialize the mutator
    if (gach->optimizer() == OptimizerAlg::GA)
    {
        mutator_ = new alexandria::PercentMutator(sii, gach->percent());
    }
    else
    {
        auto mut = new alexandria::MCMCMutator(logFile, verbose, bch, fitComp_, sii);
        mut->openParamConvFiles(oenv);
        mut->openChi2ConvFile(oenv, bch->evaluateTestset());
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
    int cont = 0;
    do
    {
        cont = cr->recv_int(0);
        if (cont)
        {
            int imsi = cr->recv_int(0);
            // Get the dataset
            iMolSelect ims;
            switch (imsi)
            {
            case 0:
                ims = iMolSelect::Train;
                break;
            case 1:
                ims = iMolSelect::Test;
                break;
            case 2:
                ims = iMolSelect::Ignore;
                break;
            }
            // Now get the parameters
            cr->recv_double_vector(0, ind_->genomePtr()->basesPtr());
            fitComp_->compute(ind_->genomePtr(), ims);
            cr->send_double(0, ind_->genome().fitness(ims));
        }
    }
    while (cont);
    // Close our files or whaterver we need to do, then we're done!
    mutator_->finalize();
}

} // namespace alexandria

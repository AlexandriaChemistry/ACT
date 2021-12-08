#include "mcmcmutator.h"


namespace alexandria
{


void MCMCMutator::changeParam(ACMIndividual *ind,
                              size_t         j)
{
    std::vector<double> *param = ind->paramPtr();
    
    GMX_RELEASE_ASSERT(j < param->size(), "Parameter out of range");
    real delta = (2*randNum()-1) * bch_->step() * (sii_->upperBound()[j] - sii_->lowerBound()[j]);
    (*param)[j] += delta;
    if (sii_->mutability()[j] == Mutability::Bounded)
    {
        if ((*param)[j] < sii_->lowerBound()[j])
        {
            (*param)[j] = sii_->lowerBound()[j];
        }
        else if ((*param)[j] > sii_->upperBound()[j])
        {
            (*param)[j] = uii_->upperBound()[j];
        }
    }
}


void MCMCMutator::printMonteCarloStatistics(ACMIndividual  *ind,
                                            FILE           *fp)
{
    if (!fp)
    {
        return;
    }

    const auto bestParam = ind->bestParam();
    const auto pmean = ind->pMean();
    const auto psigma = ind->pSigma();
    const auto paramNames = sii_->paramNames();
    const auto attemptedMoves = ind->attemptedMoves();
    const auto acceptedMoves = ind->acceptedMoves();
    const auto ntrain = sii_->nTrain();
    const auto initialParam = ind->initialParam();
    const auto weightedTemperature = sii_->weightedTemperature();

    fprintf(fp, "\nMonte Carlo statistics of parameters after optimization\n");
    fprintf(fp, "#best %zu #mean %zu #sigma %zu #param %zu\n",
            bestParam.size(), pmean.size(), psigma.size(), paramNames.size());
    if (bestParam.size() == ind->nParam())
    {
        fprintf(fp, "Parameter                     Ncopies Initial   Best    Mean    Sigma Attempt  Acceptance  T-Weight\n");
        for (size_t k = 0; k < ind->nParam(); k++)
        {
            double acceptance_ratio = 0;
            if (attemptedMoves[k] > 0)
            {
                acceptance_ratio = 100*(double(acceptedMoves[k])/attemptedMoves[k]);
            }
            fprintf(fp, "%-30s  %5d  %6.3f  %6.3f  %6.3f  %6.3f    %4d %5.1f%%  %10.5f\n",
                    paramNames[k].c_str(), ntrain[k],
                    initialParam[k], bestParam[k], pmean[k], psigma[k],
                    attemptedMoves[k], acceptance_ratio, weightedTemperature[k]);
        }
    }
}

void MCMCMutator::fprintNewMinimum(      ACMIndividual *ind,
                                   const bool           bEvaluate_testset,
                                   const double         xiter,
                                   const double         currEval,
                                   const double         currEval_testset)
{
    if (bEvaluate_testset)
    {
        fprintf(logfile_, "iter %10g. Found new minimum at %10g. Corresponding energy on the test set: %g\n",
                xiter, currEval, currEval_testset);
    }
    else
    {
        fprintf(logfile_, "iter %10g. Found new minimum at %10g\n",
                xiter, currEval);
    }
    if (debug)
    {
        ind->printParameters(debug);
    }
}             

void MCMCMutator::fprintParameterStep(      ACMIndividual   *ind,
                                      const double           xiter)
{

    const auto fpc = ind->fpc();
    const auto param = ind->param();
    const auto paramClassIndex = sii_->paramClassIndex();

    for(FILE *fp: fpc)  // Write iteration number to each parameter convergence surveillance file
    {
        fprintf(fp, "%8f", xiter);
    }
    for (size_t k = 0; k < param.size(); k++)  // Write value of each parameter to its respective surveillance file
    {
        fprintf(fpc[paramClassIndex[k]], "  %10g", param[k]);
    }
    for(FILE *fp: fpc)  // If verbose = True, flush the file to be able to add new data to surveillance plots
    {
        fprintf(fp, "\n");
        if (verbose_)
        {
            fflush(fp);
        }
    }

}                             

}
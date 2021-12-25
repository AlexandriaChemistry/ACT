/*! \internal \brief
 * Implements part of the alexandria program.
 * \author Julian Ramon Marrades Furquet <julianramon.marradesfurquet.8049@student.uu.se>
 * \author Oskar Tegby <oskar.tegby@it.uu.se>
 */


#include "mcmcmutator.h"

#include "bayes.h"
#include "memory_check.h"


namespace alexandria
{


void MCMCMutator::mutate(                 ga::Individual   *ind,
                         gmx_unused const double            prMut)
{
    
    ACMIndividual *tmpInd = static_cast<ACMIndividual*>(ind);

    double prevEval = 0;

    const size_t nParam = tmpInd->nParam();

    std::vector<bool> changed;
    // Initialize to true to make sure the parameters are 
    // all spread around the processors.
    changed.resize(nParam, true);
    tmpInd->toPoldata(changed);
    // Now set them back to false, further spreading is done
    // one parameter at a time.
    std::fill(changed.begin(), changed.end(), false);

    const std::vector<FILE*> fpc = tmpInd->fpc();
    FILE *fpe = tmpInd->fpe();
    const auto paramClassIndex = sii_->paramClassIndex();
    std::vector<double> *param = tmpInd->paramPtr();

    prevEval = fitComp_->calcDeviation(tmpInd, false, CalcDev::Parallel, iMolSelect::Train);

    print_memory_usage(debug);

    double beta0 = 1 / (BOLTZ * bch_->temperature());

    // Optimization loop
    for (int iter = 0; iter < bch_->maxIter(); iter++)
    {
        for (size_t pp = 0; pp < nParam; pp++)
        {
            // Do the step!
            stepMutation(tmpInd, param, &changed, &prevEval, pp, iter, &beta0, nParam, paramClassIndex);
        }
    }

}

void MCMCMutator::stepMutation(      ACMIndividual          *ind,
                                     std::vector<double>    *param,
                                     std::vector<bool>      *changed,
                                     double                 *prevEval,
                               const size_t                  pp,
                               const int                     iter,
                                     double                 *beta0,
                               const size_t                  nParam,
                               const std::vector<size_t>    &paramClassIndex)
{

    // Pick a random parameter index
    const size_t paramIndex = randIndex();

    // Store the original value of the parameter
    const double storeParam = (*param)[paramIndex];

    // Change the parameter
    changeParam(ind, paramIndex);
    (*changed)[paramIndex] = true;
    ind->toPoldata(*changed);

    // Evaluate energy
    const double currEval = fitComp_->calcDeviation(ind, false, CalcDev::Parallel, iMolSelect::Train);
    const double deltaEval = currEval - (*prevEval);

    // Accept any downhill move
    bool accept = (deltaEval < 0);
    // For an uphill move apply the Metropolis Criteria
    // to decide whether to accept or reject the new parameter
    if (!accept)
    {
        // Only anneal if the simulation reached a certain number of steps
        if (bch_->anneal(iter)) (*beta0) = bch_->computeBeta(iter);
        const double randProbability = randNum();
        const double mcProbability   = exp( - ( (*beta0) / (sii_->weightedTemperature())[paramIndex] ) * deltaEval );
        accept = (mcProbability > randProbability);
    }

    // Fractional iteration taking into account the inner loop with <pp> over <nParam>
    const double xiter = iter + (1.0*pp)/nParam;
    if (accept)
    {  // If the parameter change is accepted
        (*prevEval) = currEval;
    }
    else
    {  // If the parameter change is not accepted
        (*param)[paramIndex] = storeParam;  // Set the old value of the parameter back
        // poldata needs to change back as well!
        ind->toPoldata(*changed);
    }
    (*changed)[paramIndex] = false;  // Set changed[j] back to false for upcoming iterations

    fprintParameterStep(ind, xiter);
    fprintChi2Step(ind, false, xiter, *prevEval, 0);

}

bool MCMCMutator::MCMC(      ACMIndividual *ind,
                       const bool           evaluate_testset)
{

    int nsum = 0;
    const size_t nParam = ind->nParam();
    double minEval = 0;
    double prevEval = 0;
    double prevEval_testset = 0;

    std::vector<double> sum, sum_of_sq;
    sum.resize(nParam, 0);
    sum_of_sq.resize(nParam, 0);

    std::vector<bool> changed;
    // Initialize to true to make sure the parameters are 
    // all spread around the processors.
    changed.resize(nParam, true);
    ind->toPoldata(changed);
    // Now set them back to false, further spreading is done
    // one parameter at a time.
    std::fill(changed.begin(), changed.end(), false);

    const std::vector<FILE*> fpc = ind->fpc();
    FILE *fpe = ind->fpe();
    const auto paramClassIndex = sii_->paramClassIndex();
    std::vector<double> *param = ind->paramPtr();
    
    if (sii_->xvgConv().empty() || sii_->xvgEpot().empty())
    {
        gmx_fatal(FARGS, "You forgot to call setOutputFiles. Back to the drawing board.");
    }
    if (param->empty())
    {
        fprintf(stderr, "No parameters to optimize.\n");
        return 0;
    }

    // Gather initial chi2 evaluations for training and test sets
    prevEval = fitComp_->calcDeviation(ind, true, CalcDev::Parallel, iMolSelect::Train);  // This does not fill the fitness attribute in individual!
    minEval = prevEval;
    ind->setFitnessTrain(prevEval);
    if (evaluate_testset)
    {
        prevEval_testset = fitComp_->calcDeviation(ind, true, CalcDev::Parallel, iMolSelect::Test);
        ind->setFitnessTest(prevEval_testset);
    }

    print_memory_usage(debug);

    double beta0 = 1 / (BOLTZ * bch_->temperature());

    // Optimization loop
    for (int iter = 0; iter < bch_->maxIter(); iter++)
    {
        for (size_t pp = 0; pp < nParam; pp++)
        {
            // Do the step!
            stepMCMC(ind, param, &changed, &prevEval, &prevEval_testset,
                     evaluate_testset, pp, iter, &beta0, nParam, &minEval, paramClassIndex);

            // For the second half of the optimization, collect data to find the mean and standard deviation of each
            // parameter
            if (iter >= bch_->maxIter()/2)
            {
                for (size_t k = 0; k < nParam; k++)
                {
                    sum[k]       += (*param)[k];
                    sum_of_sq[k] += gmx::square((*param)[k]);
                }
                nsum++;
            }
        }
    }

    // OPTIMIZATION IS COMPLETE

    computeMeanSigma(ind->pMeanPtr(), ind->pSigmaPtr(), nParam, sum, nsum, &sum_of_sq);

    bool bMinimum = false; // Assume no better minimum was found
    if (minEval < ind->fitnessTrain())  // If better minimum was found, update the value in <*chi2> and return true
    {
        ind->setFitnessTrain(minEval);
        bMinimum = true;
    }
    return bMinimum;

}

void MCMCMutator::stepMCMC(      ACMIndividual          *ind,
                                 std::vector<double>    *param,
                                 std::vector<bool>      *changed,
                                 double                 *prevEval,
                                 double                 *prevEval_testset,
                           const bool                    evaluate_testset,
                           const size_t                  pp,
                           const int                     iter,
                                 double                 *beta0,
                           const size_t                  nParam,
                                 double                 *minEval,
                           const std::vector<size_t>    &paramClassIndex)
{

    // Get pointers for attempted and accepted moves
    std::vector<int> *attemptedMoves = ind->attemptedMovesPtr();
    std::vector<int> *acceptedMoves  = ind->acceptedMovesPtr();

    // Pick a random parameter index
    const size_t paramIndex = randIndex();

    // Store the original value of the parameter
    const double storeParam = (*param)[paramIndex];

    // Change the parameter
    changeParam(ind, paramIndex);

    (*attemptedMoves)[paramIndex] += 1;
    (*changed)[paramIndex]         = true;

    // Update FF parameter data structure with
    // the new value of parameter j
    ind->toPoldata(*changed);

    // Evaluate the energy on training set
    const double currEval = fitComp_->calcDeviation(ind, false, CalcDev::Parallel, iMolSelect::Train);
    const double deltaEval = currEval - (*prevEval);
    // Evaluate the energy on the test set only on whole steps!
    double currEval_testset = (*prevEval_testset);
    if (evaluate_testset && pp == 0)
    {
        currEval_testset = fitComp_->calcDeviation(ind, false, CalcDev::Parallel, iMolSelect::Test);
    }

    // Accept any downhill move
    bool accept = (deltaEval < 0);

    // For an uphill move apply the Metropolis Criteria
    // to decide whether to accept or reject the new parameter
    if (!accept)
    {
        // Only anneal if the simulation reached a certain number of steps
        if (bch_->anneal(iter)) (*beta0) = bch_->computeBeta(iter);
        const double randProbability = randNum();
        const double mcProbability   = exp( - ( (*beta0) / (sii_->weightedTemperature())[paramIndex] ) * deltaEval );
        accept = (mcProbability > randProbability);
    }

    // Fractional iteration taking into account the inner loop with <pp> over <nParam>
    const double xiter = iter + (1.0*pp)/nParam;
    if (accept)
    {  // If the parameter change is accepted
        if (currEval < (*minEval))
        {
            // If pointer to log file exists, write information about new minimum
            if (logfile_) fprintNewMinimum(ind, evaluate_testset, xiter, currEval, currEval_testset);
            ind->setBestParam(*param);
            (*minEval) = currEval;
            ind->saveState();
        }
        (*prevEval) = currEval;
        if (evaluate_testset)
        {
            (*prevEval_testset) = currEval_testset;
            ind->setFitnessTest(currEval_testset);
        }
        (*acceptedMoves)[paramIndex] += 1;
    }
    else
    {  // If the parameter change is not accepted
        (*param)[paramIndex] = storeParam;  // Set the old value of the parameter back
        // poldata needs to change back as well!
        ind->toPoldata(*changed);
    }
    (*changed)[paramIndex] = false;  // Set changed[j] back to false for upcoming iterations

    fprintParameterStep(ind, xiter);
    fprintChi2Step(ind, evaluate_testset, xiter, *prevEval, *prevEval_testset);

}                  

void MCMCMutator::computeMeanSigma(      std::vector<double>    *pmean,
                                         std::vector<double>    *psigma,
                                   const size_t                  nParam,
                                   const std::vector<double>    &sum,
                                   const int                     nsum,
                                         std::vector<double>    *sum_of_sq)
{
    if (nsum <= 0) return;
    double ps2 = 0.0;
    for (size_t k = 0; k < nParam; k++)
    {
        (*pmean)[k]        = (sum[k]/nsum);
        (*sum_of_sq)[k]   /= nsum;
        ps2                = std::max(0.0, (*sum_of_sq)[k]-gmx::square((*pmean)[k]));
        (*psigma)[k]       = sqrt(ps2);
    }
}                                         

void MCMCMutator::changeParam(ACMIndividual *ind,
                              size_t         j)
{
    std::vector<double> *param = ind->paramPtr();
    
    GMX_RELEASE_ASSERT(j < param->size(), "Parameter out of range");
    double rnd = randNum();
    while (rnd == 0.5) rnd = randNum();  // Make sure the parameter changes!
    real delta = (2*rnd-1) * bch_->step() * (sii_->upperBound()[j] - sii_->lowerBound()[j]);
    (*param)[j] += delta;
    if (sii_->mutability()[j] == Mutability::Bounded)
    {
        if ((*param)[j] < sii_->lowerBound()[j])
        {
            (*param)[j] = sii_->lowerBound()[j];
        }
        else if ((*param)[j] > sii_->upperBound()[j])
        {
            (*param)[j] = sii_->upperBound()[j];
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

    ind->printHeader(fp);
    fprintf(fp, "Monte Carlo statistics of parameters after optimization\n");
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
    ind->printHeader(logfile_);
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
        ind->printHeader(debug);
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

void MCMCMutator::fprintChi2Step(      ACMIndividual    *ind,
                                 const bool              bEvaluate_testset,
                                 const double            xiter,
                                 const double            prevEval,
                                 const double            prevEval_testset)
{
    auto fpe = ind->fpe();

    if (fpe == nullptr) return;  // If fpe is a null pointer, return

    if (bEvaluate_testset)
    {
        fprintf(fpe, "%8f  %10g  %10g\n", xiter, prevEval, prevEval_testset);
    }
    else
    {
        fprintf(fpe, "%8f  %10g\n", xiter, prevEval);
    }
    if (verbose_)
    {
        fflush(fpe);
    }
}                    

void MCMCMutator::sensitivityAnalysis(ACMIndividual  *ind,
                                      iMolSelect      ims)
{
    
    std::vector<double> *param_ = ind->paramPtr();
    const auto upperBound = sii_->upperBound();
    const auto lowerBound = sii_->lowerBound();
    const auto paramNames = sii_->paramNames();

    if (param_->size() == 0)
    {
        return;
    }
    std::vector<bool> changed;
    changed.resize(param_->size(), true);
    ind->toPoldata(changed);
    std::fill(changed.begin(), changed.end(), false);
    double chi2_0 = fitComp_->calcDeviation(ind, false, CalcDev::Parallel, ims);
    if (logfile_)
    {
        fprintf(logfile_, "\nStarting sensitivity analysis. chi2_0 = %g nParam = %d\n",
                chi2_0, static_cast<int>(param_->size()));
        fflush(logfile_);
    }
    for (size_t i = 0; i < param_->size(); ++i)
    {
        Sensitivity s;
        double pstore = (*param_)[i];
        double deltap = (upperBound[i]-lowerBound[i])/200;
        double pmin   = std::max((*param_)[i]-deltap, lowerBound[i]);
        double pmax   = std::min((*param_)[i]+deltap, upperBound[i]);
        double p_0    = 0.5*(pmin+pmax);
        changed[i]    = true;
        (*param_)[i]     = pmin;
        ind->toPoldata(changed);
        s.add((*param_)[i], fitComp_->calcDeviation(ind, false, CalcDev::Parallel, ims));
        (*param_)[i]     = p_0;
        ind->toPoldata(changed);
        s.add((*param_)[i], fitComp_->calcDeviation(ind, false, CalcDev::Parallel, ims));
        (*param_)[i]     = pmax;
        ind->toPoldata(changed);
        s.add((*param_)[i],  fitComp_->calcDeviation(ind, false, CalcDev::Parallel, ims));
        (*param_)[i]     = pstore;
        ind->toPoldata(changed);
        changed[i]    = false;
        s.computeForceConstants(logfile_);
        s.print(logfile_, paramNames[i]);
    }
    if (logfile_)
    {
        fprintf(logfile_, "Sensitivity analysis done.\n");
    }

}                                      


} //namespace alexandria
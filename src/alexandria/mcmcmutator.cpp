/*! \internal \brief
 * Implements part of the alexandria program.
 * \author Mohammad Mehdi Ghahremanpour <mohammad.ghahremanpour@icm.uu.se>
 * \author David van der Spoel <david.vanderspoel@icm.uu.se>
 * \author Julian Ramon Marrades Furquet <julian.marrades@hotmail.es>
 * \author Oskar Tegby <oskar.tegby@it.uu.se>
 */
#include "mcmcmutator.h"

#include <sys/types.h>
#include <sys/stat.h>

#include "gromacs/utility/basedefinitions.h"

#include "bayes.h"
#include "act/utility/memory_check.h"

namespace alexandria
{

MCMCMutator::MCMCMutator(FILE                 *logfile,
                         bool                  verbose,
                         int                   seed,
                         BayesConfigHandler   *bch,
                         ACMFitnessComputer   *fitComp,
                         StaticIndividualInfo *sii,
                         bool                  evaluateTestSet)
    : Mutator(seed), evaluateTestSet_(evaluateTestSet), gen(rd()), dis(std::uniform_int_distribution<size_t>(0, sii->nParam()-1))
{
    gen.seed(seed);
    logfile_ = logfile;
    verbose_ = verbose;
    bch_     = bch;
    fitComp_ = fitComp;
    sii_     = sii;
    size_t nParam = sii_->nParam();
    pMean_.resize(nParam, 0.0);
    pSigma_.resize(nParam, 0.0);
    attemptedMoves_.resize(nParam, 0);
    acceptedMoves_.resize(nParam, 0);
}

void MCMCMutator::mutate(ga::Genome        *genome,
                         ga::Genome        *bestGenome,
                         gmx_unused double  prMut)
{
    int    nsum             = 0;
    size_t nParam           = genome->nBase();

    std::vector<double> sum, sum_of_sq;
    sum.resize(nParam, 0);
    sum_of_sq.resize(nParam, 0);

    // changed.resize(nParam, true);
    // Update all parameters in poldata
    sii_->updatePoldata(genome);
    // Create and fill the changed vector
    // Now set them back to false, further spreading is done
    // one parameter at a time.
    // std::fill(changed.begin(), changed.end(), false);
    std::vector<bool> changed(nParam, false);
    // changed.resize(nParam, false);

    if (sii_->xvgConv().empty() || sii_->xvgEpot().empty())
    {
        gmx_fatal(FARGS, "You forgot to call setOutputFiles. Back to the drawing board.");
    }
    if (genome->nBase() == 0)
    {
        fprintf(stderr, "No parameters to optimize.\n");
        return;
    }

    // Gather initial chi2 evaluations for training and test sets
    std::map<iMolSelect, double> prevEval;
    prevEval.insert({iMolSelect::Train, fitComp_->calcDeviation(genome->basesPtr(),
                                                                CalcDev::Parallel,
                                                                iMolSelect::Train)});  
    // This does not fill the fitness attribute in individual!
    // TODO Print results
    genome->setFitness(iMolSelect::Train, prevEval[iMolSelect::Train]);
    if (evaluateTestSet_)
    {
        prevEval.insert({iMolSelect::Test, fitComp_->calcDeviation(genome->basesPtr(),
                                                                   CalcDev::Parallel,
                                                                   iMolSelect::Test)});
        // TODO printing
        genome->setFitness(iMolSelect::Test, prevEval[iMolSelect::Test]);
    }
    // Save initial evaluation and initialize a structure for the minimum evaluation
    auto initEval = prevEval;
    *bestGenome = *genome;

    print_memory_usage(debug);

    double beta0      = 1 / (BOLTZ * bch_->temperature());
    // Optimization loop
    int    iterOffset = myGeneration_*bch_->maxIter();
    for (int iter = 0; iter < bch_->maxIter(); iter++)
    {
        for (size_t pp = 0; pp < nParam; pp++)
        {
            // Do the step!
            stepMCMC(genome, bestGenome, &changed, &prevEval,
                     pp, iter, iterOffset, &beta0);

            // For the second half of the optimization, collect data 
            // to find the mean and standard deviation of each
            // parameter
            // TODO: make optional
            if (iter >= bch_->maxIter()/2)
            {
                size_t k = 0;
                for (auto &p : genome->bases())
                {
                    sum[k]       += p;
                    sum_of_sq[k] += gmx::square(p);
                    k            += 1;
                }
                nsum++;
            }
        }
    }
    myGeneration_ += 1;
    // OPTIMIZATION IS COMPLETE
    // TODO Make optional
    computeMeanSigma(sum, nsum, &sum_of_sq);

    // Assume no better minimum was found
    bMinimum_ = false;
    auto ims  = iMolSelect::Train;
    if (bestGenome->hasFitness(ims) &&
        bestGenome->fitness(ims) < initEval[ims])
    {
        // If better minimum was found, update the value in <*chi2> and return true
        // genome->setFitness(ims, minEval[ims]);
        bMinimum_ = true;
    }
}

void MCMCMutator::stepMCMC(ga::Genome                   *genome,
                           ga::Genome                   *bestGenome,
                           std::vector<bool>            *changed,
                           std::map<iMolSelect, double> *prevEval,
                           size_t                        pp,
                           int                           iter,
                           int                           iterOffset,
                           double                       *beta0)
{
    // Pick a random parameter index
    const size_t paramIndex = randIndex();

    // Shorthand for the parameters here
    auto param = genome->basesPtr();
    
    // Store the original value of the parameter
    const double storeParam = (*param)[paramIndex];

    // Change the parameter
    changeParam(genome, paramIndex);

    attemptedMoves_[paramIndex] += 1;
    (*changed)[paramIndex]         = true;

    // Update FF parameter data structure with
    // the new value of parameter j
    sii_->updatePoldata(*changed, genome);

    std::map<iMolSelect, double> currEval;
    // Evaluate the energy on training set
    auto ims = iMolSelect::Train;
    currEval.insert({ims,
            fitComp_->calcDeviation(param, CalcDev::Parallel, ims) });
    double deltaEval = currEval[ims] - prevEval->find(ims)->second;
    // Evaluate the energy on the test set only on whole steps!
    if (evaluateTestSet_ && pp == 0)
    {
        ims            = iMolSelect::Test;
        // Is the first time we insert a test value?
        if (currEval.find(ims) == currEval.end())
        {
            double oldTest = currEval[iMolSelect::Train];
            auto pef = prevEval->find(ims);
            if (prevEval->end() != pef)
            {
                oldTest = pef->second;
            }
            currEval.insert({ims, oldTest});
        }
        currEval[ims] = fitComp_->calcDeviation(param, CalcDev::Parallel, ims);
    }

    // Accept any downhill move
    bool accept = (deltaEval < 0);

    // For an uphill move apply the Metropolis Criteria
    // to decide whether to accept or reject the new parameter
    if (!accept)
    {
        // Only anneal if the simulation reached a certain number of steps
        if (bch_->anneal(iter))
        {
            *beta0 = bch_->computeBeta(iter);
        }
        const double randProbability = randNum();
        const double mcProbability   = exp( - ( (*beta0) / (sii_->weightedTemperature())[paramIndex] ) * deltaEval );
        accept = (mcProbability > randProbability);
    }

    // Fractional iteration taking into account the inner loop with <pp> over <nParam>
    const double xiter = iterOffset + iter + (1.0*pp)/genome->nBase();
    if (accept)
    {  // If the parameter change is accepted
        auto imstr = iMolSelect::Train;
        if (!bestGenome->hasFitness(imstr) || 
            currEval[imstr] < bestGenome->fitness(imstr))  // If a new minimim was found
        {
            // If pointer to log file exists, write information about new minimum
            if (logfile_)
            {
                printNewMinimum(currEval, xiter);
            }
            *bestGenome = *genome;
            bestGenome->setFitness(imstr, currEval[imstr]);  // Pass the fitness for training set to the best genome
            sii_->saveState(false);
        }
        prevEval->find(imstr)->second = currEval[imstr];
        genome->setFitness(imstr, currEval[imstr]);
        if (evaluateTestSet_)
        {
            auto imste = iMolSelect::Test;
            prevEval->find(imste)->second = currEval[imste];
            bestGenome->setFitness(imste, currEval[imste]);  // Pass the fitness for the test set to the best genome
            genome->setFitness(imste, currEval[imste]);
        }
        acceptedMoves_[paramIndex] += 1;
    }
    else
    {  
        // If the parameter change is not accepted
        // Set the old value of the parameter back
        genome->setBase(paramIndex, storeParam);
        // poldata needs to change back as well!
        sii_->updatePoldata(*changed, genome);
    }
    // Set changed[j] back to false for upcoming iterations
    (*changed)[paramIndex] = false;

    printParameterStep(genome, xiter);
    printChi2Step(*prevEval, xiter);
}                  

void MCMCMutator::computeMeanSigma(const std::vector<double>    &sum,
                                   const int                     nsum,
                                         std::vector<double>    *sum_of_sq)
{
    if (nsum <= 0) return;
    double ps2 = 0.0;
    for (size_t k = 0; k < sii_->nParam(); k++)
    {
        pMean_[k]        = (sum[k]/nsum);
        (*sum_of_sq)[k] /= nsum;
        ps2              = std::max(0.0, (*sum_of_sq)[k]-gmx::square(pMean_[k]));
        pSigma_[k]       = std::sqrt(ps2);
    }
}                                         

void MCMCMutator::changeParam(ga::Genome *genome,
                              size_t      j)
{
    std::vector<double> *param = genome->basesPtr();
    
    GMX_RELEASE_ASSERT(j < genome->nBase(), "Parameter out of range");
    double rnd = randNum();
    while (rnd == 0.5)
    {
        // Make sure the parameter changes!
        rnd = randNum();
    }
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


void MCMCMutator::printMonteCarloStatistics(FILE             *fp,
                                            const ga::Genome &initialGenome,
                                            const ga::Genome &bestGenome)
{
    if (!fp)
    {
        return;
    }

    const auto paramNames = sii_->paramNames();
    const auto ntrain = sii_->nTrain();
    const auto weightedTemperature = sii_->weightedTemperature();

    fprintf(fp, "Monte Carlo statistics of parameters after optimization\n");
    fprintf(fp, "#best %zu #mean %zu #sigma %zu #param %zu\n",
            bestGenome.nBase(), pMean_.size(), 
            pSigma_.size(), paramNames.size());
    fprintf(fp, "Parameter                     Ncopies Initial   Best    Mean    Sigma Attempt  Acceptance  T-Weight\n");
    for (size_t k = 0; k < bestGenome.nBase(); k++)
    {
        double acceptance_ratio = 0;
        if (attemptedMoves_[k] > 0)
        {
            acceptance_ratio = 100*(double(acceptedMoves_[k])/attemptedMoves_[k]);
        }
        fprintf(fp, "%-30s  %5d  %6.3f  %6.3f  %6.3f  %6.3f    %4d %5.1f%%  %10.5f\n",
                paramNames[k].c_str(), ntrain[k],
                initialGenome.base(k), bestGenome.base(k), 
                pMean_[k], pSigma_[k],
                attemptedMoves_[k], acceptance_ratio, weightedTemperature[k]);
    }
}

void MCMCMutator::printNewMinimum(const std::map<iMolSelect, double> &chi2,
                                  double                              xiter)
{
    if (!logfile_)
    {
        return;
    }
    fprintf(logfile_, "Middleman %i\n", sii_->id());
    if (evaluateTestSet_)
    {
        fprintf(logfile_, "iter %10g. Found new minimum at %10g. Corresponding energy on the test set: %g\n",
                xiter, chi2.find(iMolSelect::Train)->second, chi2.find(iMolSelect::Test)->second);
    }
    else
    {
        fprintf(logfile_, "iter %10g. Found new minimum at %10g\n",
                xiter, chi2.find(iMolSelect::Train)->second);
    }
    // if (debug)
    // {
    //     genome->print(debug);
    // }
}             

void MCMCMutator::printParameterStep(ga::Genome *genome,
                                     double      xiter)
{
    const auto &paramClassIndex = sii_->paramClassIndex();

    auto bases = genome->bases();
    if (fpc_.size() == sii_->paramClass().size())
    {
        // Write iteration number to each parameter convergence 
        // surveillance file
        for (auto &fp: fpc_) 
        {
            if (fp.get())
            {
                fprintf(fp.get(), "%8f", xiter);
            }
        }
        // Write value of each parameter to its respective surveillance file
        for (size_t k = 0; k < genome->nBase(); k++)  
        {
            if (fpc_[paramClassIndex[k]].get())
            {
                fprintf(fpc_[paramClassIndex[k]].get(), "  %10g", bases[k]);
            }
        }
        for (auto &fp: fpc_)
        {
            if (fp.get())
            {
                fprintf(fp.get(), "\n");
                // If verbose, flush the file to be able to add 
                // new data to surveillance plots
                if (verbose_)
                {
                    fflush(fp.get());
                }
            }
        }
    }
}                             

void MCMCMutator::printChi2Step(const std::map<iMolSelect, double> &chi2,
                                double                              xiter)
{
    if (fpe_.get() == nullptr)
    {
        // If fpe is a null pointer, return
        return;
    }

    if (evaluateTestSet_)
    {
        fprintf(fpe_.get(), "%8f  %10g  %10g\n",
                xiter, chi2.find(iMolSelect::Train)->second,
                chi2.find(iMolSelect::Train)->second);
    }
    else
    {
        fprintf(fpe_.get(), "%8f  %10g\n",
                xiter, chi2.find(iMolSelect::Train)->second);
    }
    if (verbose_)
    {
        fflush(fpe_.get());
    }
}                    

void MCMCMutator::sensitivityAnalysis(ga::Genome *genome,
                                      iMolSelect  ims)
{
    
    std::vector<double> *param = genome->basesPtr();
    const auto upperBound      = sii_->upperBound();
    const auto lowerBound      = sii_->lowerBound();
    const auto paramNames      = sii_->paramNames();

    if (param->size() == 0)
    {
        return;
    }
    sii_->updatePoldata(genome);
    // std::fill(changed.begin(), changed.end(), false);
    std::vector<bool> changed(param->size(), false);
    // changed.resize(param->size(), false);
    double chi2_0 = fitComp_->calcDeviation(param, CalcDev::Parallel, ims);
    if (logfile_)
    {
        fprintf(logfile_, "\nStarting sensitivity analysis. chi2_0 = %g nParam = %d\n",
                chi2_0, static_cast<int>(param->size()));
        fflush(logfile_);
    }
    for (size_t i = 0; i < param->size(); ++i)
    {
        Sensitivity s;
        double pstore = (*param)[i];
        double deltap = (upperBound[i]-lowerBound[i])/200;
        double pmin   = std::max((*param)[i]-deltap, lowerBound[i]);
        double pmax   = std::min((*param)[i]+deltap, upperBound[i]);
        double p_0    = 0.5*(pmin+pmax);
        changed[i]    = true;
        (*param)[i]     = pmin;
        sii_->updatePoldata(changed, genome);
        s.add((*param)[i], fitComp_->calcDeviation(param, CalcDev::Parallel, ims));
        (*param)[i]     = p_0;
        sii_->updatePoldata(changed, genome);
        s.add((*param)[i], fitComp_->calcDeviation(param, CalcDev::Parallel, ims));
        (*param)[i]     = pmax;
        sii_->updatePoldata(changed, genome);
        s.add((*param)[i],  fitComp_->calcDeviation(param, CalcDev::Parallel, ims));
        (*param)[i]     = pstore;
        sii_->updatePoldata(changed, genome);
        changed[i]    = false;
        s.computeForceConstants(logfile_);
        s.print(logfile_, paramNames[i]);
    }
    if (logfile_)
    {
        fprintf(logfile_, "Sensitivity analysis done.\n");
    }

}                                      

void MCMCMutator::openParamConvFiles(const gmx_output_env_t *oenv)
{
    const std::vector<std::string> &pClass = sii_->paramClass();
    for (size_t i = 0; i < pClass.size(); i++)
    {
        auto fileName = gmx::formatString("%sind%d-%s-%s",
                                          sii_->prefix().c_str(), sii_->id(),
                                          pClass[i].c_str(),
                                          sii_->xvgConv().c_str());
        gmx::FilePtr fp(xvgropen(fileName.c_str(), "Parameter convergence",
                                 "iteration", "", oenv));
        fpc_.push_back(std::move(fp));

        std::vector<const char*> tmpParamNames;
        for (size_t j = 0; j < sii_->paramNames().size(); j++)
        {
            if ( (sii_->paramClassIndex())[j] == i )
            {
                tmpParamNames.push_back( (sii_->paramNames())[j].c_str() );
            }
        }
        xvgr_legend(fpc_[i].get(), tmpParamNames.size(),
                    tmpParamNames.data(), oenv);
    }
}

void MCMCMutator::openChi2ConvFile(const gmx_output_env_t *oenv)
{
    auto fileName = gmx::formatString("%sind%d-%s", sii_->prefix().c_str(),
                                      sii_->id(), sii_->xvgEpot().c_str());
    fpe_.reset(xvgropen(fileName.c_str(),
                        "Chi squared",
                        "Iteration",
                        "Unknown units",
                        oenv));
    if (evaluateTestSet_)
    {
        std::vector<std::string> legend;
        legend.push_back(iMolSelectName(iMolSelect::Train));
        legend.push_back(iMolSelectName(iMolSelect::Test));
        xvgrLegend(fpe_.get(), legend, oenv);
    }
}

} //namespace alexandria

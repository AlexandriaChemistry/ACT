/*
 * This source file is part of the Alexandria Chemistry Toolkit.
 *
 * Copyright (C) 2014-2021
 *
 * Developers:
 *             Mohammad Mehdi Ghahremanpour,
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

/*! \internal \brief
 * Implements part of the alexandria program.
 * \author Mohammad Mehdi Ghahremanpour <mohammad.ghahremanpour@icm.uu.se>
 * \author David van der Spoel <david.vanderspoel@icm.uu.se>
 */

#include "optparam.h"

#include <assert.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include <functional>
#include <string>
#include <vector>

#include "gromacs/random.h"
#include "gromacs/commandline/pargs.h"
#include "gromacs/fileio/confio.h"
#include "gromacs/fileio/xvgr.h"
#include "gromacs/math/units.h"
#include "gromacs/math/vec.h"
#include "gromacs/utility/arraysize.h"
#include "gromacs/utility/fatalerror.h"
#include "gromacs/utility/smalloc.h"

#include "memory_check.h"
#include "regression.h"


namespace alexandria
{

void OptParam::add_pargs(std::vector<t_pargs> *pargs)
{
    t_pargs pa[] = {
        { "-maxiter", FALSE, etINT, {&maxiter_},
          "Max number of iterations for optimization" },
        { "-temp",    FALSE, etREAL, {&temperature_},
          "'Temperature' for the Monte Carlo simulation" },
        { "-tweight", FALSE, etBOOL, {&tempWeight_},
          "Weight the temperature in the MC/MC algorithm according to the square root of the number of data points. This is in order to get a lower probability ofaccepting a step in the wrong direction for parameters of which there are few copies." },
        { "-anneal", FALSE, etREAL, {&anneal_},
          "Use annealing in Monte Carlo simulation, starting from this fraction of the simulation. Value should be between 0 and 1." },
        { "-seed",   FALSE, etINT,  {&seed_},
          "Random number seed. If zero, a seed will be generated." },
        { "-step",  FALSE, etREAL, {&step_},
          "Step size for the parameter optimization. Is used as fraction of the available range per parameter which depends on the parameter type." },
        { "-v",     FALSE, etBOOL, {&verbose_},
          "Flush output immediately rather than letting the OS buffer it. Don't use for production simulations." }
    };
    for (int i = 0; i < asize(pa); i++)
    {
        pargs->push_back(pa[i]);
    }
}

void OptParam::setOutputFiles(const char                     *xvgconv,
                              const std::vector<std::string> &paramClass,
                              const char                     *xvgepot,
                              const gmx_output_env_t         *oenv)
{
    xvgconv_.assign(xvgconv);
    paramClass_ = paramClass;
    xvgepot_.assign(xvgepot);
    oenv_       = oenv;
}

double OptParam::computeBeta(int iter)
{
    double temp = temperature_;
    if (iter >= maxiter_)
    {
        temp = 1e-6;
    }
    else
    {
        temp = temperature_*(1.0 - iter/(maxiter_ + 1.0));
    }
    return 1/(BOLTZ*temp);
}

double OptParam::computeBeta(int maxiter, int iter, int ncycle)
{
    double temp = temperature_;
    if (iter >= maxiter_)
    {
        temp = 1e-6;
    }
    else
    {
        temp = (0.5*temperature_)*((exp(-iter/(0.2*(maxiter+1)))) * (1.1 + cos((ncycle*M_PI*iter)/(maxiter+1))));
    }
    return 1/(BOLTZ*temp);
}

void Sensitivity::computeForceConstants(FILE *fp)
{
    if (p_.size() >= 3)
    {
        MatrixWrapper M(3, p_.size());
        std::vector<double> solution;
        solution.resize(3, 0.0);
        for (size_t i = 0; i < p_.size(); ++i)
        {
            M.set(0, i, p_[i]*p_[i]);
            M.set(1, i, p_[i]);
            M.set(2, i, 1.0);
        }
        auto result = M.solve(chi2_, &solution);
        if (result == 0)
        {
            a_ = solution[0];
            b_ = solution[1];
            c_ = solution[2];
        }
    }
    else if (fp)
    {
        fprintf(fp, "Not enough parameters %d to do sensitivty analysis\n",
                static_cast<int>(p_.size()));
    }
}

void Sensitivity::print(FILE *fp, const std::string &label)
{
    if (fp)
    {
        fprintf(fp, "Sensitivity %s Fit to parabola: a %10g b %10g c %10g\n",
                label.c_str(), a_, b_, c_);
        for(int i = 0; i < static_cast<int>(p_.size()); ++i)
        {
            fprintf(fp, "    p[%d] %g chi2[%d] %g\n", i, p_[i], i, chi2_[i]);
        }
        if (a_ != 0.0)
        {
            double p_min = -b_/(2.0*a_);
            double chi2_min = a_*p_min*p_min + b_*p_min + c_;
            fprintf(fp, "    pmin %g chi2min %g (estimate based on parabola)\n", 
                    p_min, chi2_min);
        }
    }
}

void Bayes::printParameters(FILE *fp) const
{
    if (nullptr == fp)
    {
        return;
    }
    for(size_t i = 0; i < param_.size(); i++)
    {
        fprintf(fp, "  %s  %e,", paramNames_[i].c_str(), param_[i]);
    }
    fprintf(fp, "\n");
}

void Bayes::addParam(const std::string &name,
                     real               val,
                     Mutability         mut,
                     real               lower,
                     real               upper,
                     int                ntrain,
                     bool               bRandom)
{
    if (bRandom)
    {
        std::random_device               rd;
        std::mt19937                     gen(rd());  
        std::uniform_real_distribution<> uniform(lower, upper);
        val                              = uniform(gen);        
    }
    
    initial_param_.push_back(val);
    param_.push_back(val);
    mutability_.push_back(mut);
    ntrain_.push_back(ntrain);
    //    prevParam_.push_back(val);
    lowerBound_.push_back(lower);
    upperBound_.push_back(upper);
    paramNames_.push_back(name);
}

void Bayes::changeParam(size_t j, real rand)
{
    GMX_RELEASE_ASSERT(j < param_.size(), "Parameter out of range");
    real delta = (2*rand-1)*step()*(upperBound_[j]-lowerBound_[j]);
    param_[j] += delta;
    if (mutability_[j] == Mutability::Bounded)
    {
        if (param_[j] < lowerBound_[j])
        {
            param_[j] = lowerBound_[j];
        }
        else if (param_[j] > upperBound_[j])
        {
            param_[j] = upperBound_[j];
        }
    }
}

bool OptParam::anneal(int iter) const
{
    if (anneal_ >= 1)
    {
        return false;
    }
    else if (anneal_ <= 0)
    {
        return true;
    }
    else
    {
        return iter >= anneal_ * maxiter_;
    }   
}

void Bayes::SensitivityAnalysis(FILE *fplog, iMolSelect ims)
{
    if (param_.size() == 0)
    {
        return;
    }
    std::vector<bool> changed;
    changed.resize(param_.size(), true);
    toPoldata(changed);
    std::fill(changed.begin(), changed.end(), false);
    double chi2_0 = calcDeviation(false, CalcDev::Parallel, ims);
    if (fplog)
    {
        fprintf(fplog, "\nStarting sensitivity analysis. chi2_0 = %g nParam = %d\n",
                chi2_0, static_cast<int>(param_.size()));
        fflush(fplog);
    }
    for (size_t i = 0; i < param_.size(); ++i)
    {
        Sensitivity s;
        double pstore = param_[i];
        double deltap = (upperBound_[i]-lowerBound_[i])/200;
        double pmin   = std::max(param_[i]-deltap, lowerBound_[i]);
        double pmax   = std::min(param_[i]+deltap, upperBound_[i]);
        double p_0    = 0.5*(pmin+pmax);
        changed[i]    = true;
        param_[i]     = pmin;
        toPoldata(changed);
        s.add(param_[i], calcDeviation(false, CalcDev::Parallel, ims));
        param_[i]     = p_0;
        toPoldata(changed);
        s.add(param_[i], calcDeviation(false, CalcDev::Parallel, ims));
        param_[i]     = pmax;
        toPoldata(changed);
        s.add(param_[i],  calcDeviation(false, CalcDev::Parallel, ims));
        param_[i]     = pstore;
        toPoldata(changed);
        changed[i]    = false;
        s.computeForceConstants(fplog);
        s.print(fplog, paramNames_[i]);
    }
    if (fplog)
    {
        fprintf(fplog, "Sensitivity analysis done.\n");
    }
}

bool Bayes::MCMC(FILE *fplog, bool bEvaluate_testset, double *chi2)
{
    int                              nsum             = 0;
    int                              nParam           = 0;
    double                           minEval          = 0;
    double                           prevEval         = 0;
    double                           prevEval_testset = 0;
    parm_t                           sum, sum_of_sq;
    //! Pointers to parameter convergence surveillance files
    std::vector<FILE *>              fpc;
    std::vector<int>                 paramClassIndex;
    //! Pointer to chi2 surveillance file
    FILE                            *fpe             = nullptr;
    
    if (xvgConv().empty() || xvgEpot().empty())
    {
        gmx_fatal(FARGS, "You forgot to call setOutputFiles. Back to the drawing board.");
    }
    if (param_.empty())
    {
        fprintf(stderr, "No parameters to optimize.\n");
        return 0;
    }

    // Allocate memory for parameter class index.
    // Set to -1 to indicate not set, and to crash the program
    // in case of bugs.
    paramClassIndex.resize(paramNames_.size(), -1);
    std::vector<std::string> pClass = paramClass();
    assignParamClasses(paramClassIndex, pClass);

    openParamSurveillanceFiles(pClass, fpc, paramClassIndex);

    // Compute temperature weights if relevant, otherwise the numbers are all 1.0
    weightedTemperature_.resize(paramNames_.size(), 1.0);
    if (temperatureWeighting()) computeWeightedTemperature();
    
    fpe = openChi2SurveillanceFile(bEvaluate_testset);

    // Initialize data structures
    nParam = param_.size();
    sum.resize(nParam, 0);
    sum_of_sq.resize(nParam, 0);
    pmean_.resize(nParam, 0);
    psigma_.resize(nParam, 0);
    attemptedMoves_.resize(nParam, 0);
    acceptedMoves_.resize(nParam, 0);
    std::vector<bool> changed;
    // Initialize to true to make sure the parameters are 
    // all spread around the processors.
    changed.resize(nParam, true);
    toPoldata(changed);
    // Now set them back to false, further spreading is done
    // one parameter at a time.
    std::fill(changed.begin(), changed.end(), false);

    // training set
    prevEval = calcDeviation(true, CalcDev::Parallel, iMolSelect::Train);
    minEval  = prevEval;
    (*chi2)  = prevEval;

    if (bEvaluate_testset)
    {
        // test set
        prevEval_testset = calcDeviation(true, CalcDev::Parallel, iMolSelect::Test);
    }

    // Random number 
    std::random_device               rd;
    std::mt19937                     gen(rd());
    std::uniform_int_distribution<>  int_uniform(0, nParam-1);
    std::uniform_real_distribution<> real_uniform(0, 1);
    
    print_memory_usage(debug);

    double beta0 = 1/(BOLTZ*temperature());

    // Optimization loop
    for (int iter = 0; iter < maxIter(); iter++)
    { 
        for (int pp = 0; pp < nParam; pp++)
        {      
            // Pick a random parameter to change
            const int paramIndex = int_uniform(gen);

            // Do the step!
            stepMCMC(paramIndex, gen, real_uniform, changed, &prevEval, &prevEval_testset,
                     bEvaluate_testset, pp, iter, &beta0, nParam, &minEval, fplog, fpc, fpe,
                     paramClassIndex);

            // For the second half of the optimization, collect data to find the mean and standard deviation of each
            // parameter
            if (iter >= maxIter()/2)
            {
                for (auto k = 0; k < nParam; k++)
                {
                    sum[k]       += param_[k];
                    sum_of_sq[k] += gmx::square(param_[k]);
                }
                nsum++;
            }
        }
    }

    // OPTIMIZATION IS COMPLETE!

    computeMeanSigma(nParam, sum, nsum, sum_of_sq);

    closeConvergenceFiles(fpc, fpe);

    bool bMinimum = false;  // Assume no better minimum was found
    if (minEval < (*chi2))  // If better minimum was found, update the value in <*chi2> and return true
    {
        (*chi2) = minEval;
        bMinimum = true;
    }
    return bMinimum;
}

void Bayes::computeWeightedTemperature()
{
    for(size_t j = 0; j < paramNames_.size(); j++)
    {
        GMX_RELEASE_ASSERT(ntrain_[j] > 0, "ntrain should be > 0 for all parameters");
        // TODO: Maybe a fast inverse square root here?
        weightedTemperature_[j] = std::sqrt(1.0/ntrain_[j]);
    }
}

void Bayes::stepMCMC(const int                                  paramIndex,
                           std::mt19937&                        gen,
                           std::uniform_real_distribution<>&    real_uniform,
                           std::vector<bool>&                   changed,
                           double*                              prevEval,
                           double*                              prevEval_testset,
                     const bool                                 bEvaluate_testset,
                     const int                                  pp,
                     const int                                  iter,
                           double*                              beta0,
                     const int                                  nParam,
                           double*                              minEval,
                           FILE*                                fplog,
                           std::vector<FILE*>&                  fpc,
                           FILE*                                fpe,
                           std::vector<int>&                    paramClassIndex)
{

    // Store the original value of the parameter
    const double storeParam = param_[paramIndex];

    // Change the parameter
    double rndNumber = real_uniform(gen);
    while (rndNumber == 0.5) rndNumber = real_uniform(gen);
    changeParam(paramIndex, rndNumber);

    attemptedMoves_[paramIndex] += 1;
    changed[paramIndex]          = true;

    // Update FF parameter data structure with
    // the new value of parameter j
    toPoldata(changed);

    // Evaluate the energy on training set
    const double currEval = calcDeviation(false, CalcDev::Parallel, iMolSelect::Train);
    const double deltaEval = currEval - (*prevEval);

    // Evaluate the energy on the test set only on whole steps!
    double currEval_testset = (*prevEval_testset);
    if (bEvaluate_testset && pp == 0)
    {
        currEval_testset = calcDeviation(false, CalcDev::Parallel, iMolSelect::Test);
    }

    // Accept any downhill move
    bool accept = (deltaEval < 0);

    // For an uphill move apply the Metropolis Criteria
    // to decide whether to accept or reject the new parameter
    if (!accept)
    {
        // Only anneal if the simulation reached a certain number of steps
        if (anneal(iter)) (*beta0) = computeBeta(iter);
        const double randProbability = real_uniform(gen);
        const double mcProbability   = exp(-((*beta0)/weightedTemperature_[paramIndex])*deltaEval);
        accept = (mcProbability > randProbability);
    }

    // Fractional iteration taking into account the inner loop with <pp> over <nParam>
    const double xiter = iter + (1.0*pp)/nParam;
    if (accept)
    {  // If the parameter change is accepted
        if (currEval < (*minEval))
        {
            // If pointer to log file exists, write information about new minimum
            if (fplog) fprintNewMinimum(fplog, bEvaluate_testset, xiter, currEval, currEval_testset);
            bestParam_ = param_;
            (*minEval) = currEval;
            saveState();
        }
        (*prevEval) = currEval;
        if (bEvaluate_testset) (*prevEval_testset) = currEval_testset;
        acceptedMoves_[paramIndex] += 1;
    }
    else
    {  // If the parameter change is not accepted
        param_[paramIndex] = storeParam;  // Set the old value of the parameter back
        // poldata needs to change back as well!
        toPoldata(changed);
    }
    changed[paramIndex] = false;  // Set changed[j] back to false for upcoming iterations

    fprintParameterStep(fpc, paramClassIndex, xiter);
    // If the chi2 surveillance file exists, write progress
    if (fpe != nullptr) fprintChi2Step(bEvaluate_testset, fpe, xiter, *prevEval, *prevEval_testset);

}

void Bayes::assignParamClasses(std::vector<int>&         paramClassIndex,
                               std::vector<std::string>& pClass)
{

    for(size_t i = 0; i < pClass.size(); i++)
    {
        for (size_t j = 0; j < paramNames_.size(); j++)
        {
            if (paramNames_[j].find(pClass[i]) != std::string::npos)
            {
                paramClassIndex[j] = i;
            }
        }
    }

    // Now check for params which were not assigned a class and give them class "Other"
    bool restClass = false;
    for(size_t i = 0; i < paramClassIndex.size(); i++)
    {
        if (paramClassIndex[i] == -1)
        {
            if (!restClass)  // If <i> is the first parameter without a class
            {
                // Append "Other" to the list of classes
                pClass.push_back("Other");
                restClass = true;
            }
            // Give class "Other" to parameter <i>
            paramClassIndex[i] = pClass.size()-1;
        }
    }

}

void Bayes::openParamSurveillanceFiles(const std::vector<std::string>&  pClass,
                                             std::vector<FILE*>&        fpc,
                                             std::vector<int>&          paramClassIndex)
{
    for(size_t i = 0; i < pClass.size(); i++)
    {
        std::string fileName = pClass[i] + "-" + xvgConv();
        fpc.push_back(xvgropen(fileName.c_str(),
                               "Parameter convergence",
                               "iteration",
                               "",
                               oenv()));
        std::vector<const char*> paramNames;
        for (size_t j = 0; j < paramNames_.size(); j++)
        {
            if (paramClassIndex[j] == static_cast<int>(i))
            {
                paramNames.push_back(paramNames_[j].c_str());
            }
        }
        xvgr_legend(fpc[i], paramNames.size(), paramNames.data(), oenv());
    }
}

FILE* Bayes::openChi2SurveillanceFile(const bool bEvaluate_testset)
{
    FILE* fpe = xvgropen(xvgEpot().c_str(),
                         "Chi squared",
                         "Iteration",
                         "Unknown units",
                         oenv());
    if (bEvaluate_testset)
    {
        std::vector<std::string> legend;
        legend.push_back(iMolSelectName(iMolSelect::Train));
        legend.push_back(iMolSelectName(iMolSelect::Test));
        xvgrLegend(fpe, legend, oenv());
    }
    return fpe;
}

void Bayes::computeMeanSigma(const int     nParam,
                             const parm_t& sum,
                             const int     nsum,
                                   parm_t& sum_of_sq)
{

    if (nsum > 0)  // Compute mean and standard deviation
    {
        double ps2 = 0.0;
        for (auto k = 0; k < nParam; k++)
        {
            pmean_[k]     = (sum[k]/nsum);
            sum_of_sq[k] /= nsum;
            ps2           = std::max(0.0, sum_of_sq[k]-gmx::square(pmean_[k]));
            psigma_[k]    = sqrt(ps2);
        }
    }

}

void Bayes::closeConvergenceFiles(std::vector<FILE*>& fpc,
                                  FILE*               fpe)
{
    for(auto fp: fpc)  // Close all parameter convergence surveillance files
    {
        xvgrclose(fp);
    }
    if (nullptr != fpe)  // Close chi2 surveillance file
    {
        xvgrclose(fpe);
    }
}

void Bayes::fprintNewMinimum(      FILE*   fplog,
                             const bool    bEvaluate_testset,
                             const double  xiter,
                             const double  currEval,
                             const double  currEval_testset)
{

    if (bEvaluate_testset)
    {
        fprintf(fplog, "iter %10g. Found new minimum at %10g. Corresponding energy on the test set: %g\n",
                xiter, currEval, currEval_testset);
    }
    else
    {
        fprintf(fplog, "iter %10g. Found new minimum at %10g\n",
                xiter, currEval);
    }
    if (debug)
    {
        printParameters(debug);
    }

}

void Bayes::fprintParameterStep(      std::vector<FILE*>&   fpc,
                                const std::vector<int>&     paramClassIndex,
                                const double                xiter)
{

    for(auto fp: fpc)  // Write iteration number to each parameter convergence surveillance file
    {
        fprintf(fp, "%8f", xiter);
    }
    for (size_t k = 0; k < param_.size(); k++)  // Write value of each parameter to its respective surveillance file
    {
        fprintf(fpc[paramClassIndex[k]], "  %10g", param_[k]);
    }
    for(auto fp: fpc)  // If verbose = True, flush the file to be able to add new data to surveillance plots
    {
        fprintf(fp, "\n");
        if (verbose())
        {
            fflush(fp);
        }
    }

}

void Bayes::fprintChi2Step(const bool   bEvaluate_testset,
                                 FILE*  fpe,
                           const double xiter,
                           const double prevEval,
                           const double prevEval_testset)
{

    if (bEvaluate_testset)
    {
        fprintf(fpe, "%8f  %10g  %10g\n", xiter, prevEval, prevEval_testset);
    }
    else
    {
        fprintf(fpe, "%8f  %10g\n", xiter, prevEval);
    }
    if (verbose())
    {
        fflush(fpe);
    }

}

void Bayes::printMonteCarloStatistics(FILE *fp)
{
    if (!fp)
    {
        return;
    }
    fprintf(fp, "\nMonte Carlo statistics of parameters after optimization\n");
    fprintf(fp, "#best %zu #mean %zu #sigma %zu #param %zu\n",
            bestParam_.size(), pmean_.size(), psigma_.size(), paramNames_.size());
    if (bestParam_.size() == nParam())
    {
        fprintf(fp, "Parameter                     Ncopies Initial   Best    Mean    Sigma Attempt  Acceptance  T-Weight\n");
        for (size_t k = 0; k < Bayes::nParam(); k++)
        {
            double acceptance_ratio = 0;
            if (attemptedMoves_[k] > 0)
            {
                acceptance_ratio = 100*(double(acceptedMoves_[k])/attemptedMoves_[k]);
            }
            fprintf(fp, "%-30s  %5d  %6.3f  %6.3f  %6.3f  %6.3f    %4d %5.1f%%  %10.5f\n",
                    paramNames_[k].c_str(), ntrain_[k],
                    initial_param_[k], bestParam_[k], pmean_[k], psigma_[k],
                    attemptedMoves_[k], acceptance_ratio, weightedTemperature_[k]);
        }
    }
}


}

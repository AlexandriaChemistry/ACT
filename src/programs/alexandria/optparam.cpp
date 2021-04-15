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

namespace alexandria
{

void OptParam::add_pargs(std::vector<t_pargs> *pargs)
{
    t_pargs pa[] = {
        { "-maxiter", FALSE, etINT, {&maxiter_},
          "Max number of iterations for optimization" },
        { "-mcr", FALSE, etINT, {&mc_ratio_},
          "Reference acceptance ratio used for running Adaptive-MCMC." },
        { "-nAdapt", FALSE, etINT, {&nAdapt_},
          "Adapt Temperature every nAdapt steps when running Adaptive-MCMC." },
        { "-temp",    FALSE, etREAL, {&temperature_},
          "'Temperature' for the Monte Carlo simulation" },
        { "-anneal", FALSE, etREAL, {&anneal_},
          "Use annealing in Monte Carlo simulation, starting from this fraction of the simulation. Value should be between 0 and 1." },
        { "-adaptive", FALSE, etBOOL, {&adaptive_},
          "Use Adaptive Monte Carlo simulation." },
        { "-seed",   FALSE, etINT,  {&seed_},
          "Random number seed. If zero, a seed will be generated." },
        { "-step",  FALSE, etREAL, {&step_},
          "Step size for the parameter optimization. Is used as fraction of the available range per parameter which depends on the parameter type" }
    };
    for (size_t i = 0; i < asize(pa); i++)
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


double OptParam::adaptBeta(real mc_ratio)
{
    temperature_ *= (mc_ratio/mc_ratio_);
    return 1/(BOLTZ*temperature_);
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
                     real               lower,
                     real               upper,
                     bool               bRandom)
{
    if (bRandom)
    {
        std::random_device               rd;
        std::mt19937                     gen(rd());  
        real delta                       = std::abs((0.5*val));      
        real a                           = val - delta;
        real b                           = val + delta;       
        std::uniform_real_distribution<> uniform(a, b);
        val                              = uniform(gen);        
    }
    
    initial_param_.push_back(val);
    param_.push_back(val);
    //    prevParam_.push_back(val);
    lowerBound_.push_back(lower);
    upperBound_.push_back(upper);
    paramNames_.push_back(name);
}

void Bayes::changeParam(size_t j, real rand)
{
    GMX_RELEASE_ASSERT(j < param_.size(), "Parameter out of range");
    //real delta = (2*rand-1)*step()*fabs(param_[j]);
    real delta = (2*rand-1)*step()*(upperBound_[j]-lowerBound_[j]);
    param_[j] += delta;
    if (boxConstraint())
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

double Bayes::MCMC(FILE *fplog)
{
    double                           storeParam;
    int                              nsum            = 0;
    int                              nParam          = 0; 
    double                           currEval        = 0;
    double                           minEval         = 0;
    double                           prevEval        = 0;
    double                           deltaEval       = 0;
    double                           randProbability = 0;
    double                           mcProbability   = 0; 
    parm_t                           sum, sum_of_sq;
    
    std::vector<FILE *>              fpc;
    std::vector<int>                 paramClassIndex;
    FILE                            *fpe             = nullptr;
    
    if (xvgConv().empty() || xvgEpot().empty())
    {
        gmx_fatal(FARGS, "You forgot to call setOutputFiles. Back to the drawing board.");
    }
    if (param_.empty())
    {
        fprintf(stderr, "No parameters to optimze.\n");
        return 0;
    }
    // Allocate memory for parameter class index.
    // Set to -1 to indicate not set, and to crash the program
    // in case of bugs.
    paramClassIndex.resize(paramNames_.size(), -1);
    std::vector<std::string> pClass = paramClass();
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
    // Now check for "unclassified parameters"
    bool restClass = false;
    for(size_t i = 0; i < paramClassIndex.size(); i++)
    {
        if (paramClassIndex[i] == -1)
        {
            if (!restClass)
            {
                pClass.push_back("Other");
                restClass = true;
            }
            paramClassIndex[i] = pClass.size()-1;
        }
    }
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
    
    // Now parameter output file.
    fpe = xvgropen(xvgEpot().c_str(), 
                   "Parameter energy", 
                   "iteration",
                   "\\f{12}c\\S2\\f{4}", 
                   oenv());

    nParam = param_.size();
    sum.resize(nParam, 0);
    sum_of_sq.resize(nParam, 0);
    pmean_.resize(nParam, 0);
    psigma_.resize(nParam, 0);
    attemptedMoves_.resize(nParam, 0);
    acceptedMoves_.resize(nParam, 0);
    std::vector<bool> changed;
    changed.resize(nParam, false);
    toPoldata(changed);
    prevEval = calcDeviation();
    minEval  = prevEval;
    if (debug)
    {
        fprintf(debug, "Initial chi2 value = %g\n", prevEval);
    }

    // Random number 
    std::random_device               rd;
    std::mt19937                     gen(rd());
    std::uniform_int_distribution<>  int_uniform(0, nParam-1);
    std::uniform_real_distribution<> real_uniform(0, 1);
    
    print_memory_usage(fplog);
    // Optmization loop
    int    j                = 0;
    bool   accept           = false;
    double xiter            = 0.0;
    double beta             = 1/(BOLTZ*temperature());
    int    total_iterations = nParam*maxIter();
    
    for (int iter = 0; iter < total_iterations; iter++)
    {       
        // Pick a random parameter to change
        j                  = int_uniform(gen);
        //prevParam_         = param_;
        storeParam         = param_[j];
        attemptedMoves_[j] = attemptedMoves_[j] + 1;
        
        // Change the picked parameter
        changeParam(j, real_uniform(gen));
        changed[j]         = true;
        
        // Evaluate the energy
        toPoldata(changed);
        currEval        = calcDeviation();
        deltaEval       = currEval-prevEval; 
        
        // Accept any downhill move       
        accept          = (deltaEval < 0);
        
        // For a uphill move apply the Metropolis Criteria
        // to decide whether to accept or reject the new parameter
        if (!accept)
        {
            // Only anneal if the simulation reached a certain number of steps
            if (anneal(iter))
            {
                beta = computeBeta(iter/nParam);
            }
            randProbability = real_uniform(gen);
            mcProbability   = exp(-beta*deltaEval);
            accept          = (mcProbability > randProbability);
        }
        
        xiter = (1.0*iter)/nParam;
        if (accept)
        {
            if (currEval < minEval)
            {
                if (fplog)
                {
                    fprintf(fplog, "iter %g. Found new minimum at %g\n",
                            xiter, currEval);
                }
                bestParam_ = param_;
                minEval    = currEval;
                if (debug && false)
                {
                    fprintf(debug, "iter %g. Found new minimum at %g\n",
                            xiter, currEval);
                    printParameters(debug);
                }
                saveState();
            }
            prevEval = currEval;
            acceptedMoves_[j] = acceptedMoves_[j] + 1;
        }
        else
        {
            param_[j] = storeParam;
            // poldata needs to change back!
            toPoldata(changed);
        }
        changed[j] = false;

        for(auto fp: fpc)
        {
            fprintf(fp, "%8f", xiter);
        }
        for (size_t k = 0; k < param_.size(); k++)
        {
            fprintf(fpc[paramClassIndex[k]], "  %10g", param_[k]);
        }
        for(auto fp: fpc)
        {
            fprintf(fp, "\n");
            fflush(fp);
        }
        if (nullptr != fpe)
        {
            fprintf(fpe, "%8f  %10g\n", xiter, prevEval);
            fflush(fpe);
        }
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
    if (nsum > 0)
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
    for(auto fp: fpc)
    {
        xvgrclose(fp);
    }
    if (nullptr != fpe)
    {
        xvgrclose(fpe);
    }
    return minEval;
}

double Bayes::Adaptive_MCMC(FILE *fplog)
{
    double                           storeParam;
    int                              nsum            = 0;
    int                              nParam          = 0; 
    double                           currEval        = 0;
    double                           minEval         = 0;
    double                           prevEval        = 0;
    double                           deltaEval       = 0;
    double                           randProbability = 0;
    double                           mcProbability   = 0; 
    double                           halfIter        = maxIter()/2;   
    parm_t                           sum, sum_of_sq;
    
    std::vector<FILE *>              fpc;
    std::vector<int>                 paramClassIndex;
    FILE                            *fpe             = nullptr;
    
    if (xvgConv().empty() || xvgEpot().empty())
    {
        gmx_fatal(FARGS, "You forgot to call setOutputFiles. Back to the drawing board.");
    }
    if (paramNames_.empty())
    {
        gmx_fatal(FARGS, "You forgot to add parameterNames. Back to the drawing board.");
    }
    // Allocate memory for parameter class index.
    // Set to -1 to indicate not set, and to crash the program
    // in case of bugs.
    paramClassIndex.resize(paramNames_.size(), -1);
    std::vector<std::string> pClass = paramClass();
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
    // Now check for "unclassified parameters"
    bool restClass = false;
    for(size_t i = 0; i < paramClassIndex.size(); i++)
    {
        if (paramClassIndex[i] == -1)
        {
            if (!restClass)
            {
                pClass.push_back("Other");
                restClass = true;
            }
            paramClassIndex[i] = pClass.size()-1;
        }
    }
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
    
    // Now parameter output file.
    fpe = xvgropen(xvgEpot().c_str(), 
                   "Parameter energy", 
                   "iteration",
                   "\\f{12}c\\S2\\f{4}", 
                   oenv());

    nParam = param_.size();
    sum.resize(nParam, 0);
    sum_of_sq.resize(nParam, 0);
    pmean_.resize(nParam, 0);
    psigma_.resize(nParam, 0);
    attemptedMoves_.resize(nParam, 0);
    acceptedMoves_.resize(nParam, 0);
    std::vector<bool> changed;
    changed.resize(nParam, false);
    toPoldata(changed);
    prevEval = calcDeviation();
    minEval  = prevEval;
    if (debug)
    {
        fprintf(debug, "Initial chi2 value = %g\n", prevEval);
    }
    if (fplog)
    {
        fprintf(fplog, "minEval %g, nParam %d\n", minEval, nParam);
    }

    // Randrom number 
    std::random_device               rd;
    std::mt19937                     gen(rd());
    std::uniform_int_distribution<>  int_uniform(0, nParam-1);
    std::uniform_real_distribution<> real_uniform(0, 1);
    
    // Optmization loop
    int    j                = 0;
    int    nMCaccept        = 0;
    bool   accept           = false;
    double xiter            = 0.0;
    double beta             = 1/(BOLTZ*temperature());
    int    total_iterations = nParam*maxIter();
    
    for (int iter = 0; iter < total_iterations; iter++)
    {       
        // Pick a random parameter to change
        j                  = int_uniform(gen);        
        //        prevParam_         = param_;
        storeParam         = param_[j];
        attemptedMoves_[j] = attemptedMoves_[j] + 1;
        
        // Change the picked parameter
        changeParam(j, real_uniform(gen));
        changed[j]         = true;
        
        // Evaluate the energy
        toPoldata(changed);
        currEval        = calcDeviation();
        deltaEval       = currEval-prevEval; 
        
        // Accept any downhill move       
        accept          = (deltaEval < 0);
        
        // For a uphill move apply the Metropolis Criteria
        // to decide whether to accept or reject the new parameter
        if (!accept)
        {
            if (iter > 0 && (iter % nAdapt()) == 0)
            {
                // Adapt temperature to keep MC acceptance ratio to ref_ratio.                
                real mc_ratio = (100.0 * nMCaccept)/iter;
                if (fplog)
                {
                    fprintf(fplog, "MC-Ratio %0.2f.\n", mc_ratio);
                }
                beta       = adaptBeta(mc_ratio);
            }
            randProbability = real_uniform(gen);
            mcProbability   = exp(-beta*deltaEval);
            if (mcProbability > randProbability)
            {
                accept = true;
                nMCaccept++; 
            }          
        }
        
        xiter = (1.0*iter)/nParam;
        if (accept)
        {
            if (currEval < minEval)
            {
                if (fplog)
                {
                    fprintf(fplog, "iter %g. Found new minimum at %g\n",
                            xiter, currEval);
                }
                bestParam_ = param_;
                minEval    = currEval;
                if (false)
                {
                    printParameters(fplog);
                }
            }
            prevEval = currEval;
            acceptedMoves_[j] = acceptedMoves_[j] + 1;
        }
        else
        {
            param_[j] = storeParam;
            // poldata needs to change back!
            toPoldata(changed);
        }
        changed[j] = false;
        
        for(auto fp: fpc)
        {
            fprintf(fp, "%8f", xiter);
        }
        for (size_t k = 0; k < param_.size(); k++)
        {
            fprintf(fpc[paramClassIndex[k]], "  %10g", param_[k]);
        }
        for(auto fp: fpc)
        {
            fprintf(fp, "\n");
            fflush(fp);
        }
        if (nullptr != fpe)
        {
            fprintf(fpe, "%8f  %10g\n", xiter, prevEval);
            fflush(fpe);
        }
        if (iter >= halfIter)
        {
            for (auto k = 0; k < nParam; k++)
            {
                sum[k]       += param_[k];
                sum_of_sq[k] += gmx::square(param_[k]);
            }
            nsum++;
        }
    }
    if (nsum > 0)
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
    for(auto fp: fpc)
    {
        xvgrclose(fp);
    }
    if (nullptr != fpe)
    {
        xvgrclose(fpe);
    }
    return minEval;
}

}

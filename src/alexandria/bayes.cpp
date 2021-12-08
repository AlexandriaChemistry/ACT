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

#include "bayes.h"

#include <assert.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include <functional>
#include <string>
#include <vector>

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

    if (bch_.xvgConv().empty() || bch_.xvgEpot().empty())
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
    std::vector<std::string> pClass = bch_.paramClass();
    assignParamClassIndex(&paramClassIndex, &pClass);

    openParamSurveillanceFiles(pClass, &fpc, paramClassIndex);

    // Compute temperature weights if relevant, otherwise the numbers are all 1.0
    weightedTemperature_.resize(paramNames_.size(), 1.0);
    if (bch_.temperatureWeighting()) computeWeightedTemperature();

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

    double beta0 = 1/(BOLTZ*bch_.temperature());

    // Optimization loop
    for (int iter = 0; iter < bch_.maxIter(); iter++)
    {
        for (int pp = 0; pp < nParam; pp++)
        {
            // Pick a random parameter to change
            const int paramIndex = int_uniform(gen);

            // Do the step!
            stepMCMC(paramIndex, gen, real_uniform, &changed, &prevEval, &prevEval_testset,
                     bEvaluate_testset, pp, iter, &beta0, nParam, &minEval, fplog, fpc, fpe,
                     paramClassIndex);

            // For the second half of the optimization, collect data to find the mean and standard deviation of each
            // parameter
            if (iter >= bch_.maxIter()/2)
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

    computeMeanSigma(nParam, sum, nsum, &sum_of_sq);

    closeConvergenceFiles(fpc, fpe);

    bool bMinimum = false;  // Assume no better minimum was found
    if (minEval < (*chi2))  // If better minimum was found, update the value in <*chi2> and return true
    {
        (*chi2) = minEval;
        bMinimum = true;
    }
    return bMinimum;
}

void Bayes::stepMCMC(const int                                  paramIndex,
                           std::mt19937                        &gen,
                           std::uniform_real_distribution<>    &real_uniform,
                           std::vector<bool>                   *changed,
                           double                              *prevEval,
                           double                              *prevEval_testset,
                     const bool                                 bEvaluate_testset,
                     const int                                  pp,
                     const int                                  iter,
                           double                              *beta0,
                     const int                                  nParam,
                           double                              *minEval,
                           FILE                                *fplog,
                     const std::vector<FILE*>                  &fpc,
                           FILE                                *fpe,
                     const std::vector<int>                    &paramClassIndex)
{

    // Store the original value of the parameter
    const double storeParam = param_[paramIndex];

    // Change the parameter
    double rndNumber = real_uniform(gen);
    while (rndNumber == 0.5) rndNumber = real_uniform(gen);
    changeParam(paramIndex, rndNumber);

    attemptedMoves_[paramIndex] += 1;
    (*changed)[paramIndex]       = true;

    // Update FF parameter data structure with
    // the new value of parameter j
    toPoldata(*changed);

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
        if (bch_.anneal(iter)) (*beta0) = bch_.computeBeta(iter);
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
        toPoldata(*changed);
    }
    (*changed)[paramIndex] = false;  // Set changed[j] back to false for upcoming iterations

    fprintParameterStep(fpc, paramClassIndex, xiter);
    // If the chi2 surveillance file exists, write progress
    if (fpe != nullptr) fprintChi2Step(bEvaluate_testset, fpe, xiter, *prevEval, *prevEval_testset);

}

void Bayes::computeMeanSigma(const int     nParam,
                             const parm_t &sum,
                             const int     nsum,
                                   parm_t *sum_of_sq)
{

    if (nsum > 0)  // Compute mean and standard deviation
    {
        double ps2 = 0.0;
        for (auto k = 0; k < nParam; k++)
        {
            pmean_[k]        = (sum[k]/nsum);
            (*sum_of_sq)[k] /= nsum;
            ps2              = std::max(0.0, (*sum_of_sq)[k]-gmx::square(pmean_[k]));
            psigma_[k]       = sqrt(ps2);
        }
    }

}


}

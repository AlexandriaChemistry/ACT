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

#ifndef ALEXANDRIA_BAYES_H
#define ALEXANDRIA_BAYES_H

#include <functional>
#include <random>
#include <vector>

#include "gromacs/commandline/pargs.h"
#include "gromacs/fileio/oenv.h"
#include "gromacs/fileio/xvgr.h"
#include "gromacs/math/functions.h"
#include "gromacs/math/units.h"
#include "gromacs/utility/fatalerror.h"
#include "gromacs/utility/real.h"

#include "molselect.h"
#include "confighandler.h"

namespace alexandria
{

//! How to perform the calculation of deviations (chi-squared)
enum class CalcDev {
    //! Do it in parallel
    Parallel = 1,
    //! Do it on the master only
    Master = 2,
    //! Do the final one only (typically on the master)
    Final = 3
};

class Sensitivity
{
private:
    //! \brief parameter values used
    std::vector<double> p_;
    //! \brief chi2 values obtains
    std::vector<double> chi2_;
    //! \brief Constants for the parabola fitting
    double a_ = 0, b_ = 0, c_ = 0;
public:
    //! \brief Constructor
    Sensitivity() {}

    /*! \brief
     * Add a point
     * \param[in] p    The parameter value
     * \param[in] chi2 The chi-squared value
     */
    void add(double p, double chi2)
    {
        p_.push_back(p);
        chi2_.push_back(chi2);
    }
    /*! \brief
     * Compute the fit to the curve
     * \param[in] fp File pointer for debugging output
     */
    void computeForceConstants(FILE *fp);

    //! Return the constants after computation
    double a() const { return a_; }
    double b() const { return b_; }
    double c() const { return c_; }

    /*! \brief Print output
     * \param[in] fp    File pointer for output
     * \param[in] label Label for identifying the parameter
     */
    void print(FILE *fp, const std::string &label);
};

/*! \brief
 * Does Bayesian Monte Carlo (BMC) simulation to find the best parameter set,
 * which has the lowest chi-squared.
 *
 * \inpublicapi
 * \ingroup module_alexandria
 */
class Bayes
{
    using func_t       = std::function<double (double v[])>;
    using parm_t       = std::vector<double>;
    using mc_t         = std::vector<int>;
    using param_name_t = std::vector<std::string>;

    private:
        //! THIS IS NOT BEING USED!
        func_t                  func_;
        parm_t                  initial_param_;
        parm_t                  param_;
        std::vector<int>        ntrain_;
        parm_t                  psigma_;
        parm_t                  pmean_;
        parm_t                  lowerBound_;
        parm_t                  upperBound_;
        parm_t                  bestParam_;
        parm_t                  weightedTemperature_;
        mc_t                    attemptedMoves_;
        mc_t                    acceptedMoves_;
        std::vector<Mutability> mutability_;
        param_name_t            paramNames_;

    public:

        Bayes() {}

        /*! \brief
         * Run the Markov chain Monte carlo (MCMC) simulation
         * \param[in]  fplog            File pointer for logging info.
         *                              May be nullptr.
         * \param[in]  evaluate_testset If true, evaluate the energy on
         *                              the test set.
         * \param[out] chi2             pointer to chi2 in runMaster, at the end it will be the minimum
         * \return True if the energy decreased during the MCMC
         */
        bool MCMC(FILE *fplog, bool evaluate_testset, double *chi2);

        /*!
        * Take a step of MCMC by attempting to alter a parameter
        * @param paramIndex        index of the parameter to alter
        * @param gen               pointer to random number generator
        * @param real_uniform      pointer to random number distribution
        * @param changed           a reference to a vector which has true for parameters that change and false otherwise
        * @param prevEval          pointer to a double storage with the previous chi2 for training set
        * @param prevEval_testset  a pointer to a double storage with the previous chi2 for test set
        * @param bEvaluate_testset true if evaluation should be done on test set, false otherwise
        * @param pp                index of inner loop over number of parameters
        * @param iter              current iteration number
        * @param beta0             pointer to beta for annealing
        * @param nParam            number of parameters in the model
        * @param minEval           pointer to the minimum chi2 found so far for the training set
        * @param fplog             pointer to log file. May be nullptr
        * @param fpc               pointers to parameter surveillance files
        * @param fpe               pointer to chi2 surveillance file
        * @param paramClassIndex   class (by index) of each parameter in the model
        */
        void stepMCMC(const int                                 paramIndex,
                            std::mt19937                       &gen,
                            std::uniform_real_distribution<>   &real_uniform,
                            std::vector<bool>                  *changed,
                            double                             *prevEval,
                            double                             *prevEval_testset,
                      const bool                                bEvaluate_testset,
                      const int                                 pp,
                      const int                                 iter,
                            double                             *beta0,
                      const int                                 nParam,
                            double                             *minEval,
                            FILE                               *fplog,
                      const std::vector<FILE*>                 &fpc,
                            FILE                               *fpe,
                      const std::vector<int>                   &paramClassIndex);

        /*!
         * Write chi2 value to surveillance file
         * @param bEvaluate_testset     true if test set is evaluated, false otherwise
         * @param fpe                   pointer to chi2 surveillance file
         * @param xiter                 fractional iteration (3.6, 3.89, ...)
         * @param prevEval              chi2 fro training set
         * @param prevEval_testset      chi2 for test set
         */
        void fprintChi2Step(const bool      bEvaluate_testset,
                                  FILE     *fpe,
                            const double    xiter,
                            const double    prevEval,
                            const double    prevEval_testset);

        /*!
         * Compute mean (pmean_) and standard deviation (psigma_) for each parameter
         * @param nParam        number of parameters in the system
         * @param sum           over "nsum" iterations, the sum of each parameter
         * @param nsum          number of iterations to compute statistics over
         * @param sum_of_sq     over "nsum" iterations, the sum of each parameter squared
         */
        void computeMeanSigma(const int     nParam,
                              const parm_t &sum,
                              const int     nsum,
                                    parm_t *sum_of_sq);

        /*! \brief
         * Perform a sensitivity analysis by systematically changing
         * all parameters and re-evaluating the chi2.
         * \param[in] fplog    File pointer to print to
         * \param[in] ims   Data set to perform sensitivity analysis on
         */
        void SensitivityAnalysis(FILE *fplog, iMolSelect ims);

};

}

#endif

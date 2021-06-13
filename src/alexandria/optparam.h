/*
 * This source file is part of the Alexandria Chemistry Toolkit.
 *
 * Copyright (C) 2014-2020 
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

#ifndef ALEXANDRIA_OPTPARAM_H
#define ALEXANDRIA_OPTPARAM_H

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

namespace alexandria
{

/*! \brief
 * Does Bayesian Monte Carlo (BMC) simulation to find the best paramater set,
 * which has the lowest chi-squared.
 *
 * \inpublicapi
 * \ingroup module_alexandria
 */


class OptParam
{
    private:
        //! Maximum number of iterations
        int                      maxiter_       = 100;
        //! Reference MC acceptance ratio
        int                      mc_ratio_       = 30;
        //! Adapt Temperature every nAdapt steps
        int                      nAdapt_         = 50;
        //! Output environment structure
        const gmx_output_env_t  *oenv_          = nullptr;
        //! Use box constraints when optimizing
        gmx_bool                 bBoxConstraint_ = false;
        //! Random number seed
        real                     seed_           = -1;
        //! Relative step when optimizing
        real                     step_           = 0.02;
        //! Temperature in chi2 units
        real                     temperature_    = 5;
        //! Use annealing in the optimization. Value < 1 means annealing will happen
        real                     anneal_         = 1;
        //! Use adaptive MCMC in the optimization
        bool                     adaptive_       = false;
        //! Base name for parameter convergence file names
        std::string              xvgconv_;
        //! File name for parameter energy (chi2)
        std::string              xvgepot_;
        //! Parameter classes for printing
        std::vector<std::string> paramClass_;     
    public:
        /*! \brief Add command line arguments
         *
         * \param[in] pargs Vector of pargs
         */
        void add_pargs(std::vector<t_pargs> *pargs);

        /*! \brief Set the output file names. 
         *
         * The parameter values are split over
         * a number of files in order to make it easier to visualize the
         * results. The parameter classes should therefore match the
         * parameter names. E.g. a class could be alpha, another zeta.
         *
         * \param[in] xvgconv    The parameter convergence base name
         * \param[in] paramClass The parameter classes (e.g. zeta, alpha)
         * \param[in] xvgepot    The filename to print the chi2 value
         * \param[in] oenv       GROMACS utility structure
         */
        void setOutputFiles(const char                     *xvgconv,
                            const std::vector<std::string> &paramClass,
                            const char                     *xvgepot,
                            const gmx_output_env_t         *oenv);

        //! Return the class of parameters registered
        const std::vector<std::string> &paramClass() { return paramClass_; }
        //! \brief Return Max # iterations
        int maxIter() const { return maxiter_; }
        
        //! \brief Return temperature
        real temperature () const { return temperature_; }
        
        /*! \brief Compute and return the Boltzmann factor
         *
         * \param[in] iter  The iteration number
         * \return The Boltzmann factor
         */
        double computeBeta(int iter);
        
        /*! \brief Compute and return the Boltzmann factor
         * it applies periodic annealing
         *
         * \param[in] maxiter The maximum number of itearion
         * \param[in] iter    The iteration number
         * \param[in] ncycle  The multiplicity of the cosine function
         * \return The Boltzmann factor
         */
        double computeBeta(int maxiter, int iter, int ncycle);
        
        double adaptBeta(real ratio);

        //! \brief Return nAdapt
        int nAdapt () const { return nAdapt_; }

        //! \brief Return the step
        real step() const { return step_; }
        
        //! \brief Addapt the step size for perturbing the parameter
        void adapt_step(real factor) {step_ *= factor ;}

        //! \brief Use box constraints or not
        void setBoxConstraint(bool bBox) { bBoxConstraint_ = bBox; }

        //! \brief Return whether or not bounds are used for parameters
        bool boxConstraint() const { return bBoxConstraint_; }
        
        /*! \brief Return whether or not to do simulated annealing
         * \param iter The iteration number
         */
        bool anneal (int iter) const;
        
        // ! \brief Return whether or not call Adaptive MCMC
        bool adaptive () const { return adaptive_; }

        //! \brief Return xvg file for convergence information
        const std::string &xvgConv() const { return xvgconv_; }

        //! \brief Return xvg file for epot information
        const std::string &xvgEpot() const { return xvgepot_; }

        //! \brief Return output environment
        const gmx_output_env_t *oenv() const { return oenv_; }
        
        /*! \brief Save the current state
         * Must be overridden by child class.
         */
        virtual void saveState() = 0;
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
    
    //! \brief Add a point
    void add(double p, double chi2)
    {
        p_.push_back(p);
        chi2_.push_back(chi2);
    }
    //! \brief Compute the fit to the curve
    void computeForceConstants();
    
    double a() const { return a_; }
    double b() const { return b_; }
    double c() const { return c_; }
    
    void print(FILE *fp, const std::string &label);
};

class Bayes : public OptParam
{
    using func_t       = std::function<double (double v[])>;
    using parm_t       = std::vector<double>;
    using mc_t         = std::vector<int>;
    using param_name_t = std::vector<std::string>;

    private:
        func_t        func_;
        parm_t        initial_param_;
        parm_t        param_;
        parm_t        psigma_;
        parm_t        pmean_;
        parm_t        lowerBound_;
        parm_t        upperBound_;
        parm_t        bestParam_;
        mc_t          attemptedMoves_;
        mc_t          acceptedMoves_;
        param_name_t  paramNames_;

    public:

        Bayes() {}

        /*! \brief
         * Change parameter j based on a random unmber
         * obtained from a uniform distribution.
         */
        void changeParam(size_t j, real rand);

        //! \brief Return the number of parameters
        size_t nParam() const { return param_.size(); }

        /*! \brief
         * Append parameter and set it to value. Add bounds
         * as specified.
         * \param[in] name  String describing the parameter
         * \param[in] val   The value
         * \param[in] lower The new lower bound value
         * \param[in] upper The new lower bound value
         * \param[in] bRandom Generate random initial value for parameters if true.
         */
        void addParam(const std::string &name,
                      real val,
                      real lower,
                      real upper,
                      bool bRandom);
        /*! \brief
         * Append random parameter within the bounds specified.
         * \param[in] name  String describing the parameter
         * \param[in] lower The new lower bound value
         * \param[in] upper The new lower bound value
         */
        void addRandomParam(const std::string &name,
                            real               lower,
                            real               upper)
        {
            addParam(name, (lower+upper)*0.5, lower, upper, true);
        }

        /*! \brief
         * Set parameter j to a new value
         * \param[j]   Index
         * \param[val] The new value
         */
        void setParam(size_t j, real val)
        {
            GMX_RELEASE_ASSERT(j < param_.size(), "Parameter out of range");
            param_[j] = val;
        }

        /*! \brief
         * Set all parameters to the array passed
         */
        void setParam(parm_t param)
        {
            GMX_RELEASE_ASSERT(param.size() == param_.size() || param_.empty(),
                               "Incorrect size of input parameters");
            param_ = param;
        }
        /*! \brief
         * Returns the current vector of parameters.
         */
        const parm_t &getInitialParam() const { return initial_param_; }
        
        /*! \brief
         * Returns the current vector of parameters.
         */
        const parm_t &getParam() const { return param_; }

        /*! \brief
         * Returns the vector of best found value for each parameter.
         */
        const parm_t &getBestParam() const { return bestParam_; }

        /*! \brief
         * Returns the vector of mean value calculated for each parameter.
         */
        const parm_t &getPmean() const { return pmean_; }

        /*! \brief
         * Returns the vector of standard deviation calculated for each parameter.
         */
        const parm_t &getPsigma() const { return psigma_; };

        /*! \brief
         * Return the vector of parameter names. 
         */
        const param_name_t &getParamNames() const { return paramNames_; };

        /*! \brief
         * Print the paramters to a file
         * \param[in] fp File pointer to open file
         */
        void printParameters(FILE *fp) const;
        /*! \brief
         * Return the vector of number of attempted moves for each parameter
         */
        const mc_t &getAttemptedMoves() const {return attemptedMoves_;};
        
        /*! \brief
         * Return the vector of number of accepted moves for each parameter
         */
        const mc_t &getAcceptedMoves() const {return acceptedMoves_;};

        /*! \brief
         * Run the Markov chain Monte carlo (MCMC) simulation
         * \param[in] fplog File pointer for logging info. May be nullptr.
         * \return minimum energy value
         */
        double MCMC(FILE *fplog);
        
        /*! \brief
         * Run adaptive Markov chain Monte carlo (MCMC) simulation
         * \param[in] fplog File pointer for logging info. May be nullptr.
         * \return minimum energy value
         */
        double Adaptive_MCMC(FILE *fplog);

        /*! \brief
         * Perform a sensitivity analysis by systematically changing
         * all parameters and re-evaluating the chi2.
         */
        void SensitivityAnalysis(FILE *fplog);

        /*! \brief
         * Copy the optimization parameters to the poldata structure
         * \param[in] changed List over the parameters that have changed.
         */
        virtual void toPoldata(const std::vector<bool> &changed) = 0;

        //! Compute the chi2 from the target function
        virtual double calcDeviation() = 0;

        /*! Return number of planned function calls 
         * Return the number of calls to the objective function
         * that will be made by the Bayes::MCMC
         */
        size_t numberObjectiveFunctionCalls() const
        {
            return 1+maxIter()*nParam();
        }
};

}

#endif

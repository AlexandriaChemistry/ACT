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

#include "molselect.h"
#include "tune_eem.h"

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
        //! Output environment structure
        const gmx_output_env_t  *oenv_          = nullptr;
        //! Random number seed
        real                     seed_           = -1;
        //! Relative step when optimizing
        real                     step_           = 0.02;
        //! Temperature in chi2 units
        real                     temperature_    = 5;
        //! Weight temperature after number of training points
        bool                     tempWeight_     = false;
        //! Weighted temperatures
        std::vector<double>      weightedTemperature_;
        //! Use annealing in the optimization. Value < 1 means annealing will happen
        real                     anneal_         = 1;
        //! Flag determining whether to be verbose printing
        bool                     verbose_        = false;
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

        //! \brief Return verbosity
        bool verbose() const { return verbose_; }
        
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
        
        //! \brief Return the step
        real step() const { return step_; }
        
        //! \brief Return whether or not temperature weighting should be considered
        bool temperatureWeighting() const { return tempWeight_; }
        
        /*! \brief Return whether or not to do simulated annealing
         * \param iter The iteration number
         */
        bool anneal (int iter) const;
        
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
        std::vector<int> ntrain_;
        parm_t        psigma_;
        parm_t        pmean_;
        parm_t        lowerBound_;
        parm_t        upperBound_;
        parm_t        bestParam_;
        parm_t        weightedTemperature_;
        mc_t          attemptedMoves_;
        mc_t          acceptedMoves_;
        std::vector<Mutability> mutability_;
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
         * \param[in] name    String describing the parameter
         * \param[in] val     The value
         * \param[in] mut     Mutability
         * \param[in] lower   The new lower bound value
         * \param[in] upper   The new lower bound value
         * \param[in] ntrain  Number of copies in the training set
         * \param[in] bRandom Generate random initial value for parameters if true.
         */
        void addParam(const std::string &name,
                      real               val,
                      Mutability         mut,
                      real               lower,
                      real               upper,
                      int                ntrain,
                      bool               bRandom);
        /*! \brief
         * Append random parameter within the bounds specified.
         * \param[in] name  String describing the parameter
         * \param[in] mut     Mutability
         * \param[in] lower The new lower bound value
         * \param[in] upper The new lower bound value
         * \param[in] ntrain  Number of copies in the training set
         */
        void addRandomParam(const std::string &name,
                            Mutability         mut,
                            real               lower,
                            real               upper,
                            int                ntrain)
        {
            // TODO: Make this random for real
            addParam(name, (lower+upper)*0.5, mut, lower, upper, ntrain, true);
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
         * Returns the current vector of lower bounds
         * @return the current vector of lower bounds
         */
        const parm_t &getLowerBound() const { return lowerBound_; }

        /*! \brief
         * Returns the current vector of upper bounds
         * @return the current vector of upper bounds
         */
        const parm_t &getUpperBound() const { return upperBound_; }

        /*! \brief
         * Returns the number of training points per parameter
         */
        const std::vector<int> &getNtrain() const { return ntrain_; }

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
         * \param[in]  fplog            File pointer for logging info. 
         *                              May be nullptr.
         * \param[in]  evaluate_testset If true, evaluate the energy on 
         *                              the test set.
         * \param[out] chi2             pointer to chi2 in runMaster, at the end it will be the minimum
         * \return True if the energy decreased during the MCMC
         */
        bool MCMC(FILE *fplog, bool evaluate_testset, double *chi2);

        /*!
         * Compute weighted temperature for each parameter
         */
        void computeWeightedTemperature();

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
                            std::mt19937&                       gen,
                            std::uniform_real_distribution<>&   real_uniform,
                            std::vector<bool>&                  changed,
                            double*                             prevEval,
                            double*                             prevEval_testset,
                      const bool                                bEvaluate_testset,
                      const int                                 pp,
                      const int                                 iter,
                            double*                             beta0,
                      const int                                 nParam,
                            double*                             minEval,
                            FILE*                               fplog,
                            std::vector<FILE*>&                 fpc,
                            FILE*                               fpe,
                            std::vector<int>&                   paramClassIndex);

        /*!
         * Assign a class (by index) to each parameter
         * @param paramClassIndex   for each parameter, will have index of the class it belongs to
         * @param pClass            class types
         */
        void assignParamClasses(std::vector<int>&           paramClassIndex,
                                std::vector<std::string>&   pClass);

        /*!
         * Open parameter convergence surveillance files
         * @param pClass            different classes of parameters
         * @param fpc               vector to append pointers to parameter convergence files
         * @param paramClassIndex   for each parameter, to which class (by index) it belongs
         */
        void openParamSurveillanceFiles(const std::vector<std::string>&  pClass,
                                              std::vector<FILE*>&        fpc,
                                              std::vector<int>&          paramClassIndex);

        /*!
         * Open a chi2 surveillance file
         * @param bEvaluate_testset     whether the test set will be evaluated
         * @return                      a pointer to the opened file
         */
        FILE* openChi2SurveillanceFile(const bool bEvaluate_testset);

        /*!
         * Close chi2 and parameter convergence files
         * @param fpc   vector of pointers to parameter convergence files
         * @param fpe   pointer to chi2 convergence file
         */
        void closeConvergenceFiles(std::vector<FILE*>& fpc,
                                   FILE*               fpe);


        /*!
         * Print new minimum to log file and, if necessary, print params to debug file
         * @param fplog                 pointer to log file
         * @param bEvaluate_testset     true if test set is evaluated, false otherwise
         * @param xiter                 fractional iteration. E.g, if we are halfway through iteration 3 it is 3.5
         * @param currEval              current chi2 in training set
         * @param currEval_testset      current chi2 in test set
         */
        void fprintNewMinimum(      FILE*   fplog,
                              const bool    bEvaluate_testset,
                              const double  xiter,
                              const double  currEval,
                              const double  currEval_testset);

        /*!
         * Print parameter values to their respective surveillance files
         * @param fpc                   pointer to each parameter surveillance file
         * @param paramClassIndex       class index of each parameter
         * @param xiter                 fractional iteration (e.g. 3.5, 2.89, ...)
         */
        void fprintParameterStep(      std::vector<FILE*>&   fpc,
                                 const std::vector<int>&     paramClassIndex,
                                 const double                xiter);

        /*!
         * Write chi2 value to surveillance file
         * @param bEvaluate_testset     true if test set is evaluated, false otherwise
         * @param fpe                   pointer to chi2 surveillance file
         * @param xiter                 fractional iteration (3.6, 3.89, ...)
         * @param prevEval              chi2 fro training set
         * @param prevEval_testset      chi2 for test set
         */
        void fprintChi2Step(const bool      bEvaluate_testset,
                                  FILE*     fpe,
                            const double    xiter,
                            const double    prevEval,
                            const double    prevEval_testset);

        /*!
         * Compute mean (pmean_) and standard deviation (psigma_) for each parameter
         * @param nParam        number of parameters in the system
         * @param sum           over <nsum> iterations, the sum of each parameter
         * @param nsum          number of iterations to compute statistics over
         * @param sum_of_sq     over <nsum> iterations, the sum of each parameter squared
         */
        void computeMeanSigma(const int     nParam,
                              const parm_t& sum,
                              const int     nsum,
                              const parm_t& sum_of_sq);
        
        /*! \brief
         * Perform a sensitivity analysis by systematically changing
         * all parameters and re-evaluating the chi2.
         * \param[in] fplog    File pointer to print to
         * \param[in] ims   Data set to perform sensitivity analysis on
         */
        void SensitivityAnalysis(FILE *fplog, iMolSelect ims);

        /*! \brief
         * Copy the optimization parameters to the poldata structure
         * \param[in] changed List over the parameters that have changed.
         */
        virtual void toPoldata(const std::vector<bool> &changed) = 0;

        /*! \brief
         * Compute the chi2 from the target function
         * \param[in] verbose Whether or not to print stuff
         * \param[in] calcDev How to compute the deviation for all compounds
         * \param[in] ims     The data set to evaluate the energy upon.
         * \return the square deviation or -1 when done.
         */
        virtual double calcDeviation(bool       verbose,
                                     CalcDev    calcDev,
                                     iMolSelect ims) = 0;

        /*! Return number of planned function calls 
         * Return the number of calls to the objective function
         * that will be made by the Bayes::MCMC
         */
        size_t numberObjectiveFunctionCalls() const
        {
            return 1+maxIter()*nParam();
        }
        /* \brief
         * Print the MC statistics to a file.
         * \param[in] fp File pointer to print to
         */
        void printMonteCarloStatistics(FILE *fp);
};

}

#endif

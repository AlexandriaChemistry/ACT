#ifndef ALEXANDRIA_CONFIGHANDLER_H
#define ALEXANDRIA_CONFIGHANDLER_H

#include <vector>

#include "gromacs/commandline/pargs.h"

namespace alexandria
{

/*!
 * Abstract class to handle configuration for methods.
 * It adds its command-line arguments to a given parameter vector
 * and later checks its validity.
 */
class ConfigHandler
{

public:

    /*!
     * \brief Add command-line arguments to a vector
     * @param pargs     pointer to arguments vector
     */
    virtual void add_pargs(std::vector<t_pargs> *pargs) = 0;

    /*!
     * \brief Check the validity of the provided arguments
     * Throw an exception if an invalid combination of arguments is encountered
     */
    virtual void check_pargs() = 0;

};

class GAConfigHandler
{

private:
    int popSize_ = 8;

public:

    /*!
     * \brief Add command-line arguments to a vector
     * @param pargs     pointer to arguments vector
     */
    virtual void add_pargs(std::vector<t_pargs> *pargs);

    /*!
     * \brief Check the validity of the provided arguments
     * Throw an exception if an invalid combination of arguments is encountered
     */
    virtual void check_pargs();
};

/*!
 * Handles optimization parameters for Markov Chain Monte-Carlo method
 */
class BayesConfigHandler : public ConfigHandler
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
    /*!
     * \brief Add command-line arguments to a vector
     * @param pargs     pointer to arguments vector
     */
    virtual void add_pargs(std::vector<t_pargs> *pargs);

    /*!
     * \brief Check the validity of the provided arguments
     * Throw an exception if an invalid combination of arguments is encountered
     */
    virtual void check_pargs();

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
    real temperature() const { return temperature_; }

    /*! \brief Compute and return the Boltzmann factor
    *
    * \param[in] iter  The iteration number
    * \return The Boltzmann factor
    */
    double computeBeta(int iter);

    /*! \brief Compute and return the Boltzmann factor
    * it applies periodic annealing
    *
    * \param[in] maxiter The maximum number of iteration
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

}

#endif //ALEXANDRIA_PARAMHANDLER_H
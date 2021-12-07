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

/*!
 * Handles optimization parameters for Genetic Algorithm
 */
class GAConfigHandler : public ConfigHandler
{

private:

    // First non-NULL value indicates the default value
    // After argument parsing, first element in the array will point to the selected enum value, so optimizer_[0]
    // Static means the variable will be shared among objects (only 1 place in memory)
    //! Optimizer to use
    // static const char *optimizer_[] = {nullptr, "MCMC", "GA", "HYBRID", nullptr};
    //! Population size
    int popSize_ = 1;
    //! Amount of elites in the population
    int nElites_ = 0;
    //! Order of crossover operator
    int nCrossovers_ = 1;
    //! Sorter algorithm
    // static const char *sorter_[] = {nullptr, "QUICK", "MERGE", "NONE", nullptr};
    //! Probability computing algorithm
    // static const char *probComputer_[] = {nullptr, "RANK", "FITNESS", "BOLTZMANN", nullptr};
    //! Boltzmann probability temperature. TODO: This temperature should be lowered over time.
    real boltzTemp_ = 1;
    // TODO: Termination options???
    //! Probability of crossover
    real prCross_ = 0.35;
    //! Probability of mutation
    real prMut_ = 0.01;

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

    // TODO: Getters!

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
    //! Use annealing in the optimization. Value < 1 means annealing will happen
    real                     anneal_         = 1;
    //! Flag determining whether to be verbose printing TODO: This has to be made a global flag!
    bool                     verbose_        = false;

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

    //! \brief Return output environment
    const gmx_output_env_t *oenv() const { return oenv_; }
    
};


} //namespace alexandria


#endif //ALEXANDRIA_PARAMHANDLER_H

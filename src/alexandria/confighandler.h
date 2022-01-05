/*! \internal \brief
 * Implements part of the alexandria program.
 * \author Mohammad Mehdi Ghahremanpour <mohammad.ghahremanpour@icm.uu.se>
 * \author David van der Spoel <david.vanderspoel@icm.uu.se>
 * \author Julian Ramon Marrades Furquet <julian.marrades@hotmail.es>
 * \author Oskar Tegby <oskar.tegby@it.uu.se>
 */


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
 * FIXME: Should this be under the ga directory???
 */
class GAConfigHandler : public ConfigHandler
{

private:

    // First non-NULL value indicates the default value
    // After argument parsing, first element in the array will point to the selected enum value, so optimizer_[0]
    // Static means the variable will be shared among objects (only 1 place in memory)
    //! Optimizer to use
    const char *optimizer_[5] = {nullptr, "MCMC", "GA", "HYBRID", nullptr};
    //! Population size
    int popSize_ = 1;
    //! Amount of elites in the population
    int nElites_ = 0;
    //! Whether to initialize the individuals randomly
    bool randomInit_ = true;
    //! Order of crossover operator
    int nCrossovers_ = 1;
    //! Sorter algorithm
    const char *sorter_[5] = {nullptr, "QUICK", "MERGE", "NONE", nullptr};
    //! Probability computing algorithm
    const char *probComputer_[5] = {nullptr, "RANK", "FITNESS", "BOLTZMANN", nullptr};
    //! Boltzmann probability temperature. TODO: This temperature should be lowered over time.
    real boltzTemp_ = 1;
    //! Probability of crossover
    real prCross_ = 0.35;
    //! Probability of mutation
    real prMut_ = 0.01;
    //! For PercentMutator: Maximum allowed change in a parameter as a fraction of its allowed range
    real percent_ = 0.1;

    // TODO: Improve termination criteria
    //! Generation limit in Genetic Algorithm
    int maxGenerations_ = 10;

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

    /* * * * * * * * * * * * * * * * * * * * * *
    * BEGIN: Getters and setters               *
    * * * * * * * * * * * * * * * * * * * * * */

    //! \return the optimizer
    const char *optimizer() { return optimizer_[0]; }

    //! \return the size of the population
    int popSize() const { return popSize_; }

    /*! \brief Set the population size
     * \param[in] size The size
     */
    void setPopSize(int size) { popSize_ = size; }
    
    //! \return the amount of top individuals that pass, unchanged, to the next generation
    int nElites() const { return nElites_; }

    //! \return whether to initialize an individual randomly
    real randomInit() const { return randomInit_; }

    //! \return the order of the crossover operator
    int nCrossovers() const { return nCrossovers_; }

    /*! \brief Set the number of crossovers
     * \param[in] number The number
     */
    void setCrossovers(int number) { nCrossovers_ = number; }
    
    //! \return the sorter
    const char *sorter() const { return sorter_[0]; }

    //! \return the probability computer
    const char *probComputer() const { return probComputer_[0]; }

    //! \return the Boltzmann temperature parameter
    real boltzTemp() const { return boltzTemp_; }

    //! \return the probability of crossover
    real prCross() const { return prCross_; }

    //! \return the probability of mutation
    real prMut() const { return prMut_; }

    //! \return For PercentMutator: Maximum allowed change in a parameter as a fraction of its allowed range
    real percent() const { return percent_; }

    //! \return the generation limit
    int maxGenerations() const { return maxGenerations_; }

    /* * * * * * * * * * * * * * * * * * * * * *
    * END: Getters and setters                 *
    * * * * * * * * * * * * * * * * * * * * * */

};

/*!
 * Handles optimization parameters for Markov Chain Monte-Carlo method
 */
class BayesConfigHandler : public ConfigHandler
{

private:
    //! Maximum number of iterations
    int   maxiter_           = 100;
    //! Random number seed for the MCMCMutator random number generator
    int   seed_              = 0;
    //! Relative step when optimizing
    real  step_              = 0.02;
    //! Temperature in chi2 units
    real  temperature_       = 5;
    //! Weight temperature after number of training points
    bool  tempWeight_        = false;
    //! Use annealing in the optimization. Value < 1 means annealing will happen
    real  anneal_            = 1;
    //! Evaluate on test set during the MCMC run
    bool  evaluate_testset_  = false;

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

    //! \brief Return Max # iterations
    int maxIter() const { return maxiter_; }

    /*! \brief Set the number of iterations
     * \param[in] maxiter The max number of iterations
     */
    void setMaxIter(int maxiter) { maxiter_ = maxiter; }

    //! \brief Return temperature
    real temperature() const { return temperature_; }

    //! \return the seed
    int seed() const { return seed_; }

    /*! \brief Set the random number seed (0 is generate)
     * \param[in] seed The new seed
     */
    void setSeed(int seed) { seed_ = seed; }
    
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

    //! \return whether test set should be evaluated during the MCMC run
    bool evaluateTestset() const { return evaluate_testset_; }
    
};


} //namespace alexandria


#endif //ALEXANDRIA_PARAMHANDLER_H

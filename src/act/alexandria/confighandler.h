/*
 * This source file is part of the Alexandria Chemistry Toolkit.
 *
 * Copyright (C) 2020-2025
 *
 * Developers:
 *             Mohammad Mehdi Ghahremanpour, 
 *             Julian Marrades,
 *             Marie-Madeleine Walz,
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
 * \author Julian Ramon Marrades Furquet <julian@marrad.es>
 * \author Oskar Tegby <oskar.tegby@it.uu.se>
 */
#ifndef ALEXANDRIA_CONFIGHANDLER_H
#define ALEXANDRIA_CONFIGHANDLER_H

#include <vector>

#include "gromacs/commandline/pargs.h"

namespace alexandria
{

class MsgHandler;
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
     * @param fnms      pointer to vector of filenames
     */
    virtual void add_options(std::vector<t_pargs> *pargs,
                             std::vector<t_filenm> *fnms) = 0;

    /*!
     * \brief Check the validity of the provided arguments
     * Throw an exception if an invalid combination of arguments is encountered
     */
    virtual void check_pargs(MsgHandler *msghandler) = 0;

};

//! \brief Class to select the optimizer to use
enum class OptimizerAlg { MCMC, GA, HYBRID };

/*! \brief Convert string into OptimizerAlg
 * \param[in] str The string
 * \return the optimizer algorithm
 * \throws if there is no optimizer matching the string
 */
OptimizerAlg stringToOptimizerAlg(const std::string &str);

/*! \brief Convert OptimizerAlg to string
 * \param[in] opt The optimizer
 * \return The corresponding string
 */
const std::string &optimizerAlgToString(OptimizerAlg opt);

//! \brief Enum for probability computer
enum class ProbabilityComputerAlg { pcRANK, pcFITNESS, pcBOLTZMANN };

/*! \brief Convert string into ProbabilityComputerAlg
 * \param[in] str The string
 * \return the ProbabilityComputer algorithm
 * \throws if there is no ProbabilityComputer matching the string
 */
ProbabilityComputerAlg stringToProbabilityComputerAlg(const std::string &str);

/*! \brief Convert ProbabilityComputerAlg to string
 * \param[in] opt The ProbabilityComputer
 * \return The corresponding string
 */
const std::string &probabilityComputerAlgToString(ProbabilityComputerAlg opt);

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
    OptimizerAlg optAlg_ = OptimizerAlg::GA;
    //! Whether to evaluate loss for the test set
    bool evaluateTestset_ = false;
    //! Population size
    int popSize_ = 1;
    //! Amount of elites in the population
    int nElites_ = 0;
    //! Whether to initialize the individuals randomly
    bool randomInit_ = true;
    //! Order of crossover operator
    int nCrossovers_ = 1;
    //! Whether we sort the population or not
    bool sort_ = true;
    //! Probability computing algorithm
    ProbabilityComputerAlg pcAlg_ = ProbabilityComputerAlg::pcRANK;
    //! Boltzmann probability temperature.
    real boltzTemp_ = 1;
    //! Anneal start (as fraction of maxGenerations) for Boltzmann temperature
    real boltzAnneal_ = 1;
    //! Probability of crossover
    real prCross_ = 0.35;
    //! Probability of mutation
    real prMut_ = 0.01;
    //! For PercentMutator: Maximum allowed change in a parameter as a fraction of its allowed range
    real percent_ = 0.1;
    //! Generation limit in Genetic Algorithm
    int maxGenerations_ = 10;
    //! Generation limit for the test fitness to improve
    int maxTestGenerations_ = -1;
    //! Whether to compute the volume in logarithmic scale
    bool logVolume_ = false;
    //! For VolumeFractionPenalizer, the limit of the volume fraction
    real vfpVolFracLimit_ = -1;
    //! For VolumeFractionPenalizer, the fraction of worst genomes to reinitialize
    real vfpPopFrac_ = 0.8;
    //! For CatastrophePenalizer, the interval (in generations) between each penalty
    int cpGenInterval_ = -1;
    //! For CatastrophePenalizer, the fraction of genomes to reinitialize
    real cpPopFrac_ = 0.5;

public:

    /*!
     * \brief Add command-line arguments to a vector
     * @param pargs     pointer to arguments vector
     * @param fnms      pointer to filenames
     */
    virtual void add_options(std::vector<t_pargs>  *pargs,
                             std::vector<t_filenm> *fnms);

    /*!
     * \brief Check the validity of the provided arguments
     * \param[in] msghandler Message Handler
     * Throw an exception if an invalid combination of arguments is encountered
     */
    virtual void check_pargs(MsgHandler *msghandler);

    /* * * * * * * * * * * * * * * * * * * * * *
    * BEGIN: Getters and setters               *
    * * * * * * * * * * * * * * * * * * * * * */

    //! \return the optimizer
    OptimizerAlg optimizer() const { return optAlg_; }

    /*! \brief Set the optimizer algorithm
     * \param[in] optAlg The new algorithm
     */
    void setOptimizerAlg(OptimizerAlg optAlg) { optAlg_ = optAlg; }
    
    //! \return the size of the population
    int popSize() const { return popSize_; }

    /*! \brief Set the population size
     * \param[in] size The size
     */
    void setPopSize(int size) { popSize_ = size; }
    
    //! \return the amount of top individuals that pass, unchanged, to the next generation
    int nElites() const { return nElites_; }

    //! \return whether to initialize an individual randomly
    bool randomInit() const { return randomInit_; }

    //! \return the order of the crossover operator
    int nCrossovers() const { return nCrossovers_; }

    /*! \brief Set the number of crossovers
     * \param[in] number The number
     */
    void setCrossovers(int number) { nCrossovers_ = number; }
    
    //! \return whether we sort the population or not
    bool sort() const { return sort_; }

    //! \return the probability computer
    ProbabilityComputerAlg probabilityComputerAlg() const { return pcAlg_; }

    /*!
     * \brief Set a new probabilityComputerAlgorithm
     * \param[in] pcAlg the new algorithm
     */
    void setProbabilityComputerAlg(ProbabilityComputerAlg pcAlg) { pcAlg_ = pcAlg; }

    //! \return the Boltzmann temperature parameter
    real boltzTemp() const { return boltzTemp_; }

    //! \return the Boltzmann anneal starting fraction
    real boltzAnneal() const { return boltzAnneal_; };

    //! \return the probability of crossover
    real prCross() const { return prCross_; }

    //! \return the probability of mutation
    real prMut() const { return prMut_; }

    /*! \brief Set a new value for probability of mutation
     * \param[in] prMut the new mutation probability
     */
    void setPrMut(const real prMut) { prMut_ = prMut; }

    //! \return For PercentMutator: Maximum allowed change in a parameter as a fraction of its allowed range
    real percent() const { return percent_; }

    //! \return the generation limit
    int maxGenerations() const { return maxGenerations_; }

    //! \return the generation limit for the test fitness to improve
    int maxTestGenerations() const { return maxTestGenerations_; }

    //! \return true if we must evaluate fitness on test set, false otherwise
    bool evaluateTestset() const { return evaluateTestset_; }

    /*! \brief Turn test set evaluation on or off
     * \param[in] bOnOff On or Off
     */
    void setEvaluateTestset(bool bOnOff) {  evaluateTestset_ = bOnOff; }

    //! \return true whether the volume will be computed in log scale, false otherwise
    bool logVolume() const { return logVolume_; }

    //! \return for VolumeFractionPenalizer, the limit of the volume fraction allowed
    real vfpVolFracLimit() const { return vfpVolFracLimit_; }

    //! \return for VolumeFractionPenalizer, the fraction of the worst genomes to randomize
    real vfpPopFrac() const { return vfpPopFrac_; }

    //! \return for CatastrophePenalizer, the interval (in generations) between each penalty
    int cpGenInterval() const { return cpGenInterval_; }

    //! \return for CatastrophePenalizer, the fraction of genomes to reinitialize
    real cpPopFrac() const { return cpPopFrac_; }

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
    //! Random number seed for initializer, crossover, mutator, and selector FIXME: This should be general
    int   seed_              = 0;
    //! Relative step when optimizing
    real  step_              = 0.02;
    //! Temperature in chi2 units
    real  temperature_       = 5;
    //! Weight temperature after number of training points
    bool  tempWeight_        = false;
    //! Use annealing in the optimization. Value < 1 means annealing will happen
    real  anneal_            = 1;
    //! Whether annealing restart during each mutation round or is global over GA generations
    bool  annealGlobally_    = false;
    //! Evaluate on test set during the MCMC run
    bool  evaluate_testset_  = false;
    //! Checkpointing on or not?
    bool checkPoint_         = false;
    //! Tolerance for shell
    double shellToler_       = 1e-3;
    //! Maxiter for shell optimization
    int shellMaxIter_        = 25;
    //! Max distance between shell and core
    double shellMaxDistance_ = 0.02;
public:
    /*!
     * \brief Add command-line arguments to a vector
     * @param pargs     pointer to arguments vector
     * @param fnms      pointer to filenames
     */
    virtual void add_options(std::vector<t_pargs>  *pargs,
                             std::vector<t_filenm> *fnms);


    /*!
     * \brief Check the validity of the provided arguments
     * Throw an exception if an invalid combination of arguments is encountered
     */
    virtual void check_pargs(MsgHandler *msghandler);

    //! \brief Return Max # iterations
    int maxIter() const { return maxiter_; }

    /*! \brief Set the number of iterations
     * \param[in] maxiter The max number of iterations
     */
    void setMaxIter(int maxiter) { maxiter_ = maxiter; }

    //! \brief Return temperature
    real temperature() const { return temperature_; }

    /*! \brief Return temperature dependent on generation if global annealing is used
     * \param[in] generation      Current generation
     * \param[in] max_generations Total number of generations
     */
    real temperature(int generation,
                     int max_generations) const;

    /*! \brief set a new value for temperature
     * \param[in] temperature the new temperature
     */
    void setTemperature(const real temperature) { temperature_ = temperature; }

    //! \return whether or not to checkpoint
    bool checkPoint() const { return checkPoint_; }

    //! \return the seed
    int seed() const { return seed_; }

    /*! \brief Set the random number seed (0 is generate)
     * \param[in] seed The new seed
     */
    void setSeed(int seed) { seed_ = seed; }

    /*! \brief Compute and return the Boltzmann factor
     * Takes into account both local annealing (within a generation)
     * and global annealing (spreading over generations) depending
     * on user specified command line options. Will return 1/temperature
     * that is, not multiply with the Boltzmann constant. This is because
     * the deviation used in the MCMC algorithm is not an energy anyway.
     *
     * \param[in] generation      The generation number
     * \param[in] max_generations The total number of generations planned
     * \param[in] iteration       The iteration number
     * \return The Boltzmann factor
     */
    double computeBeta(int generation,
                       int max_generations,
                       int iteration);

    //! \brief Return the step
    real step() const { return step_; }

    //! \brief Set the step
    void setStep(real step) { step_ = step; }

    //! \brief Return whether or not temperature weighting should be considered
    bool temperatureWeighting() const { return tempWeight_; }

    /*! \brief Return whether or not to do simulated annealing
     * \param[in] generation The generation number for HYBRID
     * \param[in] iteration  The iteration number
     */
    bool anneal (int generation,
                 int iteration) const;

    //! \return start of annealing in fraction of iterations
    double annealStart() const { return anneal_; }

    //! \return true if annealing was requested
    bool annealing() const { return anneal_ < 1; }

    //! \return true if global annealing was requested
    bool globalAnnealing() const { return annealing() && annealGlobally_; }

    /*! \brief Set a new value for annealing start
     * \param[in] anneal the new starting point for simulated annealing
     */
    void setAnneal(const real anneal) { anneal_ = anneal; };

    //! \return whether test set should be evaluated during the MCMC run
    bool evaluateTestset() const { return evaluate_testset_; }

    /*! \brief Turn test set evaluation on or off
     * \param[in] bOnOff On or Off
     */
    void setEvaluateTestset(bool bOnOff) {  evaluate_testset_ = bOnOff; }

    //!\return shell tolerance
    double shellToler() const { return shellToler_; }

    //! \return shell iterations
    int shellMaxIter() const { return shellMaxIter_; }

    //! \return max shell distance from core
    double shellMaxDistance() const { return shellMaxDistance_; }
};

//! \brief Enum for selecting the coordinate minimizer in simulate
enum class eMinimizeAlgorithm {
    LBFGS, Steep, Newton
};

//! \brief Map algorithm to string
const std::string &eMinimizeAlgorithmToString(eMinimizeAlgorithm e);

/*! \brief Map string to algorithm
 * \param[in] str The string
 * \return the algorithm corresponding to the string or the default if string does not match any known minimizer
 */
eMinimizeAlgorithm stringToEMinimizeAlgorithm(const std::string &str);

class SimulationConfigHandler : ConfigHandler
{
private:
    //! Number of integration steps
    int                nsteps_               = 0;
    //! The integration time step
    double             deltat_               = 0.0002;
    //! Initial simulation temperature
    double             temperature_          = 0;
    //! Random number seed for generating velocities
    int                seed_                 = 0;
    //! How often to write coordinates
    int                nstxout_              = 1;
    //! How often to write velocities
    int                nstvout_              = 0;
    //! How often to write energies
    int                nstener_              = 1;
    //! Write shells to trajectory and coordinates
    bool               writeShells_          = false;
    //! Minmize (before MD)
    bool               minimize_             = false;
    //! Compute energies and quit
    bool               singlePoint_          = false;
    //! Minimization algorithm
    eMinimizeAlgorithm minAlg_               = eMinimizeAlgorithm::LBFGS;
    //! Tolerance on mean square force for minimizer.
    double             forceToler_           = 1e-6;
    //! Number of retries for minimizing
    int                minimizeRetries_      = 1;
    //! Max random displacement (nm) before re-trying to minize a structure
    double             minimizeDisplacement_ = 0.002;
    //! Whether to force using displacement and reminimize even if converged
    bool               forceReminimize_      = false;
    //! Apply overrelaxation (if > 1) to speed up minimization. Can be dangerous for poor energy functions.
    double             overRelax_            = 1.0;
    //! Maximum number of iterations for the energy minimizer, 0 is until convergence.
    int                maxIter_              = 100;
public:
    /*!
     * \brief Add command-line arguments to a vector
     * @param pargs     pointer to arguments vector
     */
    void add_MD_options(std::vector<t_pargs> *pargs);

    /*!
     * \brief Add command-line arguments to a vector
     * @param pargs     pointer to arguments vector
     * @param fnms      pointer to filenames
     */
    virtual void add_options(std::vector<t_pargs>  *pargs,
                             std::vector<t_filenm> *fnms);


    /*!
     * \brief Check the validity of the provided arguments
     * Throw an exception if an invalid combination of arguments is encountered
     */
    virtual void check_pargs(MsgHandler *msghandler);

    //! \return number of MD steps
    int nsteps() const { return nsteps_; }
    
    //! \return MD integration time step (ps)
    double deltat() const { return deltat_; }
    
    //! \return number of steps between coordinate output
    int nstxout() const { return nstxout_; }
    
    //! \return number of steps between energy output
    int nstener() const { return nstener_; }

    //! \return initial simulation temperature
    double temperature() const { return temperature_; }
    
    //! \return user provided random number seed
    int seed() const { return seed_; }

    //! \return whether or not to write shells to trajectory
    bool writeShells() const { return writeShells_; }
    
    //! \return whether or not to minimize the energy with respect to input coordinates
    bool minimize() const { return minimize_; }
    
    //! \return whether or not to do a single point energy calculation
    bool singlePoint() const { return singlePoint_; }

    //! Set the minimize option
    void setMinimize(bool minimize) { minimize_ = minimize; }

    //! \return the minimization algorithm
    eMinimizeAlgorithm minAlg() const { return minAlg_; }

    //! \return maximum number of retries for minimization
    int minimizeRetries() const { return minimizeRetries_; }

    /*! \brief Set number of retries
     * \param[in] maxretries The number number
     */
    void setRetries(int maxretries) { minimizeRetries_ = maxretries; }

    //! \return max displacement when retrying to minimize
    double minimizeDisplacement() const { return minimizeDisplacement_; }

    //! \return whether to force multiple minimizations
    bool forceReminimize() const { return forceReminimize_; }

    //! \brief Set the minimization algorithm
    void setMinimizeAlgorithm(eMinimizeAlgorithm minAlg) { minAlg_ = minAlg; }
    
    //! \return the convergence tolerance (RMS force) for the minimizer
    double forceTolerance() const { return forceToler_; }

    /*! \brief Set the force tolerance
     * \param[in] toler The new value
     */
    void setForceTolerance(double toler) { forceToler_ = toler; }

    //! \return overrelaxation factor for  minimization.
    double overRelax() const { return overRelax_; }
    
    //! \return the max number of iterations for the minimizer. 0 is until convergence.
    int maxIter() const { return maxIter_; }
    
    /*! Set the value of maxIter
     * \param[in] maxIter The new max number of em iterations
     */
    void setMaxIter(int maxIter) { maxIter_ = maxIter; }
};

} //namespace alexandria


#endif //ALEXANDRIA_PARAMHANDLER_H

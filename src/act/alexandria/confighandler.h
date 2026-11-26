/*
 * This source file is part of the Alexandria Chemistry Toolkit.
 *
 * Copyright (C) 2020-2026
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

#include <array>
#include <string>
#include <vector>

#include "act/alexandria/loss_function.h"
#include "act/alexandria/sensitivity.h"
#include "act/constraints/shake.h"
#include "act/properties/velocityhandler.h"
#include "act/utility/time_limit.h"
#include "gromacs/commandline/filenm.h"
#include "gromacs/commandline/pargs.h"
#include "gromacs/utility/real.h"

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
     * \param[in] msghandler Message handler
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
 *! \todo: Should this be under the ga directory???
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
    //! The number of individuals that a middleman computes at once
    int lockstep_ = 1;
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
    //! Loss function to use
    LossFunction lossFunction_ = LossFunction::MSE;
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
    //! Generation limit for the best genome to stay in a local minimum
    int maxLocalMinimumGenerations_ = -1;
    //! Sensitivity analysis after training and to test for a local minimum
    SensitivityMode sensitivity_ = SensitivityMode::OneD;
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
    //! Number of fitness evaluations per individual between writing the fitness
    int nstFitness_ = 1;
    //! Number of generations between writing the gene pool, see -gpout_interval
    int gpoutInterval_ = 5;
    //! Wall time in hours after which the training stops, no limit if zero, see -maxh
    real maxHours_ = 0;

public:

    /*!
     * \brief Add command-line arguments to a vector
     * @param pargs     pointer to arguments vector
     * @param fnms      pointer to filenames
     */
    void add_options(std::vector<t_pargs>  *pargs,
                     std::vector<t_filenm> *fnms) override;

    /*!
     * \brief Check the validity of the provided arguments
     * \param[in] msghandler Message Handler
     * Throw an exception if an invalid combination of arguments is encountered
     */
    void check_pargs(MsgHandler *msghandler) override;

    /* * * * * * * * * * * * * * * * * * * * * *
    * BEGIN: Getters and setters               *
    * * * * * * * * * * * * * * * * * * * * * */

    //! \return the optimizer
    OptimizerAlg optimizer() const { return optAlg_; }

    /*! \brief Set the optimizer algorithm
     * \param[in] optAlg The new algorithm
     */
    void setOptimizerAlg(OptimizerAlg optAlg) { optAlg_ = optAlg; }

    //! \return the LossFunction selection
    LossFunction lossFunction() const { return lossFunction_; }

    //! \return the size of the population
    int popSize() const { return popSize_; }

    //! \return the number of individuals that a middleman computes at once, see -lockstep
    int lockstep() const { return lockstep_; }

    /*! \brief Set the number of individuals that a middleman computes at once
     * \param[in] n The number
     */
    void setLockstep(int n) { lockstep_ = n; }

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

    //! \return the generation limit for the best genome to be in locla minimum
    int maxLocalMinimumGenerations() const { return maxLocalMinimumGenerations_; }

    //! \return the sensitivity analysis to do
    SensitivityMode sensitivity() const { return sensitivity_; }

    /*! \brief Set the sensitivity analysis to do
     * \param[in] mode The sensitivity mode
     */
    void setSensitivity(SensitivityMode mode) { sensitivity_ = mode; }

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

    //! \return the number of generations between writing the gene pool, see -gpout_interval
    int gpoutInterval() const { return gpoutInterval_; }

    //! \return the wall time in hours after which the training stops, zero for no limit, see -maxh
    double maxHours() const { return maxHours_; }

    /*! \brief The limit on the wall time of the training, see -maxh
     * The time counts from the start of the program, such that all
     * processes have the same limit.
     * \return The limit, which never expires if there is none
     */
    TimeLimit timeLimit() const { return TimeLimit(maxHours_); }

    //! \return for CatastrophePenalizer, the fraction of genomes to reinitialize
    real cpPopFrac() const { return cpPopFrac_; }

    //! \return the number of fitness evaluations per individual between writing the fitness
    int nstFitness() const { return nstFitness_; }

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
    //! Random number seed for initializer, crossover, mutator, and selector
    //! \todo: This should be general
    int   seed_              = 0;
    //! Relative step when optimizing
    real  step_              = 0.02;
    //! Whether to use adaptive steps
    bool adaptiveStep_       = false;
    //! Temperature in chi2 units
    real  temperature_       = 5;
    //! Weight temperature after number of training points
    bool  tempWeight_        = false;
    //! Use annealing in the MCMC optimization. Value < 1 means annealing will start at this fraction of the iterations, either per generation or in total
    real anneal_begin_       = 1;
    //! Where to stop annealing as a fraction of the number of iterations, see above
    real anneal_end_         = 1;
    //! Control annealing in HYBRID optimizationa. Value < 1 means annealing will start at this fraction of the generations.
    real anneal_global_begin_ = 1;
    //! Where to stop annealing as a fraction of the number of generations, see manual for more details.
    real anneal_global_end_   = 1;
    //! Evaluate on test set during the MCMC run
    bool  evaluate_testset_  = false;
    //! Checkpointing on or not?
    bool checkPoint_         = false;
    //! Tolerance for shell
    real   shellToler_       = 1e-3;
    //! Maxiter for shell optimization
    int shellMaxIter_        = 25;
    //! Whether to relax the shells starting from their previous positions
    bool shellKeep_          = true;
    //! Whether to compute the energies on the GPU, keeping the molecules there
    bool gpuResident_        = true;
    //! Whether the master process uses the GPU too
    bool gpuMaster_          = true;

public:
    /*!
     * \brief Add command-line arguments to a vector
     * @param pargs     pointer to arguments vector
     * @param fnms      pointer to filenames
     */
    void add_options(std::vector<t_pargs>  *pargs,
                     std::vector<t_filenm> *fnms) override;


    /*!
     * \brief Check the validity of the provided arguments
     * Throw an exception if an invalid combination of arguments is encountered
     * \param[in] msghandler Message handler
     */
    void check_pargs(MsgHandler *msghandler) override;

    //! \return Max # iterations
    int maxIter() const { return maxiter_; }

    /*! \brief Set the number of iterations
     * \param[in] maxiter The max number of iterations
     */
    void setMaxIter(int maxiter) { maxiter_ = maxiter; }

    //! \return temperature for the Monte Carlo simulation
    real temperature() const { return temperature_; }

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

    //! \return the MCMC step size, as fraction of the allowed range of each parameter
    real step() const { return step_; }

    /*! \brief Set the step
     * \param[in] step The new step size
     */
    void setStep(real step) { step_ = step; }

    //! \return whether adaptive steps should be used
    bool adaptiveStep() const { return adaptiveStep_; }

    //! \return whether or not temperature weighting should be considered
    bool temperatureWeighting() const { return tempWeight_; }

    //! \return true if any annealing was requested
    bool annealing() const { return anneal_begin_ < 1; }

    //! \return fraction of generations where global annealing will start
    double globalAnnealBegin() const { return anneal_global_begin_; }

    /*! Set global annealing
     * \param[in] anneal_begin Whether to start global annealing
     * \param[in] anneal_end   the new ending point for simulated annealing
     */
    void setGlobalAnneal(real anneal_begin,
                         real anneal_end)
    {
        anneal_global_begin_ = anneal_begin;
        anneal_global_end_   = anneal_end;
    }

    /*! \brief Set a new value for annealing start
     * \param[in] anneal_begin the new starting point for simulated annealing
     * \param[in] anneal_end   the new ending point for simulated annealing
     */
    void setAnneal(real anneal_begin,
                   real anneal_end)
    { 
        anneal_begin_ = anneal_begin;
        anneal_end_   = anneal_end;
    };

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

    //! \return whether shells are relaxed from their previous positions, see ForceComputer::setKeepShells
    bool shellKeep() const { return shellKeep_; }

    //! \return whether to compute on the GPU, see ForceComputer::setGpuResident
    bool gpuResident() const { return gpuResident_; }

    //! \return whether the master process uses the GPU too, see train_ff -gpu_master
    bool gpuMaster() const { return gpuMaster_; }
};

//! \brief Enum for selecting the coordinate minimizer in simulate
enum class eMinimizeAlgorithm {
    LBFGS, Steep, Newton
};

/*! \brief Map algorithm to string
 * \param[in] e The algorithm
 * \return the name of the algorithm
 */
const std::string &eMinimizeAlgorithmToString(eMinimizeAlgorithm e);

/*! \brief Map string to algorithm
 * \param[in] str The string
 * \return the algorithm corresponding to the string or the default if string does not match any known minimizer
 */
eMinimizeAlgorithm stringToEMinimizeAlgorithm(const std::string &str);

//! \brief Class to handle command line options for MD simulations and energy minimization
class SimulationConfigHandler : ConfigHandler
{
private:
    //! Number of integration steps
    int                nsteps_               = 0;
    //! The integration time step
    real               deltat_               = 0.0002;
    //! Initial simulation temperature
    real               temperature_          = 0;
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
    real               forceToler_           = 1e-6;
    //! Number of retries for minimizing
    int                minimizeRetries_      = 1;
    //! Max random displacement (nm) before re-trying to minize a structure
    real               minimizeDisplacement_ = 0.002;
    //! Whether to force using displacement and reminimize even if converged
    bool               forceReminimize_      = false;
    //! Apply overrelaxation (if > 1) to speed up minimization. Can be dangerous for poor energy functions.
    real               overRelax_            = 1.0;
    //! Maximum number of iterations for the energy minimizer, 0 is until convergence.
    int                maxIter_              = 0;
    //! Maximum number of iterations for relaxing the shells
    int                shellMaxIter_         = 100;
    //! Force constant for flat-bottom position restraints, zero means no restraints
    real               kFBPR_                = 0;
    //! Radius for flat-bottom position restraints (from origin)
    real               r0FBPR_               = 0;
    //! What to constrain in simulations
    ConstraintType     constrain_            = ConstraintType::No;
    //! Maxiter for SHAKE and RATTLE
    int                constraintMaxIter_    = 1000;
    //! Relative tolerance for constraining bond lengths
    real               constraintTolerance_  = 1e-6;
    //! What motion of the center of mass to remove
    COMMRemoval        commRemoval_          = COMMRemoval::Angular;
    //! Number of steps between removing the motion of the center of mass
    int                nstcomm_              = 100;
    //! Number of beads of the ring polymers, 0 for molecular dynamics at constant energy
    int                nbeads_               = 0;
    //! Time constant (ps) of the thermostat of the centroids of the ring polymers
    real               tau_                  = 0.05;
    //! Number of steps left out of the averages, -1 for a tenth of the steps
    int                nequil_               = -1;
    //! Radius (nm) of the wall that keeps the fragments together, zero for no wall
    real               rwall_                = 0;
    //! Force constant (kJ/mol nm^2) of the wall
    real               kwall_                = 1000;
    //! Whether to simulate the fragments as well, for the binding energy
    bool               monomers_             = false;
    //! Whether to compute the thermochemistry of the compound from simulations with beads
    bool               thermo_               = false;
    //! Number of simulations for each of the two integrals of the entropy, zero for no entropy
    int                nquad_                = 4;
    //! Whether the warning that the force tolerance is too strict was given, see toleranceWarning
    mutable bool       toleranceWarned_      = false;
    /*! \brief The choices of the enumerated options, as the command line
     * parser keeps them: the first element points to the choice, nullptr
     * before a command line was parsed. Each handler has its own, such that
     * the choices of one command line do not reach another handler.
     */
    std::array<const char *, 5> minAlgChoice_    = { nullptr, "LBFGS", "Newton", "Steep", nullptr };
    //! See minAlgChoice_, for -constrain
    std::array<const char *, 5> constrainChoice_ = { nullptr, "no", "bonds", "angles", nullptr };
    //! See minAlgChoice_, for -comm
    std::array<const char *, 5> commChoice_      = { nullptr, "angular", "linear", "none", nullptr };
public:
    /*!
     * \brief Add command-line arguments to a vector
     * @param pargs     pointer to arguments vector
     */
    void add_MD_options(std::vector<t_pargs> *pargs);

    /*! \brief Add the command-line arguments for simulations with ring
     * polymers at constant temperature, see RingPolymerSampler
     * \param[inout] pargs Pointer to arguments vector
     */
    void add_bead_options(std::vector<t_pargs> *pargs);

    /*! \brief Add the command-line arguments for energy minimization
     * \param[inout] pargs       Pointer to arguments vector
     * \param[in]    writeShells Whether to add the option for writing shells,
     *                           that is part of the MD options otherwise
     */
    void add_minimize_options(std::vector<t_pargs> *pargs,
                              bool                  writeShells);

    /*!
     * \brief Add command-line arguments to a vector
     * @param pargs     pointer to arguments vector
     * @param fnms      pointer to filenames
     */
    void add_options(std::vector<t_pargs>  *pargs,
                     std::vector<t_filenm> *fnms) override;


    /*!
     * \brief Take the choices of the command line and check the settings
     * \param[in] msghandler Message handler
     * \throws gmx::InvalidInputError for a setting out of range, see validate
     */
    void check_pargs(MsgHandler *msghandler) override;

    /*! \brief Check that the settings make sense
     *
     * For settings made with the setters, by a program that has no
     * command line, and by check_pargs for those of a command line.
     * \throws gmx::InvalidInputError that names the option of the first
     * setting that is out of range
     */
    void validate() const;

    /*! \brief The tolerance for relaxing the shells that goes with the force
     * tolerance of the minimization
     *
     * The shells have to be relaxed more precisely than the atoms are
     * minimized: a tolerance that is not below the force tolerance is
     * replaced by a tenth of it.
     * \param[in] requested  The tolerance asked for, as with -shelltoler
     * \param[in] msghandler Gets a warning if the tolerance is changed, may be nullptr
     * \return the tolerance to use
     */
    double shellTolerance(double requested, MsgHandler *msghandler) const;

    //! \return number of MD steps
    int nsteps() const { return nsteps_; }

    /*! \brief Set what a simulation needs, for a program without a command line
     * \param[in] nsteps      The number of steps
     * \param[in] deltat      The time step in ps
     * \param[in] temperature The temperature in K for the initial velocities, zero for none
     * \param[in] seed        The seed for the velocities, zero to make one
     * \param[in] nstener     The number of steps between energies, zero for none
     * \param[in] nstxout     The number of steps between coordinates, zero for none
     */
    void setSimulation(int nsteps, double deltat, double temperature, int seed, int nstener, int nstxout)
    {
        nsteps_      = nsteps;
        deltat_      = deltat;
        temperature_ = temperature;
        seed_        = seed;
        nstener_     = nstener;
        nstxout_     = nstxout;
    }
    
    //! \return MD integration time step (ps)
    double deltat() const { return deltat_; }

    //! \return the number of beads of the ring polymers, 0 for molecular dynamics at constant energy
    int nbeads() const { return nbeads_; }

    //! \return the time constant (ps) of the thermostat of the centroids of the ring polymers
    double tau() const { return tau_; }

    //! \return the number of steps at the start that are left out of the averages
    int nequil() const { return nequil_ < 0 ? nsteps_/10 : nequil_; }

    /*! \brief Set what a simulation with ring polymers needs, for a program without a command line
     * \param[in] nbeads The number of beads, 0 for molecular dynamics at constant
     *                   energy, 1 for classical dynamics at constant temperature
     * \param[in] tau    The time constant in ps of the thermostat of the centroids
     * \param[in] nequil The number of steps left out of the averages, -1 for a tenth of the steps
     */
    void setBeads(int nbeads, double tau, int nequil)
    {
        nbeads_ = nbeads;
        tau_    = tau;
        nequil_ = nequil;
    }

    //! \return the radius (nm) of the wall that keeps the fragments together, zero for no wall
    double rwall() const { return rwall_; }

    //! \return the force constant (kJ/mol nm^2) of the wall
    double kwall() const { return kwall_; }

    /*! \brief Set the wall that keeps the fragments together in a simulation with ring polymers
     * \param[in] radius        The radius in nm, zero for no wall
     * \param[in] forceConstant The force constant in kJ/mol nm^2
     */
    void setWall(double radius, double forceConstant)
    {
        rwall_ = radius;
        kwall_ = forceConstant;
    }

    //! \return whether to simulate the fragments as well, for the binding energy
    bool monomers() const { return monomers_; }

    /*! \brief Set whether to simulate the fragments as well
     * \param[in] monomers Whether to do so
     */
    void setMonomers(bool monomers) { monomers_ = monomers; }

    //! \return whether to compute the thermochemistry from simulations with beads
    bool thermochemistry() const { return thermo_; }

    //! \return the number of simulations for each of the two integrals of the entropy, zero for no entropy
    int nquad() const { return nquad_; }

    /*! \brief Set whether to compute the thermochemistry from simulations with beads
     * \param[in] thermo Whether to do so
     * \param[in] nquad  The number of simulations for each of the two
     *                   integrals of the entropy, zero for no entropy
     */
    void setThermochemistry(bool thermo, int nquad)
    {
        thermo_ = thermo;
        nquad_  = nquad;
    }
    
    //! \return number of steps between coordinate output
    int nstxout() const { return nstxout_; }
    
    //! \return number of steps between energy output
    int nstener() const { return nstener_; }
    //! \return number of steps between writing velocities
    int nstvout() const { return nstvout_; }

    //! \return initial simulation temperature
    double temperature() const { return temperature_; }

    //! \return force constant (kJ/mol nm) for flat-bottom position restraints, zero means no restraints
    double kFBPR() const { return kFBPR_; }

    //! \return radius (nm) for flat-bottom position restraints (from origin)
    double r0FBPR() const { return r0FBPR_; }

    //! \return user provided random number seed
    int seed() const { return seed_; }

    //! \return whether or not to write shells to trajectory
    bool writeShells() const { return writeShells_; }
    
    //! \return whether or not to minimize the energy with respect to input coordinates
    bool minimize() const { return minimize_; }
    
    //! \return whether or not to do a single point energy calculation
    bool singlePoint() const { return singlePoint_; }

    //! \return what to constrain
    ConstraintType constraintType() const { return constrain_; }

    //! \return whether to constrain anything
    bool constrain() const { return ConstraintType::No != constrain_; }

    //! \return Max number of SHAKE/RATTLE iterations
    int constraintMaxIter() const { return constraintMaxIter_; }

    //! \return The relative tolerance for SHAKE/RATTLE
    double constraintTolerance() const { return constraintTolerance_; }
    //! \return what motion of the center of mass to remove
    COMMRemoval commRemoval() const { return commRemoval_; }
    //! \return number of steps between removing the motion of the center of mass, 0 is only at the start
    int nstcomm() const { return nstcomm_; }

    /*! \brief Set what motion of the center of mass to remove
     * \param[in] commRemoval The motion to remove
     */
    void setCommRemoval(COMMRemoval commRemoval) { commRemoval_ = commRemoval; }

    /*! \brief Set the minimize option
     * \param[in] minimize Whether or not to minimize
     */
    void setMinimize(bool minimize) { minimize_ = minimize; }

    //! \return the minimization algorithm
    eMinimizeAlgorithm minAlg() const { return minAlg_; }

    //! \return maximum number of retries for minimization
    int minimizeRetries() const { return minimizeRetries_; }

    //! \return max displacement when retrying to minimize
    double minimizeDisplacement() const { return minimizeDisplacement_; }

    //! \return whether to force multiple minimizations
    bool forceReminimize() const { return forceReminimize_; }

    /*! \brief Set the minimization algorithm
     * \param[in] minAlg The new algorithm
     */
    void setMinimizeAlgorithm(eMinimizeAlgorithm minAlg) { minAlg_ = minAlg; }
    
    //! \return the convergence tolerance (RMS force) for the minimizer
    double forceTolerance() const { return forceToler_; }

    /*! \brief Set the force tolerance
     * \param[in] toler The new value
     */
    void setForceTolerance(double toler) { forceToler_ = toler; }

    /*! \brief Whether to warn that the force tolerance is too strict
     * This is only the case the first time, such that repeated
     * minimizations with the same configuration warn once.
     * \return true the first time
     */
    bool toleranceWarning() const
    {
        bool first       = !toleranceWarned_;
        toleranceWarned_ = true;
        return first;
    }

    //! \return overrelaxation factor for  minimization.
    double overRelax() const { return overRelax_; }
    
    //! \return the max number of iterations for the minimizer. 0 is until convergence.
    int maxIter() const { return maxIter_; }

    //! \return the maximum number of iterations for relaxing the shells
    int shellMaxIter() const { return shellMaxIter_; }
    
    /*! Set the value of maxIter
     * \param[in] maxIter The new max number of em iterations
     */
    void setMaxIter(int maxIter) { maxIter_ = maxIter; }
};

} //namespace alexandria


#endif //ALEXANDRIA_PARAMHANDLER_H

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

#include "confighandler.h"

#include <cstring>
#include <map>

#include "act/basics/msg_handler.h"
#include "gromacs/math/units.h"
#include "gromacs/utility/arraysize.h"
#include "gromacs/utility/exceptions.h"
#include "gromacs/utility/gmxassert.h"
#include "gromacs/utility/stringutil.h"

namespace alexandria
{

std::map<OptimizerAlg, std::string> opt2str =
    {
        { OptimizerAlg::MCMC,   "MCMC"   },
        { OptimizerAlg::GA,     "GA"     },
        { OptimizerAlg::HYBRID, "HYBRID" }
    };
    
OptimizerAlg stringToOptimizerAlg(const std::string &str)
{
    for(auto &o2s : opt2str)
    {
        if (str == o2s.second)
        {
            return o2s.first;
        }
    }
    GMX_THROW(gmx::InvalidInputError(gmx::formatString("No such Optimizer %s",
                                                       str.c_str()).c_str()));
    // To satisfy the compiler
    return OptimizerAlg::GA;
}

const std::string &optimizerAlgToString(OptimizerAlg opt)
{
    return opt2str[opt];
}

std::map<ProbabilityComputerAlg, std::string> prob2str =
    {
        { ProbabilityComputerAlg::pcRANK,      "RANK"      },
        { ProbabilityComputerAlg::pcFITNESS,   "FITNESS"   },
        { ProbabilityComputerAlg::pcBOLTZMANN, "BOLTZMANN" }
    };
    
ProbabilityComputerAlg stringToProbabilityComputerAlg(const std::string &str)
{
    for(auto &o2s : prob2str)
    {
        if (str == o2s.second)
        {
            return o2s.first;
        }
    }
    GMX_THROW(gmx::InvalidInputError(gmx::formatString("No such ProbabilityComputer %s",
                                                       str.c_str()).c_str()));
    // To satisfy the compiler
    return ProbabilityComputerAlg::pcRANK;
}

const std::string &probabilityComputerAlgToString(ProbabilityComputerAlg opt)
{
    return prob2str[opt];
}

/* * * * * * * * * * * * * * * * * * * * * *
* BEGIN: BayesConfigHandler                *
* * * * * * * * * * * * * * * * * * * * * */

void BayesConfigHandler::add_options(std::vector<t_pargs>             *pargs,
                                     gmx_unused std::vector<t_filenm> *fnms)
{
    std::vector<t_pargs> pa = {
        { "-max_iter", FALSE, etINT, {&maxiter_},
          "Max number of iterations for MCMC optimization. Also applies for the mutation step in the HYBRID optimization algoroithm." },
        { "-temp",    FALSE, etREAL, {&temperature_},
          "'Temperature' for the Monte Carlo simulation." },
        { "-tweight", FALSE, etBOOL, {&tempWeight_},
          "Weight the temperature in the MC/MC algorithm according to the square root of the number of data points. This is in order to get a lower probability of accepting a step in the wrong direction for parameters of which there are few copies." },
        { "-anneal", FALSE, etREAL, {&anneal_},
          "Use annealing in Monte Carlo simulation, starting from this fraction of the simulation. Value should be between 0 and 1." },
        { "-anneal_globally", FALSE, etBOOL, {&annealGlobally_},
          "Whether annealing restart during each mutation round or is global over GA generations" },
        { "-seed",   FALSE, etINT,  {&seed_},
          "Random number seed. If zero, a seed will be generated." },
        { "-step",  FALSE, etREAL, {&step_},
          "Step size for the MCMC parameter optimization. Is used as fraction of the available range per parameter which depends on the parameter type." },
	{ "-checkpoint", FALSE, etBOOL, {&checkPoint_},
		 "Turn on regular checkpointing, i.e. write out intermediate force field files. Expensive." },
        { "-mcmc_evaltest", FALSE, etBOOL, {&evaluate_testset_},
          "Evaluate the parameters on the test set during MCMC." },
        { "-shellTolerance", FALSE, etREAL, {&shellToler_},
          "Tolerance (RMS force) for minimizing shell positions" },
        { "-shellMaxIter", FALSE, etINT, {&shellMaxIter_},
          "Max number of iterations for minimizing shell positions" },
        { "-shellMaxDistance", FALSE, etREAL, {&shellMaxDistance_},
          "Max distance between shell and core (nm)" }
    };
    for (size_t i = 0; i < pa.size(); i++)
    {
        pargs->push_back(pa[i]);
    }
}

void BayesConfigHandler::check_pargs(MsgHandler *)
{
    // Check seed
    GMX_RELEASE_ASSERT(seed_ >= 0, "-seed must be nonnegative.");
    if (seed_ == 0)  // Randomize seed if 0 is provided
    {
      seed_ = ::time(NULL);
    }
    // maxiter_
    GMX_RELEASE_ASSERT(maxiter_ >= 0, "-maxiter must be nonnegative.");
    // step_
    GMX_RELEASE_ASSERT(step_ >= 0, "-step must be nonnegative.");
    // temperature_
    GMX_RELEASE_ASSERT(temperature_ >= 0, "-temp must be nonnegative.");
    // anneal_
    GMX_RELEASE_ASSERT(anneal_ >= 0 && anneal_ <= 1, "-anneal must be in range [0, 1].");
}

real BayesConfigHandler::temperature(int generation,
                                     int max_generations) const
{
    double temp0 = temperature_;
    if (annealGlobally_ && generation > 0)
    {
        temp0 *= (max_generations - generation)/(1.0 * max_generations);
    }
    return temp0;
}

double BayesConfigHandler::computeBeta(int generation,
                                       int max_generations,
                                       int iter)
{
    double temp0 = temperature(generation, max_generations);
    double temp  = temp0;
    if (iter >= maxiter_)
    {
        temp = 1e-6;
    }
    else if (maxiter_ > 0)
    {
        // temp = temperature_*(1.0 - iter/(1.0*maxiter_));
        // Line: temp = m * iter + b
        temp = ( temp0 / ( anneal_ * maxiter_ - maxiter_ ) ) * iter + ( ( temp0 / ( maxiter_ - anneal_ * maxiter_ ) ) * maxiter_ );
    }
    return 1/temp;
}

bool BayesConfigHandler::anneal(int generation,
                                int iter) const
{
    if (annealGlobally_ && generation > 0)
    {
        return true;
    }
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

/* * * * * * * * * * * * * * * * * * * * * *
* END: BayesConfigHandler                  *
* * * * * * * * * * * * * * * * * * * * * */

/* * * * * * * * * * * * * * * * * * * * * *
* BEGIN: GAConfigHandler                   *
* * * * * * * * * * * * * * * * * * * * * */

static const char *optimizerStr[5] = {nullptr, "MCMC", "GA", "HYBRID", nullptr};
static const char *probStr[5] = {nullptr, "RANK", "FITNESS", "BOLTZMANN", nullptr};

void GAConfigHandler::add_options(std::vector<t_pargs>             *pargs,
                                  gmx_unused std::vector<t_filenm> *fnms)
{

    t_pargs pa[] = {
        { "-ga_evaltest", FALSE, etBOOL, {&evaluateTestset_},
          "Evaluate the paramters on the test set during GA." },
        { "-optimizer", FALSE, etENUM, {optimizerStr},
          "Optimization method (see above)." },
        { "-pop_size", FALSE, etINT, {&popSize_},
          "Population size for the GA and HYBRID algorithms (must be even) alternatively number of parallel MCMC calculations that will be performed (can be any integer number > 0 in this case)." },
        { "-n_elites", FALSE, etINT, {&nElites_},
          "Amount of top individuals to be moved, unchanged, to the next generation." },
        { "-random_init", FALSE, etBOOL, {&randomInit_},
          "Initialize the individuals randomly, within the given bounds." },  
        { "-n_crossovers", FALSE, etINT, {&nCrossovers_},
          "Order of the crossover operator. That is, amount of crossover points." },
        { "-sort", FALSE, etBOOL, {&sort_},
          "Whether we sort the genomes in the population based on their fitness. Will be automatically enabled if needed." },
        { "-prob_computer", FALSE, etENUM, {probStr},
          "Probability computation algorithm, see also the two flags below." },
        { "-ga_boltz_temp", FALSE, etREAL, {&boltzTemp_},
          "Initial temperature for Boltzmann probability computing in the genetic algorithm (applies to both the GA and HYBRID optimizers)." },
        { "-ga_boltz_anneal", FALSE, etREAL, {&boltzAnneal_},
          "Starting from this fraction of the total number of generations, the Boltzmann temperature will be linearly lowered until it reaches 0 at the last generation." },
        { "-pr_cross", FALSE, etREAL, {&prCross_},
          "Probability of a crossover event." },
        { "-pr_mut", FALSE, etREAL, {&prMut_},
          "Probability of a mutation to happen in the GA algorithm." },
        { "-percent", FALSE, etREAL, {&percent_},
          "When GA optimizer is selected, -percent denotes the maximum allowed change in a parameter as a fraction of its allowed range." },
        { "-max_generations", FALSE, etINT, {&maxGenerations_},
          "Maximum number of generations before terminating the Genetic Algorithm." },
        { "-max_test_generations", FALSE, etINT, {&maxTestGenerations_},
          "Generation limit for the test fitness to improve in Genetic Algorithm. -1 stands for disabled, otherwise if the test fitness does not improve within this number of generations, the optimization will stop to prevent over-fitting." },
        { "-log_volume", FALSE, etBOOL, {&logVolume_},
          "Compute the (hyper)volume in logarithmic scale." },
        { "-vfp_vol_frac_limit", FALSE, etREAL, {&vfpVolFracLimit_},
          "Limit [0, 1] of the population_volume/total_volume to trigger the VolumeFractionPenalizer. -1 stands for disabled." },
        { "-vfp_pop_frac", FALSE, etREAL, {&vfpPopFrac_},
          "Fraction [0, 1] of the worst genomes to randomize when VolumeFractionPenalizer is triggered." },
        { "-cp_gen_interval", FALSE, etINT, {&cpGenInterval_},
          "Interval (in number of generations) between triggering the CatastrophePenalizer." },
        { "-cp_pop_frac", FALSE, etREAL, {&cpPopFrac_},
          "Fraction [0, 1] of the genomes to randomize when the CatastrophePenalizer is triggered." }
    };
    for (int i = 0; i < asize(pa); i++)
    {
        pargs->push_back(pa[i]);
    }

}

void GAConfigHandler::check_pargs(MsgHandler *msghandler)
{
    optAlg_ = stringToOptimizerAlg(optimizerStr[0]);
    pcAlg_  = stringToProbabilityComputerAlg(probStr[0]);
    
    GMX_RELEASE_ASSERT(popSize_ > 0, "-pop_size must be positive.");
    if (popSize_ % 2 != 0)  // If popSize is odd
    {
        GMX_RELEASE_ASSERT(OptimizerAlg::MCMC == optAlg_, "With odd population sizes, only the MCMC optimizer will work.");
    }

    // If MCMC is selected, change the probability of mutation to 1
    if (optAlg_ == OptimizerAlg::MCMC)
    {
        prMut_ = 1;
    }
    
    GMX_RELEASE_ASSERT(nElites_ >= 0 && nElites_ % 2 == 0, "-n_elites must be nonnegative and even.");
    
    GMX_RELEASE_ASSERT(nCrossovers_ > 0, "-n_crossovers must be nonnegative.");
    
    GMX_RELEASE_ASSERT(boltzTemp_ >= 0, "-boltz_temp must be nonnegative.");
    
    GMX_RELEASE_ASSERT(boltzAnneal_ >= 0 && boltzAnneal_ <= 1, "-boltz_anneal must be in range [0, 1].");

    GMX_RELEASE_ASSERT(prCross_ >= 0 && prCross_ <= 1, "-pr_cross must be in [0,1].");
    
    GMX_RELEASE_ASSERT(prMut_ >= 0 && prMut_ <= 1, "-pr_mut must be in [0,1].");
    
    GMX_RELEASE_ASSERT(percent_ >= 0 && percent_ <= 1, "-percent must be in [0,1].");
    
    GMX_RELEASE_ASSERT(maxGenerations_ > 0, "-max_generations must be positive.");
    
    if (maxTestGenerations_ != -1)
    {
        GMX_RELEASE_ASSERT(maxTestGenerations_ > 0, "-max_test_generations must be positive or -1 (disabled).");
    }

    if (vfpVolFracLimit_ != -1)
    {
        GMX_RELEASE_ASSERT(vfpVolFracLimit_ >= 0 && vfpVolFracLimit_ <= 1,
                           "-vfp_vol_frac_limit must be in [0, 1] when enabled.");
    }
    GMX_RELEASE_ASSERT(vfpPopFrac_ >= 0 && vfpPopFrac_ <= 1,
                       "-vfp_pop_frac must be in [0, 1].");

    if (cpGenInterval_ != -1)
    {
        GMX_RELEASE_ASSERT(cpGenInterval_ > 0,
                           "When enabled (!= -1), -cp_gen_interval must be positive.");
    }
    GMX_RELEASE_ASSERT(cpPopFrac_ >= 0 && cpPopFrac_ <= 1,
                       "-cp_pop_frac must be in [0, 1].");

    // Enable test loss if needed
    if (maxTestGenerations_ > 0)
    {
        if (!evaluateTestset_)
        {
            msghandler->write("Enabling test fitness computation in GA...");
        }
        evaluateTestset_ = true;
    }

    // Enable sorting if needed
    if (nElites_ > 0 || pcAlg_ == ProbabilityComputerAlg::pcRANK || vfpVolFracLimit_ != -1)
    {
        if (!sort_)
        {
            msghandler->write("Enabling genome sorting...");
        }
        sort_ = true;
    }

}

/* * * * * * * * * * * * * * * * * * * * * *
* END: GAConfigHandler                     *
* * * * * * * * * * * * * * * * * * * * * */

/* * * * * * * * * * * * * * * * * * * * * *
* BEGIN: SimulationConfigHandler           *
* * * * * * * * * * * * * * * * * * * * * */

static const char *eminAlgs[5] = {nullptr, "LBFGS", "Newton", "Steep", nullptr};

std::map<eMinimizeAlgorithm, std::string> eMinAlg2String = {
    { eMinimizeAlgorithm::LBFGS,  "LBFGS"  },
    { eMinimizeAlgorithm::Steep,  "Steep"  },
    { eMinimizeAlgorithm::Newton, "Newton" }
};

const std::string &eMinimizeAlgorithmToString(eMinimizeAlgorithm e)
{
    return eMinAlg2String[e];
}

eMinimizeAlgorithm stringToEMinimizeAlgorithm(const std::string &str)
{
    auto eMinA = eMinimizeAlgorithm::Newton;
    for(const auto &m : eMinAlg2String)
    {
        if (str == m.second)
        {
            eMinA = m.first;
        }
    }
    return eMinA;
}

void SimulationConfigHandler::add_options(std::vector<t_pargs>             *pargs,
                                          gmx_unused std::vector<t_filenm> *fnms)
{
    std::vector<t_pargs> extra = {
        { "-sp", FALSE, etBOOL, {&singlePoint_},
          "Do a single point energy calculation and call it a day." },
        { "-minimize", FALSE, etBOOL, {&minimize_},
          "Minimize the energy with respect to input coordinates." },
        { "-minalg",   FALSE, etENUM, {&eminAlgs}, 
          "Algorithm to use for minimization." },
        { "-maxiter",FALSE, etINT,  {&maxIter_},
          "Maximum number of iterations for the energy minimizer, 0 is until convergence (see next option)." },
        { "-maxretries", FALSE, etINT, {&minimizeRetries_},
          "Number of retries for minimizing (LBFGS only)" },
        { "-maxdisplacement", FALSE, etREAL, {&minimizeDisplacement_},
          "Max random displacement (nm) before re-trying to minize a structure (LBFGS only)" },
        { "-forcereminimize", FALSE, etBOOL, {&forceReminimize_},
          "Force using displacement and reminimize even if converged in an attempt to overcome local minima (LBFGS only)" },
        { "-toler",  FALSE, etREAL, {&forceToler_},
          "Convergence tolerance on the mean square atom force for the energy minimizer. If too small, energy minimization may not converge." },
        { "-overrelax", FALSE, etREAL, {&overRelax_},
          "Apply overrelaxation (if > 1) to speed up minimization. Can be dangerous for poor energy functions." }
    };
    for(auto &i : extra)
    {
        pargs->push_back(i);
    }
}

void SimulationConfigHandler::add_MD_options(std::vector<t_pargs> *pargs)
{
    std::vector<t_pargs> MDextra = {
        { "-nsteps", FALSE, etINT, {&nsteps_},
          "Number of integration steps." },
        { "-deltat", FALSE, etREAL, {&deltat_},
          "Integration time step (ps)." },
        { "-ws",     FALSE, etBOOL, {&writeShells_},
          "Write coordinates of shell particles to trajectory and final coordinates as well." },
        { "-temp",    FALSE, etREAL, {&temperature_},
          "Initial simulation temperature (K). If larger than zero, initial velocites will be generated from a Maxwellian distribution at the given temperature." },
        { "-seed",    FALSE, etINT, {&seed_},
          "Random number seed for generating velocities." },
        { "-nstxout", FALSE, etINT, {&nstxout_},
          "Number of steps between writing coordinates." },
        { "-nstvout", FALSE, etINT, {&nstvout_},
          "Number of steps between writing velocities." },
        { "-nstener", FALSE, etINT, {&nstener_},
          "Number of steps between writing energies." }
    };
    for(auto &i : MDextra)
    {
        pargs->push_back(i);
    }
}

void SimulationConfigHandler::check_pargs(MsgHandler *)
{
    GMX_RELEASE_ASSERT(nsteps_ >= 0, "Number of steps must be larger than zero");
    GMX_RELEASE_ASSERT(temperature_ >= 0, "Temperature must be larger than zero");
    GMX_RELEASE_ASSERT(deltat_ > 0, "Integration time step must be larger than zero");
    minAlg_ = stringToEMinimizeAlgorithm(eminAlgs[0]);
}

/* * * * * * * * * * * * * * * * * * * * * *
* END: SimulationConfigHandler             *
* * * * * * * * * * * * * * * * * * * * * */

} //namespace alexandria

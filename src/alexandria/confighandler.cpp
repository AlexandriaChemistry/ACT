/*! \internal \brief
 * Implements part of the alexandria program.
 * \author Mohammad Mehdi Ghahremanpour <mohammad.ghahremanpour@icm.uu.se>
 * \author David van der Spoel <david.vanderspoel@icm.uu.se>
 * \author Julian Ramon Marrades Furquet <julian.marrades@hotmail.es>
 * \author Oskar Tegby <oskar.tegby@it.uu.se>
 */

#include "confighandler.h"

#include <cstring>
#include <map>

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

void BayesConfigHandler::add_pargs(std::vector<t_pargs> *pargs)
{
    t_pargs pa[] = {
        { "-max_iter", FALSE, etINT, {&maxiter_},
          "Max number of iterations for MCMC optimization. Also applies for the mutation step in the HYBRID optimization algoroithm." },
        { "-temp",    FALSE, etREAL, {&temperature_},
          "'Temperature' for the Monte Carlo simulation." },
        { "-tweight", FALSE, etBOOL, {&tempWeight_},
          "Weight the temperature in the MC/MC algorithm according to the square root of the number of data points. This is in order to get a lower probability of accepting a step in the wrong direction for parameters of which there are few copies." },
        { "-anneal", FALSE, etREAL, {&anneal_},
          "Use annealing in Monte Carlo simulation, starting from this fraction of the simulation. Value should be between 0 and 1." },
        { "-seed",   FALSE, etINT,  {&seed_},
          "Random number seed. If zero, a seed will be generated." },
        { "-step",  FALSE, etREAL, {&step_},
          "Step size for the MCMC parameter optimization. Is used as fraction of the available range per parameter which depends on the parameter type." },
        { "-mcmc_evaltest", FALSE, etBOOL, {&evaluate_testset_},
          "Evaluate the parameters on the test set during MCMC." }
    };
    for (int i = 0; i < asize(pa); i++)
    {
        pargs->push_back(pa[i]);
    }
}

void BayesConfigHandler::check_pargs()
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

double BayesConfigHandler::computeBeta(int iter)
{
    double temp = temperature_;
    if (iter >= maxiter_)
    {
        temp = 1e-6;
    }
    else if (maxiter_ > 0)
    {
        // temp = temperature_*(1.0 - iter/(1.0*maxiter_));
        // Line: temp = m * iter + b
        temp = ( temperature_ / ( anneal_ * maxiter_ - maxiter_ ) ) * iter + ( ( temperature_ / ( maxiter_ - anneal_ * maxiter_ ) ) * maxiter_ );
    }
    return 1/(BOLTZ*temp);
}

double BayesConfigHandler::computeBeta(int maxiter, int iter, int ncycle)
{
    double temp = temperature_;
    if (iter >= maxiter_)
    {
        temp = 1e-6;
    }
    else
    {
        temp = (0.5*temperature_)*((exp(-iter/(0.2*(maxiter+1)))) * (1.1 + cos((ncycle*M_PI*iter)/(maxiter+1))));
    }
    return 1/(BOLTZ*temp);
}

bool BayesConfigHandler::anneal(int iter) const
{
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

void GAConfigHandler::add_pargs(std::vector<t_pargs> *pargs)
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
          "Probability computation algorithm" },
        { "-boltz_temp", FALSE, etREAL, {&boltzTemp_},
          "Initial temperature for Boltzmann probability computing." },
        { "-boltz_anneal", FALSE, etREAL, {&boltzAnneal_},
          "Starting from this fraction of generations, the Boltzmann temperature will be linearly lowered until it reaches 0 at the last generation." },
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

void GAConfigHandler::check_pargs()
{
    optAlg_ = stringToOptimizerAlg(optimizerStr[0]);
    pcAlg_  = stringToProbabilityComputerAlg(probStr[0]);
    
    GMX_RELEASE_ASSERT(popSize_ > 0, "-pop_size must be positive.");
    if (popSize_ % 2 != 0)  // If popSize is odd
    {
        GMX_RELEASE_ASSERT(OptimizerAlg::MCMC == optAlg_, "With odd population sizes, only the MCMC optimizer will do.");
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
      GMX_RELEASE_ASSERT(
        vfpVolFracLimit_ >= 0 && vfpVolFracLimit_ <= 1,
        "-vfp_vol_frac_limit must be in [0, 1] when enabled."
      );
    }
    GMX_RELEASE_ASSERT(
      vfpPopFrac_ >= 0 && vfpPopFrac_ <= 1,
      "-vfp_pop_frac must be in [0, 1]."
    );

    if (cpGenInterval_ != -1)
    {
      GMX_RELEASE_ASSERT(
        cpGenInterval_ > 0,
        "When enabled (!= -1), -cp_gen_interval must be positive."
      );
    }
    GMX_RELEASE_ASSERT(
      cpPopFrac_ >= 0 && cpPopFrac_ <= 1,
      "-cp_pop_frac must be in [0, 1]."
    );

    // Enable test loss if needed
    if (maxTestGenerations_ > 0)
    {
      if (!evaluateTestset_)
      {
        printf("Enabling test loss computation in GA...\n");
      }
      evaluateTestset_ = true;
    }

    // Enable sorting if needed
    if (nElites_ > 0 || pcAlg_ == ProbabilityComputerAlg::pcRANK || vfpVolFracLimit_ != -1)
    {
      if (!sort_)
      {
        printf("Enabling sorting...\n");
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

void SimulationConfigHandler::add_pargs(std::vector<t_pargs> *pargs)
{
    std::vector<t_pargs> extra = {
        { "-nsteps", FALSE, etINT, {&nsteps_},
          "Number of integration steps." },
        { "-deltat", FALSE, etREAL, {&deltat_},
          "Integration time step (ps)." },
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
    for(auto &i : extra)
    {
        pargs->push_back(i);
    }
}

void SimulationConfigHandler::check_pargs()
{
    GMX_RELEASE_ASSERT(nsteps_ > 0, "Number of steps must be larger than zero");
    GMX_RELEASE_ASSERT(temperature_ > 0, "Temperature must be larger than zero");
    GMX_RELEASE_ASSERT(deltat_ > 0, "Integration time step must be larger than zero");
}

/* * * * * * * * * * * * * * * * * * * * * *
* END: SimulationConfigHandler             *
* * * * * * * * * * * * * * * * * * * * * */

} //namespace alexandria

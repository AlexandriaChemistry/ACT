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
        { "-maxiter", FALSE, etINT, {&maxiter_},
          "Max number of iterations for MCMC optimization." },
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
        { "-bEvaluate_testset", FALSE, etBOOL, {&evaluate_testset_},
          "Evaluate the MCMC energy on the test set. Only used in pure MCMC optmization." }
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

void GAConfigHandler::add_pargs(std::vector<t_pargs> *pargs)
{

    t_pargs pa[] = {
        { "-optimizer", FALSE, etENUM, {optimizerStr},
          "Optimization method" },
        { "-popSize", FALSE, etINT, {&popSize_},
          "Population size." },
        { "-nElites", FALSE, etINT, {&nElites_},
          "Amount of top individuals to be moved, unchanged, to the next generation." },
        { "-randomInit", FALSE, etBOOL, {&randomInit_},
          "Initialize the individuals randomly, within the given bounds." },  
        { "-nCrossovers", FALSE, etINT, {&nCrossovers_},
          "Order of the crossover operator. That is, amount of crossover points." },
        { "-sort", FALSE, etBOOL, {&sort_},
          "Whether we sort the genomes in the population based on their fitness." },
        { "-probComputer", FALSE, etENUM, {probComputer_},
          "Probability computation algorithm" },
        { "-boltzTemp", FALSE, etREAL, {&boltzTemp_},
          "Initial temperature for Boltzmann probability computing." },
        { "-prCross", FALSE, etREAL, {&prCross_},
          "Probability of crossover." },
        { "-prMut", FALSE, etREAL, {&prMut_},
          "Probability of mutation" },
        { "-percent", FALSE, etREAL, {&percent_},
          "When GA optimizer is selected, -percent denotes the maximum allowed change in a parameter as a fraction of its allowed range." },
        { "-maxGenerations", FALSE, etINT, {&maxGenerations_},
          "Generation limit for Genetic Algorithm." },
        { "-maxTestGenerations", FALSE, etINT, {&maxTestGenerations_},
          "Generation limit for the test fitness to improve in Genetic Algorithm. -1 stands for disabled." },
        { "-vfpVolFracLimit", FALSE, etREAL, {&vfpVolFracLimit_},
          "Limit [0, 1] of the populationVolume/totalVolume to trigger the VolumeFractionPenalizer. -1 stands for disabled." },
        { "-vfpPopFrac", FALSE, etREAL, {&vfpPopFrac_},
          "Fraction [0, 1] of the worst genomes to randomize when VolumeFractionPenalizer is triggered." }
    };
    for (int i = 0; i < asize(pa); i++)
    {
        pargs->push_back(pa[i]);
    }

}

void GAConfigHandler::check_pargs()
{
    alg_ = stringToOptimizerAlg(optimizerStr[0]);
    
    GMX_RELEASE_ASSERT(popSize_ > 0, "-popSize must be positive.");
    if (popSize_ % 2 != 0)  // If popSize is odd
    {
        GMX_RELEASE_ASSERT(OptimizerAlg::MCMC == alg_, "With odd population sizes, only the MCMC optimizer will do.");
    }

    // If MCMC is selected, change the probability of mutation to 1
    if (alg_ == OptimizerAlg::MCMC)
    {
      prMut_ = 1;
    }
    
    GMX_RELEASE_ASSERT(nElites_ >= 0 && nElites_ % 2 == 0, "-nElites must be nonnegative and even.");
    if (nElites_ > 0)  // Make sure a sorter has been selected
    {
        GMX_RELEASE_ASSERT(sort_ == true,
                           "When -nElites > 0, -sort should be used.");
    }
    
    GMX_RELEASE_ASSERT(nCrossovers_ > 0, "-nCrossovers must be nonnegative.");
    
    if (strcmp(probComputer_[0], "RANK") == 0)  // If rank-based probability is requested
    {
        GMX_RELEASE_ASSERT(sort_ == true, "You must enable -sort if you want rank-based probability computing.");
    }
    
    GMX_RELEASE_ASSERT(boltzTemp_ >= 0, "-boltzTemp must be nonnegative.");
    
    GMX_RELEASE_ASSERT(prCross_ >= 0 && prCross_ <= 1, "-prCross must be in [0,1].");
    
    GMX_RELEASE_ASSERT(prMut_ >= 0 && prMut_ <= 1, "-prMut must be in [0,1].");
    
    GMX_RELEASE_ASSERT(percent_ >= 0 && percent_ <= 1, "-percent must be in [0,1].");
    
    GMX_RELEASE_ASSERT(maxGenerations_ > 0, "-maxGenerations must be positive.");
    
    if (maxTestGenerations_ != -1)
    {
      GMX_RELEASE_ASSERT(maxGenerations_ > 0, "-maxTestGenerations must be positive or -1 (disabled).");
    }

    if (vfpVolFracLimit_ != -1)
    {
      GMX_RELEASE_ASSERT(
        vfpVolFracLimit_ >= 0 && vfpVolFracLimit_ <= 1,
        "-vfpVolFracLimit must be in [0, 1] when enabled."
      );
    }
    GMX_RELEASE_ASSERT(
      vfpPopFrac_ >= 0 && vfpPopFrac_ <= 1,
      "-vfpPopFrac must be in [0, 1]."
    );

}

/* * * * * * * * * * * * * * * * * * * * * *
* END: GAConfigHandler                     *
* * * * * * * * * * * * * * * * * * * * * */

} //namespace alexandria

#include "confighandler.h"

#include "gromacs/utility/arraysize.h"
#include "gromacs/math/units.h"
#include "gromacs/utility/gmxassert.h"


namespace alexandria
{


void BayesConfigHandler::add_pargs(std::vector<t_pargs> *pargs)
{
    t_pargs pa[] = {
        { "-maxiter", FALSE, etINT, {&maxiter_},
          "Max number of iterations for optimization. Default 100." },
        { "-temp",    FALSE, etREAL, {&temperature_},
          "'Temperature' for the Monte Carlo simulation. Default 5." },
        { "-tweight", FALSE, etBOOL, {&tempWeight_},
          "Weight the temperature in the MC/MC algorithm according to the square root of the number of data points. This is in order to get a lower probability of accepting a step in the wrong direction for parameters of which there are few copies. Default false." },
        { "-anneal", FALSE, etREAL, {&anneal_},
          "Use annealing in Monte Carlo simulation, starting from this fraction of the simulation. Value should be between 0 and 1. Default 1." },
        { "-seed",   FALSE, etINT,  {&seed_},
          "Random number seed. If zero, a seed will be generated. Default -1." },
        { "-step",  FALSE, etREAL, {&step_},
          "Step size for the parameter optimization. Is used as fraction of the available range per parameter which depends on the parameter type. Default 0.02." },
        { "-v",     FALSE, etBOOL, {&verbose_},
          "Flush output immediately rather than letting the OS buffer it. Don't use for production simulations. Default false." }
    };
    for (int i = 0; i < asize(pa); i++)
    {
        pargs->push_back(pa[i]);
    }
}

void BayesConfigHandler::check_pargs()
{
    // maxiter_
    GMX_RELEASE_ASSERT(maxiter_ < 0, "-maxiter must be nonnegative.");
    // step_
    GMX_RELEASE_ASSERT(step_ < 0, "-step must be nonnegative.");
    // temperature_
    GMX_RELEASE_ASSERT(temperature_ < 0, "-temp must be nonnegative.");
    // anneal_
    GMX_RELEASE_ASSERT(anneal_ < 0 || anneal_ > 1, "-anneal must be in range [0, 1].");
}

void BayesConfigHandler::setOutputFiles(const char                     *xvgconv,
                                       const std::vector<std::string> &paramClass,
                                       const char                     *xvgepot,
                                       const gmx_output_env_t         *oenv)
{
    xvgconv_.assign(xvgconv);
    paramClass_ = paramClass;
    xvgepot_.assign(xvgepot);
    oenv_       = oenv;
}

double BayesConfigHandler::computeBeta(int iter)
{
    double temp = temperature_;
    if (iter >= maxiter_)
    {
        temp = 1e-6;
    }
    else
    {
        temp = temperature_*(1.0 - iter/(maxiter_ + 1.0));
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


}
/*! \internal \brief
 * Implements part of the alexandria program.
 * \author Mohammad Mehdi Ghahremanpour <mohammad.ghahremanpour@icm.uu.se>
 * \author David van der Spoel <david.vanderspoel@icm.uu.se>
 * \author Oskar Tegby <oskar.tegby@it.uu.se>
 * \author Julian Ramon Marrades Furquet <julian.marrades@hotmail.es>
 */


#include "acminitializer.h"

#include "act/ga//Individual.h"
#include "acmindividual.h"


namespace alexandria
{


/* * * * * * * * * * * * * * * * * * * *
* BEGIN: ACMInitializer                *
* * * * * * * * * * * * * * * * * * * */
ACMInitializer::ACMInitializer(StaticIndividualInfo *sii,
                               bool                  randInit,
                               int                   seed)
    : gen_(rd_()), dis_(std::uniform_real_distribution<double>(0.0, 1.0)),
      sii_(sii), randInit_(randInit)
{
    gen_.seed(seed);
}

ga::Individual *ACMInitializer::initialize()
{
    int id   = sii_->commRec()->middleManOrdinal();
    auto ind = new ACMIndividual(id, sii_);
    if (randInit_)
    // Insert random value in range
    {
        for (size_t i = 0; i < sii_->nParam(); i++)
        {
            ind->addParam(dis_(gen_)*(sii_->upperBoundAtIndex(i) - sii_->lowerBoundAtIndex(i))
                          + sii_->lowerBoundAtIndex(i));
        }
    }
    else  // Insert default values in StaticIndividualInfo
    {
        for (const double val : sii_->defaultParam())
        {
            ind->addParam(val);
        }
    }
    ind->setBestGenome(ind->genome());
    return ind;
}

/* * * * * * * * * * * * * * * * * * * *
* END: ACMInitializer                  *
* * * * * * * * * * * * * * * * * * * */


} // namespace alexandria

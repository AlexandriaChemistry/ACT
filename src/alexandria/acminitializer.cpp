/*! \internal \brief
 * Implements part of the alexandria program.
 * \author Mohammad Mehdi Ghahremanpour <mohammad.ghahremanpour@icm.uu.se>
 * \author David van der Spoel <david.vanderspoel@icm.uu.se>
 * \author Oskar Tegby <oskar.tegby@it.uu.se>
 * \author Julian Ramon Marrades Furquet <julian.marrades@hotmail.es>
 */


#include "acminitializer.h"

#include "ga/Individual.h"
#include "acmindividual.h"


namespace alexandria
{


/* * * * * * * * * * * * * * * * * * * *
* BEGIN: ACMInitializer                *
* * * * * * * * * * * * * * * * * * * */

void ACMInitializer::initialize(ga::Individual **ind)
{
    nCreated_++;
    ACMIndividual *tmpInd = new ACMIndividual(nCreated_, sii_, outputFile_);
    (*ind) = tmpInd;
    if (randInit_)  // Insert random value in range
    {
        for (size_t i = 0; i < sii_->nParam(); i++)
        {
            tmpInd->addParam(dis(gen)*(sii_->upperBoundAtIndex(i) - sii_->lowerBoundAtIndex(i))
                             + sii_->lowerBoundAtIndex(i));
        }
    }
    else  // Insert default values in SharedIndividualInfo
    {
        for (const double val : sii_->defaultParam()) tmpInd->addParam(val);
    }
}

/* * * * * * * * * * * * * * * * * * * *
* END: ACMInitializer                  *
* * * * * * * * * * * * * * * * * * * */


} // namespace alexandria
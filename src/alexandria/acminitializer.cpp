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
    Poldata* pd = sii_->poldata();
    // FIXME: We could store the file-in value of parameters in SharedIndividualInfo and avoid doing all this again here.
    // Then, we could just either use that value from sii_ or a randomly generated one. We also remove dependency in mindata_
    for(auto &optIndex : sii_->optIndex())
    {
        auto                iType = optIndex.iType();
        ForceFieldParameter p;
        if (iType == InteractionType::CHARGE)
        {
            if (pd->hasParticleType(optIndex.particleType()))
            {
                p = pd->findParticleType(optIndex.particleType())->parameterConst(optIndex.parameterType());
            }
        }
        else if (pd->interactionPresent(iType))
        {
            p = pd->findForcesConst(iType).findParameterTypeConst(optIndex.id(), optIndex.parameterType());
        }
        if (p.ntrain() >= mindata_)
        {
            if (randInit_)
            {
                tmpInd->addParam(dis(gen)*(p.maximum()-p.minimum()) + p.minimum());
            }
            else
            {
                tmpInd->addParam(p.value());
            }
        }
    }
}

/* * * * * * * * * * * * * * * * * * * *
* END: ACMInitializer                  *
* * * * * * * * * * * * * * * * * * * */


} // namespace alexandria
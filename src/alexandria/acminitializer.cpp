/*! \internal \brief
 * Implements part of the alexandria program.
 * \author Mohammad Mehdi Ghahremanpour <mohammad.ghahremanpour@icm.uu.se>
 * \author David van der Spoel <david.vanderspoel@icm.uu.se>
 * \author Oskar Tegby <oskar.tegby@it.uu.se>
 * \author Julian Ramon Marrades Furquet <julianramon.marradesfurquet.8049@student.uu.se>
 */


#include "acminitializer.h"

#include "ga/Individual.h"
#include "acmindividual.h"


namespace alexandria
{


/* * * * * * * * * * * * * * * * * * * *
* BEGIN: ACMInitializer                *
* * * * * * * * * * * * * * * * * * * */

void ACMInitializer::initialize(ga::Individual **individual)
{
    nCreated_++;
    ACMIndividual *tmpInd = new ACMIndividual(nCreated_, sii_, outputFile_);
    (*individual) = tmpInd;
    Poldata* pd = sii_->poldata();
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
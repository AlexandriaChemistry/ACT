#include "acminitializer.h"
#include "acmindividual.h"
#include "tune_eem.h"

#include "ga/Individual.h"

namespace alexandria
{


/* * * * * * * * * * * * * * * * * * * *
* BEGIN: ACMInitializer                *
* * * * * * * * * * * * * * * * * * * */

void ACMInitializer::initialize(ga::Individual *individual)
{
    ACMIndividual *tmpInd = static_cast<ACMIndividual*>(individual);
    Poldata* pd = tmpInd->poldata();
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
            // FIXME: This information has been moved to the individual.
            // Bayes::addParam(optIndex.name(),
            //                 p.value(), p.mutability(),
            //                 p.minimum(), p.maximum(),
            //                 p.ntrain(), bRandom);
        }
    }
}

/* * * * * * * * * * * * * * * * * * * *
* END: ACMInitializer                  *
* * * * * * * * * * * * * * * * * * * */


} // namespace alexandria
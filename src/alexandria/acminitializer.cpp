#include "acminitializer.h"
#include "tune_eem.h"

namespace alexandria
{


/* * * * * * * * * * * * * * * * * * * *
* BEGIN: ACMInitializer                *
* * * * * * * * * * * * * * * * * * * */

void ACMInitializer::initialize(Individual *individual)
{
    for(auto &optIndex : optIndex_)
    {
        auto                iType = optIndex.iType();
        ForceFieldParameter p;
        if (iType == InteractionType::CHARGE)
        {
            if (poldata()->hasParticleType(optIndex.particleType()))
            {
                p = poldata()->findParticleType(optIndex.particleType())->parameterConst(optIndex.parameterType());
            }
        }
        else if (poldata()->interactionPresent(iType))
        {
            p = poldata()->findForcesConst(iType).findParameterTypeConst(optIndex.id(), optIndex.parameterType());
        }
        if (p.ntrain() >= mindata())
        {
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
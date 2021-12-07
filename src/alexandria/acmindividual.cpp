#include "acmindividual.h"

#include "poldata_xml.h"


namespace alexandria
{

/* * * * * * * * * * * * * * * * * * * * * *
* BEGIN: Output stuff                      *
* * * * * * * * * * * * * * * * * * * * * */

void ACMIndividual::saveState()
{
    writePoldata(outputFile_, pd_, false);
}

/* * * * * * * * * * * * * * * * * * * * * *
* END: Output stuff                        *
* * * * * * * * * * * * * * * * * * * * * */

/* * * * * * * * * * * * * * * * * * * * * *
* BEGIN: Poldata stuff                     *
* * * * * * * * * * * * * * * * * * * * * */

void ACMIndividual::toPoldata(const std::vector<bool> &changed)
{

    size_t n = 0;
    if (psigma_.empty())
    {
        psigma_.resize(param_.size(), 0);
    }
    printParameters(debug);
    for (const auto &optIndex : sii_.optIndex())
    {
        if (changed[n])
        {
            auto                 iType = optIndex.iType();
            ForceFieldParameter *p = nullptr;
            if (iType != InteractionType::CHARGE)
            {
                p = pd_->findForces(iType)->findParameterType(optIndex.id(), optIndex.parameterType());
            }
            else if (pd_->hasParticleType(optIndex.particleType()))
            {
                p = pd_->findParticleType(optIndex.particleType())->parameter(optIndex.parameterType());
            }
            GMX_RELEASE_ASSERT(p, gmx::formatString("Could not find parameter %s", optIndex.id().id().c_str()).c_str());
            if (p)
            {
                p->setValue(param_[n]);
                p->setUncertainty(psigma_[n]);
            }
        }
        n++;
    }
    GMX_RELEASE_ASSERT(n == changed.size(),
                       gmx::formatString("n = %zu changed.size() = %zu",
                                         n, changed.size()).c_str());

}

/* * * * * * * * * * * * * * * * * * * * * *
* END: Poldata stuff                       *
* * * * * * * * * * * * * * * * * * * * * */

} //namespace alexandria
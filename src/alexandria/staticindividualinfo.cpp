/*! \internal \brief
 * Implements part of the alexandria program.
 * \author Mohammad Mehdi Ghahremanpour <mohammad.ghahremanpour@icm.uu.se>
 * \author David van der Spoel <david.vanderspoel@icm.uu.se>
 * \author Julian Ramon Marrades Furquet <julian.marrades@hotmail.es>
 */


#include "staticindividualinfo.h"

#include "ga/Genome.h"
#include "memory_check.h"
#include "poldata_xml.h"

namespace alexandria
{

StaticIndividualInfo::StaticIndividualInfo(const CommunicationRecord *cr) : cr_(cr)
{
    fillFittingTargets();
    id_ = cr_->middleManOrdinal();
    if (id_ >= 0)
    {
        prefix_ = gmx::formatString("ind%d/", id_);
    }
}

void StaticIndividualInfo::setOutputFile(const std::string &outputFile)
{
    outputFile_ = prefix_ + "ind" + std::to_string(id_) + "-" + outputFile;
}

/* * * * * * * * * * * * * * * * * * * * * *
* BEGIN: Poldata stuff                     *
* * * * * * * * * * * * * * * * * * * * * */

void StaticIndividualInfo::fillPoldata(      FILE *fp,
                                       const char *pd_fn)
{
    if (!cr_->isHelper())
    {
        GMX_RELEASE_ASSERT(nullptr != pd_fn, "Give me a poldata file name");
        if (debug)
        {
            fprintf(debug, "On node %d, will try to read %s\n", cr_->rank(), pd_fn);
        }
        try
        {
            alexandria::readPoldata(pd_fn, &pd_);
        }
        GMX_CATCH_ALL_AND_EXIT_WITH_FATAL_ERROR;
        print_memory_usage(debug);
    }
    /* Broadcasting Force Field Data from Master to Helper nodes */
    if (cr_->isParallel())
    {
        pd_.sendToHelpers(cr_);
    }
    if (nullptr != fp)
    {
        fprintf(fp, "There are %d atom types in the input file %s.\n\n",
                static_cast<int>(pd_.getNatypes()), pd_fn);
    }
}

void StaticIndividualInfo::updatePoldata(const std::vector<bool> &changed,
                                         const ga::Genome        *genome)
{
    size_t n = 0;
    for (const auto &optIndex : optIndex())
    {
        if (changed[n])
        {
            auto                 iType = optIndex.iType();
            ForceFieldParameter *p  = nullptr;
            if (iType != InteractionType::CHARGE)
            {
                p = pd_.findForces(iType)->findParameterType(optIndex.id(), optIndex.parameterType());
            }
            else if (pd_.hasParticleType(optIndex.particleType()))
            {
                p = pd_.findParticleType(optIndex.particleType())->parameter(optIndex.parameterType());
            }
            GMX_RELEASE_ASSERT(p, gmx::formatString("Could not find parameter %s", optIndex.id().id().c_str()).c_str());
            if (p)
            {
                p->setValue(genome->base(n));
                // TODO fix the uncertainty
                // p->setUncertainty(psigma_[n]);
            }
        }
        n++;
    }
    GMX_RELEASE_ASSERT(n == changed.size(),
                       gmx::formatString("n = %zu changed.size() = %zu",
                                         n, changed.size()).c_str());
}

void StaticIndividualInfo::saveState(bool updateCheckSum)
{
    pd_.updateTimeStamp();
    if (updateCheckSum)
    {
        pd_.updateCheckSum();
    }
    writePoldata(outputFile_, &pd_, false);
}

/* * * * * * * * * * * * * * * * * * * * * *
* END: Poldata stuff                       *
* * * * * * * * * * * * * * * * * * * * * */

/* * * * * * * * * * * * * * * * * * * * * *
* BEGIN: FittingTarget stuff               *
* * * * * * * * * * * * * * * * * * * * * */

void StaticIndividualInfo::sumChiSquared(bool             parallel,
                                         iMolSelect       ims)
{
    // Now sum over processors, except the master!
    if (cr_->isParallel() && parallel)
    {
        for (auto &ft : targets_[ims])
        {
            auto chi2 = ft.second.chiSquared();
            cr_->sumd_helpers(1, &chi2);
            ft.second.setChiSquared(chi2);
            auto ndp = ft.second.numberOfDatapoints();
            cr_->sumi_helpers(1, &ndp);
            ft.second.setNumberOfDatapoints(ndp);
        }
    }
    auto myft = targets_.find(ims);
    if (myft !=  targets_.end())
    {
        auto etot = myft->second.find(eRMS::TOT);
        GMX_RELEASE_ASSERT(etot != myft->second.end(), "Cannot find etot");
        etot->second.reset();
        auto tc = targets_.find(ims);
        if (tc != targets_.end())
        {
            for (const auto &ft : tc->second)
            {
                if (ft.first != eRMS::TOT)
                { 
                    etot->second.increase(1.0, ft.second.chiSquaredWeighted());
                }
            }
        }
        // Weighting is already included.
        etot->second.setNumberOfDatapoints(1);
    }
}

void StaticIndividualInfo::resetChiSquared(iMolSelect ims)
{
    auto fts = targets_.find(ims);
    if (fts != targets_.end())
    {
        for (auto &ft : fts->second)
        {
            ft.second.reset();
        }
    }
}

void StaticIndividualInfo::fillFittingTargets()
{
    for (const auto &ims : iMolSelectNames()) 
    {
        std::map<eRMS, FittingTarget> ft;
        for ( auto &rms : geteRMSNames() )
        {
            FittingTarget    fft(rms.first, ims.first);
            ft.insert(std::pair<eRMS, FittingTarget>(rms.first, std::move(fft)));
        }
        targets_.insert(std::pair<iMolSelect, std::map<eRMS, FittingTarget>>(ims.first, std::move(ft)));
        auto etot = target(ims.first, eRMS::TOT);
        GMX_RELEASE_ASSERT(etot != nullptr, "Could not find etot");
        etot->setWeight(1);
    }
}

void StaticIndividualInfo::propagateWeightFittingTargets()
{
    for(auto &ims : iMolSelectNames())
    {
        if (ims.first != iMolSelect::Train)
        {
            auto ft = fittingTargets(ims.first);
            if (ft != nullptr)
            {
                for(auto &rms : geteRMSNames())
                {
                    auto w = target(iMolSelect::Train, rms.first)->weight();
                    target(ims.first, rms.first)->setWeight(w);
                }
            }
        }
    }
}

/* * * * * * * * * * * * * * * * * * * * * *
* END: FittingTarget stuff                 *
* * * * * * * * * * * * * * * * * * * * * */

/* * * * * * * * * * * * * * * * * * * * * *
* BEGIN: Weighted temperature stuff        *
* * * * * * * * * * * * * * * * * * * * * */

void StaticIndividualInfo::computeWeightedTemperature(const bool tempWeight)
{
    if (tempWeight)
    {
        for (size_t j = 0; j < paramNames_.size(); j++)
        {
            GMX_RELEASE_ASSERT(ntrain_[j] > 0, "ntrain should be > 0 for all parameters");
            // TODO: Maybe a fast inverse square root here?
            weightedTemperature_.push_back(std::sqrt(1.0/ntrain_[j]));
        }
    }
    else
    {  
        weightedTemperature_.resize(paramNames_.size(), 1.0);
    }
}

/* * * * * * * * * * * * * * * * * * * * * *
* END: Weighted temperature stuff        *
* * * * * * * * * * * * * * * * * * * * * */

/* * * * * * * * * * * * * * * * * * * * * *
* BEGIN: OptimizationIndex stuff           *
* * * * * * * * * * * * * * * * * * * * * */

void StaticIndividualInfo::generateOptimizationIndex(FILE   *fp,
                                                     MolGen *mg)
{
    for(auto &fs : pd_.forcesConst())
    {
        if (mg->optimize(fs.first))
        {
            for(auto &fpl : fs.second.parametersConst())
            {
                for(auto &param : fpl.second)
                {
                    if (mg->fit(param.first))
                    {
                        if (param.second.isMutable() && param.second.ntrain() >= mg->mindata())
                        {
                            optIndex_.push_back(OptimizationIndex(fs.first, fpl.first, param.first));
                        }
                    }
                }
            }
        }
    }
    for(auto &pt : pd_.particleTypesConst())
    {
        for(auto &p : pt.parametersConst())
        {
            if (mg->fit(p.first) && p.second.ntrain() > 0)
            {
                optIndex_.push_back(OptimizationIndex(pt.id().id(), p.first));
            }
        }
    }
    if (fp)
    {
        fprintf(fp, "There are %zu variables to optimize.\n", optIndex_.size());
        for(auto &i : optIndex_)
        {
            fprintf(fp, "Will optimize %s\n", i.name().c_str());
        }
    }
}

/* * * * * * * * * * * * * * * * * * * * * *
* END: OptimizationIndex stuff             *
* * * * * * * * * * * * * * * * * * * * * */

/* * * * * * * * * * * * * * * * * * * * * *
* BEGIN: Vector stuff                      *
* * * * * * * * * * * * * * * * * * * * * */

void StaticIndividualInfo::fillVectors(const int mindata)
{

    for(auto &optIndex : optIndex_)
    {
        auto                iType = optIndex.iType();
        ForceFieldParameter p;
        if (iType == InteractionType::CHARGE)
        {
            if (pd_.hasParticleType(optIndex.particleType()))
            {
                p = pd_.findParticleType(optIndex.particleType())->parameterConst(optIndex.parameterType());
            }
        }
        else if (pd_.interactionPresent(iType))
        {
            p = pd_.findForcesConst(iType).findParameterTypeConst(optIndex.id(), optIndex.parameterType());
        }
        if (p.ntrain() >= mindata)
        {
            defaultParam_.push_back(p.value());
            paramNames_.push_back(optIndex.name());
            mutability_.push_back(p.mutability());
            lowerBound_.push_back(p.minimum());
            upperBound_.push_back(p.maximum());
            ntrain_.push_back(p.ntrain());
        }
    }

}

/* * * * * * * * * * * * * * * * * * * * * *
* END: Vector stuff                        *
* * * * * * * * * * * * * * * * * * * * * */

/* * * * * * * * * * * * * * * * * * * * * *
* BEGIN: ParamClassIndex stuff             *
* * * * * * * * * * * * * * * * * * * * * */

void StaticIndividualInfo::assignParamClassIndex()
{
    size_t notFound = ~0;
    paramClassIndex_.resize(paramNames_.size(), notFound);

    for (size_t i = 0; i < paramClass_.size(); i++)
    {
        for (size_t j = 0; j < paramNames_.size(); j++)
        {
            if (paramNames_[j].find(paramClass_[i]) != std::string::npos)
            {
                paramClassIndex_[j] = i;
            }
        }
    }

    // Now check for params which were not assigned a class and give them class "Other"
    bool restClass = false;
    for (size_t i = 0; i < paramClassIndex_.size(); i++)
    {
        if (paramClassIndex_[i] == notFound)
        {
            if (!restClass)  // If <i> is the first parameter without a class
            {
                // Append "Other" to the list of classes
                paramClass_.push_back("Other");
                restClass = true;
            }
            // Give class "Other" to parameter <i>
            paramClassIndex_[i] = paramClass_.size() - 1;
        }
    }
}

/* * * * * * * * * * * * * * * * * * * * * *
* END: ParamClassIndex stuff               *
* * * * * * * * * * * * * * * * * * * * * */

/* * * * * * * * * * * * * * * * * * * * * *
* BEGIN: File stuff                        *
* * * * * * * * * * * * * * * * * * * * * */

void StaticIndividualInfo::setOutputFiles(const char                     *xvgconv,
                                          const std::vector<std::string> &paramClass,
                                          const char                     *xvgepot)
{
    xvgconv_.assign(xvgconv);
    paramClass_ = paramClass;
    xvgepot_.assign(xvgepot);
}

/* * * * * * * * * * * * * * * * * * * * * *
* END: File stuff                          *
* * * * * * * * * * * * * * * * * * * * * */


} //namespace alexandria

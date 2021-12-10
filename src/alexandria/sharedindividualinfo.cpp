#include "sharedindividualinfo.h"

#include "poldata_xml.h"
#include "memory_check.h"


namespace alexandria
{


/* * * * * * * * * * * * * * * * * * * * * *
* BEGIN: Poldata stuff                     *
* * * * * * * * * * * * * * * * * * * * * */

void SharedIndividualInfo::fillPoldata(      FILE *fp,
                                       const char *pd_fn)
{
    if (MASTER(cr_))
    {
        GMX_RELEASE_ASSERT(nullptr != pd_fn, "Give me a poldata file name");
        try
        {
            alexandria::readPoldata(pd_fn, &pd_);
        }
        GMX_CATCH_ALL_AND_EXIT_WITH_FATAL_ERROR;
        print_memory_usage(debug);
    }
    /* Broadcasting Force Field Data from Master to Helper nodes */
    if (PAR(cr_))
    {
        pd_.broadcast(cr_);
    }
    if (nullptr != fp)
    {
        fprintf(fp, "There are %d atom types in the input file %s.\n\n",
                static_cast<int>(pd_.getNatypes()), pd_fn);
    }
}

/* * * * * * * * * * * * * * * * * * * * * *
* END: Poldata stuff                       *
* * * * * * * * * * * * * * * * * * * * * */

/* * * * * * * * * * * * * * * * * * * * * *
* BEGIN: FittingTarget stuff               *
* * * * * * * * * * * * * * * * * * * * * */

void SharedIndividualInfo::fillFittingTargets()
{
    for (const auto &ims : iMolSelectNames()) 
    {
        std::map<eRMS, FittingTarget> ft;
        for ( auto &rms : ermsNames )
        {
            FittingTarget    fft(rms.first, ims.first);
            ft.insert(std::pair<eRMS, FittingTarget>(rms.first, std::move(fft)));
        }
        targets_.insert(std::pair<iMolSelect, RmsFittingTarget>(ims.first, std::move(ft)));
        auto etot = target(ims.first, eRMS::TOT);
        GMX_RELEASE_ASSERT(etot != nullptr, "Could not find etot");
        etot->setWeight(1);
    }
}

void SharedIndividualInfo::propagateWeightFittingTargets()
{
    for(auto &ims : iMolSelectNames())
    {
        if (ims.first != iMolSelect::Train)
        {
            auto ft = fittingTargets(ims.first);
            if (ft != nullptr)
            {
                for(auto &rms : ermsNames)
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

void SharedIndividualInfo::computeWeightedTemperature(const bool tempWeight)
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

void SharedIndividualInfo::generateOptimizationIndex(FILE   *fp,
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
            if (fit(p.first) && p.second.ntrain() > 0)
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

void SharedIndividualInfo::fillVectors(const int mindata)
{

    for(auto &optIndex : optIndex_)
    {
        auto                iType = optIndex.iType();
        ForceFieldParameter p;
        if (iType == InteractionType::CHARGE)
        {
            if (pd_->hasParticleType(optIndex.particleType()))
            {
                p = pd_->findParticleType(optIndex.particleType())->parameterConst(optIndex.parameterType());
            }
        }
        else if (pd_->interactionPresent(iType))
        {
            p = pd_->findForcesConst(iType).findParameterTypeConst(optIndex.id(), optIndex.parameterType());
        }
        if (p.ntrain() >= mindata)
        {
            paramNames_.push_back(optIndex_.name())
            mutability_.push_back(p.mutability())
            lowerBound_.push_back(p.minimum())
            upperBound_.push_back(p.maximum())
            ntrain_.push_back(p.ntrain())
        }
    }

}

/* * * * * * * * * * * * * * * * * * * * * *
* END: Vector stuff                        *
* * * * * * * * * * * * * * * * * * * * * */

/* * * * * * * * * * * * * * * * * * * * * *
* BEGIN: ParamClassIndex stuff             *
* * * * * * * * * * * * * * * * * * * * * */

void SharedIndividualInfo::assignParamClassIndex()
{

    paramClassIndex_.resize(paramNames_.size(), -1);

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
        if (paramClassIndex_[i] == -1)
        {
            if (!restClass)  // If <i> is the first parameter without a class
            {
                // Append "Other" to the list of classes
                paramClass_->push_back("Other");
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

void SharedIndividualInfo::setOutputFiles(const char                     *xvgconv,
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
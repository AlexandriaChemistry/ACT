/*
 * This source file is part of the Alexandria Chemistry Toolkit.
 *
 * Copyright (C) 2014-2025
 *
 * Developers:
 *             Mohammad Mehdi Ghahremanpour, 
 *             Julian Marrades,
 *             Marie-Madeleine Walz,
 *             Paul J. van Maaren, 
 *             David van der Spoel (Project leader)
 *
 * This program is free software; you can redistribute it and/or
 * modify it under the terms of the GNU General Public License
 * as published by the Free Software Foundation; either version 2
 * of the License, or (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this program; if not, write to the Free Software
 * Foundation, Inc., 51 Franklin Street, Fifth Floor, 
 * Boston, MA  02110-1301, USA.
 */
/*! \internal \brief
 * Implements part of the alexandria program.
 * \author Mohammad Mehdi Ghahremanpour <mohammad.ghahremanpour@icm.uu.se>
 * \author David van der Spoel <david.vanderspoel@icm.uu.se>
 * \author Julian Ramon Marrades Furquet <julian@marrad.es>
 */


#include "staticindividualinfo.h"

#include <algorithm>
#include <filesystem>

#include "act/forces/combinationrules.h"
#include "act/ga/genome.h"
#include "act/utility/memory_check.h"
#include "act/forcefield/forcefield_xml.h"

#include "gromacs/utility/stringutil.h"

#include <sys/types.h>
#include <sys/stat.h>

namespace alexandria
{

StaticIndividualInfo::StaticIndividualInfo(CommunicationRecord *cr) : cr_(cr)
{
    fillFittingTargets(iMolSelect::Train);
}

void StaticIndividualInfo::fillIdAndPrefix()
{
    id_ = cr_->middleManOrdinal();
    if (id_ >= 0)
    {
        prefix_ = "inds/";
    }
}

void StaticIndividualInfo::setOutputFile(const std::string &outputFile)
{
    if (cr_->isMasterOrMiddleMan())
    {   
        outputFile_     = prefix_ + "ind" + std::to_string(id_) + "-" + outputFile;
        outputFileLast_ = prefix_ + "ind" + std::to_string(id_) + "-last-" + outputFile;
    }
}

void StaticIndividualInfo::setFinalOutputFile(const std::string &outputFile)
{
    outputFile_ = outputFile;
}

/* * * * * * * * * * * * * * * * * * * * * *
* BEGIN: ForceField stuff                     *
* * * * * * * * * * * * * * * * * * * * * */

void StaticIndividualInfo::fillForceField(      FILE *fp,
                                       const char *pd_fn)
{
    int root = 0;
    if (!cr_->isHelper())
    {
        GMX_RELEASE_ASSERT(nullptr != pd_fn, "Give me a forcefield file name");
        if (debug)
        {
            fprintf(debug, "On node %d, will try to read %s\n", cr_->rank(), pd_fn);
        }
        try
        {
            alexandria::readForceField(pd_fn, &pd_);
        }
        GMX_CATCH_ALL_AND_EXIT_WITH_FATAL_ERROR;
        print_memory_usage(debug);
        root = cr_->rank();
    }
    else
    {
        root = cr_->superior();
    }
    /* Broadcasting Force Field Data from Middlemen to Helper nodes */
    if (cr_->isParallelCalc())
    {
        pd_.sendToHelpers(cr_, root);
    }
    if (nullptr != fp)
    {
        fprintf(fp, "There are %d atom types in the input file %s.\n\n",
                static_cast<int>(pd_.getNatypes()), pd_fn);
    }
}

void StaticIndividualInfo::updateForceField(const std::set<int>       &changed,
                                            const std::vector<double> &bases)
{
    if (bases.size() == 0)
    {
        return;
    }
    if (bases.size() != optIndex().size())
    {
        GMX_THROW(gmx::InternalError(gmx::formatString("Number of parameters to update in forcefield (%zu) does not match the known number of parameters (%zu)",
                                                       bases.size(),
                                                       optIndex().size()).c_str()));
    }
    std::set<int> mychanged = changed;
    if (mychanged.empty())
    {
        for(size_t i = 0; i < optIndex().size(); i++)
        {
            mychanged.insert(i);
        }
    }
    for(int n : mychanged)
    {
        auto p = optIndex_[n].forceFieldParameter();
        if (!p)
        {
            optIndex_[n].findForceFieldParameter(&pd_);
        }
        p = optIndex_[n].forceFieldParameter();
        if (p)
        {
            auto iType = optIndex_[n].iType();
            if (debug)
            {
                fprintf(debug, "Updating %s parameter %d (set size %lu) to %g\n",
                        interactionTypeToString(iType).c_str(), n, mychanged.size(), bases[n]);
            }
            p->setValue(bases[n]);
            p->setUpdated(true);
            // TODO fix the uncertainty
            // p->setUncertainty(psigma_[n]);
        }
    }
    // TODO: This will generate the whole matrix of parameters from scratch.
    // It would be more efficient to only generate the atom type pair parameters for
    // which one of the atomic Van der Waals or Coulomb parameters changed.
    generateDependentParameter(&pd_);
    for(int n : mychanged)
    {
        optIndex_[n].forceFieldParameter()->setUpdated(false);
    }
}

void StaticIndividualInfo::saveState(bool updateCheckSum)
{
    pd_.updateTimeStamp();
    if (updateCheckSum)
    {
        pd_.updateCheckSum();
    }
    if (!outputFile_.empty())
    {
        writeForceField(outputFile_, &pd_, false);
    }
}

void StaticIndividualInfo::saveState(bool updateCheckSum, const std::string &fname)
{
    pd_.updateTimeStamp();
    if (updateCheckSum)
    {
        pd_.updateCheckSum();
    }
    if (!fname.empty())
    {
        writeForceField(fname, &pd_, false);
    }
}

/* * * * * * * * * * * * * * * * * * * * * *
* END: ForceField stuff                       *
* * * * * * * * * * * * * * * * * * * * * */

/* * * * * * * * * * * * * * * * * * * * * *
* BEGIN: FittingTarget stuff               *
* * * * * * * * * * * * * * * * * * * * * */

void StaticIndividualInfo::sumChiSquared(bool             parallel,
                                         iMolSelect       ims)
{
    // Now sum over processors, except the master/middleman!
    if (cr_->isParallelCalc() && parallel)
    {
        for (auto &ft : targets_[ims])
        {
            auto chi2 = ft.second.chiSquared();
            auto tw = ft.second.totalWeight();
            if (tw > 0 && debug)
            {
                fprintf(debug, "node %d %s before sum chi2 %g tw %g weighted %g\n",
                        cr_->rank(), rmsName(ft.first), chi2, tw, ft.second.chiSquaredWeighted());
            }
            cr_->sumd_helpers(1, &chi2);
            ft.second.setChiSquared(chi2);
            cr_->sumd_helpers(1, &tw);
            ft.second.setTotalWeight(tw);
            if (tw > 0 && debug)
            {
                fprintf(debug, "node %d %s after sum chi2 %g tw %g weighted %g\n",
                        cr_->rank(), rmsName(ft.first), chi2, tw, ft.second.chiSquaredWeighted());
            }
        }
    }
    auto myft = targets_.find(ims);
    if (myft != targets_.end())
    {
        auto etot = myft->second.find(eRMS::TOT);
        GMX_RELEASE_ASSERT(etot != myft->second.end(),
                           gmx::formatString("Cannot find etot for %s",
                                             iMolSelectName(ims)).c_str());
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
        etot->second.setTotalWeight(1);
        if (debug)
        {
            fprintf(debug, "node %d %s after sum chi2 %g tw %g\n", cr_->rank(), rmsName(eRMS::TOT),
                    etot->second.chiSquared(), etot->second.totalWeight());
        }
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

void StaticIndividualInfo::fillFittingTargets(iMolSelect ims)
{
    std::map<eRMS, FittingTarget> ft;
    for ( auto &rms : geteRMSNames() )
    {
        FittingTarget    fft(rms.first, ims);
        ft.insert(std::pair<eRMS, FittingTarget>(rms.first, std::move(fft)));
    }
    targets_.insert(std::pair<iMolSelect, std::map<eRMS, FittingTarget>>(ims, std::move(ft)));
    auto etot = target(ims, eRMS::TOT);
    GMX_RELEASE_ASSERT(etot != nullptr, "Could not find etot");
    etot->setWeight(1);
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

void StaticIndividualInfo::generateOptimizationIndex(FILE                      *fp,
                                                     const MolGen              *mg,
                                                     const CommunicationRecord *cr)
{
    if (!cr->isHelper())
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
                            else if (debug)
                            {
                                fprintf(debug, "WARNING: Not enough data to train %s-%s\n",
                                        fpl.first.id().c_str(), param.first.c_str());
                            }
                        }
                    }
                }
            }
        }
        for(auto &pt : pd_.particleTypesConst())
        {
            for(auto &p : pt.second.parametersConst())
            {
                if (p.second.isMutable() && mg->fit(p.first) && p.second.ntrain() > 0)
                {
                    optIndex_.push_back(OptimizationIndex(pt.second.id().id(), p.first));
                }
            }
        }
        if (fp)
        {
            fprintf(fp, "There are %zu parameters to train.\n", optIndex_.size());
            fprintf(fp, "Identifier        Parameter     Minimum     Maximum\n");
            auto fcs = pd_.forcesConst();
            for(auto &i : optIndex_)
            {
                auto ff = fcs.find(i.iType());
                if (fcs.end() != ff)
                {
                    auto pp = ff->second.parametersConst();
                    auto gg = pp.find(i.id());
                    if (pp.end() != gg)
                    {
                        const auto &fs = gg->second.find(i.parameterType());
                        if (gg->second.end() != fs)
                        {
                            fprintf(fp, "%-15s %11s  %10g  %10g\n",
                                    i.id().id().c_str(), i.parameterType().c_str(),
                                    fs->second.minimum(), fs->second.maximum());
                        }
                    }
                }
                else if (InteractionType::CHARGE == i.iType())
                {
                    auto pt = pd_.findParticleType(i.particleType());
                    if (!pt->hasParameter(i.parameterType()))
                    {
                        GMX_THROW(gmx::InvalidInputError(gmx::formatString("No such parameter type '%s' for particle '%s'", i.parameterType().c_str(), i.particleType().c_str()).c_str()));
                    }
                    else
                    {
                        auto &pp = pt->parameterConst(i.parameterType());
                        fprintf(fp, "%-15s %11s  %10g  %10g\n",
                                i.particleType().c_str(), i.parameterType().c_str(),
                                pp.minimum(), pp.maximum());
                    }
                }
                else
                {
                    GMX_THROW(gmx::InternalError(gmx::formatString("Incomprehensible parameter '%s' to train", i.name().c_str()).c_str()));
                }
            }
        }
        // Now send the data over to my helpers
        for(const int dst : cr->helpers())
        {
            cr->send(dst, static_cast<int>(optIndex_.size()));
            for(size_t i = 0; i < optIndex_.size(); i++)
            {
                optIndex_[i].send(cr, dst);
            }
        }
    }
    else
    {
        // Receive the data from my superior
        int src = cr->superior();
        GMX_RELEASE_ASSERT(src >= 0, "No superior");
        int nopt;
        cr->recv(src, &nopt);
        optIndex_.resize(nopt);
        for(int i = 0; i < nopt; i++)
        {
            optIndex_[i].receive(cr, src);
        }
    }
}

/* * * * * * * * * * * * * * * * * * * * * *
* END: OptimizationIndex stuff             *
* * * * * * * * * * * * * * * * * * * * * */

/* * * * * * * * * * * * * * * * * * * * * *
* BEGIN: Vector stuff                      *
* * * * * * * * * * * * * * * * * * * * * */

void StaticIndividualInfo::fillVectors(unsigned int mindata)
{
    if (cr_->isMasterOrMiddleMan())
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
                paramNamesWOClass_.push_back(
                                             optIndex.name().substr(0, optIndex.name().rfind("-"))
                                             );
                mutability_.push_back(p.mutability());
                lowerBound_.push_back(p.minimum());
                upperBound_.push_back(p.maximum());
                ntrain_.push_back(p.ntrain());
            }
        }
        for(int dest : cr_->helpers())
        {
            cr_->send(dest, defaultParam_);
            cr_->send(dest, lowerBound_);
            cr_->send(dest, upperBound_);
            cr_->send(dest, static_cast<int>(paramNames_.size()));
            for(size_t i = 0; i < paramNames_.size(); i++)
            {
                cr_->send(dest, paramNames_[i]);
                cr_->send(dest, paramNamesWOClass_[i]);
                cr_->send(dest, ntrain_[i]);
                std::string mutab = mutabilityName(mutability_[i]);
                cr_->send(dest, mutab);
            }
        }
    }
    else
    {
        int src = cr_->superior();
        cr_->recv(src, &defaultParam_);
        cr_->recv(src, &lowerBound_);
        cr_->recv(src, &upperBound_);
        int nparam;
        cr_->recv(src, &nparam);
        for(int i = 0; i < nparam; i++)
        {
            std::string paramName, paramNameWOClass, mutab;
            cr_->recv(src, &paramName);
            cr_->recv(src, &paramNameWOClass);
            paramNames_.push_back(paramName);
            paramNamesWOClass_.push_back(paramNameWOClass);
            int nt;
            cr_->recv(src, &nt);
            ntrain_.push_back(nt);
            cr_->recv(src, &mutab);
            Mutability mmm;
            if (nameToMutability(mutab, &mmm))
            {
                mutability_.push_back(mmm);
            }
            else
            {
                GMX_THROW(gmx::InternalError(gmx::formatString("Received unintelligable string %s instead of a mutability", mutab.c_str()).c_str()));
            }
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
            const auto tokenizedName = gmx::splitDelimitedString(paramNames_[j], '-');
            if (std::find(tokenizedName.begin(), tokenizedName.end(), paramClass_[i]) != tokenizedName.end())
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

void StaticIndividualInfo::makeIndividualDir()
{
    if (!prefix_.empty())
    {
        (void) std::filesystem::create_directory(prefix_);
    }
}

/* * * * * * * * * * * * * * * * * * * * * *
* END: File stuff                          *
* * * * * * * * * * * * * * * * * * * * * */

/* * * * * * * * * * * * * * * * * * * * * *
* BEGIN: Getters and setters               *
* * * * * * * * * * * * * * * * * * * * * */

double StaticIndividualInfo::getParamSpaceVolume(const bool logScale) const
{
    double aggregate = 1;
    for (size_t i = 0; i < nParam(); i++)
    {
        aggregate *= upperBound_[i] - lowerBound_[i];
    }
    return logScale ? log(aggregate) : aggregate;
}

/* * * * * * * * * * * * * * * * * * * * * *
* END: Getters and setters                 *
* * * * * * * * * * * * * * * * * * * * * */

} //namespace alexandria

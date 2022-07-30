/*
 * This source file is part of the Alexandria Chemistry Toolkit.
 *
 * Copyright (C) 2014-2022
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
 * \author Julian Ramon Marrades Furquet <julian.marrades@hotmail.es>
 * \author Oskar Tegby <oskar.tegby@it.uu.se>
 */
#include "devcomputer.h"

#include <numeric>
#include <string>
#include <vector>

#include "act/basics/identifier.h"
#include "act/utility/communicationrecord.h"
#include "act/utility/units.h"
#include "gromacs/math/vecdump.h"

namespace alexandria
{

/*! \brief Compute penalty for variables that are out of bounds
 * \param[in] x       The actual value
 * \param[in] min     The minimum allowed value
 * \param[in] max     The maximum allowed value
 * \return 0 when in bounds, square deviation from bounds otherwise.
 */
static double l2_regularizer(double x, double min, double max)
{
    double p = 0;
    if (x < min)
    {
        p = (0.5 * gmx::square(x-min));
    }
    else if (x > max)
    {
        p = (0.5 * gmx::square(x-max));
    }
    return p;
}


/* * * * * * * * * * * * * * * * * * * * * *
* BEGIN: BoundsDevComputer                 *
* * * * * * * * * * * * * * * * * * * * * */

void BoundsDevComputer::calcDeviation(gmx_unused const ForceComputer       *forceComputer,
                                      gmx_unused MyMol                     *mymol,
                                      std::map<eRMS, FittingTarget>        *targets,
                                      Poldata                              *poldata,
                                      const std::vector<double>            &param,
                                      gmx_unused const CommunicationRecord *commrec)
{
    auto   mytarget = targets->find(eRMS::BOUNDS);
    if (targets->end() != mytarget)
    {
        double bound = 0;
        size_t n     = 0;
        for (auto &optIndex : *optIndex_)
        {
            InteractionType iType = optIndex.iType();
            ForceFieldParameter p;
            if (iType == InteractionType::CHARGE)
            {
                p = poldata->findParticleType(optIndex.particleType())->parameterConst(optIndex.parameterType());
            }
            else if (poldata->interactionPresent(iType))
            {
                auto &fs = poldata->findForcesConst(iType);
                p = fs.findParameterTypeConst(optIndex.id(), optIndex.parameterType());
            }
            if (p.mutability() == Mutability::Bounded)
            {
                real db = l2_regularizer(param[n], p.minimum(), p.maximum());
                bound += db;
                if (verbose_ && logfile_ && db != 0.0)
                {
                    fprintf(logfile_, "Variable %s is %g, should be within %g and %g\n",
                            optIndex.name().c_str(), param[n], p.minimum(), p.maximum());
                }
            }

            n++;
        }
        mytarget->second.increase(1, bound);
        GMX_RELEASE_ASSERT(n == param.size(),
                           gmx::formatString("Death horror error. n=%zu param.size()=%zu", n, param.size()).c_str());
    }
    mytarget = targets->find(eRMS::UNPHYSICAL);
    if (targets->end() != mytarget)
    {
        double bound = 0;
        // Check whether shell zeta > core zeta. Only for polarizable models.
        auto   itype = InteractionType::COULOMB;
        if (poldata->polarizable() && poldata->interactionPresent(itype))
        {
            auto &fs             = poldata->findForcesConst(itype);
            std::string poltype  = "poltype";
            std::string zetatype = "zetatype";
            for(const auto &p : poldata->particleTypesConst())
            {
                if (p.hasOption(poltype))
                {
                    auto coreID  = Identifier(p.optionValue(zetatype));
                    auto shell   = poldata->findParticleType(p.optionValue(poltype));
                    auto shellID = Identifier(shell->optionValue(zetatype));
                    auto fpshell = fs.findParameterTypeConst(shellID, "zeta");
                    auto fpcore  = fs.findParameterTypeConst(coreID, "zeta");
                    double deltaZeta = fpshell.value() - fpcore.value();
                    if (deltaZeta > 0)
                    {
                        bound += gmx::square(deltaZeta);
                    }
                }
            }
        }
        mytarget->second.increase(1, bound);
    }
}

/* * * * * * * * * * * * * * * * * * * * * *
* END: BoundsDevComputer                   *
* * * * * * * * * * * * * * * * * * * * * */

/* * * * * * * * * * * * * * * * * * * * * *
* BEGIN: ChargeCM5DevComputer              *
* * * * * * * * * * * * * * * * * * * * * */

void ChargeCM5DevComputer::calcDeviation(gmx_unused const ForceComputer       *forceComputer,
                                         MyMol                                *mymol,
                                         std::map<eRMS, FittingTarget>        *targets,
                                         Poldata                              *poldata,
                                         gmx_unused const std::vector<double> &param,
                                         gmx_unused const CommunicationRecord *commrec)
{
    double qtot = 0;
    int i = 0;
    const auto &myatoms = mymol->atomsConst();
    std::vector<double> qcm5;
    QtypeProps *qp = mymol->qTypeProps(qType::CM5);
    if (qp)
    {
        qcm5 = qp->charge();
        if (debug)
        {
            for (size_t j = 0; j < myatoms.size(); j++)
            {
                fprintf(debug, "Charge %lu. CM5 = %g ACM = %g\n", j, qcm5[j], myatoms[j].charge());
            }
        }
    }
    // Iterate over the atoms
    for (size_t j = 0; j < myatoms.size(); j++)
    {
        if (myatoms[j].pType() == eptShell)
        {
            continue;
        }
        ParticleTypeIterator       atype = poldata->findParticleType(myatoms[j].ffType());
        const ForceFieldParameter &qparm = atype->parameterConst("charge");
        double qj  = myatoms[j].charge();
        double qjj = qj;
        // TODO: only count in real shells
        if (mymol->haveShells() && j < myatoms.size()-1 &&
            myatoms[j+1].pType() == eptShell)
        {
            qjj += myatoms[j+1].charge();
        }
        qtot += qjj;
        switch (qparm.mutability())
        {
        case Mutability::Fixed:
            if (qparm.value() != qj)
            {
                GMX_THROW(gmx::InternalError(gmx::formatString("Fixed charge for atom %s in %s was changed from %g to %g",
                                                               myatoms[j].name().c_str(),
                                                               mymol->getMolname().c_str(), qparm.value(), qj).c_str()));
            }
            break;
        case Mutability::ACM:
        case Mutability::Bounded:
            {
                if ((*targets).find(eRMS::CHARGE)->second.weight() > 0)
                {
                    if (qparm.maximum() > qparm.minimum())
                    {
                        (*targets).find(eRMS::CHARGE)->second.increase(1, l2_regularizer(qjj, qparm.minimum(),
                                                                                         qparm.maximum()));
                    }
                }
            }
            break;
        default:
            break;
        }
        if (qp &&
            qparm.mutability() != Mutability::Fixed &&
            (*targets).find(eRMS::CM5)->second.weight() > 0)
        {
            real dq2 = gmx::square(qjj - qcm5[i]);
            (*targets).find(eRMS::CM5)->second.increase(1, dq2);
        }
        i += 1;
    }
    (*targets).find(eRMS::CHARGE)->second.increase(1, gmx::square(qtot - mymol->totalCharge()));

}

/* * * * * * * * * * * * * * * * * * * * * *
* END: ChargeCM5DevComputer                *
* * * * * * * * * * * * * * * * * * * * * */

/* * * * * * * * * * * * * * * * * * * * * *
* BEGIN: EspDevComputer                    *
* * * * * * * * * * * * * * * * * * * * * */

void EspDevComputer::calcDeviation(gmx_unused const ForceComputer       *forceComputer,
                                   MyMol                                *mymol,
                                   std::map<eRMS, FittingTarget>        *targets,
                                   Poldata                              *poldata,
                                   gmx_unused const std::vector<double> &param,
                                   gmx_unused const CommunicationRecord *commrec)
{

    real rrms     = 0;
    real cosangle = 0;
    QgenResp *qgr = mymol->qTypeProps(qType::Calc)->qgenResp();
    if (mymol->haveShells())
    {
        qgr->updateAtomCoords(mymol->x());
    }
    if (fit_)
    {
        qgr->updateZeta(mymol->atomsConst(), poldata);
    }
    dumpQX(logfile_, mymol, "ESP");
    qgr->updateAtomCharges(mymol->atomsConst());
    qgr->calcPot(poldata->getEpsilonR());
    real mae, mse;
    real rms = qgr->getStatistics(&rrms, &cosangle, &mae, &mse);
    double myRms = convertToGromacs(rms, "Hartree/e");
    size_t nEsp = qgr->nEsp();
    (*targets).find(eRMS::ESP)->second.increase(nEsp, gmx::square(myRms)*nEsp);
    if (debug)
    {
        fprintf(debug, "%s ESPrms = %g cosangle = %g\n",
                mymol->getMolname().c_str(),
                myRms, cosangle);
    }

}

void EspDevComputer::dumpQX(FILE *fp, MyMol *mol, const std::string &info)
{
    if (false && fp)
    {
        std::string label = mol->getMolname() + "-" + info;
        fprintf(fp, "%s q:", label.c_str());
        auto myatoms = mol->topology()->atoms();
        for (size_t i = 0; i < myatoms.size(); i++)
        {
            fprintf(fp, " %g", myatoms[i].charge());
        }
        fprintf(fp, "\n");
        auto top = mol->topology();
        if (top->hasEntry(InteractionType::POLARIZATION))
        {
            fprintf(fp, "%s alpha", label.c_str());
            for(auto &topentry : top->entry(InteractionType::POLARIZATION))
            {
                //fprintf(fp, " %g", mol->ltop_->idef.iparams[topentry->gromacsType()].polarize.alpha);
            }
            fprintf(fp, "\n");
        }
        pr_rvecs(fp, 0, label.c_str(), mol->x().rvec_array(), myatoms.size());
    }
}

/* * * * * * * * * * * * * * * * * * * * * *
* END: EspDevComputer                      *
* * * * * * * * * * * * * * * * * * * * * */


/* * * * * * * * * * * * * * * * * * * * * *
* BEGIN: PolarDevComputer                  *
* * * * * * * * * * * * * * * * * * * * * */

PolarDevComputer::PolarDevComputer(    FILE  *logfile,
                                   const bool verbose)
    : DevComputer(logfile, verbose)
{
    convert_ = convertFromGromacs(1.0, mpo_unit2(MolPropObservable::POLARIZABILITY));
}

void PolarDevComputer::calcDeviation(const ForceComputer                  *forceComputer,
                                     MyMol                                *mymol,
                                     std::map<eRMS, FittingTarget>        *targets,
                                     gmx_unused Poldata                   *poldata,
                                     gmx_unused const std::vector<double> &param,
                                     gmx_unused const CommunicationRecord *commrec)
{
    mymol->CalcPolarizability(forceComputer);
    auto aelec = mymol->qTypeProps(qType::Elec)->polarizabilityTensor();
    auto acalc = mymol->qTypeProps(qType::Calc)->polarizabilityTensor();
    double diff2 = 0;
    for(int i = 0; i < DIM; i++)
    {
        for(int j = 0; j < DIM; j++)
        {
            diff2 += gmx::square(convert_*(aelec[i][j]-acalc[i][j]));
        }
    }
    
    if (false && logfile_)
    {
        fprintf(logfile_, "DIFF %s %g\n", mymol->getMolname().c_str(), diff2);
    }
    (*targets).find(eRMS::Polar)->second.increase(1, diff2);
}

/* * * * * * * * * * * * * * * * * * * * * *
* END: PolarDevComputer                    *
* * * * * * * * * * * * * * * * * * * * * */

/* * * * * * * * * * * * * * * * * * * * * *
* BEGIN: MultiPoleDevComputer              *
* * * * * * * * * * * * * * * * * * * * * */

void MultiPoleDevComputer::calcDeviation(gmx_unused const ForceComputer       *forceComputer,
                                         MyMol                                *mymol,
                                         std::map<eRMS, FittingTarget>        *targets,
                                         gmx_unused Poldata                   *poldata,
                                         gmx_unused const std::vector<double> &param,
                                         gmx_unused const CommunicationRecord *commrec)
{
    if (!(mymol->qTypeProps(qType::Elec)->hasMultipole(mpo_) &&
          mymol->qTypeProps(qType::Calc)->hasMultipole(mpo_)))
    {
        return;
    }
    auto qelec = mymol->qTypeProps(qType::Elec)->getMultipole(mpo_);
    auto qcalc = mymol->qTypeProps(qType::Calc)->getMultipole(mpo_);
    double delta = 0;
    for (size_t mm = 0; mm < qelec.size(); mm++)
    {
        delta += gmx::square(qcalc[mm] - qelec[mm]);
    }
    eRMS rms;
    switch (mpo_)
    {
    case MolPropObservable::DIPOLE:
        rms = eRMS::MU;
        break;
    case MolPropObservable::QUADRUPOLE:
        rms = eRMS::QUAD;
        break;
    case MolPropObservable::OCTUPOLE:
        rms = eRMS::OCT;
        break;
    case MolPropObservable::HEXADECAPOLE:
        rms = eRMS::HEXADEC;
        break;
    default:
        GMX_THROW(gmx::InternalError(gmx::formatString("Not support MolPropObservable %s", mpo_name(mpo_)).c_str()));
    }
    
    (*targets).find(rms)->second.increase(1, delta);
}

/* * * * * * * * * * * * * * * * * * * * * *
* END: MultipoleDevComputer                *
* * * * * * * * * * * * * * * * * * * * * */

/* * * * * * * * * * * * * * * * * * * * * *
* BEGIN: HarmonicsDevComputer              *
* * * * * * * * * * * * * * * * * * * * * */

HarmonicsDevComputer::HarmonicsDevComputer(      FILE              *logfile,
                                           const bool               verbose,
                                                 MolPropObservable  mpo)
    : DevComputer(logfile, verbose), mpo_(mpo) 
{
}

void HarmonicsDevComputer::calcDeviation(const ForceComputer                  *forceComputer,
                                         MyMol                                *mymol,
                                         std::map<eRMS, FittingTarget>        *targets,
                                         gmx_unused Poldata                   *poldata,
                                         gmx_unused const std::vector<double> &param,
                                         gmx_unused const CommunicationRecord *commrec)
{
    // Only compute frequencies for structures that have an optimize reference
    if (JobType::OPT != mymol->jobType())
    {
        return;
    }
    std::vector<gmx::RVec> coords = mymol->xOriginal();
    auto eMin = handler_.minimizeCoordinates(mymol, forceComputer, simConfig_, &coords,
                                             nullptr, nullptr);
    if (eMinimizeStatus::OK != eMin)
    {
        // Something fishy happened, but it means we cannot use this structure
        // for computing frequencies.
        return;
    }
    // Compute frequencies
    std::vector<double> frequencies, intensities;
    handler_.nma(mymol, forceComputer, &coords, &frequencies, &intensities);

    switch (mpo_)
    {
    case MolPropObservable::FREQUENCY:
        {
            auto ref_freqs = mymol->referenceFrequencies();
            if (ref_freqs.size() != frequencies.size())
            {
                fprintf(stderr, "Reference frequencies size %zu, but calculated %zu for %s. Ignoring frequencies for this compound.\n",
                        ref_freqs.size(), frequencies.size(), mymol->getMolname().c_str());
            }
            else
            {
                double delta = 0;
                for(size_t k = 0; k < frequencies.size(); k++)
                {
                    delta += gmx::square(frequencies[k]-ref_freqs[k]);
                }
                (*targets).find(eRMS::FREQUENCY)->second.increase(1, delta);
            }
        }
        break;
    case MolPropObservable::INTENSITY:
        {
            auto ref_intens = mymol->referenceIntensities();
            if (ref_intens.size() != intensities.size())
            {
                fprintf(stderr, "Reference intensities size %zu, but calculated %zu for %s. Ignoring frequencies for this compound.\n",
                        ref_intens.size(), intensities.size(), mymol->getMolname().c_str());
            }
            else
            {
                double delta = 0;
                for(size_t k = 0; k < intensities.size(); k++)
                {
                    delta += gmx::square(intensities[k]-ref_intens[k]);
                }
                (*targets).find(eRMS::INTENSITY)->second.increase(1, delta);
            }
        }
        break;
    default:
        fprintf(stderr, "Don't know how to handle %s in this devcomputer\n",
                mpo_name(mpo_));
    }
}

/* * * * * * * * * * * * * * * * * * * * * *
* END: HarmonicsDevComputer                *
* * * * * * * * * * * * * * * * * * * * * */


/* * * * * * * * * * * * * * * * * * * * * *
* BEGIN: ForceEnergyDevComputer                 *
* * * * * * * * * * * * * * * * * * * * * */

void ForceEnergyDevComputer::calcDeviation(const ForceComputer                  *forceComputer,
                                           MyMol                                *mymol,
                                           std::map<eRMS, FittingTarget>        *targets,
                                           gmx_unused Poldata                   *poldata,
                                           gmx_unused const std::vector<double> &param,
                                           gmx_unused const CommunicationRecord *commrec)
{
    std::vector<std::pair<double, double> >                 eMap;
    std::vector<std::vector<std::pair<double, double> > >   fMap;
    std::vector<std::pair<double, std::map<InteractionType, double> > > enerAllMap;
    mymol->forceEnergyMaps(forceComputer, &fMap, &eMap, &enerAllMap);

    auto tf = targets->find(eRMS::Force2);
    if (tf != targets->end() && !fMap.empty())
    {
        for(const auto &fstruct : fMap)
        {
            for(const auto &ff : fstruct)
            {
                tf->second.increase(1, gmx::square(ff.first-ff.second));
            }
        }
    }
    auto te = targets->find(eRMS::EPOT);
    if (te != targets->end() && !eMap.empty())
    {
        for(const auto &ff : eMap)
        {
            // TODO Double check if the atomizationEnergy is needed.
            // auto enerexp = mymol->atomizationEnergy() + ff.first;
            te->second.increase(1, gmx::square(ff.first-ff.second));
        }
    }
    auto ti = targets->find(eRMS::Interaction);
    if (ti != targets->end() && !eMap.empty())
    {
        for(const auto &ff : eMap)
        {
            ti->second.increase(1, gmx::square(ff.first-ff.second));
        }
    }
}

/* * * * * * * * * * * * * * * * * * * * * *
* END: ForceEnergyDevComputer              *
* * * * * * * * * * * * * * * * * * * * * */

} // namespace alexandria

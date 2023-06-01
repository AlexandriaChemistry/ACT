/*
 * This source file is part of the Alexandria Chemistry Toolkit.
 *
 * Copyright (C) 2014-2023
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
#include "act/forcefield/forcefield_parametername.h"
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

void BoundsDevComputer::calcDeviation(gmx_unused const ForceComputer    *forceComputer,
                                      gmx_unused ACTMol                 *actmol,
                                      gmx_unused std::vector<gmx::RVec> *coords,
                                      std::map<eRMS, FittingTarget>     *targets,
                                      const ForceField                  *forcefield)
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
                p = forcefield->findParticleType(optIndex.particleType())->parameterConst(optIndex.parameterType());
            }
            else if (forcefield->interactionPresent(iType))
            {
                auto &fs = forcefield->findForcesConst(iType);
                p = fs.findParameterTypeConst(optIndex.id(), optIndex.parameterType());
            }
            if (p.mutability() == Mutability::Bounded)
            {
                real db = l2_regularizer(p.value(), p.minimum(), p.maximum());
                bound += db;
                if (verbose_ && logfile_ && db != 0.0)
                {
                    fprintf(logfile_, "Variable %s is %g, should be within %g and %g\n",
                            optIndex.name().c_str(), p.value(), p.minimum(), p.maximum());
                }
            }

            n++;
        }
        mytarget->second.increase(1, bound);
    }
    mytarget = targets->find(eRMS::UNPHYSICAL);
    if (targets->end() != mytarget)
    {
        double bound = 0;
        // Check whether shell zeta > core zeta. Only for polarizable models.
        auto   itype = InteractionType::COULOMB;
        if (forcefield->polarizable() && forcefield->interactionPresent(itype))
        {
            auto &fs             = forcefield->findForcesConst(itype);
            std::string poltype  = "poltype";
            std::string zetatype = "zetatype";
            for(const auto &p : forcefield->particleTypesConst())
            {
                if (p.hasOption(poltype))
                {
                    auto coreID  = Identifier(p.optionValue(zetatype));
                    auto shell   = forcefield->findParticleType(p.optionValue(poltype));
                    auto shellID = Identifier(shell->optionValue(zetatype));
                    auto fpshell = fs.findParameterTypeConst(shellID, "zeta");
                    auto fpcore  = fs.findParameterTypeConst(coreID, "zeta");
                    if (fpshell.ntrain() > 0 && fpcore.ntrain() > 0)
                    {
                        double deltaZeta = zetaDiff_ - (fpcore.value() - fpshell.value());
                        if (deltaZeta > 0)
                        {
                            bound += gmx::square(deltaZeta);
                        }
                    }
                }
            }
        }
        itype = InteractionType::BONDS;
        if (forcefield->interactionPresent(itype))
        {
            auto &fs = forcefield->findForcesConst(itype);
            if (fs.gromacsType() == F_CUBICBONDS)
            {
                for(const auto &ffp : fs.parametersConst())
                {
                    auto param = ffp.second;
                    auto rmax  = param[cubic_name[cubicRMAX]].value();
                    auto blen  = param[cubic_name[cubicLENGTH]].value();
                    // We want the maximum in the potential to be at least 0.1 nm further away than then minimum
                    if (rmax < blen+0.1)
                    {
                        bound += gmx::square(rmax-blen-0.1);
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
                                         ACTMol                                *actmol,
                                         gmx_unused std::vector<gmx::RVec>    *coords,
                                         std::map<eRMS, FittingTarget>        *targets,
                                         const ForceField                        *forcefield)
{
    double qtot = 0;
    int i = 0;
    const auto &myatoms = actmol->atomsConst();
    std::vector<double> qcm5;
    QtypeProps *qp = actmol->qTypeProps(qType::CM5);
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
        ParticleTypeConstIterator  atype = forcefield->findParticleType(myatoms[j].ffType());
        const ForceFieldParameter &qparm = atype->parameterConst("charge");
        double qj  = myatoms[j].charge();
        double qjj = qj;
        for(const auto &ss : myatoms[j].shells())
        {
            qjj += myatoms[ss].charge();
        }
        qtot += qjj;
        switch (qparm.mutability())
        {
        case Mutability::Fixed:
            if (qparm.value() != qj)
            {
                GMX_THROW(gmx::InternalError(gmx::formatString("Fixed charge for atom %s in %s was changed from %g to %g",
                                                               myatoms[j].name().c_str(),
                                                               actmol->getMolname().c_str(), qparm.value(), qj).c_str()));
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
    (*targets).find(eRMS::CHARGE)->second.increase(1, gmx::square(qtot - actmol->totalCharge()));

}

/* * * * * * * * * * * * * * * * * * * * * *
* END: ChargeCM5DevComputer                *
* * * * * * * * * * * * * * * * * * * * * */

/* * * * * * * * * * * * * * * * * * * * * *
* BEGIN: EspDevComputer                    *
* * * * * * * * * * * * * * * * * * * * * */

void EspDevComputer::calcDeviation(gmx_unused const ForceComputer *forceComputer,
                                   ACTMol                         *actmol,
                                   std::vector<gmx::RVec>         *coords,
                                   std::map<eRMS, FittingTarget>  *targets,
                                   const ForceField               *forcefield)
{
    real rrms     = 0;
    real cosangle = 0;
    QgenResp *qgr = actmol->qTypeProps(qType::Calc)->qgenResp();
    if (actmol->haveShells())
    {
        qgr->updateAtomCoords(*coords);
    }
    if (fit_)
    {
        qgr->updateZeta(actmol->atomsConst(), forcefield);
    }
    dumpQX(logfile_, actmol, *coords, "ESP");
    qgr->updateAtomCharges(actmol->atomsConst());
    double epsilonr;
    if (!ffOption(*forcefield, InteractionType::COULOMB, "epsilonr", &epsilonr))
    {
        epsilonr = 1;
    }
    qgr->calcPot(epsilonr);
    real mae, mse;
    real rms = qgr->getStatistics(&rrms, &cosangle, &mae, &mse);
    double myRms = convertToGromacs(rms, "Hartree/e");
    size_t nEsp = qgr->nEsp();
    (*targets).find(eRMS::ESP)->second.increase(nEsp, gmx::square(myRms)*nEsp);
    if (debug)
    {
        fprintf(debug, "%s ESPrms = %g cosangle = %g\n",
                actmol->getMolname().c_str(),
                myRms, cosangle);
    }

}

void EspDevComputer::dumpQX(FILE                         *fp,
                            const ACTMol                  *mol,
                            const std::vector<gmx::RVec> &coords,
                            const std::string            &info)
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
                fprintf(fp, " %g", topentry->params()[polALPHA]);
            }
            fprintf(fp, "\n");
        }
        pr_rvecs(fp, 0, label.c_str(), as_rvec_array(coords.data()),
                 myatoms.size());
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

void PolarDevComputer::calcDeviation(const ForceComputer               *forceComputer,
                                     ACTMol                            *actmol,
                                     gmx_unused std::vector<gmx::RVec> *coords,
                                     std::map<eRMS, FittingTarget>     *targets,
                                     const ForceField                  *forcefield)
{
    actmol->CalcPolarizability(forcefield, forceComputer);
    auto aelec = actmol->qTypeProps(qType::Elec)->polarizabilityTensor();
    auto acalc = actmol->qTypeProps(qType::Calc)->polarizabilityTensor();
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
        fprintf(logfile_, "DIFF %s %g\n", actmol->getMolname().c_str(), diff2);
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
                                         ACTMol                                *actmol,
                                         gmx_unused std::vector<gmx::RVec>    *coords,
                                         std::map<eRMS, FittingTarget>        *targets,
                                         gmx_unused const ForceField             *forcefield)
{
    if (!(actmol->qTypeProps(qType::Elec)->hasMultipole(mpo_) &&
          actmol->qTypeProps(qType::Calc)->hasMultipole(mpo_)))
    {
        return;
    }
    auto qelec = actmol->qTypeProps(qType::Elec)->getMultipole(mpo_);
    auto qcalc = actmol->qTypeProps(qType::Calc)->getMultipole(mpo_);
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

void HarmonicsDevComputer::calcDeviation(const ForceComputer           *forceComputer,
                                         ACTMol                        *actmol,
                                         std::vector<gmx::RVec>        *coords,
                                         std::map<eRMS, FittingTarget> *targets,
                                         const ForceField              *forcefield)
{
    // Only compute frequencies for structures that have an optimize reference
    if (JobType::OPT != actmol->jobType())
    {
        return;
    }
    auto eMin = handler_.minimizeCoordinates(forcefield, actmol, forceComputer, simConfig_, coords,
                                             nullptr, nullptr, {});
    if (eMinimizeStatus::OK != eMin)
    {
        // Something fishy happened, but it means we cannot use this structure
        // for computing frequencies.
        return;
    }
    // Compute frequencies
    std::vector<double> frequencies, intensities;
    handler_.nma(forcefield, actmol, forceComputer, coords, &frequencies, &intensities);

    switch (mpo_)
    {
    case MolPropObservable::FREQUENCY:
        {
            auto ref_freqs = actmol->referenceFrequencies();
            if (ref_freqs.size() != frequencies.size())
            {
                fprintf(stderr, "Reference frequencies size %zu, but calculated %zu for %s. Ignoring frequencies for this compound.\n",
                        ref_freqs.size(), frequencies.size(), actmol->getMolname().c_str());
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
            auto ref_intens = actmol->referenceIntensities();
            if (ref_intens.size() != intensities.size())
            {
                fprintf(stderr, "Reference intensities size %zu, but calculated %zu for %s. Ignoring frequencies for this compound.\n",
                        ref_intens.size(), intensities.size(), actmol->getMolname().c_str());
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

void ForceEnergyDevComputer::calcDeviation(const ForceComputer               *forceComputer,
                                           ACTMol                            *actmol,
                                           gmx_unused std::vector<gmx::RVec> *coords,
                                           std::map<eRMS, FittingTarget>     *targets,
                                           const ForceField                  *forcefield)
{
    std::vector<ACTEnergy>                                  energyMap;
    std::vector<std::vector<std::pair<double, double> > >   forceMap;
    std::vector<std::pair<double, std::map<InteractionType, double> > > enerComponentMap, interactionEnergyMap;
    actmol->forceEnergyMaps(forcefield, forceComputer, &forceMap, &energyMap,
                            &interactionEnergyMap, &enerComponentMap);

    auto tf = targets->find(eRMS::Force2);
    if (tf != targets->end() && !forceMap.empty())
    {
        for(const auto &fstruct : forceMap)
        {
            for(const auto &ff : fstruct)
            {
                if (std::isnan(ff.second))
                {
                    printf("Force for %s is NaN\n", actmol->getMolname().c_str());
                }
                tf->second.increase(1, gmx::square(ff.first-ff.second));
            }
        }
    }
    auto te = targets->find(eRMS::EPOT);
    if (te != targets->end() && !energyMap.empty())
    {
        for(const auto &ff : energyMap)
        {
            if (std::isnan(ff.eact()))
            {
                printf("Energy for %s is NaN\n", actmol->getMolname().c_str());
            }
            else
            {
                // TODO Double check if the atomizationEnergy is needed.
                // auto enerexp = actmol->atomizationEnergy() + ff.first;
                double mydev2 = gmx::square(ff.eqm()-ff.eact());
                if (mydev2 == 0)
                {
                    printf("Energy difference exactly zero for %s. Ref ener %g\n",
                           actmol->getMolname().c_str(), ff.eqm());
                }
                te->second.increase(1, mydev2);
            }
        }
    }
    auto ti = targets->find(eRMS::Interaction);
    if (ti != targets->end() && !interactionEnergyMap.empty())
    {
        for(const auto &ff : interactionEnergyMap)
        {
            ti->second.increase(1, gmx::square(ff.first-ff.second.find(InteractionType::EPOT)->second));
        }
    }
}

/* * * * * * * * * * * * * * * * * * * * * *
* END: ForceEnergyDevComputer              *
* * * * * * * * * * * * * * * * * * * * * */

} // namespace alexandria

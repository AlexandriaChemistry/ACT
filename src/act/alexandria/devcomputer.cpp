/*
 * This source file is part of the Alexandria Chemistry Toolkit.
 *
 * Copyright (C) 2020-2025
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
 * \author Oskar Tegby <oskar.tegby@it.uu.se>
 */
#include "devcomputer.h"

#include <map>
#include <numeric>
#include <string>
#include <vector>

#include "act/basics/identifier.h"
#include "act/forcefield/forcefield_parametername.h"
#include "act/utility/communicationrecord.h"
#include "act/utility/units.h"
#include "gromacs/math/vecdump.h"
#include "gromacs/topology/atoms.h"
#include "gromacs/utility/real.h"
#include "gromacs/utility/textwriter.h"

// Precompute max argument for exp function at compile time.
static const double DOUBLE_MAX_LOG = std::log(GMX_DOUBLE_MAX);

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

void BoundsDevComputer::calcDeviation(MsgHandler                        *msghandler,
                                      const ForceComputer               *,
                                      ACTMol                            *,
                                      std::vector<gmx::RVec>            *,
                                      std::map<eRMS, FittingTarget>     *targets,
                                      const ForceField                  *forcefield)
{
    auto   mytarget = targets->find(eRMS::BOUNDS);
    if (targets->end() != mytarget)
    {
        double bound = 0;
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
                if (db != 0.0)
                {
                    msghandler->msg(ACTStatus::Warning,
                                    gmx::formatString("Variable %s is %g, should be within %g and %g\n",
                                                      optIndex.name().c_str(),
                                                      p.value(), p.minimum(), p.maximum()));
                }
            }
        }
        mytarget->second.increase(1, bound);
    }
    mytarget = targets->find(eRMS::UNPHYSICAL);
    if (targets->end() != mytarget)
    {
        double bound = 0;
        // Check whether shell zeta > core zeta. Only for polarizable models.
        auto   itype = InteractionType::ELECTROSTATICS;
        if (forcefield->polarizable() && forcefield->interactionPresent(itype))
        {
            auto &fs             = forcefield->findForcesConst(itype);
            std::string poltype  = "poltype";
            std::string zetatype = "zetatype";
            for(const auto &p : forcefield->particleTypesConst())
            {
                if (p.second.hasOption(poltype))
                {
                    auto coreID  = Identifier(p.second.optionValue(zetatype));
                    auto shell   = forcefield->findParticleType(p.second.optionValue(poltype));
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
            if (fs.potential() == Potential::CUBIC_BONDS)
            {
                auto cubic_name = potentialToParameterName(fs.potential());
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

void ChargeCM5DevComputer::calcDeviation(MsgHandler                    *msghandler,
                                         const ForceComputer           *,
                                         ACTMol                        *actmol,
                                         std::vector<gmx::RVec>        *,
                                         std::map<eRMS, FittingTarget> *targets,
                                         const ForceField              *forcefield)
{
    double qtot = 0;
    int i = 0;
    const auto &myatoms = actmol->atomsConst();
    std::vector<double> qcm5;
    auto qProps = actmol->qPropsConst();
    for(auto qp = qProps.begin(); qp < qProps.end(); ++qp)
    {
        auto qqm = qp->qPqmConst();
        if (qType::CM5 == qqm.qtype())
        {
            qcm5 = qqm.charge();
            if (debug)
            {
                for (size_t j = 0; j < myatoms.size(); j++)
                {
                    msghandler->msg(ACTStatus::Debug,
                                    gmx::formatString("Charge %lu. CM5 = %g ACM = %g\n",
                                                      j, qcm5[j], myatoms[j].charge()));
                }
            }
        }
    }
    // Iterate over the atoms
    for (size_t j = 0; j < myatoms.size(); j++)
    {
        if (myatoms[j].pType() == ActParticle::Shell)
        {
            continue;
        }
        auto                       atype = forcefield->findParticleType(myatoms[j].ffType());
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
                msghandler->msg(ACTStatus::Fatal,
                                gmx::formatString("Fixed charge for atom %s in %s was changed from %g to %g",
                                                    myatoms[j].name().c_str(),
                                                    actmol->getMolname().c_str(), qparm.value(), qj).c_str());
            }
            break;
        case Mutability::ACM:
        case Mutability::Bounded:
            {
                if ((*targets).find(eRMS::CHARGE)->second.weight() > 0)
                {
                    if (qparm.maximum() > qparm.minimum())
                    {
                        (*targets).find(eRMS::CHARGE)->second.increase(1, l2_regularizer(qjj, qparm.minimum(), qparm.maximum()));
                    }
                }
            }
            break;
        case Mutability::Dependent:
        case Mutability::Free:
            break;
        }
        if (!qcm5.empty() &&
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

void EspDevComputer::calcDeviation(MsgHandler                    *msghandler,
                                   const ForceComputer           *forceComputer,
                                   ACTMol                        *actmol,
                                   std::vector<gmx::RVec>        *,
                                   std::map<eRMS, FittingTarget> *targets,
                                   const ForceField              *forcefield)
{
    real rrms     = 0;
    real cosangle = 0;
    auto qProps   = actmol->qProps();
    auto topology = actmol->topology();
    std::vector<gmx::RVec> forces(topology->atoms().size());
    std::map<InteractionType, double> energies;
    bool doForce = forcefield->polarizable() || topology->hasVsites();
    // Loop over zero or more ESP data sets.
    for(auto qc = qProps->begin(); qc < qProps->end(); qc++)
    {
        QgenResp *qgr = qc->qgenResp();
        if (qgr->nEsp() == 0)
        {
            continue;
        }
        if (fitZeta_)
        {
            qgr->updateZeta(actmol->atomsConst(), forcefield);
        }
        auto coords = qgr->coords();
        qgr->updateAtomCharges(actmol->atomsConst());
        // Need to call the force routine to update shells and/or vsites
        if (doForce)
        {
            forceComputer->compute(msghandler, forcefield, topology, &coords,
                                   &forces, &energies);
            qgr->updateAtomCoords(coords);
        }
        if (msghandler->debug())
        {
            dumpQX(msghandler->tw(), actmol, "ESP");
        }
        double epsilonr;
        if (!ffOption(*forcefield, InteractionType::ELECTROSTATICS, "epsilonr", &epsilonr))
        {
            epsilonr = 1;
        }
        qgr->calcPot(msghandler, epsilonr);
        real   mae, mse;
        real   rms   = qgr->getStatistics(msghandler, &rrms, &cosangle, &mae, &mse);
        double myRms = convertToGromacs(rms, "Hartree/e");
        size_t nEsp  = qgr->nEsp();
        (*targets).find(eRMS::ESP)->second.increase(nEsp, gmx::square(myRms)*nEsp);
        if (debug)
        {
            fprintf(debug, "%s ESPrms = %g cosangle = %g\n",
                    actmol->getMolname().c_str(),
                    myRms, cosangle);
        }
    }
}

void EspDevComputer::dumpQX(gmx::TextWriter              *tw,
                            const ACTMol                 *mol,
                            const std::string            &info)
{
    if (false && tw)
    {
        std::string label = mol->getMolname() + "-" + info;
        tw->writeStringFormatted("%s q:", label.c_str());
        auto myatoms = mol->topology()->atoms();
        for (size_t i = 0; i < myatoms.size(); i++)
        {
            tw->writeStringFormatted(" %g", myatoms[i].charge());
        }
        tw->writeStringFormatted("\n");
        auto top = mol->topology();
        if (top->hasEntry(InteractionType::POLARIZATION))
        {
            tw->writeStringFormatted("%s alpha", label.c_str());
            for(auto &topentry : top->entry(InteractionType::POLARIZATION))
            {
                tw->writeStringFormatted(" %g", topentry->params()[polALPHA]);
            }
            tw->writeStringFormatted("\n");
        }
    }
}

/* * * * * * * * * * * * * * * * * * * * * *
* END: EspDevComputer                      *
* * * * * * * * * * * * * * * * * * * * * */


/* * * * * * * * * * * * * * * * * * * * * *
* BEGIN: PolarDevComputer                  *
* * * * * * * * * * * * * * * * * * * * * */

PolarDevComputer::PolarDevComputer()
    : DevComputer(mpo_name(MolPropObservable::POLARIZABILITY))
{
    convert_ = convertFromGromacs(1.0, mpo_unit2(MolPropObservable::POLARIZABILITY));
}

void PolarDevComputer::calcDeviation(MsgHandler                    *,
                                     const ForceComputer           *forceComputer,
                                     ACTMol                        *actmol,
                                     std::vector<gmx::RVec>        *,
                                     std::map<eRMS, FittingTarget> *targets,
                                     const ForceField              *forcefield)
{
    auto qProps = actmol->qProps();
    int    ndiff = 0;
    double diff2 = 0;
    for(auto qp = qProps->begin(); qp < qProps->end(); ++qp)
    {
        auto &qref = qp->qPqmConst();
        auto *qact = qp->qPact();
        if (qref.hasPolarizability())
        {
            tensor aelec;
            copy_mat(qref.polarizabilityTensor(), aelec);
            qact->setQ(actmol->topology()->atoms());
            qact->calcPolarizability(forcefield, actmol->topology(), forceComputer);
            auto acalc = qact->polarizabilityTensor();
            for(int i = 0; i < DIM; i++)
            {
                for(int j = 0; j < DIM; j++)
                {
                    diff2 += gmx::square(convert_*(aelec[i][j]-acalc[i][j]));
                }
            }
            ndiff += 1;
        }
    }
    
    (*targets).find(eRMS::Polar)->second.increase(ndiff, diff2);
}

/* * * * * * * * * * * * * * * * * * * * * *
* END: PolarDevComputer                    *
* * * * * * * * * * * * * * * * * * * * * */

/* * * * * * * * * * * * * * * * * * * * * *
* BEGIN: MultiPoleDevComputer              *
* * * * * * * * * * * * * * * * * * * * * */

void MultiPoleDevComputer::calcDeviation(MsgHandler                    *msghandler,
                                         const ForceComputer           *forceComputer,
                                         ACTMol                        *actmol,
                                         std::vector<gmx::RVec>        *,
                                         std::map<eRMS, FittingTarget> *targets,
                                         const ForceField              *forcefield)
{
    auto   qProps   = actmol->qProps();
    int    ndiff    = 0;
    double delta    = 0;
    auto   topology = actmol->topology();
    bool doForce    = forcefield->polarizable() || topology->hasVsites();
    for(auto qp = qProps->begin(); qp < qProps->end(); ++qp)
    {
        auto qqm  = qp->qPqmConst();
        auto qact = qp->qPactConst();
        if (qqm.hasMultipole(mpo_) && qact.hasMultipole(mpo_))
        {
            auto qelec = qqm.getMultipole(mpo_);
            // TODO: Compute this only once if both dipole and quadrupole are used in fitting
            qact.setQ(*actmol->atoms());
            if (doForce)
            {
                std::vector<gmx::RVec>            forces;
                std::map<InteractionType, double> energies;
                auto                              myx = qact.x();
                forceComputer->compute(msghandler, forcefield, topology, &myx, &forces, &energies);
                qact.setX(myx);
            }
            qact.calcMoments();
            auto qcalc = qact.getMultipole(mpo_);
            for (size_t mm = 0; mm < qelec.size(); mm++)
            {
                delta += gmx::square(qcalc[mm] - qelec[mm]);
            }
            ndiff += 1;
        }
    }
    // Not correct, but just to make compiler happy.
    eRMS rms = eRMS::TOT;
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
    default: // throws
        msghandler->fatal(gmx::formatString("Not support MolPropObservable %s", mpo_name(mpo_)).c_str());
    }
    
    (*targets).find(rms)->second.increase(ndiff, delta);
}

/* * * * * * * * * * * * * * * * * * * * * *
* END: MultipoleDevComputer                *
* * * * * * * * * * * * * * * * * * * * * */

/* * * * * * * * * * * * * * * * * * * * * *
* BEGIN: HarmonicsDevComputer              *
* * * * * * * * * * * * * * * * * * * * * */

HarmonicsDevComputer::HarmonicsDevComputer(MolPropObservable  mpo)
    : DevComputer(mpo_name(mpo)), mpo_(mpo) 
{
}

void HarmonicsDevComputer::calcDeviation(MsgHandler                    *msghandler,
                                         const ForceComputer           *forceComputer,
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
    auto eMin = handler_.minimizeCoordinates(msghandler, forcefield, actmol, forceComputer,
                                             simConfig_, coords, nullptr, {});
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
    default: // throws
        GMX_THROW(gmx::InternalError(gmx::formatString("Don't know how to handle %s in this devcomputer\n", mpo_name(mpo_)).c_str()));
    }
}

/* * * * * * * * * * * * * * * * * * * * * *
* END: HarmonicsDevComputer                *
* * * * * * * * * * * * * * * * * * * * * */


/* * * * * * * * * * * * * * * * * * * * * *
* BEGIN: ForceEnergyDevComputer                 *
* * * * * * * * * * * * * * * * * * * * * */

ForceEnergyDevComputer::ForceEnergyDevComputer(std::map<eRMS, double>  boltzmannTemperature)
    : DevComputer("ForceEnergy")
{
    boltzmannTemperature_ = boltzmannTemperature;
}

double ForceEnergyDevComputer::computeBeta(eRMS ermsi)
{
    double beta = 0;
    if (boltzmannTemperature_.find(ermsi) != boltzmannTemperature_.end())
    {
        double T = boltzmannTemperature_[ermsi];
        if (T > 0)
        {
            beta = 1.0/(BOLTZ*T);
        }
    }
    return beta;
}

void ForceEnergyDevComputer::calcDeviation(MsgHandler                    *msghandler,
                                           const ForceComputer           *forceComputer,
                                           ACTMol                        *actmol,
                                           std::vector<gmx::RVec>        *,
                                           std::map<eRMS, FittingTarget> *targets,
                                           const ForceField              *forcefield)
{
    std::vector<ACTEnergy>                                              energyMap;
    std::vector<std::vector<std::pair<double, double> > >               forceMap;
    std::vector<std::pair<double, std::map<InteractionType, double> > > enerComponentMap;
    ACTEnergyMapVector                                                  interactionEnergyMap;
    // Check what we need to do first!
    auto ermse    = eRMS::EPOT;
    auto te       = targets->find(ermse);
    bool doEpot   = (te != targets->end() && te->second.weight() > 0);
    auto tf       = targets->find(eRMS::Force2);
    bool doForce2 = (tf != targets->end() && tf->second.weight() > 0);
    bool                            doInter = false;
    std::map<eRMS, InteractionType> rmsE    = {
        { eRMS::Interaction,    InteractionType::EPOT         },
        { eRMS::Electrostatics, InteractionType::ELECTROSTATICS      },
        { eRMS::Dispersion,     InteractionType::DISPERSION   },
        { eRMS::Exchange,       InteractionType::EXCHANGE     },
        { eRMS::Induction,      InteractionType::INDUCTION    },
        { eRMS::DeltaHF,        InteractionType::INDUCTIONCORRECTION },
        { eRMS::ExchInd,        InteractionType::EXCHIND      },
        { eRMS::AllElec,        InteractionType::ALLELEC      }
    };
    for (auto &rms : rmsE)
    {
        auto ttt = targets->find(rms.first);
        if (ttt != targets->end() && ttt->second.weight() > 0)
        {
            doInter = true;
        }
    }
    if (doForce2 || doEpot || doInter)
    {
        actmol->forceEnergyMaps(msghandler, forcefield,
                                forceComputer, &forceMap, &energyMap,
                                &interactionEnergyMap, &enerComponentMap, separateInductionCorrection_);
        if (doForce2 && !forceMap.empty())
        {
            for(const auto &fstruct : forceMap)
            {
                for(const auto &ff : fstruct)
                {
                    if (!std::isfinite(ff.second))
                    {
                        msghandler->msg(ACTStatus::Warning,
                                        gmx::formatString("Force for %s is NaN\n", actmol->getMolname().c_str()));
                    }
                    else
                    {
                        tf->second.increase(1, gmx::square(ff.first-ff.second));
                    }
                }
            }
        }
        if (doEpot && !energyMap.empty())
        {
            auto beta = computeBeta(ermse);
            double eqmMin = 1e8;
            if (beta > 0)
            {
                for(const auto &ff : energyMap)
                {
                    if (ff.haveQM())
                    {
                        eqmMin = std::min(eqmMin, ff.eqm());
                    }
                }
            }
            for(const auto &ff : energyMap)
            {
                // std::isfinite checks for both Inf and NaN values.
                if (ff.haveACT())
                {
                    if (!std::isfinite(ff.eact()))
                    {
                        msghandler->msg(ACTStatus::Warning,
                                        gmx::formatString("Energy for %s is NaN\n", actmol->getMolname().c_str()));
                    }
                    else
                    {
                        if (ff.haveQM())
                        {
                            double eqm    = ff.eqm();
                            double mydev2 = gmx::square(eqm-ff.eact());
                            double weight = 1;
                            if (beta > 0)
                            {
                                weight = exp(-beta*(eqm-eqmMin));
                            }
                            te->second.increase(weight, mydev2);
                        }
                    }
                }
            }
        }
        if (doInter)
        {
            // We can only use weighting on the total interaction energy not on
            // components.
            auto   rms    = eRMS::Interaction;
            double beta   = computeBeta(rms);
            double eqmMin = 1e8;
            if (beta > 0)
            {
                for(auto &iem : interactionEnergyMap)
                {
                    auto ff = iem.find(InteractionType::EPOT);
                    if (iem.end() != ff && ff->second.haveQM())
                    {
                        eqmMin = std::min(eqmMin, ff->second.eqm());
                    }
                }
            }
            // Check whether we should sum the induction terms
            bool sumInductionTerms = false;
            auto tdhf = targets->find(eRMS::DeltaHF);
            if (tdhf == targets->end() || tdhf->second.weight() == 0)
            {
                // No -fc_deltahf option was passed
                if (targets->find(eRMS::Induction) != targets->end())
                {
                    sumInductionTerms = true;
                }
                std::map<bool, const char *> bool_names = { { false, "False" }, { true, "True" } };
                msghandler->msg(ACTStatus::Debug, gmx::formatString("Sum induction terms: %s", bool_names[sumInductionTerms]));
            }
            for(auto &iem : interactionEnergyMap)
            {
                double weight = 1;
                if (beta > 0)
                {
                    auto ff = iem.find(InteractionType::EPOT);
                    if (iem.end() != ff && ff->second.haveQM())
                    {
                        auto eqm =  ff->second.eqm();
                        if (beta > 0)
                        {
                            // Argument should always be positive
                            double earg = beta*(eqm-eqmMin);
                            // Check whether exponentiating will work
                            if (earg < DOUBLE_MAX_LOG)
                            {
                                weight = exp(-earg);
                            }
                            else
                            {
                                // Otherwise the weight is (close to) zero.
                                // Totally zero will not work due to checks later.
                                weight = 1e-8;
                            }
                        }
                    }
                }
                for(const auto &rms: rmsE)
                {
                    auto ti = targets->find(rms.first);
                    if (ti == targets->end() || ti->second.weight() == 0)
                    {
                        continue;
                    }
                    auto myff = iem.find(rms.second);
                    if (myff != iem.end())
                    {
                        auto &ff = myff->second;
                        if (ff.haveQM() && ff.haveACT())
                        {
                            auto eqm  = ff.eqm();
                            auto eact = ff.eact();
                            if (std::isfinite(eact))
                            {
                                ti->second.increase(weight, gmx::square(eqm-eact));
                            }
                            else
                            {
                                // We do not want to deal with infinite numbers but it should
                                // be clear that this is a very bad parameter set.
                                eact = 1e16;
                                ti->second.increase(weight, gmx::square(eqm-eact));
                            }
                        }
                    }
                }
            }
        }
    }
}

/* * * * * * * * * * * * * * * * * * * * * *
* END: ForceEnergyDevComputer              *
* * * * * * * * * * * * * * * * * * * * * */

} // namespace alexandria

/*
 * This source file is part of the Alexandria Chemistry Toolkit.
 *
 * Copyright (C) 2020-2026
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

#include <algorithm>
#include <map>
#include <numeric>
#include <string>
#include <vector>

#include "act/basics/identifier.h"
#include "act/forcefield/forcefield_parametername.h"
#include "act/forces/forces_simd.h" 
#include "act/utility/communicationrecord.h"
#include "act/utility/units.h"
#include "gromacs/math/vecdump.h"
#include "gromacs/utility/textwriter.h"

// Precompute max argument for exp function at compile time.
static const double DOUBLE_MAX_LOG = std::log(GMX_DOUBLE_MAX);

namespace alexandria
{

/* * * * * * * * * * * * * * * * * * * * * *
* BEGIN: BoundsDevComputer                 *
* * * * * * * * * * * * * * * * * * * * * */

void BoundsDevComputer::calcDeviation(MsgHandler                    *,
                                      const ForceComputer           *,
                                      std::vector<ACTMol>           *actmol,
                                      std::map<eRMS, FittingTarget> *targets,
                                      const ForceField              *forcefield,
                                      LossFunction                   lossFunction)
{
    auto   mytarget = targets->find(eRMS::BOUNDS);
    if (targets->end() == mytarget)
    {
        return;
    }
    for(auto mol = actmol->begin(); mol < actmol->end(); ++mol)
    {
        double bound = 0;
        // Check ACM charges
        auto atoms = mol->atomsConst();
        for(const auto &atom: atoms)
        {
            auto pt = atom.ffType();
            auto p  = forcefield->findParticleType(pt)->parameterConst("charge");
            if (p.mutability() == Mutability::ACM && p.maximum() > p.minimum())
            {
                double q = atom.charge();
                // Add shell charge(s)
                for(auto ss : atom.shells())
                {
                    q += atoms[ss].charge();
                }
                // Add vsite 1 charge
                for(auto vv : atom.vsites())
                {
                    if (atoms[vv].parents().size() == 1)
                    {
                        q += atoms[vv].charge();
                }
                }
                if (q < p.minimum())
                {
                    bound += loss(lossFunction, 1, p.minimum() - q);
                }
                else if (q > p.maximum())
                {
                    bound += loss(lossFunction, 1, q - p.maximum());
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
* BEGIN: EspDevComputer                    *
* * * * * * * * * * * * * * * * * * * * * */

void EspDevComputer::calcDeviation(MsgHandler                    *msghandler,
                                   const ForceComputer           *forceComputer,
                                   std::vector<ACTMol>           *actmols,
                                   std::map<eRMS, FittingTarget> *targets,
                                   const ForceField              *forcefield,
                                   gmx_unused LossFunction        lossFunction)
{
    for(auto actmol = actmols->begin(); actmol < actmols->end(); ++actmol)
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
            qgr->updateAtomCharges(msghandler, actmol->atomsConst());
            // Need to call the force routine to update shells and/or vsites
            if (doForce)
            {
                forceComputer->compute(msghandler, forcefield, topology, &coords,
                                       &forces, &energies);
                qgr->updateAtomCoords(coords);
            }
            if (msghandler->debug())
            {
                dumpQX(msghandler->tw(), &(*actmol), "ESP");
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
            if (msghandler->debug())
            {
                msghandler->writeDebug(gmx::formatString("%s ESPrms = %g cosangle = %g\n",
                                                         actmol->getMolname().c_str(),
                                                         myRms, cosangle));
            }
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

void PolarDevComputer::calcDeviation(MsgHandler                    *msghandler,
                                     const ForceComputer           *forceComputer,
                                     std::vector<ACTMol>           *actmols,
                                     std::map<eRMS, FittingTarget> *targets,
                                     const ForceField              *forcefield,
                                     LossFunction                   lossFunction)
{
    for(auto actmol = actmols->begin(); actmol < actmols->end(); ++actmol)
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
                qact->calcPolarizability(msghandler, forcefield, actmol->topology(), forceComputer);
                auto acalc = qact->polarizabilityTensor();
                for(int i = 0; i < DIM; i++)
                {
                    for(int j = 0; j < DIM; j++)
                    {
                        diff2 += loss(lossFunction, convert_, (aelec[i][j]-acalc[i][j]));
                    }
                }
                ndiff += 1;
            }
        }
        (*targets).find(eRMS::Polar)->second.increase(ndiff, diff2);
    }
}

/* * * * * * * * * * * * * * * * * * * * * *
* END: PolarDevComputer                    *
* * * * * * * * * * * * * * * * * * * * * */

/* * * * * * * * * * * * * * * * * * * * * *
* BEGIN: MultiPoleDevComputer              *
* * * * * * * * * * * * * * * * * * * * * */

void MultiPoleDevComputer::calcDeviation(MsgHandler                    *msghandler,
                                         const ForceComputer           *forceComputer,
                                         std::vector<ACTMol>           *actmols,
                                         std::map<eRMS, FittingTarget> *targets,
                                         const ForceField              *forcefield,
                                         LossFunction                   lossFunction)
{
    for(auto actmol = actmols->begin(); actmol < actmols->end(); ++actmol)
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
                auto qelec = qqm.getMultipole(msghandler, mpo_);
                //! \todo Compute this only once if both dipole and quadrupole are used in fitting
                qact.setQ(*actmol->atoms());
                if (doForce)
                {
                    std::vector<gmx::RVec>            forces(actmol->atomsConst().size());
                    std::map<InteractionType, double> energies;
                    auto                              myx = qact.x();
                    forceComputer->compute(msghandler, forcefield, topology, &myx, &forces, &energies);
                    qact.setX(myx);
                }
                qact.calcMoments();
                auto qcalc = qact.getMultipole(msghandler, mpo_);
                for (size_t mm = 0; mm < qelec.size(); mm++)
                {
                    delta += loss(lossFunction, 1, qcalc[mm] - qelec[mm]);
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
                                         std::vector<ACTMol>           *actmols,
                                         std::map<eRMS, FittingTarget> *targets,
                                         const ForceField              *forcefield,
                                         LossFunction                   lossFunction)
{
    for(auto actmol = actmols->begin(); actmol < actmols->end(); ++actmol)
    {
        // Only compute frequencies for structures that have an optimize reference
        if (JobType::OPT != actmol->jobType())
        {
            continue;
        }
        auto coords = actmol->xOriginal();
        auto eMin   = handler_.minimizeCoordinates(msghandler, forcefield, actmol->topology(), forceComputer,
                                                   simConfig_, &coords, nullptr, {});
        if (eMinimizeStatus::OK != eMin)
        {
            // Something fishy happened, but it means we cannot use this structure
            // for computing frequencies.
            return;
        }
        // Compute frequencies
        std::vector<double> frequencies, intensities;
        handler_.nma(forcefield, &(*actmol), forceComputer, &coords, &frequencies, &intensities);

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
                        delta += loss(lossFunction, 1, frequencies[k]-ref_freqs[k]);
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
                    fprintf(stderr, "Reference intensities size %zu, but calculated %zu for %s. Ignoring intensities for this compound.\n",
                            ref_intens.size(), intensities.size(), actmol->getMolname().c_str());
                }
                else
                {
                    double delta = 0;
                    for(size_t k = 0; k < intensities.size(); k++)
                    {
                        delta += loss(lossFunction, 1, intensities[k]-ref_intens[k]);
                    }
                    (*targets).find(eRMS::INTENSITY)->second.increase(1, delta);
                }
            }
            break;
        default: // throws
            GMX_THROW(gmx::InternalError(gmx::formatString("Don't know how to handle %s in this devcomputer\n", mpo_name(mpo_)).c_str()));
        }
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

/*! \brief Compute loss due to forces
 * \param[in]  coords       The data including fcalc and fobs
 * \param[in]  lossFunction The loss function to use
 * \param[out] weight       The total weight corresponding to the loss
 * \param[out] loss         The aggregate loss
 */
static inline void calc_force_loss(const std::vector<CoordsEnergy> *coords,
                                   LossFunction                     lossFunction,
                                   double                          *weight,
                                   double                          *loss2)
{
    std::vector<double> weights(coords->size(), 0);
    std::vector<double> losses(coords->size(), 0);
#pragma omp simd
    for(size_t ix = 0; ix < coords->size(); ix++)
    {
        auto &x = (*coords)[ix];
        if (x.fcalc.size() == x.fobs.size())
        {
            double www  = 0;
            double lss = 0;
            for(size_t i = 0; i < x.fcalc.size(); i++)
            {
#pragma unroll
                for(int m = 0; m < DIM; m++)
                {
                    if (std::isfinite(x.fcalc[i][m]))
                    {
                        lss += loss(lossFunction, 1, x.fcalc[i][m] - x.fobs[i][m]);
                        www += 1;
                    }
                }
            }
            weights[ix] = www;
            losses[ix]  = lss;
        }
    }
    // Gather results, the initial value must be a double to prevent rounding
    *weight = std::accumulate(weights.begin(), weights.end(), 0.0);
    *loss2  = std::accumulate(losses.begin(), losses.end(), 0.0);
}

/*! \brief Find the QM reference energy for an interaction energy term
 * Some terms combine several QM observables. The induction correction
 * may be missing if it is not treated separately, in which case it is also
 * added to the induction.
 * \param[in]  eMap         The energies
 * \param[in]  itype        The term
 * \param[in]  sepIndCorr   Whether the induction correction is separate
 * \param[out] eqm          The QM reference
 * \return whether there is a QM reference
 */
static bool qmReference(const ACTEnergyMap &eMap,
                        InteractionType     itype,
                        bool                sepIndCorr,
                        double             *eqm)
{
    auto itICorr = InteractionType::INDUCTIONCORRECTION;
    std::vector<InteractionType> terms;
    switch (itype)
    {
    case InteractionType::ALLELEC:
        terms = { InteractionType::ELECTROSTATICS, InteractionType::INDUCTION, itICorr };
        break;
    case InteractionType::EXCHIND:
        terms = { InteractionType::EXCHANGE, InteractionType::INDUCTION, itICorr };
        break;
    case InteractionType::INDUCTION:
        if (sepIndCorr)
        {
            terms = { itype };
        }
        else
        {
            // The ACT induction includes the induction correction
            terms = { itype, itICorr };
        }
        break;
    default:
        terms = { itype };
    }
    *eqm         = 0;
    size_t found = 0;
    for(auto t : terms)
    {
        auto ee = eMap.find(t);
        if (eMap.end() != ee && ee->second.haveQM())
        {
            *eqm  += ee->second.eqm();
            found += 1;
        }
        else if (t == itICorr && !sepIndCorr && terms.size() > 1)
        {
            // Allowed to be missing
            found += 1;
        }
    }
    return found == terms.size();
}

void ForceEnergyDevComputer::calcDeviation(MsgHandler                     *msghandler,
                                           gmx_unused const ForceComputer *forceComputer,
                                           std::vector<ACTMol>            *actmols,
                                           std::map<eRMS, FittingTarget>  *targets,
                                           const ForceField               *forcefield,
                                           LossFunction                    lossFunction)
{
    // Check what we need to do first!
    auto ermse    = eRMS::EPOT;
    auto te       = targets->find(ermse);
    bool doEpot   = (te != targets->end() && te->second.weight() > 0);
    auto tf       = targets->find(eRMS::Force2);
    bool doForce2 = (tf != targets->end() && tf->second.weight() > 0);
    bool doInter  = false;
    typedef struct {
        eRMS            rms;
        InteractionType itype;
        double          weight;
        double          loss;
    } rmsInter;
    std::vector<rmsInter> rmsE    = {
        { eRMS::Interaction,    InteractionType::EPOT                },
        { eRMS::Electrostatics, InteractionType::ELECTROSTATICS      },
        { eRMS::Dispersion,     InteractionType::DISPERSION          },
        { eRMS::Exchange,       InteractionType::EXCHANGE            },
        { eRMS::Induction,      InteractionType::INDUCTION           },
        { eRMS::DeltaHF,        InteractionType::INDUCTIONCORRECTION },
        { eRMS::ExchInd,        InteractionType::EXCHIND             },
        { eRMS::AllElec,        InteractionType::ALLELEC             }
    };
    // Remove the terms without a target or with zero weight
    rmsE.erase(std::remove_if(rmsE.begin(), rmsE.end(),
                              [targets](const rmsInter &r)
                              {
                                  auto ttt = targets->find(r.rms);
                                  return ttt == targets->end() || ttt->second.weight() <= 0;
                              }), rmsE.end());
    doInter = !rmsE.empty();
    // Whether to treat the induction correction separately is set by the
    // caller, see ACMFitnessComputer::fillDevComputers.
    bool sepInductionTerms = separateInductionCorrection_;
    if (doForce2 || doEpot || doInter)
    {
        // Compute all energies at once, and the forces if needed
        forceComputer->computeRelaxedEnergies(msghandler, forcefield, actmols, true, sepInductionTerms, doForce2);
        for(auto actmol = actmols->begin(); actmol < actmols->end(); ++actmol)
        {
            if (doForce2)
            {
                // Forces can be from dimers or monomers
                double weight = 0, loss = 0;
                if (actmol->isDimer())
                {
                    calc_force_loss(actmol->dimerCoords(), lossFunction, &weight, &loss);
                }
                else
                {
                    calc_force_loss(actmol->monomerCoords(0), lossFunction, &weight, &loss);
                }
                tf->second.increase(weight, loss);
            }
            if (doInter)
            {
                // Store weights and losses per interaction term
                for (auto iter = rmsE.begin(); iter < rmsE.end(); ++iter)
                {
                    iter->weight = 0;
                    iter->loss   = 0;
                }
                double     beta   = computeBeta(eRMS::Interaction);
                const auto coords = actmol->dimerCoords();
#pragma omp simd
                for(size_t i = 0; i < coords->size(); i++)
                {
                    const auto &eMap = (*coords)[i].eMap;
                    double w         = 1;
                    auto   dd        = eMap.find(InteractionType::EPOT);
                    if (dd != eMap.end() && beta > 0)
                    {
                        double earg = beta*(dd->second.eqm() - (*coords)[i].eMin);
                        if (earg < DOUBLE_MAX_LOG)
                        {
                            w = std::exp(-earg);
                        }
                        else
                        {
                            // Otherwise the weight is (close to) zero.
                            // Totally zero will not work due to checks later.
                            w = 1e-8;
                        }
                    }
                    for (auto iter = rmsE.begin(); iter != rmsE.end(); ++iter)
                    {
                        auto   dd  = eMap.find(iter->itype);
                        double eqm = 0;
                        if (dd != eMap.end() && dd->second.haveACT() &&
                            qmReference(eMap, iter->itype, sepInductionTerms, &eqm))
                        {
                            if (std::isfinite(dd->second.eact()))
                            {
                                iter->loss   += w * loss(lossFunction, 1, eqm - dd->second.eact());
                                iter->weight += w;
                            }
                        }
                    }
                }
                // Collect losses, maybe they need to be per dimer conf to allow SIMD above
                for (auto iter = rmsE.begin(); iter != rmsE.end(); ++iter)
                {
                    targets->find(iter->rms)->second.increase(iter->weight, iter->loss);
                }    
            }
            if (doEpot)
            {
                double     beta   = computeBeta(eRMS::EPOT);
                const auto coords = actmol->monomerCoords(0);
#pragma omp simd
                for(size_t i = 0; i < coords->size(); i++)
                {
                    const auto &eMap = (*coords)[i].eMap;
                    double w         = 1;
                    auto   dd        = eMap.find(InteractionType::EPOT);
                    if (dd != eMap.end() && dd->second.haveQM())
                    {
                        double eqm  = dd->second.eqm();
                        if (beta > 0)
                        {
                            double earg = beta*(eqm - (*coords)[i].eMin);
                            if (earg < DOUBLE_MAX_LOG)
                            {
                                w = std::exp(-earg);
                            }
                            else
                            {
                                // Otherwise the weight is (close to) zero.
                                // Totally zero will not work due to checks later.
                                w = 1e-8;
                            }
                        }
                        // std::isfinite checks for both Inf and NaN values.
                        if (dd->second.haveACT() && std::isfinite(dd->second.eact()))
                        {
                            double myloss = w * loss(lossFunction, w, eqm-dd->second.eact());
                            te->second.increase(w, myloss);
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

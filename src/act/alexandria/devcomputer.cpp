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
                                      const std::vector<ACTMol *>   &actmols,
                                      std::map<eRMS, FittingTarget> *targets,
                                      const ForceField              *forcefield,
                                      LossFunction                   lossFunction)
{
    const auto mytarget = targets->find(eRMS::BOUNDS);
    if (targets->end() == mytarget)
    {
        return;
    }
    for(auto mol : actmols)
    {
        mytarget->second.increase(1, moleculeBound(mol, cachedBounds(mol, forcefield), lossFunction));
    }
}

std::vector<BoundsDevComputer::ChargeBound> BoundsDevComputer::chargeBounds(const ACTMol     *mol,
                                                                            const ForceField *forcefield)
{
    std::vector<ChargeBound> bounds;
    const auto &atoms = mol->atomsConst();
    for(size_t i = 0; i < atoms.size(); i++)
    {
        const auto &atom = atoms[i];
        const auto &p    = forcefield->findParticleType(atom.ffType())->parameterConst("charge");
        if (p.mutability() == Mutability::ACM && p.maximum() > p.minimum())
        {
            ChargeBound cb = { static_cast<int>(i), p.minimum(), p.maximum(), {} };
            // The shell charge(s)
            for(auto ss : atom.shells())
            {
                cb.extra.push_back(ss);
            }
            // The charge of the virtual sites of this atom alone
            for(auto vv : atom.vsites())
            {
                if (atoms[vv].parents().size() == 1)
                {
                    cb.extra.push_back(vv);
                }
            }
            bounds.push_back(std::move(cb));
        }
    }
    return bounds;
}

double BoundsDevComputer::moleculeBound(const ACTMol                   *mol,
                                        const std::vector<ChargeBound> &bounds,
                                        LossFunction                    lossFunction)
{
    double bound = 0;
    const auto &atoms = mol->atomsConst();
    for(const auto &cb : bounds)
    {
        // Including the second charge, if any
        double q = atoms[cb.atom].totalCharge();
        for(auto e : cb.extra)
        {
            q += atoms[e].charge();
        }
        if (q < cb.minimum)
        {
            bound += loss(lossFunction, 1, cb.minimum - q);
        }
        else if (q > cb.maximum)
        {
            bound += loss(lossFunction, 1, q - cb.maximum);
        }
    }
    return bound;
}

double BoundsDevComputer::moleculeBound(const ACTMol     *mol,
                                        const ForceField *forcefield,
                                        LossFunction      lossFunction)
{
    return moleculeBound(mol, chargeBounds(mol, forcefield), lossFunction);
}

const std::vector<BoundsDevComputer::ChargeBound> &BoundsDevComputer::cachedBounds(const ACTMol     *mol,
                                                                                   const ForceField *forcefield)
{
    if (forcefield != boundsForceField_)
    {
        bounds_.clear();
        boundsForceField_ = forcefield;
    }
    auto b = bounds_.find(mol);
    if (b == bounds_.end())
    {
        b = bounds_.emplace(mol, chargeBounds(mol, forcefield)).first;
    }
    return b->second;
}

/* * * * * * * * * * * * * * * * * * * * * *
* END: BoundsDevComputer                   *
* * * * * * * * * * * * * * * * * * * * * */

/* * * * * * * * * * * * * * * * * * * * * *
* BEGIN: EspDevComputer                    *
* * * * * * * * * * * * * * * * * * * * * */

void EspDevComputer::calcDeviation(MsgHandler                    *msghandler,
                                   const ForceComputer           *forceComputer,
                                   const std::vector<ACTMol *>   &actmols,
                                   std::map<eRMS, FittingTarget> *targets,
                                   const ForceField              *forcefield,
                                   gmx_unused LossFunction        lossFunction)
{
    for(auto actmol : actmols)
    {
        double rrms     = 0;
        double cosangle = 0;
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
            auto coords = qgr->coords();
            qgr->updateAtomCharges(msghandler, actmol->atomsConst());
            // Need to call the force routine to update shells and/or vsites
            if (doForce)
            {
                forceComputer->compute(msghandler, forcefield, topology, &coords,
                                       &forces, &energies);
                qgr->updateAtomCoords(coords);
            }
            double epsilonr;
            epsilonr = forcefield->epsilonr();
            // Refresh both ordinary and second-charge widths, independently of the fit selection.
            qgr->updateZeta(actmol->atomsConst(), forcefield);
            qgr->calcPot(msghandler, epsilonr);
            double mae, mse;
            double rms   = qgr->getStatistics(msghandler, &rrms, &cosangle, &mae, &mse);
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
                                     const std::vector<ACTMol *>   &actmols,
                                     std::map<eRMS, FittingTarget> *targets,
                                     const ForceField              *forcefield,
                                     LossFunction                   lossFunction)
{
    for(auto actmol : actmols)
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
                        // The deviation in the unit of the polarizability
                        diff2 += loss(lossFunction, 1, convert_*(aelec[i][j]-acalc[i][j]));
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
                                         const std::vector<ACTMol *>   &actmols,
                                         std::map<eRMS, FittingTarget> *targets,
                                         const ForceField              *forcefield,
                                         LossFunction                   lossFunction)
{
    for(auto actmol : actmols)
    {
        auto   qProps   = actmol->qProps();
        int    ndiff    = 0;
        double delta    = 0;
        auto   topology = actmol->topology();
        bool doForce    = forcefield->polarizable() || topology->hasVsites();
        for(auto qp = qProps->begin(); qp < qProps->end(); ++qp)
        {
            const auto qqm  = qp->qPqmConst();
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
            msghandler->fatal(gmx::formatString("Not support MolPropObservable %s", mpo_name(mpo_)));
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
                                         const std::vector<ACTMol *>   &actmols,
                                         std::map<eRMS, FittingTarget> *targets,
                                         const ForceField              *forcefield,
                                         LossFunction                   lossFunction)
{
    for(auto actmol : actmols)
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
            continue;
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
 * The reference forces are those of the atoms only, the computed forces
 * are those of all particles, therefore the particles that are not atoms,
 * such as shells and virtual sites, are skipped.
 * \param[in]  msghandler   For error messages
 * \param[in]  molname      Name of the molecule, for error messages
 * \param[in]  atoms        The particles of the topology of the conformations
 * \param[in]  cdata        The conformations, with the computed and reference forces
 * \param[in]  lossFunction The loss function to use
 * \param[out] weight       The total weight corresponding to the loss
 * \param[out] loss2        The aggregate loss
 */
static void calc_force_loss(MsgHandler                      *msghandler,
                            const std::string               &molname,
                            const std::vector<ActAtom>      &atoms,
                            const ConformationData          *cdata,
                            LossFunction                     lossFunction,
                            double                          *weight,
                            double                          *loss2)
{
    std::vector<size_t> atomIndex;
    for(size_t i = 0; i < atoms.size(); i++)
    {
        if (atoms[i].pType() == ActParticle::Atom)
        {
            atomIndex.push_back(i);
        }
    }
    std::vector<double> weights(cdata->nconf(), 0);
    std::vector<double> losses(cdata->nconf(), 0);
    for(size_t ix = 0; ix < cdata->nconf(); ix++)
    {
        auto &x  = cdata->confs[ix];
        const auto fx = cdata->f(ix);
        // Not all conformations need to have reference forces
        if (x.fobs.empty())
        {
            continue;
        }
        if (x.fobs.size() != atomIndex.size() || fx.size() != atoms.size())
        {
            msghandler->msg(ACTStatus::Fatal,
                            gmx::formatString("Inconsistency for %s: there are %zu atoms and %zu particles, but %zu reference and %zu computed forces",
                                              molname.c_str(), atomIndex.size(), atoms.size(),
                                              x.fobs.size(), fx.size()));
            return;
        }
        double www = 0;
        double lss = 0;
        for(size_t i = 0; i < atomIndex.size(); i++)
        {
            const auto fc = fx[atomIndex[i]];
#pragma unroll
            for(int m = 0; m < DIM; m++)
            {
                if (std::isfinite(fc[m]))
                {
                    lss += loss(lossFunction, 1, fc[m] - x.fobs[i][m]);
                    www += 1;
                }
            }
        }
        weights[ix] = www;
        losses[ix]  = lss;
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
        const auto ee = eMap.find(t);
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

/*! \brief Whether the reference data was made for these conformations and settings
 * \param[in] ref    The reference data
 * \param[in] cdata  The conformation data it belongs to
 * \param[in] nconf  Number of conformations
 * \param[in] itypes The interaction types
 * \param[in] beta   The inverse Boltzmann temperature
 * \param[in] sep    Whether the induction correction is separate
 */
static bool referenceValid(const DeviationReference         &ref,
                           const ConformationData           &cdata,
                           size_t                            nconf,
                           const std::vector<InteractionType> &itypes,
                           double                            beta,
                           bool                              sep)
{
    return (ref.conformationVersion == cdata.conformations.value() && ref.beta == beta &&
            ref.separate == sep && ref.itypes == itypes && ref.weight.size() == nconf);
}

/*! \brief The reference interaction energies and Boltzmann weights of dimers
 * They are computed once for the conformations, since they only depend on
 * the quantum chemistry data, not on the positions of virtual sites.
 * \param[inout] cdata  The conformations of the dimer with their data
 * \param[in]    itypes The interaction types
 * \param[in]    beta   The inverse Boltzmann temperature, no weighting if zero
 * \param[in]    sep    Whether the induction correction is separate
 * \return The reference data
 */
static const DeviationReference &interactionReference(ConformationData                   *cdata,
                                                      const std::vector<InteractionType> &itypes,
                                                      double                              beta,
                                                      bool                                sep)
{
    const auto &coords = cdata->confs;
    auto       &ref    = cdata->deviation;
    if (referenceValid(ref, *cdata, coords.size(), itypes, beta, sep))
    {
        return ref;
    }
    ref.conformationVersion = cdata->conformations.value();
    ref.beta              = beta;
    ref.separate          = sep;
    ref.itypes            = itypes;
    ref.weight.assign(coords.size(), 1.0);
    ref.energy.assign(itypes.size(), std::vector<double>(coords.size(), 0.0));
    ref.have.assign(itypes.size(), std::vector<char>(coords.size(), 0));
    for(size_t i = 0; i < coords.size(); i++)
    {
        const auto &eMap = coords[i].eMap;
        const auto  dd   = eMap.find(InteractionType::EPOT);
        if (dd != eMap.end() && beta > 0)
        {
            double earg = beta*(dd->second.eqm() - coords[i].eMin);
            // Otherwise the weight is (close to) zero.
            // Totally zero will not work due to checks later.
            ref.weight[i] = earg < DOUBLE_MAX_LOG ? std::exp(-earg) : 1e-8;
        }
        for(size_t j = 0; j < itypes.size(); j++)
        {
            ref.have[j][i] = qmReference(eMap, itypes[j], sep, &ref.energy[j][i]);
        }
    }
    return ref;
}

/*! \brief The reference potential energies and Boltzmann weights of monomers
 * \param[inout] cdata  The conformations of the monomer with their data
 * \param[in]    beta   The inverse Boltzmann temperature, no weighting if zero
 * \return The reference data, with only EPOT
 */
static const DeviationReference &epotReference(ConformationData *cdata,
                                               double            beta)
{
    const std::vector<InteractionType> itypes = { InteractionType::EPOT };
    const auto &coords = cdata->confs;
    auto       &ref    = cdata->deviation;
    if (referenceValid(ref, *cdata, coords.size(), itypes, beta, false))
    {
        return ref;
    }
    ref.conformationVersion = cdata->conformations.value();
    ref.beta              = beta;
    ref.separate          = false;
    ref.itypes            = itypes;
    ref.weight.assign(coords.size(), 1.0);
    ref.energy.assign(1, std::vector<double>(coords.size(), 0.0));
    ref.have.assign(1, std::vector<char>(coords.size(), 0));
    for(size_t i = 0; i < coords.size(); i++)
    {
        const auto &eMap = coords[i].eMap;
        const auto  dd   = eMap.find(InteractionType::EPOT);
        if (dd != eMap.end() && dd->second.haveQM())
        {
            double eqm = dd->second.eqm();
            if (beta > 0)
            {
                double earg   = beta*(eqm - coords[i].eMin);
                ref.weight[i] = earg < DOUBLE_MAX_LOG ? std::exp(-earg) : 1e-8;
            }
            ref.energy[0][i] = eqm;
            ref.have[0][i]   = 1;
        }
    }
    return ref;
}

/*! \brief Add the weighted losses of the conformations that have both energies
 * Conformations without a reference energy, without a computed one or with
 * one that is not finite are left out. Without branches, such that the
 * loop is fast, and in the order of the conformations, such that the sums
 * are the same as when adding the losses one by one.
 * \tparam     lf       The loss function
 * \param[in]  nconf    Number of conformations
 * \param[in]  w        The weights of the conformations
 * \param[in]  eqm      The reference energies
 * \param[in]  have     Whether there is a reference energy
 * \param[in]  eact     The computed energies
 * \param[in]  haveAct  Whether there is a computed energy
 * \param[inout] weight  The sum of the weights
 * \param[inout] losses  The sum of the weighted losses
 */
template<LossFunction lf>
static void addLosses(size_t        nconf,
                      const double *w,
                      const double *eqm,
                      const char   *have,
                      const double *eact,
                      const char   *haveAct,
                      double       *weight,
                      double       *losses)
{
    double sw = *weight;
    double sl = *losses;
    for(size_t i = 0; i < nconf; i++)
    {
        bool   use = haveAct[i] && have[i] && std::isfinite(eact[i]);
        double de  = use ? eqm[i] - eact[i] : 0.0;
        double wi  = use ? w[i] : 0.0;
        sl += wi*loss(lf, 1, de);
        sw += wi;
    }
    *weight = sw;
    *losses = sl;
}

//! \copydoc addLosses
static void addLosses(LossFunction lf, size_t nconf, const double *w, const double *eqm,
                      const char *have, const double *eact, const char *haveAct,
                      double *weight, double *losses)
{
    switch (lf)
    {
    case LossFunction::MSE:
        addLosses<LossFunction::MSE>(nconf, w, eqm, have, eact, haveAct, weight, losses);
        break;
    case LossFunction::MAE:
        addLosses<LossFunction::MAE>(nconf, w, eqm, have, eact, haveAct, weight, losses);
        break;
    case LossFunction::Huber:
        addLosses<LossFunction::Huber>(nconf, w, eqm, have, eact, haveAct, weight, losses);
        break;
    case LossFunction::Asinh:
        addLosses<LossFunction::Asinh>(nconf, w, eqm, have, eact, haveAct, weight, losses);
        break;
    }
}

ForceEnergyDevComputer::EnergyTask ForceEnergyDevComputer::energyTask(const std::map<eRMS, FittingTarget> &targets) const
{
    EnergyTask task;
    const auto te = targets.find(eRMS::EPOT);
    const auto tf = targets.find(eRMS::Force2);
    bool       doEpot = (te != targets.end() && te->second.weight() > 0);
    task.forces       = (tf != targets.end() && tf->second.weight() > 0);
    bool       doInter = false;
    for (auto rms : { eRMS::Interaction, eRMS::Electrostatics, eRMS::Dispersion, eRMS::Exchange,
                      eRMS::Induction, eRMS::DeltaHF, eRMS::ExchInd, eRMS::AllElec })
    {
        const auto t = targets.find(rms);
        doInter = doInter || (t != targets.end() && t->second.weight() > 0);
    }
    task.compute                     = doEpot || task.forces || doInter;
    task.separateInductionCorrection = separateInductionCorrection_;
    return task;
}

void ForceEnergyDevComputer::calcDeviation(MsgHandler                     *msghandler,
                                           gmx_unused const ForceComputer *forceComputer,
                                           const std::vector<ACTMol *>    &actmols,
                                           std::map<eRMS, FittingTarget>  *targets,
                                           const ForceField               *forcefield,
                                           LossFunction                    lossFunction)
{
    // Check what we need to do first!
    auto ermse    = eRMS::EPOT;
    const auto te       = targets->find(ermse);
    bool doEpot   = (te != targets->end() && te->second.weight() > 0);
    const auto tf       = targets->find(eRMS::Force2);
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
                                  const auto ttt = targets->find(r.rms);
                                  return ttt == targets->end() || ttt->second.weight() <= 0;
                              }), rmsE.end());
    doInter = !rmsE.empty();
    // Whether to treat the induction correction separately is set by the
    // caller, see ACMFitnessComputer::fillDevComputers.
    bool sepInductionTerms = separateInductionCorrection_;
    if (doForce2 || doEpot || doInter)
    {
        // Compute all energies at once, and the forces if needed, unless
        // the caller did, see energyTask
        if (computeEnergies_)
        {
            forceComputer->computeRelaxedEnergies(msghandler, forcefield, actmols, true, sepInductionTerms, doForce2);
        }
        for(auto actmol : actmols)
        {
            if (doForce2)
            {
                // Forces can be from dimers or monomers
                double weight = 0, loss = 0;
                if (actmol->isDimer())
                {
                    calc_force_loss(msghandler, actmol->getMolname(), actmol->topology()->atoms(),
                                    actmol->conformationData(2), lossFunction, &weight, &loss);
                }
                else
                {
                    calc_force_loss(msghandler, actmol->getMolname(),
                                    actmol->fragmentHandler()->topologies()[0].atoms(),
                                    actmol->conformationData(0), lossFunction, &weight, &loss);
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
                // The ACT energies are taken from the conformation data,
                // the QM energies and weights from the reference data.
                auto      *cdata  = actmol->conformationData(2);
                auto      &ener   = cdata->energies;
                std::vector<InteractionType> itypes;
                for (const auto &r : rmsE)
                {
                    itypes.push_back(r.itype);
                }
                const auto &ref   = interactionReference(cdata, itypes, beta, sepInductionTerms);
                size_t      nconf = std::min(cdata->confs.size(), ener.nconf());
                for (size_t j = 0; j < rmsE.size(); j++)
                {
                    auto &rms  = rmsE[j];
                    const auto eact = ener.find(rms.itype);
                    if (eact.have)
                    {
                        addLosses(lossFunction, nconf, ref.weight.data(), ref.energy[j].data(),
                                  ref.have[j].data(), eact.act, eact.have,
                                  &rms.weight, &rms.loss);
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
                auto      *cdata  = actmol->conformationData(0);
                auto      &ener   = cdata->energies;
                const auto epot   = ener.find(InteractionType::EPOT);
                const auto &ref   = epotReference(cdata, beta);
                const auto &eqm   = ref.energy[0];
                const auto &have  = ref.have[0];
                const auto &w     = ref.weight;
                size_t      nconf = std::min(cdata->confs.size(), ener.nconf());
                for(size_t i = 0; i < nconf; i++)
                {
                    // std::isfinite checks for both Inf and NaN values.
                    if (have[i] && epot.has(i) && std::isfinite(epot.act[i]))
                    {
                        double myloss = w[i] * loss(lossFunction, 1, eqm[i]-epot.act[i]);
                        te->second.increase(w[i], myloss);
                    }
                }
            }
        }
    }
}

//! The energies of one interaction type of all conformations, as ConformationEnergies::Term
struct CompactTerm
{
    //! The energies
    std::vector<double> act;
    //! The backups
    std::vector<double> before;
    //! Whether there is an energy
    std::vector<char>   have;
};

/*! \brief The interaction energies of a dimer as postProcessRelaxedEnergies makes them
 * From the energies of the topologies as the GPU gives them: the monomers
 * are subtracted from the dimer, then the induction and the sums are
 * computed and the induction correction of a potential is moved to the
 * induction, in the same order of operations as subtractMonomers,
 * calcInduction, potentialCorrection and addCorrection in forcecomputer.cpp.
 * The induction correction is not separate.
 * \param[in]  dimer    The energies of the dimer, in increasing order of type
 * \param[in]  monomers The energies of the monomers
 * \return the energies of each type
 */
static std::map<InteractionType, CompactTerm> compactInteraction(const std::vector<GpuTypeEnergies> &dimer,
                                                                 const std::vector<GpuTypeEnergies>  monomers[2])
{
    std::map<InteractionType, CompactTerm> inter;
    size_t nconf = dimer.empty() ? 0 : dimer[0].act.size();
    // subtractMonomers
    for (const auto &d : dimer)
    {
        auto &term = inter[d.type];
        term.act    = d.act;
        term.before = d.before;
        term.have.assign(nconf, 1);
        const GpuTypeEnergies *em[2] = { nullptr, nullptr };
        for (int j = 0; j < 2; j++)
        {
            for (const auto &m : monomers[j])
            {
                if (m.type == d.type)
                {
                    em[j] = &m;
                }
            }
        }
        for (size_t c = 0; c < nconf; c++)
        {
            double edim = term.act[c];
            double edb4 = term.before[c];
            for (int j = 0; j < 2; j++)
            {
                if (em[j] && c < em[j]->act.size())
                {
                    edim -= em[j]->act[c];
                }
            }
            for (int j = 0; j < 2; j++)
            {
                if (em[j] && c < em[j]->act.size())
                {
                    edb4 -= em[j]->before[c];
                }
            }
            term.act[c]    = edim;
            term.before[c] = edb4;
        }
    }
    // calcInduction, the terms in increasing order of type
    const auto isElec = [](InteractionType t)
    {
        return (InteractionType::ELECTROSTATICS == t || InteractionType::POLARIZATION == t ||
                InteractionType::CHARGETRANSFER == t);
    };
    std::vector<std::pair<InteractionType, CompactTerm *>> terms;
    for (auto &t : inter)
    {
        if (isElec(t.first) || InteractionType::EXCHANGE == t.first)
        {
            terms.emplace_back(t.first, &t.second);
        }
    }
    bool         haveInductionTerm = inter.count(InteractionType::INDUCTION) > 0;
    CompactTerm  induction;
    if (haveInductionTerm)
    {
        induction = inter[InteractionType::INDUCTION];
    }
    CompactTerm allElec, exchInd;
    for (auto *sum : { &allElec, &exchInd })
    {
        sum->act.assign(nconf, 0.0);
        sum->before.assign(nconf, 0.0);
        sum->have.assign(nconf, 1);
    }
    for (size_t c = 0; c < nconf; c++)
    {
        bool   haveInduction = false;
        double eInduction    = 0;
        double eAllElec      = 0;
        double eExchInd      = 0;
        for (auto &tt : terms)
        {
            auto &et = *tt.second;
            if (isElec(tt.first))
            {
                if (et.before[c] != et.act[c])
                {
                    eInduction   += et.act[c] - et.before[c];
                    et.act[c]     = et.before[c];
                    et.before[c]  = 0;
                    haveInduction = true;
                }
                eAllElec += et.act[c];
            }
            else
            {
                eExchInd += et.act[c];
            }
        }
        if (haveInduction)
        {
            if (!haveInductionTerm)
            {
                induction.act.assign(nconf, 0.0);
                induction.before.assign(nconf, 0.0);
                induction.have.assign(nconf, 0);
                haveInductionTerm = true;
            }
            if (!induction.have[c])
            {
                induction.before[c] = 0;
                induction.have[c]    = 1;
            }
            induction.act[c] = eInduction;
            eAllElec        += eInduction;
            eExchInd        += eInduction;
        }
        else if (haveInductionTerm && induction.have[c])
        {
            eExchInd += induction.act[c];
        }
        allElec.act[c] = eAllElec;
        exchInd.act[c] = eExchInd;
    }
    // potentialCorrection and addCorrection, to the induction
    const auto icorr = inter.find(InteractionType::INDUCTIONCORRECTION);
    if (icorr != inter.end())
    {
        std::vector<double> correction(nconf, 0.0);
        for (size_t c = 0; c < nconf; c++)
        {
            if (icorr->second.have[c])
            {
                correction[c] = icorr->second.act[c];
            }
        }
        for (size_t c = 0; c < nconf; c++)
        {
            allElec.act[c] += correction[c];
            exchInd.act[c] += correction[c];
        }
        if (!haveInductionTerm)
        {
            induction.act.assign(nconf, 0.0);
            induction.before.assign(nconf, 0.0);
            induction.have.assign(nconf, 0);
            haveInductionTerm = true;
        }
        for (size_t c = 0; c < nconf; c++)
        {
            if (icorr->second.have[c])
            {
                double e = (induction.have[c] ? induction.act[c] : 0.0) + correction[c];
                if (!induction.have[c])
                {
                    induction.before[c] = 0;
                    induction.have[c]    = 1;
                }
                induction.act[c]        = e;
                icorr->second.have[c]   = 0;
            }
        }
    }
    if (haveInductionTerm)
    {
        inter[InteractionType::INDUCTION] = std::move(induction);
    }
    inter[InteractionType::ALLELEC] = std::move(allElec);
    inter[InteractionType::EXCHIND] = std::move(exchInd);
    return inter;
}

//! The energy targets of the deviations, as ForceEnergyDevComputer::calcDeviation selects them
struct EnergyTargets
{
    //! Whether there are the potential energy and the force targets
    bool                         epot   = false;
    bool                         force2 = false;
    //! The interaction energy targets and their interaction types
    std::vector<eRMS>            rms;
    std::vector<InteractionType> itypes;
};

/*! \brief Select the energy targets with a weight
 * \param[in] targets The fitting targets
 * \return the targets
 */
static EnergyTargets energyTargets(const std::map<eRMS, FittingTarget> &targets)
{
    EnergyTargets et;
    const auto te   = targets.find(eRMS::EPOT);
    et.epot   = (te != targets.end() && te->second.weight() > 0);
    const auto tf   = targets.find(eRMS::Force2);
    et.force2 = (tf != targets.end() && tf->second.weight() > 0);
    const std::vector<std::pair<eRMS, InteractionType>> all = {
        { eRMS::Interaction,    InteractionType::EPOT                },
        { eRMS::Electrostatics, InteractionType::ELECTROSTATICS      },
        { eRMS::Dispersion,     InteractionType::DISPERSION          },
        { eRMS::Exchange,       InteractionType::EXCHANGE            },
        { eRMS::Induction,      InteractionType::INDUCTION           },
        { eRMS::DeltaHF,        InteractionType::INDUCTIONCORRECTION },
        { eRMS::ExchInd,        InteractionType::EXCHIND             },
        { eRMS::AllElec,        InteractionType::ALLELEC             }
    };
    for (const auto &r : all)
    {
        const auto t = targets.find(r.first);
        if (t != targets.end() && t->second.weight() > 0)
        {
            et.rms.push_back(r.first);
            et.itypes.push_back(r.second);
        }
    }
    return et;
}

bool ForceEnergyDevComputer::setupGpuDeviations(const ForceComputer                 *forceComputer,
                                                const std::vector<ACTMol *>         &actmols,
                                                const std::map<eRMS, FittingTarget> &targets,
                                                LossFunction                         lossFunction)
{
    const auto et = energyTargets(targets);
    if (et.force2 || separateInductionCorrection_)
    {
        return false;
    }
    GpuDeviationSetup setup;
    setup.targets      = et.itypes;
    setup.epot         = et.epot;
    setup.lossFunction = static_cast<int>(lossFunction);
    double betaInter   = computeBeta(eRMS::Interaction);
    double betaEpot    = computeBeta(eRMS::EPOT);
    for (auto actmol : actmols)
    {
        setup.mols.push_back(actmol);
        const DeviationReference *inter  = nullptr;
        size_t                    nInter = 0;
        if (!et.rms.empty() && actmol->isDimer())
        {
            auto *cdata = actmol->conformationData(2);
            inter       = &interactionReference(cdata, et.itypes, betaInter, false);
            nInter      = std::min(cdata->confs.size(), cdata->nconf());
        }
        setup.interaction.push_back(inter);
        setup.interactionConfs.push_back(nInter);
        const DeviationReference *epot  = nullptr;
        size_t                    nEpot = 0;
        if (et.epot)
        {
            auto *cdata = actmol->conformationData(0);
            epot        = &epotReference(cdata, betaEpot);
            nEpot       = std::min(cdata->confs.size(), cdata->nconf());
        }
        setup.epotReference.push_back(epot);
        setup.epotConfs.push_back(nEpot);
    }
    return forceComputer->setGpuDeviations(setup);
}

bool ForceEnergyDevComputer::calcDeviationGpu(const ForceComputer           *forceComputer,
                                              std::map<eRMS, FittingTarget> *targets) const
{
    std::vector<std::pair<double, double>> sums;
    auto et = energyTargets(*targets);
    if (!forceComputer->gpuDeviationSums(&sums) || sums.size() != et.rms.size() + 1)
    {
        return false;
    }
    for (size_t j = 0; j < et.rms.size(); j++)
    {
        targets->find(et.rms[j])->second.increase(sums[j].first, sums[j].second);
    }
    if (et.epot)
    {
        targets->find(eRMS::EPOT)->second.increase(sums.back().first, sums.back().second);
    }
    return true;
}

void ForceEnergyDevComputer::gpuTargets(const std::map<eRMS, FittingTarget> &targets,
                                        std::vector<eRMS>                   *rms,
                                        bool                                *epot)
{
    const auto et = energyTargets(targets);
    *rms    = et.rms;
    *epot   = et.epot;
}

void ForceEnergyDevComputer::addGpuDeviations(const std::vector<const GpuMoleculeDeviations *> &molecules,
                                              std::map<eRMS, FittingTarget>                    *targets) const
{
    // In the order of GpuTrainingSystem::deviationSums
    auto et = energyTargets(*targets);
    for (size_t j = 0; j < et.rms.size(); j++)
    {
        double sw = 0, sl = 0;
        for (const auto *dev : molecules)
        {
            if (dev->dimer)
            {
                sw += dev->targets[2*j];
                sl += dev->targets[2*j + 1];
            }
        }
        targets->find(et.rms[j])->second.increase(sw, sl);
    }
    if (et.epot)
    {
        double sw = 0, sl = 0;
        for (const auto *dev : molecules)
        {
            for (size_t c = 0; dev->epot && c < dev->epotTerms.size(); c += 2)
            {
                sw += dev->epotTerms[c];
                sl += dev->epotTerms[c + 1];
            }
        }
        targets->find(eRMS::EPOT)->second.increase(sw, sl);
    }
}

bool ForceEnergyDevComputer::calcDeviationCompact(const ForceComputer           *forceComputer,
                                                  const std::vector<ACTMol *>   &actmols,
                                                  std::map<eRMS, FittingTarget> *targets,
                                                  LossFunction                   lossFunction)
{
    // The same targets as calcDeviation
    const auto te       = targets->find(eRMS::EPOT);
    bool doEpot   = (te != targets->end() && te->second.weight() > 0);
    const auto tf       = targets->find(eRMS::Force2);
    bool doForce2 = (tf != targets->end() && tf->second.weight() > 0);
    struct RmsInter
    {
        eRMS            rms;
        InteractionType itype;
    };
    std::vector<RmsInter> rmsE = {
        { eRMS::Interaction,    InteractionType::EPOT                },
        { eRMS::Electrostatics, InteractionType::ELECTROSTATICS      },
        { eRMS::Dispersion,     InteractionType::DISPERSION          },
        { eRMS::Exchange,       InteractionType::EXCHANGE            },
        { eRMS::Induction,      InteractionType::INDUCTION           },
        { eRMS::DeltaHF,        InteractionType::INDUCTIONCORRECTION },
        { eRMS::ExchInd,        InteractionType::EXCHIND             },
        { eRMS::AllElec,        InteractionType::ALLELEC             }
    };
    rmsE.erase(std::remove_if(rmsE.begin(), rmsE.end(),
                              [targets](const RmsInter &r)
                              {
                                  const auto ttt = targets->find(r.rms);
                                  return ttt == targets->end() || ttt->second.weight() <= 0;
                              }), rmsE.end());
    bool doInter = !rmsE.empty();
    // The separate induction correction needs more calculations, see
    // ForceComputer::computeRelaxedEnergies, and the forces are not kept
    if (doForce2 || separateInductionCorrection_)
    {
        return false;
    }
    std::vector<InteractionType> itypes;
    for (const auto &r : rmsE)
    {
        itypes.push_back(r.itype);
    }
    for (auto actmol : actmols)
    {
        if (doInter)
        {
            std::vector<double> weight(rmsE.size(), 0.0), losses(rmsE.size(), 0.0);
            double              beta  = computeBeta(eRMS::Interaction);
            auto               *cdata = actmol->conformationData(2);
            if (actmol->isDimer())
            {
                std::vector<GpuTypeEnergies> dimer, monomers[2];
                if (!forceComputer->gpuTypeEnergies(actmol, 2, &dimer))
                {
                    return false;
                }
                for (int j = 0; j < 2; j++)
                {
                    (void) forceComputer->gpuTypeEnergies(actmol, j, &monomers[j]);
                }
                auto        inter = compactInteraction(dimer, monomers);
                const auto &ref   = interactionReference(cdata, itypes, beta, false);
                size_t      nconf = std::min(cdata->confs.size(), dimer.empty() ? 0 : dimer[0].act.size());
                for (size_t j = 0; j < rmsE.size(); j++)
                {
                    const auto eact = inter.find(rmsE[j].itype);
                    if (eact != inter.end())
                    {
                        addLosses(lossFunction, nconf, ref.weight.data(), ref.energy[j].data(),
                                  ref.have[j].data(), eact->second.act.data(), eact->second.have.data(),
                                  &weight[j], &losses[j]);
                    }
                }
            }
            for (size_t j = 0; j < rmsE.size(); j++)
            {
                targets->find(rmsE[j].rms)->second.increase(weight[j], losses[j]);
            }
        }
        if (doEpot)
        {
            std::vector<GpuTypeEnergies> monomer;
            if (!forceComputer->gpuTypeEnergies(actmol, 0, &monomer))
            {
                return false;
            }
            const GpuTypeEnergies *epot = nullptr;
            for (const auto &m : monomer)
            {
                if (m.type == InteractionType::EPOT)
                {
                    epot = &m;
                }
            }
            double      beta  = computeBeta(eRMS::EPOT);
            auto       *cdata = actmol->conformationData(0);
            const auto &ref   = epotReference(cdata, beta);
            const auto &eqm   = ref.energy[0];
            const auto &have  = ref.have[0];
            const auto &w     = ref.weight;
            size_t      nconf = std::min(cdata->confs.size(), epot ? epot->act.size() : 0);
            for (size_t i = 0; i < nconf; i++)
            {
                if (have[i] && std::isfinite(epot->act[i]))
                {
                    double myloss = w[i] * loss(lossFunction, 1, eqm[i]-epot->act[i]);
                    te->second.increase(w[i], myloss);
                }
            }
        }
    }
    return true;
}

/* * * * * * * * * * * * * * * * * * * * * *
* END: ForceEnergyDevComputer              *
* * * * * * * * * * * * * * * * * * * * * */

} // namespace alexandria

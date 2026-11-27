/*
 * This source file is part of the Alexandria Chemistry Toolkit.
 *
 * Copyright (C) 2026
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
 * \author David van der Spoel <david.vanderspoel@icm.uu.se>
 */

#include <algorithm>
#include <cmath>
#include <cstddef>
#include <map>
#include <random>
#include <utility>
#include <vector>

#include "act/basics/act_particle.h"
#include "act/basics/coordinate.h"
#include "act/basics/interactiontype.h"
#include "act/topology/actatom.h"
#include "gromacs/math/units.h"
#include "gromacs/math/vectypes.h"
#include "gromacs/utility/exceptions.h"
#include "gromacs/utility/stringutil.h"
#include "ringpolymer.h"

namespace alexandria
{

namespace
{

//! Planck's constant divided by 2 pi in kJ/mol ps
const double hbar = PLANCK/(2*M_PI);

/*! Largest product of the time step and the highest frequency of the free
 * ring polymer. Its motion is integrated exactly, but together with the
 * forces of the molecule it is not stable when this product is close to pi.
 */
const double maxSpringPhase = 2.0;

//! An integer coordinate, to assign it to a component of CoordinateVectors
struct Units
{
    //! The coordinate in units of coordinateUnit
    Coordinate u;
    //! \return the coordinate in units of coordinateUnit
    Coordinate units() const { return u; }
};

} // namespace

void RingPolymerSettings::check() const
{
    if (nbeads < 1)
    {
        GMX_THROW(gmx::InvalidInputError(gmx::formatString("The number of beads %d should be one or more", nbeads)));
    }
    if (!std::isfinite(temperature) || temperature <= 0)
    {
        GMX_THROW(gmx::InvalidInputError(gmx::formatString("The temperature %g K of the ring polymers should be more than zero", temperature)));
    }
    if (!std::isfinite(timeStep) || timeStep <= 0)
    {
        GMX_THROW(gmx::InvalidInputError(gmx::formatString("The time step %g ps of the ring polymers should be more than zero", timeStep)));
    }
    if (!std::isfinite(centroidTau) || centroidTau <= 0)
    {
        GMX_THROW(gmx::InvalidInputError(gmx::formatString("The time constant %g ps of the thermostat of the ring polymers should be more than zero", centroidTau)));
    }
    if (!std::isfinite(planckScale) || planckScale <= 0)
    {
        GMX_THROW(gmx::InvalidInputError(gmx::formatString("The factor %g for the constant of Planck of the ring polymers should be more than zero", planckScale)));
    }
}

RingPolymerSampler::RingPolymerSampler(const ForceField          *pd,
                                       const Topology            *top,
                                       const ForceComputer       *forceComp,
                                       const RingPolymerSettings &settings) :
    pd_(pd), top_(top), forceComp_(forceComp), settings_(settings), normal_(0.0, 1.0)
{
    settings_.check();
    nbeads_     = settings_.nbeads;
    auto &atoms = top_->atoms();
    for(size_t i = 0; i < atoms.size(); i++)
    {
        if (atoms[i].mass() > 0)
        {
            massive_.push_back(i);
            for(int m = 0; m < DIM; m++)
            {
                mass_.push_back(atoms[i].mass());
            }
        }
        else if (ActParticle::Shell == atoms[i].pType() && atoms[i].parents().size() == 1)
        {
            int parent = atoms[i].parents()[0];
            if (parent < 0 || static_cast<size_t>(parent) >= atoms.size())
            {
                GMX_THROW(gmx::InternalError(gmx::formatString("Shell %zu has parent %d, which is not a particle of its molecule of %zu particles",
                                                               i, parent, atoms.size()).c_str()));
            }
            shells_.emplace_back(i, static_cast<size_t>(parent));
        }
    }
    if (massive_.empty())
    {
        GMX_THROW(gmx::InvalidInputError("There are no particles with a mass for the ring polymer"));
    }
    kTP_    = BOLTZ*settings_.temperature*nbeads_;
    omegaP_ = kTP_/(settings_.planckScale*hbar);
    // The normal modes of the free ring polymer: the centroid, then
    // cosines up to half the number of beads and sines beyond
    size_t P = nbeads_;
    omega_.resize(P);
    cmat_.resize(P*P);
    for(size_t k = 0; k < P; k++)
    {
        omega_[k] = (k == 0) ? 0.0 : 2*omegaP_*std::sin(M_PI*k/P);
        for(size_t j = 0; j < P; j++)
        {
            double c;
            if (0 == k)
            {
                c = std::sqrt(1.0/P);
            }
            else if (2*k == P)
            {
                c = ((j % 2 == 0) ? 1.0 : -1.0)*std::sqrt(1.0/P);
            }
            else if (2*k < P)
            {
                c = std::sqrt(2.0/P)*std::cos(2*M_PI*j*k/P);
            }
            else
            {
                c = std::sqrt(2.0/P)*std::sin(2*M_PI*j*k/P);
            }
            cmat_[P*j + k] = c;
        }
    }
    double omegaMax = *std::max_element(omega_.begin(), omega_.end());
    if (omegaMax*settings_.timeStep > maxSpringPhase)
    {
        GMX_THROW(gmx::InvalidInputError(gmx::formatString("The time step %g ps is too large for ring polymers of %zu beads at %g K, it should be at most %g ps",
                                                           settings_.timeStep, P, settings_.temperature,
                                                           maxSpringPhase/omegaMax)));
    }
    size_t n = mass_.size();
    q_.assign(n*P, 0.0);
    p_.assign(n*P, 0.0);
    f_.assign(n*P, 0.0);
    qn_.assign(n*P, 0.0);
    pn_.assign(n*P, 0.0);
    epot_.assign(P, 0.0);
    if (settings_.seed == 0)
    {
        std::random_device rd;
        gen_.seed(rd());
    }
    else
    {
        gen_.seed(settings_.seed);
    }
}

void RingPolymerSampler::toNormalModes(const std::vector<double> &beads,
                                       std::vector<double>       *modes) const
{
    size_t P = nbeads_;
    size_t n = mass_.size();
    std::fill(modes->begin(), modes->end(), 0.0);
    for(size_t k = 0; k < P; k++)
    {
        double *mk = modes->data() + n*k;
        for(size_t j = 0; j < P; j++)
        {
            double        c  = cmat_[P*j + k];
            const double *bj = beads.data() + n*j;
            for(size_t d = 0; d < n; d++)
            {
                mk[d] += c*bj[d];
            }
        }
    }
}

void RingPolymerSampler::fromNormalModes(const std::vector<double> &modes,
                                         std::vector<double>       *beads) const
{
    size_t P = nbeads_;
    size_t n = mass_.size();
    std::fill(beads->begin(), beads->end(), 0.0);
    for(size_t j = 0; j < P; j++)
    {
        double *bj = beads->data() + n*j;
        for(size_t k = 0; k < P; k++)
        {
            double        c  = cmat_[P*j + k];
            const double *mk = modes.data() + n*k;
            for(size_t d = 0; d < n; d++)
            {
                bj[d] += c*mk[d];
            }
        }
    }
}

void RingPolymerSampler::freeMotion(double dt)
{
    size_t n = mass_.size();
    for(size_t k = 0; k < nbeads_; k++)
    {
        double *q = qn_.data() + n*k;
        double *p = pn_.data() + n*k;
        if (0 == omega_[k])
        {
            for(size_t d = 0; d < n; d++)
            {
                q[d] += dt*p[d]/mass_[d];
            }
        }
        else
        {
            double w = omega_[k];
            double c = std::cos(w*dt);
            double s = std::sin(w*dt);
            for(size_t d = 0; d < n; d++)
            {
                double mw = mass_[d]*w;
                double pd = c*p[d] - mw*s*q[d];
                q[d]      = (s/mw)*p[d] + c*q[d];
                p[d]      = pd;
            }
        }
    }
}

void RingPolymerSampler::thermostat(double dt)
{
    size_t n = mass_.size();
    for(size_t k = 0; k < nbeads_; k++)
    {
        // Critical damping of the internal modes
        double gamma = (0 == omega_[k]) ? 1.0/settings_.centroidTau : 2*omega_[k];
        double c1    = std::exp(-gamma*dt);
        double c2    = std::sqrt((1 - c1*c1)*kTP_);
        double *p    = pn_.data() + n*k;
        for(size_t d = 0; d < n; d++)
        {
            p[d] = c1*p[d] + c2*std::sqrt(mass_[d])*normal_(gen_);
        }
    }
}

void RingPolymerSampler::stopCenterOfMass()
{
    // The centroid is the first normal mode
    size_t n = mass_.size();
    for(int m = 0; m < DIM; m++)
    {
        double ptot = 0, mtot = 0;
        for(size_t d = m; d < n; d += DIM)
        {
            ptot += pn_[d];
            mtot += mass_[d];
        }
        double v = ptot/mtot;
        for(size_t d = m; d < n; d += DIM)
        {
            pn_[d] -= mass_[d]*v;
        }
    }
}

void RingPolymerSampler::setCentroidPotential(CentroidPotential potential)
{
    centroidPotential_ = std::move(potential);
}

void RingPolymerSampler::computeForces(MsgHandler *msghandler,
                                       bool        first,
                                       bool        energyTerms)
{
    size_t n     = mass_.size();
    size_t natom = top_->atoms().size();
    if (!first)
    {
        if (scratch_.size() != natom)
        {
            scratch_ = Coordinates(natom);
        }
        for(size_t j = 0; j < nbeads_; j++)
        {
            const auto xb = set_.data.x(j);
            for(size_t a = 0; a < natom; a++)
            {
                for(int m = 0; m < DIM; m++)
                {
                    scratch_[a][m] = Units{ xb[a].units(m) };
                }
            }
            const double *q = q_.data() + n*j;
            for(size_t i = 0; i < massive_.size(); i++)
            {
                for(int m = 0; m < DIM; m++)
                {
                    scratch_[massive_[i]][m] = q[DIM*i + m];
                }
            }
            // The shells move with their parents, which may be virtual
            // sites, and are relaxed from there
            forceComp_->constructVsiteCoordinates(top_, &scratch_);
            for(const auto &sp : shells_)
            {
                for(int m = 0; m < DIM; m++)
                {
                    Coordinate offset = coordinateDifference(xb[sp.first].units(m), xb[sp.second].units(m));
                    scratch_[sp.first][m] = Units{ coordinateSum(scratch_[sp.second].units(m), offset) };
                }
            }
            set_.data.setX(j, std::as_const(scratch_).view());
        }
    }
    set_.invalidate();
    // The induction energy is the change of the energy from the shells on
    // their parents. For the energy terms the shells are put there first,
    // and are relaxed from where they were.
    ShellOptions options;
    options.resetShells     = first || energyTerms;
    options.shellsFromInput = energyTerms && !first;
    options.energyTerms     = energyTerms;
    std::vector<std::map<InteractionType, double>> beadTerms;
    forceComp_->computeTotal(msghandler, pd_, &set_, options, energyTerms ? &beadTerms : nullptr);
    terms_.clear();
    for(const auto &bt : beadTerms)
    {
        for(const auto &t : bt)
        {
            terms_[t.first] += t.second/nbeads_;
        }
    }
    const auto epot = set_.data.energies.find(InteractionType::EPOT);
    for(size_t j = 0; j < nbeads_; j++)
    {
        if (!epot.has(j))
        {
            GMX_THROW(gmx::InternalError("No potential energy for a bead of the ring polymer"));
        }
        epot_[j] = epot.act[j];
        const auto fb = set_.data.f(j);
        double    *f  = f_.data() + n*j;
        for(size_t i = 0; i < massive_.size(); i++)
        {
            for(int m = 0; m < DIM; m++)
            {
                f[DIM*i + m] = fb[massive_[i]][m];
            }
        }
    }
    ecentroid_ = 0;
    if (centroidPotential_)
    {
        // Each bead of an atom gets the force on its centroid
        const gmx::DVec zero = { 0, 0, 0 };
        xcentroid_.assign(natom, zero);
        fcentroid_.assign(natom, zero);
        for(size_t i = 0; i < massive_.size(); i++)
        {
            for(int m = 0; m < DIM; m++)
            {
                double qc = 0;
                for(size_t j = 0; j < nbeads_; j++)
                {
                    qc += q_[n*j + DIM*i + m];
                }
                xcentroid_[massive_[i]][m] = qc/nbeads_;
            }
        }
        ecentroid_ = centroidPotential_(xcentroid_, &fcentroid_);
        for(size_t j = 0; j < nbeads_; j++)
        {
            double *f = f_.data() + n*j;
            for(size_t i = 0; i < massive_.size(); i++)
            {
                for(int m = 0; m < DIM; m++)
                {
                    f[DIM*i + m] += fcentroid_[massive_[i]][m];
                }
            }
        }
    }
}

void RingPolymerSampler::centroids(std::vector<gmx::DVec> *x) const
{
    size_t n     = mass_.size();
    size_t natom = top_->atoms().size();
    const gmx::DVec zero = { 0, 0, 0 };
    x->assign(natom, zero);
    // The shells and virtual sites
    for(size_t j = 0; j < nbeads_; j++)
    {
        const auto xb = set_.data.x(j);
        for(size_t a = 0; a < natom; a++)
        {
            for(int m = 0; m < DIM; m++)
            {
                (*x)[a][m] += static_cast<double>(xb[a][m])/nbeads_;
            }
        }
    }
    // The particles with a mass, as they are integrated
    for(size_t i = 0; i < massive_.size(); i++)
    {
        for(int m = 0; m < DIM; m++)
        {
            double qc = 0;
            for(size_t j = 0; j < nbeads_; j++)
            {
                qc += q_[n*j + DIM*i + m];
            }
            (*x)[massive_[i]][m] = qc/nbeads_;
        }
    }
}

void RingPolymerSampler::start(MsgHandler                   *msghandler,
                               const std::vector<gmx::RVec> &x,
                               bool                          energyTerms)
{
    if (x.size() != top_->atoms().size())
    {
        GMX_THROW(gmx::InvalidInputError(gmx::formatString("There are %zu coordinates for a molecule of %zu particles",
                                                           x.size(), top_->atoms().size())));
    }
    size_t n = mass_.size();
    set_.top = top_;
    set_.setCoordinates(std::vector<std::vector<gmx::RVec>>(nbeads_, x));
    for(size_t j = 0; j < nbeads_; j++)
    {
        for(size_t i = 0; i < massive_.size(); i++)
        {
            for(int m = 0; m < DIM; m++)
            {
                size_t d = n*j + DIM*i + m;
                q_[d] = x[massive_[i]][m];
                p_[d] = std::sqrt(mass_[DIM*i + m]*kTP_)*normal_(gen_);
            }
        }
    }
    toNormalModes(p_, &pn_);
    stopCenterOfMass();
    fromNormalModes(pn_, &p_);
    computeForces(msghandler, true, energyTerms);
    started_ = true;
}

void RingPolymerSampler::step(MsgHandler *msghandler,
                              bool        energyTerms)
{
    double dt = settings_.timeStep;
    for(size_t d = 0; d < p_.size(); d++)
    {
        p_[d] += 0.5*dt*f_[d];
    }
    toNormalModes(q_, &qn_);
    toNormalModes(p_, &pn_);
    freeMotion(0.5*dt);
    if (settings_.thermostat)
    {
        thermostat(dt);
        stopCenterOfMass();
    }
    freeMotion(0.5*dt);
    fromNormalModes(qn_, &q_);
    fromNormalModes(pn_, &p_);
    computeForces(msghandler, false, energyTerms);
    for(size_t d = 0; d < p_.size(); d++)
    {
        p_[d] += 0.5*dt*f_[d];
    }
}

void RingPolymerSampler::run(MsgHandler *msghandler,
                             int         nsteps,
                             bool        energyTerms)
{
    if (!started_)
    {
        GMX_THROW(gmx::InternalError("The ring polymer has no coordinates, call start first"));
    }
    for(int i = 0; i < nsteps; i++)
    {
        step(msghandler, energyTerms && i + 1 == nsteps);
    }
}

double RingPolymerSampler::potentialEnergy() const
{
    double e = 0;
    for(auto ej : epot_)
    {
        e += ej;
    }
    return e/nbeads_;
}

double RingPolymerSampler::kineticEnergy() const
{
    size_t n      = mass_.size();
    double virial = 0;
    for(size_t d = 0; d < n; d++)
    {
        double qc = 0;
        for(size_t j = 0; j < nbeads_; j++)
        {
            qc += q_[n*j + d];
        }
        qc /= nbeads_;
        for(size_t j = 0; j < nbeads_; j++)
        {
            virial += (q_[n*j + d] - qc)*f_[n*j + d];
        }
    }
    return 0.5*n*kTP_/nbeads_ - 0.5*virial/nbeads_;
}

double RingPolymerSampler::primitiveKineticEnergy() const
{
    size_t n      = mass_.size();
    double spring = 0;
    for(size_t j = 0; j < nbeads_; j++)
    {
        size_t jn = (j + 1) % nbeads_;
        for(size_t d = 0; d < n; d++)
        {
            double dq = q_[n*j + d] - q_[n*jn + d];
            spring   += 0.5*mass_[d]*omegaP_*omegaP_*dq*dq;
        }
    }
    // The springs are those of the ring polymer at P times the temperature
    return 0.5*n*kTP_ - spring/nbeads_;
}

double RingPolymerSampler::gyration2() const
{
    size_t n  = mass_.size();
    double r2 = 0;
    for(size_t d = 0; d < n; d++)
    {
        double qc = 0;
        for(size_t j = 0; j < nbeads_; j++)
        {
            qc += q_[n*j + d];
        }
        qc /= nbeads_;
        for(size_t j = 0; j < nbeads_; j++)
        {
            double dq = q_[n*j + d] - qc;
            r2 += dq*dq;
        }
    }
    return r2/(nbeads_*massive_.size());
}

double RingPolymerSampler::hamiltonian() const
{
    size_t n = mass_.size();
    double h = 0;
    for(size_t j = 0; j < nbeads_; j++)
    {
        size_t jn = (j + 1) % nbeads_;
        for(size_t d = 0; d < n; d++)
        {
            double p  = p_[n*j + d];
            double dq = q_[n*j + d] - q_[n*jn + d];
            h += 0.5*p*p/mass_[d] + 0.5*mass_[d]*omegaP_*omegaP_*dq*dq;
        }
        h += epot_[j];
    }
    return h + nbeads_*ecentroid_;
}

double RingPolymerSampler::beadTemperature() const
{
    size_t n    = mass_.size();
    double ekin = 0;
    for(size_t j = 0; j < nbeads_; j++)
    {
        for(size_t d = 0; d < n; d++)
        {
            double p = p_[n*j + d];
            ekin += 0.5*p*p/mass_[d];
        }
    }
    // The center of mass of the centroid does not move
    double ndof = static_cast<double>(n*nbeads_) - (settings_.thermostat ? DIM : 0);
    return ndof > 0 ? 2*ekin/(ndof*BOLTZ) : 0.0;
}

RingPolymer RingPolymerSampler::polymer() const
{
    size_t n = mass_.size();
    // The center of mass of the centroid
    double com[DIM] = { 0, 0, 0 };
    double mtot     = 0;
    for(size_t i = 0; i < massive_.size(); i++)
    {
        mtot += mass_[DIM*i];
        for(size_t j = 0; j < nbeads_; j++)
        {
            for(int m = 0; m < DIM; m++)
            {
                com[m] += mass_[DIM*i]*q_[n*j + DIM*i + m];
            }
        }
    }
    for(int m = 0; m < DIM; m++)
    {
        com[m] /= mtot*nbeads_;
    }
    RingPolymer rp;
    rp.x.resize(nbeads_);
    rp.epot = epot_;
    for(size_t j = 0; j < nbeads_; j++)
    {
        const auto xb = set_.data.x(j);
        rp.x[j].resize(xb.size());
        for(size_t a = 0; a < xb.size(); a++)
        {
            for(int m = 0; m < DIM; m++)
            {
                rp.x[j][a][m] = static_cast<double>(xb[a][m]) - com[m];
            }
        }
    }
    return rp;
}

RingPolymerAverages RingPolymerSampler::sample(MsgHandler               *msghandler,
                                               int                       nequil,
                                               int                       stride,
                                               size_t                    npolymers,
                                               std::vector<RingPolymer> *pool)
{
    RingPolymerAverages aver;
    run(msghandler, nequil);
    pool->reserve(pool->size() + npolymers);
    for(size_t i = 0; i < npolymers; i++)
    {
        if (i > 0)
        {
            run(msghandler, stride);
        }
        pool->push_back(polymer());
        double epot = potentialEnergy();
        double ekin = kineticEnergy();
        aver.n         += 1;
        aver.epot      += epot;
        aver.epot2     += epot*epot;
        aver.ekin      += ekin;
        aver.ekin2     += ekin*ekin;
        aver.gyration2 += gyration2();
    }
    return aver;
}

} // namespace alexandria

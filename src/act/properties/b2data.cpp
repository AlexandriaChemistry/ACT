/*
 * This source file is part of the Alexandria Chemistry Toolkit.
 *
 * Copyright (C) 2024-2026
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
#include "b2data.h"

#include <algorithm>
#include <cmath>
#include <cstdio>
#include <sstream>

#include "act/forces/forces_explicit_simd.h"
#include "act/forces/simd_information.h"
#include "act/utility/simpson.h"
#include "act/utility/timings.h"
#include "gromacs/fileio/oenv.h"
#include "gromacs/fileio/xvgr.h"
#include "gromacs/math/units.h"
#include "gromacs/math/vec.h"
#include "gromacs/utility/exceptions.h"
#include "gromacs/utility/stringutil.h"

namespace alexandria
{

B2Data::B2Data(int                        nbins,
               const std::vector<double> &temperatures,
               double                     eMin) : eMin_(eMin)
{
    temperatures_ = temperatures;
    nbins_        = nbins;
    size_t width  = explicitSimdWidth<double>();
    tstride_      = ((temperatures.size() + width - 1)/width)*width;
    sums_.resize(nbins_*tstride_*nSum, 0.0);
    mayer_.resize(temperatures.size());
    for(auto T : temperatures)
    {
        beta_.push_back(T == 0 ? 0.0 : 1.0/(BOLTZ*T));
        allFinite_ = allFinite_ && T != 0;
    }
    g_.resize(temperatures.size());
    betaPadded_ = beta_;
    betaPadded_.resize(tstride_, 0.0);
    dist_.resize(nbins, -1);
    for(auto T : temperatures)
    {
        if (T == 0)
        {
            fprintf(stderr, "Please provide a finite temperature to compute second virial.\n");
        }
    }
}

void B2Data::dump(MsgHandler *msghandler)
{
    flush();
    if (msghandler)
    {
        std::ostringstream oss;
        for(size_t i = 0; i < nbins_; i++)
        {
            oss << " " << sum(i, 0, sumN);
        }
        msghandler->writeDebug(oss.str());
    }
}

template<typename Func>
void B2Data::forEachSum(Func func)
{
    flush();
    for(auto &x : sums_)
    {
        func(x);
    }
    func(nDimer_);
    func(nOverlap_);
}

void B2Data::aggregate(CommunicationRecord *cr)
{
    // All sums of all processes are added with a reduction on the master,
    // in parallel on the processes rather than by the master alone, which
    // took most of the time with many processes.
    std::vector<double> buf;
    forEachSum([&buf](double &x) { buf.push_back(x); });
    cr->sumd_master(static_cast<int>(buf.size()), buf.data());
    if (cr->isMaster())
    {
        size_t k = 0;
        forEachSum([&buf, &k](double &x) { x = buf[k++]; });
    }
}

void B2Data::addData(size_t iTemp, size_t index, double dist,
                     double exp_U12, double exp_F0, double exp_F1,
                     const gmx::RVec exp_tau0, const gmx::RVec exp_tau1,
                     double weight)
{
    // After the orientations added before
    flush();
    if (iTemp >= temperatures_.size())
    {
        GMX_THROW(gmx::InternalError(gmx::formatString("iTemp = %zu, should be less than %zu", iTemp, temperatures_.size()).c_str()));
    }
    if (index >= nbins_)
    {
        GMX_THROW(gmx::InternalError(gmx::formatString("index = %zu, should be less than %zu", index, nbins_).c_str()));
    }
    dist_[index]               = dist;
    sum(index, iTemp, sumU12) += weight*exp_U12;
    sum(index, iTemp, sumF0)  += weight*exp_F0;
    sum(index, iTemp, sumF1)  += weight*exp_F1;
    for(int m = 0; m < DIM; m++)
    {
        sum(index, iTemp, sumT0+m) += weight*exp_tau0[m];
        sum(index, iTemp, sumT1+m) += weight*exp_tau1[m];
    }
    sum(index, iTemp, sumN)   += weight;
}

void B2Data::addOrientation(const std::vector<size_t>    &index,
                            const std::vector<double>    &dist,
                            const std::vector<double>    &energy,
                            const std::vector<double>     force2[2],
                            const std::vector<gmx::RVec>  torque[2],
                            const std::vector<double>    &weight)
{
    size_t ndimer = index.size();
    if (dist.size() != ndimer || energy.size() != ndimer || weight.size() != ndimer ||
        force2[0].size() != ndimer || force2[1].size() != ndimer ||
        torque[0].size() != ndimer || torque[1].size() != ndimer)
    {
        GMX_THROW(gmx::InternalError("Inconsistent number of dimers"));
    }
    // Whether the dimers are in the order of their distances
    bool inOrder = true;
    for(size_t ii = 0; ii < ndimer; ii++)
    {
        if (index[ii] >= dist_.size())
        {
            GMX_THROW(gmx::InternalError(gmx::formatString("index = %zu, should be less than %zu", index[ii], dist_.size()).c_str()));
        }
        dist_[index[ii]] = dist[ii];
        inOrder          = inOrder && index[ii] == ii;
    }
    // Dimers with an energy below eMin_ are overlapping atoms
    for(size_t ii = 0; ii < ndimer; ii++)
    {
        if (energy[ii] < eMin_)
        {
            nOverlap_ += 1;
        }
    }
    nDimer_ += ndimer;
    if (!inOrder)
    {
        // Added at once, after the orientations before it
        flush();
        const Orientation one = { energy, weight, { force2[0], force2[1] }, { torque[0], torque[1] } };
        addRange(one, index.data(), 0, ndimer);
        return;
    }
    // Kept until maxPending_ orientations are there, see flush
    if (pending_.size() <= npending_)
    {
        pending_.resize(npending_ + 1);
    }
    auto &o = pending_[npending_++];
    o.energy.assign(energy.begin(), energy.end());
    o.weight.assign(weight.begin(), weight.end());
    for(int kk = 0; kk < 2; kk++)
    {
        o.force2[kk].assign(force2[kk].begin(), force2[kk].end());
        o.torque[kk].assign(torque[kk].begin(), torque[kk].end());
    }
    if (npending_ == maxPending_)
    {
        flush();
    }
}

void B2Data::addTerms(const Orientation &o,
                      size_t             ii,
                      double             w,
                      const double      *g,
                      double            *s,
                      size_t             stride,
                      size_t             nT)
{
    // Terms with a Boltzmann factor of zero add zero, the sums are
    // never negative zero, so adding zero does not change them.
    // The loops over the temperatures can be vectorized.
    double f0 = o.force2[0][ii];
    double f1 = o.force2[1][ii];
    for(size_t t = 0; t < nT; t++)
    {
        s[sumU12*stride + t] += w*(g[t]-1);
    }
    // Gray and Gubbins Eqn. 3.272
    for(size_t t = 0; t < nT; t++)
    {
        s[sumF0*stride + t] += w*(g[t] > 0 ? g[t]*f0 : 0.0);
    }
    for(size_t t = 0; t < nT; t++)
    {
        s[sumF1*stride + t] += w*(g[t] > 0 ? g[t]*f1 : 0.0);
    }
    // Gray and Gubbins Eqn. 3.282
    for(int kk = 0; kk < 2; kk++)
    {
        for(int m = 0; m < DIM; m++)
        {
            double  tau = o.torque[kk][ii][m];
            double *st  = s + ((kk == 0 ? sumT0 : sumT1) + m)*stride;
            for(size_t t = 0; t < nT; t++)
            {
                st[t] += w*(g[t] > 0 ? g[t]*tau*tau : 0.0);
            }
        }
    }
    for(size_t t = 0; t < nT; t++)
    {
        s[sumN*stride + t] += w;
    }
}

void B2Data::addRange(const Orientation &o,
                      const size_t      *index,
                      size_t             i0,
                      size_t             i1)
{
    size_t ntemp = 0;
    for(auto beta : beta_)
    {
        ntemp += (beta > 0) ? 1 : 0;
    }
    // The dimers in the outer loop, such that the sums of a distance, of
    // all temperatures, are contiguous. Each sum gets one term.
    size_t  nT       = beta_.size();
    double *g        = g_.data();
    bool    simd     = useExplicitSimd();
    size_t  nsimd    = 0;
    for(size_t ii = i0; ii < i1; ii++)
    {
        size_t  j       = index ? index[ii] : ii;
        double  w       = o.weight[ii];
        double  energy  = o.energy[ii];
        bool    overlap = energy < eMin_;
        if (simd)
        {
            // The weighted squared forces and torques, as the sums
            // sumF0 to sumT1+ZZ, which follow each other
            double wx[8] = { w*o.force2[0][ii], w*o.force2[1][ii] };
            bool   finite = std::isfinite(wx[0]) && std::isfinite(wx[1]);
            for(int kk = 0; kk < 2; kk++)
            {
                for(int m = 0; m < DIM; m++)
                {
                    double tau       = o.torque[kk][ii][m];
                    wx[2 + DIM*kk + m] = w*(tau*tau);
                    finite           = finite && std::isfinite(wx[2 + DIM*kk + m]);
                }
            }
            // Otherwise the terms at temperatures with a Boltzmann factor
            // of zero would not be zero
            if (finite &&
                addB2TermsExplicitSimd(tstride_, betaPadded_.data(), energy, overlap, w, wx,
                                       &sum(j, 0, 0), tstride_))
            {
                nsimd += 1;
                continue;
            }
        }
        // The Boltzmann factors at all temperatures. Not vectorized, since
        // a vector exp is not faster here, but can differ from exp in the
        // last bit. Gray and Gubbins Eqn. 3.261
#if defined(__clang__)
#pragma clang loop vectorize(disable)
#endif
        for(size_t iTemp = 0; iTemp < nT; iTemp++)
        {
            double yb = energy*beta_[iTemp];
            g[iTemp]  = (!overlap && yb < 700) ? std::exp(-yb) : 0;
        }
        double *s = &sum(j, 0, 0);
        if (!allFinite_)
        {
            // Temperatures of zero are left out
            for(size_t iTemp = 0; iTemp < nT; iTemp++)
            {
                if (beta_[iTemp] != 0)
                {
                    addTerms(o, ii, w, g + iTemp, s + iTemp, tstride_, 1);
                }
            }
        }
        else
        {
            addTerms(o, ii, w, g, s, tstride_, nT);
        }
    }
    // Per dimer and temperature: the Boltzmann factor with the exponential
    // 31, the Mayer function 3, the forces 6, the torques 24 and the weight 1,
    // see Timings::addFlops. With explicit SIMD the forces and torques are
    // multiplied by the weight once per dimer, 14, and added with 2 per
    // temperature.
    Timings::instance().addFlops(65.0*(i1 - i0 - nsimd)*ntemp + nsimd*(14.0 + 51.0*ntemp));
}

void B2Data::flush()
{
    if (npending_ == 0)
    {
        return;
    }
    // A block of distances of all pending orientations at all temperatures,
    // such that the sums of the block stay in the cache while the terms of
    // all orientations are added. Each sum gets the terms in the order of
    // the orientations, as without the blocks, such that the sums are the
    // same. Without the blocks, the sums of all distances and temperatures
    // are read and written for each orientation, which is limited by the
    // memory bandwidth.
    constexpr size_t block = 16;
    size_t nmax = 0;
    for(size_t b = 0; b < npending_; b++)
    {
        nmax = std::max(nmax, pending_[b].energy.size());
    }
    for(size_t i0 = 0; i0 < nmax; i0 += block)
    {
        for(size_t b = 0; b < npending_; b++)
        {
            const auto &o  = pending_[b];
            size_t      i1 = std::min(i0 + block, o.energy.size());
            if (i0 < i1)
            {
                addRange(o, nullptr, i0, i1);
            }
        }
    }
    npending_ = 0;
}

void B2Data::integrate(MsgHandler *msghandler,
                       int iTemp, double beta,
                       const std::vector<double>    &mass,
                       const std::vector<gmx::RVec> &inertia,
                       double *Bclass, double *BqmForce,
                       double *BqmTorque1, double *BqmTorque2)
{
    flush();
    // Constant for quantum corrections
    double    hbarfac  = beta*gmx::square(beta*PLANCK/(2*M_PI))/24;
    // Initiate arrays
    mayer_[iTemp].resize(nbins_, -1);
    std::vector<double> myBforce(nbins_, 0);
    std::vector<double> myBtorq[2];
    for(int kk = 0; kk < 2; kk++)
    {
        myBtorq[kk].resize(nbins_, 0);
    }
    for(size_t ii = 0; ii < nbins_; ii++)
    {
        double n = sum(ii, iTemp, sumN);
        if (n > 0)
        {
            double Unew = sum(ii, iTemp, sumU12)/n;
            mayer_[iTemp][ii] = Unew;
            for(int kk = 0; kk < 2; kk++)
            {
                // Weighted square force
                // We follow Eqn. 9 in Schenter, JCP 117 (2002) 6573
                myBforce[ii] += sum(ii, iTemp, sumF0+kk)/(mass[kk]*n);
                // Contributions from torque
                for(int m = 0; m < DIM; m++)
                {
                    if (inertia[kk][m] > 0)
                    {
                        myBtorq[kk][ii] += sum(ii, iTemp, (kk == 0 ? sumT0 : sumT1)+m)/(n*inertia[kk][m]);
                    }
                }
            }
        }
        else
        {
            if (msghandler)
            {
                msghandler->msg(ACTStatus::Error,
                                gmx::formatString("No data for point %zu (dist %g) to integrate Mayer function", ii, dist_[ii]));
            }
        }
    }
    *Bclass     = -0.5*simpsonIntegrate(msghandler, true, dist_, mayer_[iTemp]);
    // Below the shortest distance the Mayer function is -1, which gives
    // the volume of a sphere with that radius times one half
    double xmin = dist_.empty() ? 0.0 : dist_[0];
    *Bclass    += (2.0*M_PI/3.0)*xmin*xmin*xmin;
    *BqmForce   =  0.5*hbarfac*simpsonIntegrate(msghandler, true, dist_, myBforce);
    *BqmTorque1 = hbarfac*simpsonIntegrate(msghandler, true, dist_, myBtorq[0]);
    *BqmTorque2 = hbarfac*simpsonIntegrate(msghandler, true, dist_, myBtorq[1]);
}

void B2Data::plotMayer(const char                *fmayer,
                       gmx_output_env_t          *oenv)
{
    if (!fmayer || strlen(fmayer) == 0)
    {
        return;
    }
    GMX_RELEASE_ASSERT(mayer_.size() == temperatures_.size(), "Mismatch between Mayer and temperatures");
    FILE *fp = xvgropen(fmayer, "Mayer function", "r (nm)",
                        "< exp[-U12/kBT]-1 >", oenv);
    std::vector<std::string> label;
    for(auto T : temperatures_)
    {
        if (T != 0)
        {
            label.push_back(gmx::formatString("T = %g K", T));
        }
    }
    xvgrLegend(fp, label, oenv);
    for(size_t i = 0; i < mayer_[0].size(); i++)
    {
        fprintf(fp, "%10g", dist_[i]);
        for(size_t j = 0; j < mayer_.size(); j++)
        {
            fprintf(fp, "  %10g", mayer_[j][i]);
        }
        fprintf(fp, "\n");
    }
    xvgrclose(fp);
}

} // namespace


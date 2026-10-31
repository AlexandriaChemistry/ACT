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

#include <cmath>
#include <cstdio>
#include <sstream>

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
    exp_U12_.resize(temperatures.size());
    for(int kk = 0; kk < 2; kk++)
    {
        exp_F2_[kk].resize(temperatures.size());
        exp_tau_[kk].resize(temperatures.size());
    }
    n_U12_.resize(temperatures.size());
    mayer_.resize(temperatures.size());
    for(size_t i = 0; i < temperatures.size(); i++)
    {
        exp_U12_[i].resize(nbins, 0);
        for(int kk = 0; kk < 2; kk++)
        {
            exp_F2_[kk][i].resize(nbins, 0);
            exp_tau_[kk][i].resize(nbins, { 0.0, 0.0, 0.0 });
        }
        n_U12_[i].resize(nbins, 0);
    }
    dist_.resize(nbins, -1);
}

void B2Data::dump(MsgHandler *msghandler) const
{
    if (msghandler)
    {
        std::ostringstream oss;
        for(size_t i = 0; i < n_U12_[0].size(); i++)
        {
            oss << " " << n_U12_[0][i];
        }
        msghandler->writeDebug(oss.str());
    }
}

template<typename Func>
void B2Data::forEachSum(Func func)
{
    for(auto &t : exp_U12_)
    {
        for(auto &x : t)
        {
            func(x);
        }
    }
    for(auto &t : n_U12_)
    {
        for(auto &x : t)
        {
            func(x);
        }
    }
    for(int kk = 0; kk < 2; kk++)
    {
        for(auto &t : exp_F2_[kk])
        {
            for(auto &x : t)
            {
                func(x);
            }
        }
        for(auto &t : exp_tau_[kk])
        {
            for(auto &x : t)
            {
                for(int m = 0; m < DIM; m++)
                {
                    func(x[m]);
                }
            }
        }
    }
    func(nDimer_);
    func(nOverlap_);
}

void B2Data::aggregate(CommunicationRecord *cr)
{
    // All sums are sent in one message per process, the master adds them
    // in the order of the processes.
    size_t nsum = 0;
    forEachSum([&nsum](double &) { nsum += 1; });
    std::vector<double> buf;
    if (cr->isMaster())
    {
        for(int src = 1; src < cr->size(); src++)
        {
            if (CommunicationStatus::RECV_DATA == cr->recv_data(src))
            {
                cr->recv(src, &buf);
                if (buf.size() != nsum)
                {
                    GMX_THROW(gmx::InternalError(gmx::formatString("Received %zu sums from process %d, expected %zu", buf.size(), src, nsum).c_str()));
                }
                size_t k = 0;
                forEachSum([&buf, &k](double &x) { x += buf[k++]; });
            }
        }
    }
    else
    {
        int dest = cr->superior();
        if (CommunicationStatus::SEND_DATA == cr->send_data(dest))
        {
            buf.reserve(nsum);
            forEachSum([&buf](double &x) { buf.push_back(x); });
            cr->send(dest, buf);
        }
    }
}

void B2Data::addData(size_t iTemp, size_t index, double dist,
                     double exp_U12, double exp_F0, double exp_F1,
                     const gmx::RVec exp_tau0, const gmx::RVec exp_tau1,
                     double weight)
{
    if (iTemp >= exp_U12_.size())
    {
        GMX_THROW(gmx::InternalError(gmx::formatString("iTemp = %zu, should be less than %zu", iTemp, exp_U12_.size()).c_str()));
    }
    if (index >= exp_U12_[iTemp].size())
    {
        GMX_THROW(gmx::InternalError(gmx::formatString("index = %zu, should be less than %zu", index, exp_U12_[iTemp].size()).c_str()));
    }
    dist_[index]               = dist;
    exp_U12_[iTemp][index]    += weight*exp_U12;
    exp_F2_[0][iTemp][index]  += weight*exp_F0;
    exp_F2_[1][iTemp][index]  += weight*exp_F1;
    for(int m = 0; m < DIM; m++)
    {
        exp_tau_[0][iTemp][index][m] += weight*exp_tau0[m];
        exp_tau_[1][iTemp][index][m] += weight*exp_tau1[m];
    }
    n_U12_[iTemp][index]      += weight;
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
    std::vector<bool> overlap(ndimer);
    for(size_t ii = 0; ii < ndimer; ii++)
    {
        overlap[ii] = energy[ii] < eMin_;
        if (overlap[ii])
        {
            nOverlap_ += 1;
        }
    }
    nDimer_ += ndimer;
    // The Boltzmann factors of the dimers at one temperature
    std::vector<double> g0_12(ndimer);
    size_t              ntemp = 0;
    for(size_t iTemp = 0; iTemp < temperatures_.size(); iTemp++)
    {
        auto T = temperatures_[iTemp];
        if (T == 0)
        {
            fprintf(stderr, "Please provide a finite temperature to compute second virial.\n");
            continue;
        }
        double beta = 1.0/(BOLTZ*T);
        ntemp      += 1;
        auto  &eU12 = exp_U12_[iTemp];
        auto  &eF0  = exp_F2_[0][iTemp];
        auto  &eF1  = exp_F2_[1][iTemp];
        auto  &eT0  = exp_tau_[0][iTemp];
        auto  &eT1  = exp_tau_[1][iTemp];
        auto  &nU12 = n_U12_[iTemp];
        for(size_t ii = 0; ii < ndimer; ii++)
        {
            // Gray and Gubbins Eqn. 3.261
            double yb = energy[ii]*beta;
            g0_12[ii] = (!overlap[ii] && yb < 700) ? std::exp(-yb) : 0;
        }
        // Terms with a Boltzmann factor of zero add zero, the sums are
        // never negative zero, so adding zero does not change them.
        auto add = [&](size_t ii, size_t j)
        {
            double g    = g0_12[ii];
            bool   some = g > 0;
            double w    = weight[ii];
            eU12[j] += w*(g-1);
            // Gray and Gubbins Eqn. 3.272
            double ipf0 = some ? g*force2[0][ii] : 0.0;
            double ipf1 = some ? g*force2[1][ii] : 0.0;
            eF0[j] += w*ipf0;
            eF1[j] += w*ipf1;
            for(int m = 0; m < DIM; m++)
            {
                // Gray and Gubbins Eqn. 3.282
                double tau0 = some ? g*torque[0][ii][m]*torque[0][ii][m] : 0.0;
                double tau1 = some ? g*torque[1][ii][m]*torque[1][ii][m] : 0.0;
                eT0[j][m] += w*tau0;
                eT1[j][m] += w*tau1;
            }
            nU12[j] += w;
        };
        if (inOrder)
        {
            // Without indirect indexing the loop can be vectorized
            for(size_t ii = 0; ii < ndimer; ii++)
            {
                add(ii, ii);
            }
        }
        else
        {
            for(size_t ii = 0; ii < ndimer; ii++)
            {
                add(ii, index[ii]);
            }
        }
    }
    // Per dimer and temperature: the Boltzmann factor with the exponential
    // 27, the Mayer function 3, the forces 6, the torques 24 and the weight 1,
    // see Timings::addFlops
    Timings::instance().addFlops(61.0*ndimer*ntemp);
}

void B2Data::integrate(MsgHandler *msghandler,
                       int iTemp, double beta,
                       const std::vector<double>    &mass,
                       const std::vector<gmx::RVec> &inertia,
                       double *Bclass, double *BqmForce,
                       double *BqmTorque1, double *BqmTorque2)
{
    // Constant for quantum corrections
    double    hbarfac  = beta*gmx::square(beta*PLANCK/(2*M_PI))/24;
    // Initiate arrays
    mayer_[iTemp].resize(n_U12_[iTemp].size(), -1);
    std::vector<double> myBforce(n_U12_[iTemp].size(), 0);
    std::vector<double> myBtorq[2];
    for(int kk = 0; kk < 2; kk++)
    {
        myBtorq[kk].resize(n_U12_[iTemp].size(), 0);
    }
    std::vector<double> myBtorq2(n_U12_[iTemp].size(), 0);
    for(size_t ii = 0; ii < n_U12_[iTemp].size(); ii++)
    {
        if (n_U12_[iTemp][ii] > 0)
        {
            double Unew = exp_U12_[iTemp][ii]/n_U12_[iTemp][ii];
            mayer_[iTemp][ii] = Unew;
            for(int kk = 0; kk < 2; kk++)
            {
                // Weighted square force
                // We follow Eqn. 9 in Schenter, JCP 117 (2002) 6573
                myBforce[ii] += exp_F2_[kk][iTemp][ii]/(mass[kk]*n_U12_[iTemp][ii]);
                // Contributions from torque
                for(int m = 0; m < DIM; m++)
                {
                    if (inertia[kk][m] > 0)
                    {
                        myBtorq[kk][ii] += exp_tau_[kk][iTemp][ii][m]/(n_U12_[iTemp][ii]*inertia[kk][m]);
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


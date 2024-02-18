/*
 * This source file is part of the Alexandria Chemistry Toolkit.
 *
 * Copyright (C) 2024
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
 
#include "gromacs/fileio/oenv.h"
#include "gromacs/fileio/xvgr.h"
#include "gromacs/math/units.h"
#include "gromacs/math/vec.h"
#include "gromacs/utility/stringutil.h"

namespace alexandria
{

double sphereIntegrator(double r1, double r2, double val1, double val2)
{
    // Approximate trapezium by y = ax + b (a == slope)
    // then integrate that multiplied by x^2 to get
    // a/4 x^4 + b/3 x^3
    // insert old and new point
    double a        = (val2-val1)/(r2-r1);
    double b        = val1 - a*r1;
    double integral = ((a/4)*(r2*r2*r2*r2 - r1*r1*r1*r1) + 
                       (b/3)*(r2*r2*r2 - r1*r1*r1));
    return 4*M_PI*integral;
}

B2Data::B2Data(int                        nbins,
               double                     binwidth,
               const std::vector<double> &temperatures)
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
    for(int i = 0; i < nbins; i++)
    {
        dist_.push_back(i*binwidth);
    }
}

void B2Data::dump(FILE *fp) const
{
    for(size_t i = 0; i < n_U12_[0].size(); i++)
    {
        fprintf(fp, " %d", n_U12_[0][i]);
    }
    fprintf(fp, "\n");
}

void B2Data::aggregate(CommunicationRecord *cr)
{
    if (cr->isMaster())
    {
        for(int src = 1; src < cr->size(); src++)
        {
            if (CommunicationStatus::RECV_DATA == cr->recv_data(src))
            {
                for(size_t i = 0; i < exp_U12_.size(); i++)
                {
                    std::vector<double> d;
                    cr->recv_double_vector(src, &d);
                    for(size_t j = 0; j < d.size(); j++)
                    {
                        exp_U12_[i][j] += d[j];
                    }
                }
                for(size_t i = 0; i < n_U12_.size(); i++)
                {
                    for(size_t j = 0; j < n_U12_[i].size(); j++)
                    {
                        n_U12_[i][j] += cr->recv_int(src);
                    }
                }
                for(int kk = 0; kk < 2; kk++)
                {
                    for(size_t i = 0; i < exp_F2_[kk].size(); i++)
                    {
                        std::vector<double> d;
                        cr->recv_double_vector(src, &d);
                        for(size_t j = 0; j < d.size(); j++)
                        {
                            exp_F2_[kk][i][j] += d[j];
                        }
                    }
                    for(size_t i = 0; i < exp_tau_[kk].size(); i++)
                    {
                        for(size_t j = 0; j < exp_tau_[kk][i].size(); j++)
                        {
                            std::vector<double> d;
                            cr->recv_double_vector(src, &d);
                            for(size_t m = 0; m < d.size(); m++)
                            {
                                exp_tau_[kk][i][j][m] += d[m];
                            }
                        }
                    }
                }
            }
        }
    }
    else
    {
        int dest = cr->superior();
        if (CommunicationStatus::SEND_DATA == cr->send_data(dest))
        {
            for(size_t i = 0; i < exp_U12_.size(); i++)
            {
                std::vector<double> d = exp_U12_[i];
                cr->send_double_vector(dest, &d);
            }
            for(size_t i = 0; i < n_U12_.size(); i++)
            {
                for(size_t j = 0; j < n_U12_[i].size(); j++)
                {
                    cr->send_int(dest, n_U12_[i][j]);
                }
                }
            for(int kk = 0; kk < 2; kk++)
            {
                for(size_t i = 0; i < exp_F2_[kk].size(); i++)
                {
                    std::vector<double> d = exp_F2_[kk][i];
                    cr->send_double_vector(dest, &d);
                }
                for(size_t i = 0; i < exp_tau_[kk].size(); i++)
                {
                    for(size_t j = 0; j < exp_tau_[kk][i].size(); j++)
                    {
                        std::vector<double> d(DIM);
                        for(int m = 0; m < DIM; m++)
                        {
                            d[m] = exp_tau_[kk][i][j][m];
                        }
                        cr->send_double_vector(dest, &d);
                    }
                }
            }
        }
    }
}
    
void B2Data::addData(int iTemp, int index,
                     double exp_U12, double exp_F0, double exp_F1,
                     gmx::RVec exp_tau0, gmx::RVec exp_tau1)
{
    exp_U12_[iTemp][index]    += exp_U12;
    exp_F2_[0][iTemp][index]  += exp_F0;
    exp_F2_[1][iTemp][index]  += exp_F1;
    rvec_inc(exp_tau_[0][iTemp][index], exp_tau0);
    rvec_inc(exp_tau_[1][iTemp][index], exp_tau1);
    n_U12_[iTemp][index]      += 1;
}

void B2Data::fillToXmin(int iTemp, double xmin, double binWidth)
{
    size_t jj = 0;
    while(jj*binWidth < xmin && jj < exp_U12_[iTemp].size())
    {
        exp_U12_[iTemp][jj] = -1;
        n_U12_[iTemp][jj]   = 1;
        jj += 1;
    }
}

void B2Data::integrate(int iTemp, double binWidth, double beta,
                       const std::vector<double>    &mass,
                       const std::vector<gmx::RVec> &inertia,
                       double *Bclass, double *BqmForce,
                       double *BqmTorque1, double *BqmTorque2)
{
    // Starting force
    double    Fprev[2] =  { 0, 0 };
    // Starting torque
    gmx::RVec Tprev[2] = { { 0, 0, 0 }, { 0, 0, 0 } };
    double    hbarfac  = beta*gmx::square(beta*PLANCK/(2*M_PI))/24;
    double    r1       = 0;
    double    Uprev    = -1;
    *Bclass = *BqmForce = *BqmTorque1 = *BqmTorque2 = 0;
    mayer_[iTemp].push_back(Uprev);
    for(size_t ii = 1; ii < n_U12_[iTemp].size(); ii++)
    {
        double r2 = ii*binWidth;
        if (n_U12_[iTemp][ii] > 0)
        {
            double Unew = exp_U12_[iTemp][ii]/n_U12_[iTemp][ii];
            mayer_[iTemp].push_back(Unew);
            auto dB       = sphereIntegrator(r1, r2, Uprev, Unew);
            // TODO: There is factor 0.5 here
            *Bclass       -= 0.5*dB;
            Uprev         = Unew;
            for(int kk = 0; kk < 2; kk++)
            {
                // Weighted square force
                // We follow Eqn. 9 in Schenter, JCP 117 (2002) 6573
                double Fnew  = exp_F2_[kk][iTemp][ii]/(mass[kk]*n_U12_[iTemp][ii]);
                *BqmForce    += 0.5*hbarfac*sphereIntegrator(r1, r2, Fprev[kk], Fnew);
                Fprev[kk]    = Fnew;
                // Contributions from torque
                double bt[2] = { 0, 0 };
                for(int m = 0; m < DIM; m++)
                {
                    if (inertia[kk][m] > 0)
                    {
                        double Tnew   = exp_tau_[kk][iTemp][ii][m]/(n_U12_[iTemp][ii]*inertia[kk][m]);
                        bt[kk]       += hbarfac*sphereIntegrator(r1, r2, Tprev[kk][m], Tnew);
                        Tprev[kk][m]  = Tnew;
                    }
                }
                *BqmTorque1 += bt[0];
                *BqmTorque2 += bt[1];
            }
        }
        r1 = r2;
    }
}

void B2Data::plotMayer(const char                *fmayer,
                       gmx_output_env_t          *oenv)
{
    if (!fmayer || strlen(fmayer) == 0)
    {
        return;
    }
    GMX_RELEASE_ASSERT(mayer_.size() == temperatures_.size()+1, "Mismatch between Mayer and temperatures");
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


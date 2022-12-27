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
#include <map>

#include "forcecomputerimpl.h"

#include "act/coulombintegrals/coulombintegrals.h"
#include "act/poldata/forcefieldparametername.h"
#include "act/forces/forcecomputerutils.h"

#include "gromacs/gmxlib/nonbonded/nb_generic.h"
#include "gromacs/topology/ifunc.h"
#include "gromacs/math/vec.h"
#include "gromacs/math/units.h"

namespace alexandria
{

static void computeLJ(const std::vector<TopologyEntry *>    &pairs,
                      gmx_unused const std::vector<ActAtom> &atoms,
                      const std::vector<gmx::RVec>          *coordinates,
                      std::vector<gmx::RVec>                *forces,
                      std::map<InteractionType, double>     *energies)
{
    double ebond = 0;

    auto   x     = *coordinates;
    auto  &f     = *forces;
    for (const auto &b : pairs)
    {
        // Get the parameters. We have to know their names to do this.
        auto &params    = b->params();
        auto c6         = params[ljC6_IJ];
        auto c12        = params[ljC12_IJ];
        // Get the atom indices
        auto &indices   = b->atomIndices();
        auto ai         = indices[0];
        auto aj         = indices[1];
        rvec dx;
        rvec_sub(x[ai], x[aj], dx);
        auto dr2        = iprod(dx, dx);
        auto rinv       = gmx::invsqrt(dr2);
        auto rinv2      = rinv*rinv;
        auto rinv6      = rinv2*rinv2*rinv2;
        auto vvdw_disp  = c6*rinv6;
        auto vvdw_rep   = c12*rinv6*rinv6;
        auto elj        = vvdw_rep - vvdw_disp;
        auto flj        = (12*vvdw_rep - 6*vvdw_disp)*rinv2;
        if (debug)
        {
            fprintf(debug, "ACT ai %d aj %d vvdw: %10g c6: %10g c12: %10g\n", ai, aj, elj, c6, c12);
        }
        
        ebond      += elj;
        for (int m = 0; (m < DIM); m++)
        {
            auto fij          = flj*dx[m];
            f[indices[0]][m] += fij;
            f[indices[1]][m] -= fij;
        }
    }
    if (debug)
    {
        fprintf(debug, "ACT vvdwtot: %10g \n", ebond);
    }
    energies->insert({InteractionType::VDW, ebond});
}

static void computeWBH(const std::vector<TopologyEntry *>    &pairs,
                       gmx_unused const std::vector<ActAtom> &atoms,
                       const std::vector<gmx::RVec>          *coordinates,
                       std::vector<gmx::RVec>                *forces,
                       std::map<InteractionType, double>                 *energies)
{
    double erep  = 0;
    double edisp = 0;
    auto   x     = *coordinates;
    auto  &f     = *forces;
    for (const auto &b : pairs)
    {
        // Get the parameters. We have to know their names to do this.
        auto &params    = b->params();
        auto sigma      = params[wbhSIGMA_IJ];
        auto epsilon    = params[wbhEPSILON_IJ];
        auto gamma      = params[wbhGAMMA_IJ];
        if (epsilon > 0 && gamma > 0 && sigma > 0)
        {
            // Get the atom indices
            auto &indices   = b->atomIndices();
            rvec dx;
            rvec_sub(x[indices[0]], x[indices[1]], dx);
            auto dr2        = iprod(dx, dx);
            auto rinv       = gmx::invsqrt(dr2);
            real eerep = 0, eedisp = 0, fwbh = 0;
            wang_buckingham(sigma, epsilon, gamma, dr2, rinv, &eerep, &eedisp, &fwbh);
            erep     += eerep;
            edisp    += eedisp;
            real fbond  = fwbh*rinv;
            for (int m = 0; (m < DIM); m++)
            {
                auto fij          = fbond*dx[m];
                f[indices[0]][m] += fij;
                f[indices[1]][m] -= fij;
            }
        }
    }
    //    energies->insert({InteractionType::VDW, erep+edisp});
    energies->insert({InteractionType::REPULSION, erep});
    energies->insert({InteractionType::DISPERSION, edisp});
}

static void computeNonBonded(const std::vector<TopologyEntry *>    &pairs,
                             gmx_unused const std::vector<ActAtom> &atoms,
                             const std::vector<gmx::RVec>          *coordinates,
                             std::vector<gmx::RVec>                *forces,
                             std::map<InteractionType, double>     *energies)
{
    double erep  = 0;
    double edisp = 0;
    auto   x     = *coordinates;
    auto  &f     = *forces;
    for (const auto &b : pairs)
    {
        // Get the parameters. We have to know their names to do this.
        auto &params    = b->params();
        auto rmin       = params[gbhRMIN_IJ];
        auto epsilon    = params[gbhEPSILON_IJ];
        auto gamma      = params[gbhGAMMA_IJ];
        auto delta      = params[gbhGAMMA_IJ];
        if (epsilon > 0 && gamma > 0 && rmin > 0 && delta > 0)
        {
            // Get the atom indices
            auto &indices   = b->atomIndices();
            rvec dx;
            rvec_sub(x[indices[0]], x[indices[1]], dx);
            auto dr2        = iprod(dx, dx);
            auto rinv       = gmx::invsqrt(dr2);

            real delta6     = 6+delta;
            real delta6gam2 = delta6 + 2*gamma;
            real rstar      = dr2*rinv/rmin;
            real sixterm    = 1 + std::pow(rstar,6);
            real delterm    = 1 + std::pow(rstar,delta);
            real expterm    = std::exp(gamma*(1 - rstar));
            real sixdenom   = 1/(2*gamma*sixterm);
            real eerep      = epsilon*delta6*expterm*sixdenom;
            real eedisp     = -epsilon*(delta6gam2*sixdenom + 1/delterm);
            real fgbham     = (epsilon*((-6*(6 + delta - (6 + delta)*std::exp(gamma - gamma*rstar) + 2*gamma)*std::pow(rstar,5))/(gamma*std::pow(1 + std::pow(rstar,6),2)) + 
                                        ((6 + delta)*std::exp(gamma - gamma*rstar))/(1 + std::pow(rstar,6)) - (2*delta*std::pow(rstar,-1 + delta))/std::pow(1 + std::pow(rstar,delta),2)))/2.;

            erep     += eerep;
            edisp    += eedisp;
            real fbond  = fgbham*rinv;
            for (int m = 0; (m < DIM); m++)
            {
                auto fij          = fbond*dx[m];
                f[indices[0]][m] += fij;
                f[indices[1]][m] -= fij;
            }
        }
    }
    energies->insert({InteractionType::REPULSION, erep});
    energies->insert({InteractionType::DISPERSION, edisp});
}

static void computeCoulomb(const std::vector<TopologyEntry *> &pairs,
                           const std::vector<ActAtom>         &atoms,
                           const std::vector<gmx::RVec>       *coordinates,
                           std::vector<gmx::RVec>             *forces,
                           std::map<InteractionType, double>              *energies)
{
    double ebond = 0;
    auto   x     = *coordinates;
    auto  &f     = *forces;
    for (const auto &b : pairs)
    {
        // Get the parameters. We have to know their names to do this.
        auto ai     = b->atomIndices()[0];
        auto aj     = b->atomIndices()[1];
        auto &params= b->params();
        auto izeta  = params[coulZETAI];
        auto jzeta  = params[coulZETAJ];
        // Get the atom indices
        real qq         = ONE_4PI_EPS0*atoms[ai].charge()*atoms[aj].charge();
        rvec dx;
        rvec_sub(x[ai], x[aj], dx);
        auto dr2        = iprod(dx, dx);
        real velec, felec;
        coulomb_gaussian(qq, izeta, jzeta, std::sqrt(dr2), &velec, &felec);
        ebond += velec;
        if (dr2 > 0)
        {
            felec *= gmx::invsqrt(dr2);
            for (int m = 0; (m < DIM); m++)
            {
                auto fij  = felec*dx[m];
                f[ai][m] += fij;
                f[aj][m] -= fij;
            }
        }
    }
    energies->insert({InteractionType::COULOMB, ebond});
}

#ifdef FUTURE
static void computePartridge(const std::vector<TopologyEntry *>    &angles,
                             gmx_unused const std::vector<ActAtom> &atoms,
                             const std::vector<gmx::RVec>          *coordinates,
                             std::vector<gmx::RVec>                *forces,
                             std::map<InteractionType, double>                 *energies)
{
    // Energy function according to 
    // The determination of an accurate isotope dependent potential energy 
    // surface for water from extensive ab initio calculations and 
    // experimental data
    // Harry Partridge and David W. Schwenke
    // J Chem Phys 106 (1997) p. 4618
    double  energy = 0, costh = 0;
    auto    x     = *coordinates;
    auto   &f     = *forces;
    const  real half = 0.5;
    for (const auto &a : angles)
    {
        // Get the parameters. We have to know their names to do this.
        auto &params    = a->params();
        auto theta0     = params[psANGLE];
        auto rij0       = params[psRIJ0];
        auto rjk0       = params[psRJK0];
        // Get the atom indices
        auto &indices   = a->atomIndices();

        rvec r_ij, r_kj;
        auto theta = bond_angle(x[indices[0]], x[indices[1]], x[indices[2]], 
                                r_ij, r_kj, &costh);
        // Compute the components of the potential
        auto costh0 = std::cos(theta0);
        auto drij   = rij0-norm(r_ij);
        auto drkj   = rjk0-norm(r_kj);
        
        auto costh2 = gmx::square(costh);
        auto st     = theta0*gmx::invsqrt(1 - costh2);   
        auto sth    = st*costh; 
        
        auto nrij2 = iprod(r_ij, r_ij);                 
        auto nrkj2 = iprod(r_kj, r_kj);                 
        
        auto nrij_1 = gmx::invsqrt(nrij2);              
        auto nrkj_1 = gmx::invsqrt(nrkj2);              
        
        auto cik = st*nrij_1*nrkj_1;                    
        auto cii = sth*nrij_1*nrij_1;                   
        auto ckk = sth*nrkj_1*nrkj_1;                  
        
        for (auto m = 0; m < DIM; m++)
        {
            auto f_im    = -(cik*r_kj[m] - cii*r_ij[m]);
            auto f_km    = -(cik*r_ij[m] - ckk*r_kj[m]);
            auto f_jm    = -f_im - f_km;
            f[indices[0]][m] += f_im;
            f[indices[1]][m] += f_km;
            f[indices[2]][m] += f_jm;
        }                                          
    }
    energies->insert({InteractionType::ANGLES, , energy});
}
#endif

static void harmonic(real k, real x0, real x, real *V, real *F)
{
    const real half = 0.5;

    real dx  = x-x0;
    real dx2 = dx*dx;

    *F  = -k*dx;
    *V  = half*k*dx2;
}

static void computeDummy(gmx_unused const std::vector<TopologyEntry *>    &bonds,
                         gmx_unused const std::vector<ActAtom> &atoms,
                         gmx_unused const std::vector<gmx::RVec>          *coordinates,
                         gmx_unused std::vector<gmx::RVec>                *forces,
                         gmx_unused std::map<InteractionType, double>                 *energies)
{
}

static void computeBonds(const std::vector<TopologyEntry *>    &bonds,
                         gmx_unused const std::vector<ActAtom> &atoms,
                         const std::vector<gmx::RVec>          *coordinates,
                         std::vector<gmx::RVec>                *forces,
                         std::map<InteractionType, double>     *energies)
{
    if (nullptr == coordinates || nullptr == forces)
    {
        GMX_THROW(gmx::InternalError("nullptr"));
    }
    double ebond = 0;
    auto   x     = *coordinates;
    auto  &f     = *forces;
    
    for (const auto &b : bonds)
    {
        // Get the parameters. We have to know their names to do this.
        auto &params    = b->params();
        auto bondlength = params[bondLENGTH];
        auto kb         = params[bondKB];
        auto D0         = params[bondENERGY];
        // Get the atom indices
        auto &indices   = b->atomIndices();

        rvec dx;
        rvec_sub(x[indices[0]], x[indices[1]], dx);
        auto dr2        = iprod(dx, dx);
        auto r_1        = gmx::invsqrt(dr2);
        auto r1         = dr2*r_1;

        real vB, fbond;
        harmonic(kb, bondlength, r1, &vB, &fbond);
        
        fbond *= r_1;
        ebond += vB+D0;
        
        for (int m = 0; (m < DIM); m++)
        {
            auto fij          = fbond*dx[m];
            f[indices[0]][m] += fij;
            f[indices[1]][m] -= fij;
        }
    }
    energies->insert({InteractionType::BONDS, ebond});
}

static void computeCubic(const std::vector<TopologyEntry *>    &bonds,
                         gmx_unused const std::vector<ActAtom> &atoms,
                         const std::vector<gmx::RVec>          *coordinates,
                         std::vector<gmx::RVec>                *forces,
                         std::map<InteractionType, double>     *energies)
{
    if (nullptr == coordinates || nullptr == forces)
    {
        GMX_THROW(gmx::InternalError("nullptr"));
    }
    double ebond = 0;
    auto   x     = *coordinates;
    auto  &f     = *forces;
    
    for (const auto &b : bonds)
    {
        // Get the parameters. We have to know their names to do this.
        auto &params    = b->params();
        auto bondlength = params[cubicLENGTH];
        auto rmax       = params[cubicRMAX];
        auto kb         = params[cubicKB];
        auto De         = params[cubicDE];
        // Get the atom indices
        auto &indices   = b->atomIndices();

        rvec dx;
        rvec_sub(x[indices[0]], x[indices[1]], dx);
        auto r2  = iprod(dx, dx);
        auto r_1 = gmx::invsqrt(r2);
        auto r1  = r2*r_1;

        real vB, fbond;
        vB    = kb*gmx::square(r1-bondlength)*(rmax-r1) - De;
        fbond = r_1*kb*(bondlength+2*rmax-3*r1)*(bondlength-r1);
        
        ebond += vB;
        
        for (int m = 0; (m < DIM); m++)
        {
            auto fij          = fbond*dx[m];
            f[indices[0]][m] += fij;
            f[indices[1]][m] -= fij;
        }
    }
    energies->insert({InteractionType::BONDS, ebond});
}

static void computeMorse(const std::vector<TopologyEntry *>    &bonds,
                         gmx_unused const std::vector<ActAtom> &atoms,
                         const std::vector<gmx::RVec>          *coordinates,
                         std::vector<gmx::RVec>                *forces,
                         std::map<InteractionType, double>     *energies)
{
    double  ebond = 0;
    auto    x     = *coordinates;
    auto   &f     = *forces;
    for (const auto &b : bonds)
    {
        // Get the parameters. We have to know their names to do this.
        auto &params    = b->params();
        auto bondlength = params[morseLENGTH];
        auto beta       = params[morseBETA];
        auto De         = params[morseDE];
        auto D0         = params[morseD0];
        // Get the atom indices
        auto &indices   = b->atomIndices();
        rvec dx;
        rvec_sub(x[indices[0]], x[indices[1]], dx);
        auto dr2        = iprod(dx, dx);
        auto dr         = std::sqrt(dr2);
        auto temp       = std::exp(-beta*(dr-bondlength));
        auto omtemp     = 1.0-temp;
        auto cbomtemp   = De*omtemp;
        auto vbond      = D0+cbomtemp*omtemp;
        auto fbond      = -2.0*beta*temp*cbomtemp/dr;
        ebond          += vbond;

        for (int m = 0; (m < DIM); m++)
        {
            auto fij          = fbond*dx[m];
            f[indices[0]][m] += fij;
            f[indices[1]][m] -= fij;
        }
    }
    energies->insert({InteractionType::BONDS, ebond});
}

static void computeLinearAngles(const std::vector<TopologyEntry *>    &angles,
                                gmx_unused const std::vector<ActAtom> &atoms,
                                const std::vector<gmx::RVec>          *coordinates,
                                std::vector<gmx::RVec>                *forces,
                                std::map<InteractionType, double>                 *energies)
{
    double  ebond = 0;
    auto    x     = *coordinates;
    auto   &f     = *forces;
    for (const auto &aaa : angles)
    {
        // Get the parameters. We have to know their names to do this.
        auto &params= aaa->params();
        auto a      = params[linangA];
        auto klin   = params[linangKLIN];
        // Get the atom indices
        auto &indices= aaa->atomIndices();
        
        rvec r_ij, r_kj, r_ik;
        rvec_sub(x[indices[0]], x[indices[1]], r_ij);
        rvec_sub(x[indices[2]], x[indices[1]], r_kj);
        rvec_sub(r_ij, r_kj, r_ik);

        real dr2 = 0;
        real b   = 1-a;
        for (int m = 0; (m < DIM); m++)
        {
            // Minus sign due to F = -dV/dr
            real dr   = -a * r_ij[m] - b * r_kj[m];
            dr2      += dr*dr;
            real f_i  = a*klin*dr;
            real f_k  = b*klin*dr;
            real f_j  = -(f_i+f_k);
            f[indices[0]][m] += f_i;
            f[indices[1]][m] += f_j;
            f[indices[2]][m] += f_k;
        }
        ebond += 0.5*klin*dr2;
    }
    energies->insert({InteractionType::LINEAR_ANGLES, ebond});
}

// It is not finished yes, needs to be double checked. 
static void computeAngles(const std::vector<TopologyEntry *>    &angles,
                          gmx_unused const std::vector<ActAtom> &atoms,
                          const std::vector<gmx::RVec>          *coordinates,
                          std::vector<gmx::RVec>                *forces,
                          std::map<InteractionType, double>                 *energies)
{
    double  energy = 0, costh = 0;
    auto    x     = *coordinates;
    auto   &f     = *forces;

    for (const auto &a : angles)
    {
        // Get the parameters. We have to know their names to do this.
        auto &params    = a->params();
        auto theta0     = params[angleANGLE];
        auto ka         = params[angleKT];
        // Get the atom indices
        auto &indices   = a->atomIndices();

        rvec r_ij, r_kj;
        auto theta = bond_angle(x[indices[0]], x[indices[1]], x[indices[2]], 
                                r_ij, r_kj, &costh);

        real vA, fangle;
        harmonic(ka, theta0, theta, &vA, &fangle);
        
        energy += vA;
        
        auto costh2 = gmx::square(costh);
        if (costh2 < 1)
        {
            real st, sth;
            real cik, cii, ckk;
            real nrkj2, nrij2;
            real nrkj_1, nrij_1;

            st    = fangle*gmx::invsqrt(1 - costh2);   
            sth   = st*costh;                      

            nrij2 = iprod(r_ij, r_ij);                 
            nrkj2 = iprod(r_kj, r_kj);                 

            nrij_1 = gmx::invsqrt(nrij2);              
            nrkj_1 = gmx::invsqrt(nrkj2);              

            cik = st*nrij_1*nrkj_1;                    
            cii = sth*nrij_1*nrij_1;                   
            ckk = sth*nrkj_1*nrkj_1;                  

            for (auto m = 0; m < DIM; m++)
            {
                auto f_im    = -(cik*r_kj[m] - cii*r_ij[m]);
                auto f_km    = -(cik*r_ij[m] - ckk*r_kj[m]);
                auto f_jm    = -f_im - f_km;
                f[indices[0]][m] += f_im;
                f[indices[1]][m] += f_jm;
                f[indices[2]][m] += f_km;
            }
        }                                          
    }
    energies->insert({InteractionType::ANGLES, energy});
}

// It is not finished yes, needs to be double checked. 
static void computeUreyBradley(const std::vector<TopologyEntry *>    &angles,
                               gmx_unused const std::vector<ActAtom> &atoms,
                               const std::vector<gmx::RVec>          *coordinates,
                               std::vector<gmx::RVec>                *forces,
                               std::map<InteractionType, double>                 *energies)
{
    double  energy = 0, costh = 0;
    auto    x     = *coordinates;
    auto   &f     = *forces;
    const  real half = 0.5;
    for (const auto &a : angles)
    {
        // Get the parameters. We have to know their names to do this.
        auto &params    = a->params();
        auto theta0     = params[ubANGLE];
        auto ka         = params[ubKT];
        auto r13        = params[ubR13];
        auto kUB        = params[ubKUB];
        // Get the atom indices
        auto &indices   = a->atomIndices();
        auto ai = indices[0];
        auto aj = indices[1];
        auto ak = indices[2];
        
        rvec r_ij, r_kj;
        auto theta = bond_angle(x[ai], x[aj], x[ak], r_ij, r_kj, &costh);

        // Compute deviation from the reference angle
        auto da  = theta - theta0;
        auto da2 = da*da;
        
        auto fangle      = -ka*da;
        energy          += half*ka*da2;
        
        // Now UB part
        rvec r_ik;
        rvec_sub(x[ai], x[ak], r_ik);
        auto rik2        = iprod(r_ik, r_ik);
        auto norm_rik    = std::sqrt(rik2);
        
        real vUB, fUB;
        harmonic(kUB, r13, norm_rik, &vUB, &fUB);

        energy += vUB;
        fUB    *= gmx::invsqrt(rik2);

        auto costh2 = gmx::square(costh);
        if (costh2 < 1)
        {
            real st, sth;
            real cik, cii, ckk;
            real nrkj2, nrij2;
            real nrkj_1, nrij_1;

            st    = fangle*gmx::invsqrt(1 - costh2);   
            sth   = st*costh;                      

            nrij2 = iprod(r_ij, r_ij);                 
            nrkj2 = iprod(r_kj, r_kj);                 

            nrij_1 = gmx::invsqrt(nrij2);              
            nrkj_1 = gmx::invsqrt(nrkj2);              

            cik = st*nrij_1*nrkj_1;                    
            cii = sth*nrij_1*nrij_1;                   
            ckk = sth*nrkj_1*nrkj_1;                  

            for (auto m = 0; m < DIM; m++)
            {
                real f_i = -(cik*r_kj[m] - cii*r_ij[m]);
                real f_k = -(cik*r_ij[m] - ckk*r_kj[m]);
                real f_j = -f_i - f_k;
                f[ai][m] += f_i + fUB*r_ik[m];
                f[aj][m] += f_j;
                f[ak][m] += f_k - fUB*r_ik[m];
            }
        }                                          
    }
    energies->insert({InteractionType::ANGLES, energy});
}

static void computePolarization(const std::vector<TopologyEntry *>    &bonds,
                                gmx_unused const std::vector<ActAtom> &atoms,
                                const std::vector<gmx::RVec>          *coordinates,
                                std::vector<gmx::RVec>                *forces,
                                std::map<InteractionType, double>                 *energies)
{
    double ebond = 0;
    auto   x     = *coordinates;
    auto  &f     = *forces;
    const  real half = 0.5;
    for (const auto &b : bonds)
    {
        // Get the parameters. We have to know their names to do this.
        auto &params   = b->params();
        auto ksh       = params[polKSH];
        // Get the atom indices
        auto &indices  = b->atomIndices();
        rvec dx;
        rvec_sub(x[indices[0]], x[indices[1]], dx);
        auto dr2        = iprod(dx, dx);
        
        auto fbond      = -ksh;
        ebond          += half*ksh*dr2;
        
        for (int m = 0; m < DIM; m++)
        {
            auto fij          = fbond*dx[m];
            f[indices[0]][m] += fij;
            f[indices[1]][m] -= fij;
        }
    }
    energies->insert({InteractionType::POLARIZATION, ebond});

}

static void do_dih_fup_noshiftf(int i, int j, int k, int l, real ddphi,
                                rvec r_ij, rvec r_kj, rvec r_kl,
                                rvec m, rvec n,
                                std::vector<gmx::RVec> *forces)
{
    rvec f_i, f_j, f_k, f_l;
    rvec uvec, vvec, svec;
    real iprm, iprn, nrkj, nrkj2, nrkj_1, nrkj_2;
    real a, b, p, q, toler;

    iprm  = iprod(m, m);       /*  5    */
    iprn  = iprod(n, n);       /*  5	*/
    nrkj2 = iprod(r_kj, r_kj); /*  5	*/
    toler = nrkj2*GMX_REAL_EPS;
    if ((iprm > toler) && (iprn > toler))
    {
        nrkj_1 = gmx::invsqrt(nrkj2); /* 10	*/
        nrkj_2 = nrkj_1*nrkj_1;       /*  1	*/
        nrkj   = nrkj2*nrkj_1;        /*  1	*/
        a      = -ddphi*nrkj/iprm;    /* 11	*/
        svmul(a, m, f_i);             /*  3	*/
        b     = ddphi*nrkj/iprn;      /* 11	*/
        svmul(b, n, f_l);             /*  3  */
        p     = iprod(r_ij, r_kj);    /*  5	*/
        p    *= nrkj_2;               /*  1	*/
        q     = iprod(r_kl, r_kj);    /*  5	*/
        q    *= nrkj_2;               /*  1	*/
        svmul(p, f_i, uvec);          /*  3	*/
        svmul(q, f_l, vvec);          /*  3	*/
        rvec_sub(uvec, vvec, svec);   /*  3	*/
        rvec_sub(f_i, svec, f_j);     /*  3	*/
        rvec_add(f_l, svec, f_k);     /*  3	*/
        rvec_inc((*forces)[i], f_i);          /*  3	*/
        rvec_dec((*forces)[j], f_j);          /*  3	*/
        rvec_dec((*forces)[k], f_k);          /*  3	*/
        rvec_inc((*forces)[l], f_l);          /*  3	*/
    }
}

static real dih_angle(const rvec xi, const rvec xj, const rvec xk, const rvec xl,
                      rvec r_ij, rvec r_kj, rvec r_kl, rvec m, rvec n)
{
    rvec_sub(xi, xj, r_ij); /*  3        */
    rvec_sub(xk, xj, r_kj); /*  3		*/
    rvec_sub(xk, xl, r_kl); /*  3		*/

    cprod(r_ij, r_kj, m);                  /*  9        */
    cprod(r_kj, r_kl, n);                  /*  9		*/
    real phi  = gmx_angle(m, n);           /* 49 (assuming 25 for atan2) */
    real ipr  = iprod(r_ij, n);            /*  5        */
    real sign = (ipr < 0.0) ? -1.0 : 1.0;
    phi       = sign*phi;                  /*  1		*/
    /* 82 TOTAL	*/
    return phi;
}

static void computeFourDihs(const std::vector<TopologyEntry *>    &propers,
                            gmx_unused const std::vector<ActAtom> &atoms,
                            const std::vector<gmx::RVec>          *coordinates,
                            std::vector<gmx::RVec>                *forces,
                            std::map<InteractionType, double>                 *energies)
{
    double energy = 0;
    auto   x      = *coordinates;
    for (const auto &a : propers)
    {
        // Get the parameters. We have to know their names to do this.
        auto &params  = a->params();
        // Assume we have the parameters in the correct order 0-nparam-1
        size_t nparam = params.size();
        // Get the atom indices
        auto &indices = a->atomIndices();
        auto ai       = indices[0];
        auto aj       = indices[1];
        auto ak       = indices[2];
        auto al       = indices[3];

        rvec r_ij, r_kj, r_kl, m, n;
        auto phi = dih_angle(x[ai], x[aj], x[ak], x[al],
                             r_ij, r_kj, r_kl, m, n);

        real cos_phi = std::cos(phi);
        
        /* Calculate cosine powers */
        /* Calculate the energy */
        /* Calculate the derivative */

        // Energy is
        // \sum_{j=0}^{nparam} params[j] cos^j(phi)
        // Force then becomes
        // - \sum_{j=1}^{nparam} j * params[j] cos^{j-1}(phi) sin(phi)
        real cosfac = 1;
        real v      = params[0];
        real dvdphi = 0;
        for (size_t j = 1; (j < nparam); j++)
        {
            dvdphi += j*params[j]*cosfac;
            cosfac *= cos_phi;
            v      += params[j]*cosfac;
        }

        /* Beware of accuracy loss, cannot use 1-sqrt(cos^2) ! */
        dvdphi *= -std::sin(phi);

        do_dih_fup_noshiftf(ai, aj, ak, al, dvdphi, r_ij, r_kj, r_kl, m, n,
                            forces);
        energy += v;
    }
    energies->insert({InteractionType::PROPER_DIHEDRALS, energy});
}

static void computeImpropers(const std::vector<TopologyEntry *>    &impropers,
                             gmx_unused const std::vector<ActAtom> &atoms,
                             const std::vector<gmx::RVec>          *coordinates,
                             std::vector<gmx::RVec>                *forces,
                             std::map<InteractionType, double>                 *energies)
{
    double  energy = 0;
    auto    x     = *coordinates;
    const  real half = 0.5;
    for (const auto &a : impropers)
    {
        // Get the parameters. We have to know their names to do this.
        auto &params  = a->params();
        auto kA       = params[idihKPHI];
        // Get the atom indices
        auto &indices = a->atomIndices();
        auto ai       = indices[0];
        auto aj       = indices[1];
        auto ak       = indices[2];
        auto al       = indices[3];

        rvec r_ij, r_kj, r_kl, m, n;
        auto phi = dih_angle(x[ai], x[aj], x[ak], x[al],
                             r_ij, r_kj, r_kl, m, n);

        auto dp2 = phi*phi;

        energy     += half*kA*dp2;
        auto ddphi  = kA*phi;

        do_dih_fup_noshiftf(ai, aj, ak, al, ddphi, r_ij, r_kj, r_kl, m, n,
                            forces);
    }
    energies->insert({InteractionType::IMPROPER_DIHEDRALS, energy});
}

static void computePropers(const std::vector<TopologyEntry *>    &propers,
                           gmx_unused const std::vector<ActAtom> &atoms,
                           const std::vector<gmx::RVec>          *coordinates,
                           std::vector<gmx::RVec>                *forces,
                           std::map<InteractionType, double>                 *energies)
{
    double energy = 0;
    auto   x      = *coordinates;
    for (const auto &a : propers)
    {
        // Get the parameters. We have to know their names to do this.
        auto &params  = a->params();
        auto phi0     = params[pdihANGLE];
        auto kphi     = params[pdihKP];
        auto mult     = params[pdihMULT];
        // Get the atom indices
        auto &indices = a->atomIndices();
        auto ai       = indices[0];
        auto aj       = indices[1];
        auto ak       = indices[2];
        auto al       = indices[3];

        rvec r_ij, r_kj, r_kl, m, n;
        auto phi = dih_angle(x[ai], x[aj], x[ak], x[al],
                             r_ij, r_kj, r_kl, m, n);

        auto mdphi  =  mult*phi - phi0;
        auto sdphi  =  std::sin(mdphi);
        auto ddphi  = -kphi*mult*sdphi;
        auto v      = kphi*(1.0 + std::cos(mdphi));
        energy     += v;
#ifdef DEBUG
        printf("ai %d aj %d ak %d al %d phi %g phi0 %g kphi %g, mult %g, mdphi %g, v %g ddphi %g\n",
               ai, aj, ak, al, phi, phi0, kphi, mult, mdphi, v, ddphi);
#endif
        do_dih_fup_noshiftf(ai, aj, ak, al, ddphi, r_ij, r_kj, r_kl, m, n,
                            forces);
    }
    energies->insert({InteractionType::PROPER_DIHEDRALS, energy});
}

std::map<int, bondForceComputer> bondForceComputerMap = {
    { F_DUMMY,         computeDummy        },
    { F_BONDS,         computeBonds        },
    { F_MORSE,         computeMorse        },
    { F_CUBICBONDS,    computeCubic        },
    { F_ANGLES,        computeAngles       },
    { F_LINEAR_ANGLES, computeLinearAngles },
    { F_LJ,            computeLJ           },
    { F_BHAM,          computeWBH          },
    { F_GBHAM,         computeNonBonded    },
    { F_COUL_SR,       computeCoulomb      },
    { F_POLARIZATION,  computePolarization },
    { F_IDIHS,         computeImpropers    },
    { F_PDIHS,         computePropers      },
    { F_FOURDIHS,      computeFourDihs     },
    { F_UREY_BRADLEY,  computeUreyBradley  }
};

bondForceComputer getBondForceComputer(int gromacs_index)
{
    auto bfc = bondForceComputerMap.find(gromacs_index);
    if (bondForceComputerMap.end() != bfc)
    {
        return *bfc->second;
    }

    // Keep the compiler happy
    return nullptr;
}

} // namespace

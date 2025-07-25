/*
 * This source file is part of the Alexandria Chemistry Toolkit.
 *
 * Copyright (C) 2014-2025
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

#include "act/coulombintegrals/gaussian_integrals.h"
#include "act/coulombintegrals/slater_integrals.h"
#include "act/forcefield/forcefield_parametername.h"
#include "act/forces/forcecomputerutils.h"

#include "gromacs/gmxlib/nonbonded/nb_generic.h"
#include "gromacs/math/vec.h"
#include "gromacs/math/units.h"
#include "gromacs/topology/ifunc.h"

namespace alexandria
{

static inline void pairforces(double                  fscalar,
                              const gmx::RVec        &dx,
                              const std::vector<int> &indices,
                              std::vector<gmx::RVec> *forces)
{
    auto  &f = *forces;

    auto fijX          = fscalar*dx[XX];
    f[indices[0]][XX] += fijX;
    f[indices[1]][XX] -= fijX;
    auto fijY          = fscalar*dx[YY];
    f[indices[0]][YY] += fijY;
    f[indices[1]][YY] -= fijY;
    auto fijZ          = fscalar*dx[ZZ];
    f[indices[0]][ZZ] += fijZ;
    f[indices[1]][ZZ] -= fijZ;
}

static double computeLJ12_6(const TopologyEntryVector             &pairs,
                            gmx_unused const std::vector<ActAtom> &atoms,
                            const std::vector<gmx::RVec>          *coordinates,
                            std::vector<gmx::RVec>                *forces,
                            std::map<InteractionType, double>     *energies)
{
    double erep  = 0;
    double edisp = 0;
    auto   &x    = *coordinates;
    for (const auto &b : pairs)
    {
        // Get the parameters. We have to know their names to do this.
        auto &params    = b->params();
        auto sig_ij     = params[lj12_6SIGMA];
        auto eps_ij     = params[lj12_6EPSILON];
        auto sig6       = gmx::square(sig_ij*sig_ij*sig_ij);
        auto c6         = 4*eps_ij*sig6;
        auto c12        = c6*sig6;
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
        auto vvdw_disp  = -c6*rinv6;
        auto vvdw_rep   = c12*rinv6*rinv6;
        auto flj        = (12*vvdw_rep + 6*vvdw_disp)*rinv2;
        if (debug)
        {
            fprintf(debug, "ACT ai %d aj %d vvdw_rep: %10g vvdw_disp: %10g c6: %10g c12: %10g\n",
                    ai, aj, vvdw_rep, vvdw_disp, c6, c12);
        }
        pairforces(flj, dx, indices, forces);
        erep      += vvdw_rep;
        edisp     += vvdw_disp;
    }
    if (debug)
    {
        fprintf(debug, "ACT erep: %10g edisp %10g\n", erep, edisp);
    }
    energies->insert({InteractionType::EXCHANGE, erep});
    energies->insert({InteractionType::DISPERSION, edisp});

    return erep + edisp;
}

static double computeLJ8_6(const TopologyEntryVector             &pairs,
                           gmx_unused const std::vector<ActAtom> &atoms,
                           const std::vector<gmx::RVec>          *coordinates,
                           std::vector<gmx::RVec>                *forces,
                           std::map<InteractionType, double>     *energies)
{   
    double erep  = 0;
    double edisp = 0;
    auto   x     = *coordinates;
    for (const auto &b : pairs)
    {   
        // Get the parameters. We have to know their names to do this.
        auto &params    = b->params();
        auto sig_ij     = params[lj8_6SIGMA];
        auto eps_ij     = params[lj8_6EPSILON];
        auto sig6       = gmx::square(sig_ij*sig_ij*sig_ij);
        auto c6         = 4*eps_ij*sig6;
        auto c8         = 0.75*c6*sig_ij*sig_ij;
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
        auto vvdw_disp  = -c6*rinv6;     
        auto vvdw_rep   = c8*rinv6*rinv2;
        auto flj        = (8*vvdw_rep + 6*vvdw_disp)*rinv2;
        if (debug)
        {       
            fprintf(debug, "ACT ai %d aj %d vrep: %10g vdisp: %10g c6: %10g c8: %10g\n", ai, aj, vvdw_rep, vvdw_disp, c6, c8);
        }

        erep  += vvdw_rep;
        edisp += vvdw_disp;
        pairforces(flj, dx, indices, forces);
    }
    if (debug)               
    {                        
        fprintf(debug, "ACT vvdw_rep: %10g vvdw_disp: %10g\n", erep, edisp); 
    }
    energies->insert({InteractionType::EXCHANGE, erep});
    energies->insert({InteractionType::DISPERSION, edisp});

    return erep + edisp;
} 

static double computeLJ14_7(const TopologyEntryVector             &pairs,
                            gmx_unused const std::vector<ActAtom> &atoms,
                            const std::vector<gmx::RVec>          *coordinates,
                            std::vector<gmx::RVec>                *forces,
                            std::map<InteractionType, double>     *energies)
{
    double erep  = 0;
    double edisp = 0;
    auto   x     = *coordinates;
    for (const auto &b : pairs)
    {
        // Get the parameters. We have to know their names to do this.
        auto &params    = b->params();
        auto sigma      = params[lj14_7SIGMA];
        auto epsilon    = params[lj14_7EPSILON];
        real f147       = 0;
        // Get the atom indices
        auto &indices   = b->atomIndices();
        auto ai         = indices[0];
        auto aj         = indices[1];
        rvec dx; 
        rvec_sub(x[ai], x[aj], dx);
        auto dr2        = iprod(dx, dx);
        auto rinv       = gmx::invsqrt(dr2);
        real eerep = 0, eedisp = 0;
        if (epsilon > 0)
        {
            auto gamma      = params[lj14_7GAMMA];
            auto delta      = params[lj14_7DELTA];
            real rstar      = dr2*rinv/sigma;
            real delta1     = delta + 1;
            real gamma1     = gamma + 1;
            real deltars    = delta + rstar;
            real repfac     = epsilon * std::pow( (delta1/deltars), 7);
            eerep           = repfac * (gamma1/(std::pow((rstar), 7) + gamma ));
            eedisp          = -2 * repfac;
            real gamrstar7  = gamma + std::pow(rstar, 7);
            f147            = 7*epsilon*std::pow(delta1, 7)*(gamma1*std::pow(rstar, 6)*(delta+rstar)/gmx::square(gamrstar7) + gamma1/(gamrstar7) - 2)/(sigma*std::pow(delta+rstar, 8));

            if (debug)
            {
                fprintf(debug, "ACT ai %d aj %d vvdw: %10g epsilon: %10g gamma: %10g sigma: %10g delta: %10g\n",
                        ai, aj, eerep + eedisp, epsilon, gamma, sigma, delta);
            }
            erep     += eerep;
            edisp    += eedisp;
        }
        real fbond  = f147*rinv;
        pairforces(fbond, dx, indices, forces);
    }
    energies->insert({InteractionType::EXCHANGE, erep});
    energies->insert({InteractionType::DISPERSION, edisp});

    return erep + edisp;
}

static double computeBornMayer(const TopologyEntryVector             &pairs,
                               gmx_unused const std::vector<ActAtom> &atoms,
                               const std::vector<gmx::RVec>          *coordinates,
                               std::vector<gmx::RVec>                *forces,
                               std::map<InteractionType, double>     *energies)
{
    double eexp  = 0;
    auto   x     = *coordinates;
    for (const auto &b : pairs)
    {
        // Get the parameters. We have to know their names to do this.
        auto &params    = b->params();
        // Get the atom indices
        auto &indices   = b->atomIndices();
        auto ai         = indices[0];
        auto aj         = indices[1];
        rvec dx;
        rvec_sub(x[ai], x[aj], dx);
        auto dr2        = iprod(dx, dx);
        auto rinv       = gmx::invsqrt(dr2);
        // Charge transfer correction, according to Eqn. 20, Walker et al.
        // https://doi.org/10.1002/jcc.26954
        real aexp       = params[expA];
        if (aexp > 0)
        {
            real bexp   = params[expB];
            auto eeexp  = -aexp*std::exp(-bexp*dr2*rinv);
            if (debug)
            {
                fprintf(debug, "vdwcorr r  %g  eeexp %g aexp %g bexp %g\n", dr2*rinv, eeexp, aexp, bexp);
            }
            real fexp  = bexp*eeexp;
            eexp      += eeexp;

            real fbond  = fexp*rinv;
            pairforces(fbond, dx, indices, forces);
        }
    }
    energies->insert({InteractionType::VDWCORRECTION, eexp});

    return eexp;
}

static double computeDoubleExponential(const TopologyEntryVector             &pairs,
                                       gmx_unused const std::vector<ActAtom> &atoms,
                                       const std::vector<gmx::RVec>          *coordinates,
                                       std::vector<gmx::RVec>                *forces,
                                       std::map<InteractionType, double>     *energies)
{
    double eexp  = 0;
    auto   x     = *coordinates;
    for (const auto &b : pairs)
    {
        // Get the parameters. We have to know their names to do this.
        auto &params    = b->params();
        real aexp       = params[dexpA1] - params[dexpA2];
        if (aexp == 0)
        {
            continue;
        }
        real bexp       = params[dexpB];
        // Get the atom indices
        auto &indices   = b->atomIndices();
        auto ai         = indices[0];
        auto aj         = indices[1];
        rvec dx;
        rvec_sub(x[ai], x[aj], dx);
        auto dr2        = iprod(dx, dx);
        auto rinv       = gmx::invsqrt(dr2);
        auto eeexp      = -aexp*std::exp(-bexp*dr2*rinv);
        if (debug)
        {
            fprintf(debug, "r  %g  dexp %g aexp %g bexp %g\n", dr2*rinv, eeexp, aexp, bexp);
        }
        real fexp  = bexp*eeexp;
        eexp      += eeexp;
        
        real fbond  = fexp*rinv;
        pairforces(fbond, dx, indices, forces);
    }
    energies->insert({InteractionType::INDUCTIONCORRECTION, eexp});

    return eexp;
}

static double computeWBH(const TopologyEntryVector             &pairs,
                         gmx_unused const std::vector<ActAtom> &atoms,
                         const std::vector<gmx::RVec>          *coordinates,
                         std::vector<gmx::RVec>                *forces,
                         std::map<InteractionType, double>     *energies)
{
    double erep  = 0;
    double edisp = 0;
    auto   x     = *coordinates;
    for (const auto &b : pairs)
    {
        // Get the parameters. We have to know their names to do this.
        auto &params    = b->params();
        auto sigma      = params[wbhSIGMA];
        auto epsilon    = params[wbhEPSILON];
        auto gamma      = params[wbhGAMMA];
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
            if (debug)
            {
                fprintf(debug, "WBHAM ai %d aj %d sigma %g epsilon %g gamma %g erep: %g edisp: %g\n",
                        indices[0], indices[1], sigma, epsilon, gamma, eerep, eedisp);
            }
            real fbond  = fwbh*rinv;
            pairforces(fbond, dx, indices, forces);
        }
    }
    energies->insert({InteractionType::EXCHANGE, erep});
    energies->insert({InteractionType::DISPERSION, edisp});

    return erep + edisp;
}

static double computeBuckingham(const TopologyEntryVector             &pairs,
                                gmx_unused const std::vector<ActAtom> &atoms,
                                const std::vector<gmx::RVec>          *coordinates,
                                std::vector<gmx::RVec>                *forces,
                                std::map<InteractionType, double>     *energies)
{
    double erep  = 0;
    double edisp = 0;
    auto   x     = *coordinates;
    for (const auto &b : pairs)
    {
        // Get the parameters. We have to know their names to do this.
        auto &params    = b->params();
        auto Abh  = params[bhA];
        auto bbh  = params[bhB];
        auto c6bh = params[bhC6];
        if (Abh > 0 && bbh > 0 && c6bh > 0)
        {
            // Get the atom indices
            auto &indices   = b->atomIndices();
            rvec dx;
            rvec_sub(x[indices[0]], x[indices[1]], dx);
            auto dr2    = iprod(dx, dx);
            auto rinv   = gmx::invsqrt(dr2);
            auto rinv2  = rinv*rinv;
            real eerep  = Abh*std::exp(-bbh*dr2*rinv);
            real eedisp = -c6bh*rinv2*rinv2*rinv2;
            erep     += eerep;
            edisp    += eedisp;
            if (debug)
            {
                fprintf(debug, "BHAM ai %d aj %d A %g b %g c6 %g erep: %g edisp: %g\n",
                        indices[0], indices[1], Abh, bbh, c6bh, eerep, eedisp);
            }
            real fbond  = (bbh*eerep + 6*edisp*rinv)*rinv;
            pairforces(fbond, dx, indices, forces);
        }
    }
    energies->insert({InteractionType::EXCHANGE, erep});
    energies->insert({InteractionType::DISPERSION, edisp});

    return erep + edisp;
}

static void lowTT(const std::vector<int>       &indices,
                  const std::vector<gmx::RVec> &x,
                  double                        Att,
                  double                        bDisp,
                  double                        bExch,
                  const std::vector<double>    &ctt,
                  std::vector<gmx::RVec>       *forces,
                  double                       *erep,
                  double                       *edisp)
{
    static const int fac[11] = { 1, 1, 2, 6, 24, 120, 720, 5040, 40320, 362880, 3628800 };
    std::vector<double> br(11);
    rvec dx;
    rvec_sub(x[indices[0]], x[indices[1]], dx);
    auto dr2     = iprod(dx, dx);
    auto rinv1   = gmx::invsqrt(dr2);
    auto rinv2   = rinv1*rinv1;
    real eerep   = Att*std::exp(-bExch*dr2*rinv1);
    real frep    = bExch*eerep;
    real bDispr  = bDisp*dr2*rinv1;
    real ebrDisp = std::exp(-bDispr);
    real eedisp  = 0;
    real fdisp   = 0;
    real rinvn   = rinv2*rinv2*rinv2;
    br[0] = 1;
    for(int k = 1; k < 11; k++)
    {
        br[k] = br[k-1]*bDispr;
    }
    for (int m = 3; m <= 5; m++)
    {
        real fk  = 0;
        real dfk = 0;
        int  nn  = 2*m;
        for (int k = 0; k <= nn; k++)
        {
            fk  += br[k]/fac[k];
            if (k > 0)
            {
                dfk += k*bDisp*br[k-1]/fac[k];
            }
        }
        real ed = -(1-ebrDisp*fk)*ctt[m-3]*rinvn;
        eedisp += ed;
        fdisp  += ed*nn*rinv1;
        fdisp  += ctt[m-3]*rinvn*(bDisp*ebrDisp*fk - ebrDisp*dfk);
        rinvn  *= rinv2;
    }
    *erep     += eerep;
    *edisp    += eedisp;
    if (debug)
    {
        fprintf(debug, "lowTT ai %d aj %d dr %g A %g bExch %g bDisp %g c6 %g c8 %g c10 %g erep: %g edisp: %g frep: %g fdisp: %g\n",
                indices[0], indices[1], 1/rinv1, Att, bExch, bDisp, ctt[0], ctt[1], ctt[2], eerep, eedisp, frep, fdisp);
    }
    real fbond = (frep+fdisp)*rinv1;
    pairforces(fbond, dx, indices, forces);
}

static double computeTangToennies(const TopologyEntryVector             &pairs,
                                  gmx_unused const std::vector<ActAtom> &atoms,
                                  const std::vector<gmx::RVec>          *coordinates,
                                  std::vector<gmx::RVec>                *forces,
                                  std::map<InteractionType, double>     *energies)
{
    double erep  = 0;
    double edisp = 0;
    for (const auto &b : pairs)
    {
        // Get the parameters. We have to know their names to do this.
        auto &params    = b->params();
        // Call low level routine
        lowTT(b->atomIndices(), *coordinates, params[ttA], params[ttB], params[ttB],
              { params[ttC6], params[ttC8], params[ttC10] },
              forces, &erep, &edisp);
    }
    energies->insert({InteractionType::EXCHANGE, erep});
    energies->insert({InteractionType::DISPERSION, edisp});

    return erep + edisp;
}

static double computeTT2b(const TopologyEntryVector             &pairs,
                          gmx_unused const std::vector<ActAtom> &atoms,
                          const std::vector<gmx::RVec>          *coordinates,
                          std::vector<gmx::RVec>                *forces,
                          std::map<InteractionType, double>     *energies)
{
    double erep  = 0;
    double edisp = 0;
    for (const auto &b : pairs)
    {
        // Get the parameters. We have to know their names to do this.
        auto &params    = b->params();
        // Call low level routine
        lowTT(b->atomIndices(), *coordinates, params[tt2bA],
              params[tt2bBdisp], params[tt2bBexch],
              { params[tt2bC6], params[tt2bC8], params[tt2bC10] },
              forces, &erep, &edisp);
    }
    energies->insert({InteractionType::EXCHANGE, erep});
    energies->insert({InteractionType::DISPERSION, edisp});

    return erep + edisp;
}

static double computeNonBonded(const TopologyEntryVector             &pairs,
                               gmx_unused const std::vector<ActAtom> &atoms,
                               const std::vector<gmx::RVec>          *coordinates,
                               std::vector<gmx::RVec>                *forces,
                               std::map<InteractionType, double>     *energies)
{
    double erep  = 0;
    double edisp = 0;
    auto   x     = *coordinates;
    for (const auto &b : pairs)
    {
        // Get the parameters. We have to know their names to do this.
        auto &params    = b->params();
        auto rmin       = params[gbhRMIN];
        auto epsilon    = params[gbhEPSILON];
        auto gamma      = params[gbhGAMMA];
        auto delta      = params[gbhDELTA];
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
                                        ((6 + delta)*std::exp(gamma - gamma*rstar))/(1 + std::pow(rstar,6)) - (2*delta*std::pow(rstar,-1 + delta))/std::pow(1 + std::pow(rstar,delta),2)))/(2.*rmin);
            if (debug)
             {
                 fprintf(debug, "vrep: ai %d aj %d %g vdisp: %g epsilon: %10g gamma: %10g sigma: %10g delta: %10g\n",
                         indices[0], indices[1], eerep, eedisp, epsilon, gamma, rmin, delta);
             }

            erep     += eerep;
            edisp    += eedisp;
            real fbond  = fgbham*rinv;
            pairforces(fbond, dx, indices, forces);
        }
    }
    energies->insert({InteractionType::EXCHANGE, erep});
    energies->insert({InteractionType::DISPERSION, edisp});

    return erep + edisp;
}

static double gmx_unused computeNonBondedTest(const TopologyEntryVector             &pairs,
                                              gmx_unused const std::vector<ActAtom> &atoms,
                                              const std::vector<gmx::RVec>          *coordinates,
                                              std::vector<gmx::RVec>                *forces,
                                              std::map<InteractionType, double>     *energies)
{
    double erep  = 0;
    double edisp = 0;
    auto   x     = *coordinates;
    for (const auto &b : pairs)
    {
        // Get the parameters. We have to know their names to do this.
        auto &params    = b->params();
        auto rmin       = params[gbhRMIN];
        auto epsilon    = params[gbhEPSILON];
        auto gamma      = params[gbhGAMMA];
        auto delta      = params[gbhDELTA];
        if (epsilon > 0 && gamma > 0 && rmin > 0 && delta > 0)
        {
            // Get the atom indices
            auto &indices   = b->atomIndices();
            rvec dx;
            rvec_sub(x[indices[0]], x[indices[1]], dx);
            auto dr2        = iprod(dx, dx);
            auto rinv       = gmx::invsqrt(dr2);
            real rstar      = dr2*rinv/rmin;
            real sixterm    = 1 + std::pow(rstar,6);

            real r6term     = 1.0/sixterm;
            real delterm    = 1 + std::pow(rstar,delta);
            real delta6gam2 = 6 + delta + 2*gamma;
            real expterm    = std::exp(gamma*(1 - rstar));
            
            real eerep      = delta6gam2*r6term*(6+delta)*expterm/(2*gamma*delta6gam2);
            real eedisp     = -delta6gam2*r6term/(2*gamma) - 1/delterm;
            //real sixdenom   = 1/(2*gamma*sixterm);
            //real eerep      = epsilon*delta6*expterm*sixdenom;
            //real eedisp     = -epsilon*(delta6gam2*sixdenom + 1/delterm);
            real fgbham     = (epsilon*((-6*(6 + delta - (6 + delta)*std::exp(gamma - gamma*rstar) + 2*gamma)*std::pow(rstar,5))/(gamma*std::pow(1 + std::pow(rstar,6),2)) + 
                                        ((6 + delta)*std::exp(gamma - gamma*rstar))/(1 + std::pow(rstar,6)) - (2*delta*std::pow(rstar,-1 + delta))/std::pow(1 + std::pow(rstar,delta),2)))/2.;

            erep     += eerep;
            edisp    += eedisp;
            real fbond  = fgbham*rinv;
            pairforces(fbond, dx, indices, forces);
        }
    }
    energies->insert({InteractionType::EXCHANGE, erep});
    energies->insert({InteractionType::DISPERSION, edisp});

    return erep + edisp;
}

static double computeCoulombGaussian(const TopologyEntryVector         &pairs,
                                     const std::vector<ActAtom>        &atoms,
                                     const std::vector<gmx::RVec>      *coordinates,
                                     std::vector<gmx::RVec>            *forces,
                                     std::map<InteractionType, double> *energies)
{
    double ebond = 0;
    auto   x     = *coordinates;
    double vc_pt = 0;
    for (const auto &b : pairs)
    {
        // Get the parameters.
        auto &indices = b->atomIndices();
        auto ai       = indices[0];
        auto aj       = indices[1];
        auto &params= b->params();
        auto izeta  = params[coulZETA];
        auto jzeta  = params[coulZETA2];
        // Get the atom indices
        real qq         = ONE_4PI_EPS0*atoms[ai].charge()*atoms[aj].charge();
        rvec dx;
        rvec_sub(x[ai], x[aj], dx);
        auto dr2        = iprod(dx, dx);
        auto r          = std::sqrt(dr2);
        real velec =  qq*Coulomb_GG(r, izeta, jzeta);
        real felec = -qq*DCoulomb_GG(r, izeta, jzeta);

        if (debug)
        {
            auto r1 = std::sqrt(dr2);
            fprintf(debug, "vcoul ai %d aj %d %g izeta %g jzeta %g qi %g qj %g vcoul_pc %g dist %g felec %g\n",
                    ai, aj, velec, izeta, jzeta, 
                    atoms[ai].charge(), atoms[aj].charge(), qq/r1, r1, felec/r1);
            vc_pt += qq/r1;
        }
        ebond += velec;
        if (dr2 > 0)
        {
            felec *= gmx::invsqrt(dr2);
            pairforces(felec, dx, indices, forces);
        }
    }
    if (debug)
    {
        fprintf(debug, "vcoul %g vcoul_pc %g\n", ebond, vc_pt);
    }
    energies->insert({InteractionType::ELECTROSTATICS, ebond});

    return ebond;
}

static double computeCoulombSlater(const TopologyEntryVector         &pairs,
                                   const std::vector<ActAtom>        &atoms,
                                   const std::vector<gmx::RVec>      *coordinates,
                                   std::vector<gmx::RVec>            *forces,
                                   std::map<InteractionType, double> *energies)
{
    double ebond = 0;
    auto   x     = *coordinates;
    for (const auto &b : pairs)
    {
        // Get the parameters.
        auto &indices = b->atomIndices();
        auto ai       = indices[0];
        auto aj       = indices[1];
        // Get the parameters. We have to know their names to do this.
        auto &params= b->params();
        auto izeta  = params[coulZETA];
        auto jzeta  = params[coulZETA2];
        auto irow   = atoms[ai].row();
        auto jrow   = atoms[aj].row();
        real qq     = ONE_4PI_EPS0*atoms[ai].charge()*atoms[aj].charge();
        rvec dx;
        rvec_sub(x[ai], x[aj], dx);
        auto dr2   = iprod(dx, dx);
        real r1    = std::sqrt(dr2);
        real velec =  qq*Coulomb_SS(r1, irow, jrow, izeta, jzeta);
        // The DCoulomb_SS code returns the derivative of the energy
        // so we still need a minus sign.
        // https://github.com/dspoel/ACT/issues/369
        real felec = -qq*DCoulomb_SS(r1, irow, jrow, izeta, jzeta);
        if (debug)
        {
            auto r1 = std::sqrt(dr2);
            fprintf(debug, "vcoul ai %d aj %d %g fcoul %g izeta %g jzeta %g qi %g qj %g vcoul_pc %g fcoul_pc %g dist %g\n",
                    ai, aj, velec, felec, izeta, jzeta, atoms[ai].charge(),
                    atoms[aj].charge(), qq/r1, qq/dr2, r1);
        }
        ebond += velec;
        if (dr2 > 0)
        {
            felec *= gmx::invsqrt(dr2);
            pairforces(felec, dx, indices, forces);
        }
    }
    energies->insert({InteractionType::ELECTROSTATICS, ebond});

    return ebond;
}

#ifdef FUTURE
static double computePartridge(const TopologyEntryVector             &angles,
                               gmx_unused const std::vector<ActAtom> &atoms,
                               const std::vector<gmx::RVec>          *coordinates,
                               std::vector<gmx::RVec>                *forces,
                               std::map<InteractionType, double>     *energies)
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

    return energy;
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

static double computeDummy(gmx_unused const TopologyEntryVector         &bonds,
                           gmx_unused const std::vector<ActAtom>        &atoms,
                           gmx_unused const std::vector<gmx::RVec>      *coordinates,
                           gmx_unused std::vector<gmx::RVec>            *forces,
                           gmx_unused std::map<InteractionType, double> *energies)
{
    return 0;
}

static double computeBonds(const TopologyEntryVector             &bonds,
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
        
        pairforces(fbond, dx, indices, forces);
    }
    energies->insert({InteractionType::BONDS, ebond});

    return ebond;
}

static double computeCubic(const TopologyEntryVector             &bonds,
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
        double rhi     = (2*rmax+bondlength)/3; 
        double r_1, r1 = rhi;
        {
            // Do not allow the potential to go to -infinity
            auto r2 = std::min(rhi*rhi, iprod(dx, dx));
            r_1     = gmx::invsqrt(r2);
            r1      = r2*r_1;
        }
        real vB, fbond;
        vB    = kb*gmx::square(r1-bondlength)*(rmax-r1) - De;
        fbond = r_1*kb*(bondlength+2*rmax-3*r1)*(bondlength-r1);
        
        ebond += vB;
        
        pairforces(fbond, dx, indices, forces);
    }
    energies->insert({InteractionType::BONDS, ebond});

    return ebond;
}

static double computeMorse(const TopologyEntryVector             &bonds,
                           gmx_unused const std::vector<ActAtom> &atoms,
                           const std::vector<gmx::RVec>          *coordinates,
                           std::vector<gmx::RVec>                *forces,
                           std::map<InteractionType, double>     *energies)
{
    double  ebond = 0;
    auto    x     = *coordinates;
    for (const auto &b : bonds)
    {
        // Get the atom indices
        auto &indices   = b->atomIndices();
        // Get the parameters. We have to know their names to do this.
        auto &params    = b->params();
        auto bondlength = params[morseLENGTH];
        auto beta       = params[morseBETA];
        auto De         = params[morseDE];
        auto D0         = params[morseD0];
        rvec dx;
        rvec_sub(x[indices[0]], x[indices[1]], dx);
        auto dr2        = iprod(dx, dx);
        auto dr         = std::sqrt(dr2);
        auto expterm    = std::exp(-beta*(dr-bondlength));
        auto vbond      = De*expterm*(expterm - 2) + D0;
        auto fbond      = 2*De*beta*expterm*(expterm-1)/dr;
        ebond          += vbond;
        if (debug)
        {
            fprintf(debug, "Morse i %d j %d dr %g vbond %g fbond %g\n",
                    indices[0], indices[1], dr, vbond, fbond);
        }
        pairforces(fbond, dx, indices, forces);
    }
    energies->insert({InteractionType::BONDS, ebond});

    return ebond;
}

static double computeHua(const TopologyEntryVector             &bonds,
                         gmx_unused const std::vector<ActAtom> &atoms,
                         const std::vector<gmx::RVec>          *coordinates,
                         std::vector<gmx::RVec>                *forces,
                         std::map<InteractionType, double>     *energies)
{
    double  ebond = 0;
    auto    x     = *coordinates;
    for (const auto &b : bonds)
    {
        // Get the parameters. We have to know their names to do this.
        auto &params = b->params();
        auto re      = params[huaLENGTH];
        auto hb      = params[huaB];
        auto De      = params[huaDE];
        auto c       = params[huaC];
        // Get the atom indices
        auto &indices   = b->atomIndices();
        rvec dx;
        rvec_sub(x[indices[0]], x[indices[1]], dx);
        auto dr2        = iprod(dx, dx);
        auto dr         = std::sqrt(dr2);
        auto expterm    = std::exp(-hb*(dr-re));
        auto exp2       = (1-expterm)/(1-c*expterm);
        auto vbond      = De*(exp2*exp2-1);
        auto denom      = (c*expterm - 1);
        auto fbond      = 2*De*hb*expterm*(expterm-1)*(c-1)/(dr*denom*denom*denom);
        ebond          += vbond;

        pairforces(fbond, dx, indices, forces);
    }
    energies->insert({InteractionType::BONDS, ebond});

    return ebond;
}

static double computeLinearAngles(const TopologyEntryVector             &angles,
                                  gmx_unused const std::vector<ActAtom> &atoms,
                                  const std::vector<gmx::RVec>          *coordinates,
                                  std::vector<gmx::RVec>                *forces,
                                  std::map<InteractionType, double>     *energies)
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

    return ebond;
}

static double computeAngles(const TopologyEntryVector             &angles,
                            gmx_unused const std::vector<ActAtom> &atoms,
                            const std::vector<gmx::RVec>          *coordinates,
                            std::vector<gmx::RVec>                *forces,
                            std::map<InteractionType, double>     *energies)
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
        if (debug)
        {
            fprintf(debug, "Angle %g refangle %g, energy %g\n", theta, theta0, vA);
        }
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

    return energy;
}

static double computeUreyBradley(const TopologyEntryVector             &angles,
                                 gmx_unused const std::vector<ActAtom> &atoms,
                                 const std::vector<gmx::RVec>          *coordinates,
                                 std::vector<gmx::RVec>                *forces,
                                 std::map<InteractionType, double>     *energies)
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

    return energy;
}

static double computePolarization(const TopologyEntryVector             &bonds,
                                  gmx_unused const std::vector<ActAtom> &atoms,
                                  const std::vector<gmx::RVec>          *coordinates,
                                  std::vector<gmx::RVec>                *forces,
                                  std::map<InteractionType, double>     *energies)
{
    const real half = 0.5;
    double ebond    = 0;
    auto   x        = *coordinates;
    for (const auto &b : bonds)
    {
        // Get the parameters. We have to know their names to do this.
        auto &params   = b->params();
        // Get polarizability
        auto alpha     = params[polALPHA];
        if (alpha == 0)
        {
            continue;
        }
        auto rhyper    = params[polRHYPER];
        auto fchyper   = params[polFCHYPER];
        // Get the atom indices
        auto &indices  = b->atomIndices();
        // Get shell charge
        auto q         = atoms[indices[1]].charge();
        double ksh     = ONE_4PI_EPS0*q*q/alpha;
        rvec dx;
        rvec_sub(x[indices[0]], x[indices[1]], dx);
        auto dr2  = iprod(dx, dx);
        
        if (dr2 == 0.0)
        {
            continue;
        }

        auto dr         = dr2 * gmx::invsqrt(dr2);
        auto fbond      = -ksh;
        ebond          += half*ksh*dr2;
        if (debug)
        {
            fprintf(debug, "vpol %g ai %d aj %d alpha %g ksh %g dist %g\n", half*ksh*dr2, indices[0],
                    indices[1], alpha, ksh, dr);
        }
        if (dr > rhyper && fchyper > 0)
        {
            auto ddr  = dr - rhyper;
            auto ddr2 = ddr * ddr;
            ebond += fchyper * ddr2 * ddr2;
            // We will multiply by the distance vector components below
            fbond -= 4 * fchyper * ddr2 * ddr/dr;
            if (debug)
            {
                fprintf(debug, "Adding hyperpolarization energy correction %g kJ/mol core pos %g %g %g ddr = %g.\n",
                        fchyper * ddr2 * ddr2, x[indices[0]][XX],
                        x[indices[0]][YY],x[indices[0]][ZZ], ddr);
            }
        }
        pairforces(fbond, dx, indices, forces);
    }
    energies->insert({InteractionType::POLARIZATION, ebond});

    return ebond;
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

static double computeFourDihs(const TopologyEntryVector             &propers,
                              gmx_unused const std::vector<ActAtom> &atoms,
                              const std::vector<gmx::RVec>          *coordinates,
                              std::vector<gmx::RVec>                *forces,
                              std::map<InteractionType, double>     *energies)
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
        auto phi = M_PI + dih_angle(x[ai], x[aj], x[ak], x[al],
                                    r_ij, r_kj, r_kl, m, n);

        real cos_phi = std::cos(phi);
        
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

    return energy;
}

static double computeImpropers(const TopologyEntryVector             &impropers,
                               gmx_unused const std::vector<ActAtom> &atoms,
                               const std::vector<gmx::RVec>          *coordinates,
                               std::vector<gmx::RVec>                *forces,
                               std::map<InteractionType, double>     *energies)
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

    return energy;
}

static double computePropers(const TopologyEntryVector             &propers,
                             gmx_unused const std::vector<ActAtom> &atoms,
                             const std::vector<gmx::RVec>          *coordinates,
                             std::vector<gmx::RVec>                *forces,
                             std::map<InteractionType, double>     *energies)
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
        if (debug)
        {
            fprintf(debug, "ai %d aj %d ak %d al %d phi %g phi0 %g kphi %g, mult %g, mdphi %g, v %g ddphi %g\n",
                    ai, aj, ak, al, phi, phi0, kphi, mult, mdphi, v, ddphi);
        }
        do_dih_fup_noshiftf(ai, aj, ak, al, ddphi, r_ij, r_kj, r_kl, m, n,
                            forces);
    }
    energies->insert({InteractionType::PROPER_DIHEDRALS, energy});

    return energy;
}

std::map<Potential, bondForceComputer> bondForceComputerMap = {
    { Potential::NONE,                   computeDummy           },
    { Potential::HARMONIC_BONDS,         computeBonds           },
    { Potential::MORSE_BONDS,            computeMorse           },
    { Potential::HUA_BONDS,              computeHua             },
    { Potential::CUBIC_BONDS,            computeCubic           },
    { Potential::HARMONIC_ANGLES,        computeAngles          },
    { Potential::LINEAR_ANGLES,          computeLinearAngles    },
    { Potential::LJ12_6,                 computeLJ12_6          },
    { Potential::LJ8_6,                  computeLJ8_6           },
    { Potential::LJ14_7,                 computeLJ14_7          },
    { Potential::BUCKINGHAM,             computeBuckingham      },
    { Potential::WANG_BUCKINGHAM,        computeWBH             },
    { Potential::TANG_TOENNIES,          computeTangToennies    },
    { Potential::TT2b,                   computeTT2b            },
    { Potential::GENERALIZED_BUCKINGHAM, computeNonBonded       },
    { Potential::BORN_MAYER,             computeBornMayer       },
    { Potential::MACDANIEL_SCHMIDT,      computeDoubleExponential },
    { Potential::COULOMB_POINT,          computeCoulombGaussian },
    { Potential::COULOMB_GAUSSIAN,       computeCoulombGaussian },
    { Potential::COULOMB_SLATER,         computeCoulombSlater   },
    { Potential::POLARIZATION,           computePolarization    },
    { Potential::HARMONIC_DIHEDRALS,     computeImpropers       },
    { Potential::PROPER_DIHEDRALS,       computePropers         },
    { Potential::FOURIER_DIHEDRALS,      computeFourDihs        },
    { Potential::UREY_BRADLEY_ANGLES,    computeUreyBradley     }
};

bondForceComputer getBondForceComputer(Potential pot)
{
    auto bfc = bondForceComputerMap.find(pot);
    if (bondForceComputerMap.end() != bfc)
    {
        return *bfc->second;
    }
    // Keep the compiler happy
    return nullptr;
}

} // namespace

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
 * Pair potentials for many conformations at once with explicit SIMD.
 *
 * The SIMD type and math functions come from simd_math.h, such that this
 * code is the same for all processors.
 * The SIMD lanes hold different conformations. The parameters of each
 * interaction are prepared once, before looping over the conformations.
 * The loop over interactions is the outer loop. The other order, with
 * the energies kept in registers, was tried as well, but it is much
 * slower for many conformations, since the coordinates of one block of
 * conformations are then spread over few cache sets.
 *
 * Each potential is a struct with two functions, prepare and compute,
 * see LJ12_6 for the documentation. The formulas are those of the
 * corresponding xxxOne functions in forces_simd.cpp, which compute the
 * same potentials for one conformation and are the reference.
 * \author David van der Spoel <david.vanderspoel@icm.uu.se>
 */
#include "actpre.h"

#include "forces_explicit_simd.h"

#include "act/forces/simd_math.h"

#ifdef ACT_HAVE_SIMD

#include "act/coulombintegrals/gaussian_integrals.h"
#include "act/forcefield/forcefield_parametername.h"
#include "gromacs/math/functions.h"
#include "gromacs/math/units.h"
#include "gromacs/math/vectypes.h"

namespace alexandria
{

namespace
{

using simd::SimdDouble;

//! Maximum number of prepared parameters per interaction
constexpr int maxPairParams = 8;

//! A pair interaction with prepared parameters
struct PairParams
{
    //! Offset of the coordinates of the first atom
    size_t i;
    //! Offset of the coordinates of the second atom
    size_t j;
    //! The parameters
    double p[maxPairParams];
};

//! Lennard-Jones 12-6
struct LJ12_6
{
    /*! \brief Prepare the parameters of one interaction
     * \param[in]  b        The interaction
     * \param[in]  atoms    The atoms
     * \param[in]  epsilonr Relative dielectric constant
     * \param[out] p        The parameters
     * \return false if there is no interaction
     */
    static bool prepare(const TopologyEntry &b, const std::vector<ActAtom> &, double, double *p)
    {
        auto sig6 = gmx::power6(b.params()[lj12_6SIGMA]);
        p[0]      = 4*b.params()[lj12_6EPSILON]*sig6;
        p[1]      = p[0]*sig6;
        return true;
    }
    /*! \brief Compute the interaction
     * \param[in]    p     The parameters
     * \param[in]    rsq   The distance squared
     * \param[inout] e     The energy terms to add to
     * \param[out]   fscal The scalar force divided by the distance
     */
    static void compute(const double *p, SimdDouble rsq, SimdDouble *e, SimdDouble *fscal)
    {
        SimdDouble rinv2 = 1.0/rsq;
        SimdDouble rinv6 = rinv2*rinv2*rinv2;
        SimdDouble disp  = -p[0]*rinv6;
        SimdDouble rep   = p[1]*rinv6*rinv6;
        e[0]  += rep;
        e[1]  += disp;
        *fscal = (12.0*rep + 6.0*disp)*rinv2;
    }
};

//! Lennard-Jones 12-6-4
struct LJ12_6_4
{
    //! \copydoc LJ12_6::prepare
    static bool prepare(const TopologyEntry &b, const std::vector<ActAtom> &, double, double *p)
    {
        auto sig6 = gmx::power6(b.params()[lj12_6_4SIGMA]);
        p[0]      = 4*b.params()[lj12_6_4EPSILON]*sig6;
        p[1]      = p[0]*sig6;
        p[2]      = b.params()[lj12_6_4GAMMA];
        return true;
    }
    //! \copydoc LJ12_6::compute
    static void compute(const double *p, SimdDouble rsq, SimdDouble *e, SimdDouble *fscal)
    {
        SimdDouble rinv2 = 1.0/rsq;
        SimdDouble rinv4 = rinv2*rinv2;
        SimdDouble rinv6 = rinv4*rinv2;
        SimdDouble mid   = -p[2]*rinv4;
        SimdDouble disp  = -p[0]*rinv6;
        SimdDouble rep   = p[1]*rinv6*rinv6;
        e[0]  += rep;
        e[1]  += disp;
        e[2]  += mid;
        *fscal = (12.0*rep + 6.0*disp + 4.0*mid)*rinv2;
    }
};

//! Lennard-Jones 8-6
struct LJ8_6
{
    //! \copydoc LJ12_6::prepare
    static bool prepare(const TopologyEntry &b, const std::vector<ActAtom> &, double, double *p)
    {
        auto sigma = b.params()[lj8_6SIGMA];
        p[0]       = 4*b.params()[lj8_6EPSILON]*gmx::power6(sigma);
        p[1]       = 0.75*p[0]*sigma*sigma;
        return true;
    }
    //! \copydoc LJ12_6::compute
    static void compute(const double *p, SimdDouble rsq, SimdDouble *e, SimdDouble *fscal)
    {
        SimdDouble rinv2 = 1.0/rsq;
        SimdDouble rinv6 = rinv2*rinv2*rinv2;
        SimdDouble disp  = -p[0]*rinv6;
        SimdDouble rep   = p[1]*rinv6*rinv2;
        e[0]  += rep;
        e[1]  += disp;
        *fscal = (8.0*rep + 6.0*disp)*rinv2;
    }
};

//! Buffered 14-7
struct LJ14_7
{
    //! \copydoc LJ12_6::prepare
    static bool prepare(const TopologyEntry &b, const std::vector<ActAtom> &, double, double *p)
    {
        auto &params  = b.params();
        auto  sigma   = params[lj14_7SIGMA];
        auto  epsilon = params[lj14_7EPSILON];
        p[0] = 1/sigma;
        p[1] = epsilon;
        p[2] = params[lj14_7GAMMA];
        p[3] = params[lj14_7DELTA];
        return epsilon > 0 && sigma > 0;
    }
    //! \copydoc LJ12_6::compute
    static void compute(const double *p, SimdDouble rsq, SimdDouble *e, SimdDouble *fscal)
    {
        double     gamma1    = p[2] + 1;
        SimdDouble rinv      = simd::invsqrt(rsq);
        SimdDouble rstar     = rsq*rinv*p[0];
        SimdDouble rstar2    = rstar*rstar;
        SimdDouble rstar6    = rstar2*rstar2*rstar2;
        SimdDouble deltars   = p[3] + rstar;
        SimdDouble invdrs    = 1.0/deltars;
        SimdDouble t         = (p[3] + 1)*invdrs;
        SimdDouble t2        = t*t;
        SimdDouble repfac    = p[1]*t2*t2*t2*t;
        SimdDouble invgr7    = 1.0/(p[2] + rstar6*rstar);
        e[0]  += repfac*gamma1*invgr7;
        e[1]  += -2.0*repfac;
        // Since epsilon*delta1^7/deltars^8 = repfac/deltars
        SimdDouble f147      = 7.0*repfac*(gamma1*rstar6*deltars*invgr7*invgr7 + gamma1*invgr7 - 2.0)*p[0]*invdrs;
        *fscal = f147*rinv;
    }
};

//! Born-Mayer, used also in Tang-Toennies
static SimdDouble bornMayer(double aexp, double bexp, SimdDouble r, SimdDouble rinv, SimdDouble *fscal)
{
    SimdDouble eexp = aexp*simd::exp(-bexp*r);
    *fscal += bexp*eexp*rinv;
    return eexp;
}

//! Slater-ISA, used also in Slater-ISA-Tang-Toennies
static SimdDouble slaterISA(double aexp, double bexp, SimdDouble r, SimdDouble rinv, SimdDouble *fscal)
{
    const double third = 1.0/3.0;
    SimdDouble   br    = bexp*r;
    SimdDouble   aterm = aexp*simd::exp(-br);
    SimdDouble   eexp  = aterm*(br*br*third + br + 1.0);
    *fscal += (bexp*eexp - bexp*aterm*(2.0*third*br + 1.0))*rinv;
    return eexp;
}

/*! \brief Tang-Toennies damped dispersion
 * \param[in]    bDisp Exponent of the damping
 * \param[in]    ctt   C6, C8 and C10
 * \param[in]    r     Distance
 * \param[in]    rinv  Inverse distance
 * \param[inout] fscal Scalar force divided by the distance to add to
 * \return The dispersion energy
 */
static SimdDouble ttDispersion(double bDisp, const double *ctt, SimdDouble r, SimdDouble rinv, SimdDouble *fscal)
{
    SimdDouble rinv2  = rinv*rinv;
    SimdDouble bDispr = bDisp*r;
    SimdDouble ebr    = simd::exp(-bDispr);
    SimdDouble rinvn  = rinv2*rinv2*rinv2;
    SimdDouble edisp  = 0.0;
    SimdDouble fdisp  = 0.0;
    // The damping function for C_2m is f_2m = 1 - exp(-br) sum_k=0^2m (br)^k/k!
    // The sum is accumulated in fk, with term holding (br)^k/k!
    SimdDouble fk     = 1.0;
    SimdDouble term   = 1.0;
    int        k      = 0;
    for (int m = 3; m <= 5; m++)
    {
        int nn = 2*m;
        for (; k < nn; k++)
        {
            term *= bDispr*(1.0/(k+1));
            fk   += term;
        }
        SimdDouble ed = -(1.0 - ebr*fk)*ctt[m-3]*rinvn;
        edisp += ed;
        fdisp += ed*static_cast<double>(nn)*rinv;
        // The derivative of the sum is b times the sum up to 2m-1,
        // the terms up to 2m-1 cancel against the derivative of exp(-br).
        fdisp += ctt[m-3]*rinvn*bDisp*ebr*term;
        rinvn *= rinv2;
    }
    *fscal += fdisp*rinv;
    return edisp;
}

//! Born-Mayer
struct BornMayer
{
    //! \copydoc LJ12_6::prepare
    static bool prepare(const TopologyEntry &b, const std::vector<ActAtom> &, double, double *p)
    {
        // Charge transfer correction, according to Eqn. 20, Walker et al.
        // https://doi.org/10.1002/jcc.26954
        p[0] = -b.params()[expA];
        p[1] = b.params()[expB];
        return true;
    }
    //! \copydoc LJ12_6::compute
    static void compute(const double *p, SimdDouble rsq, SimdDouble *e, SimdDouble *fscal)
    {
        SimdDouble rinv = simd::invsqrt(rsq);
        *fscal = 0.0;
        e[0]  += bornMayer(p[0], p[1], rsq*rinv, rinv, fscal);
    }
};

//! Slater-ISA
struct SlaterISA
{
    //! \copydoc LJ12_6::prepare
    static bool prepare(const TopologyEntry &b, const std::vector<ActAtom> &, double, double *p)
    {
        p[0] = b.params()[expA];
        p[1] = b.params()[expB];
        return true;
    }
    //! \copydoc LJ12_6::compute
    static void compute(const double *p, SimdDouble rsq, SimdDouble *e, SimdDouble *fscal)
    {
        SimdDouble rinv = simd::invsqrt(rsq);
        *fscal = 0.0;
        e[0]  += slaterISA(p[0], p[1], rsq*rinv, rinv, fscal);
    }
};

//! Double exponential (MacDaniel-Schmidt)
struct DoubleExponential
{
    //! \copydoc LJ12_6::prepare
    static bool prepare(const TopologyEntry &b, const std::vector<ActAtom> &, double, double *p)
    {
        p[0] = b.params()[dexpA2] - b.params()[dexpA1];
        p[1] = b.params()[dexpB];
        return true;
    }
    //! \copydoc LJ12_6::compute
    static void compute(const double *p, SimdDouble rsq, SimdDouble *e, SimdDouble *fscal)
    {
        SimdDouble rinv = simd::invsqrt(rsq);
        *fscal = 0.0;
        e[0]  += bornMayer(p[0], p[1], rsq*rinv, rinv, fscal);
    }
};

//! Wang-Buckingham
struct WangBuckingham
{
    //! \copydoc LJ12_6::prepare
    static bool prepare(const TopologyEntry &b, const std::vector<ActAtom> &, double, double *p)
    {
        auto &params  = b.params();
        auto  sigma   = params[wbhSIGMA];
        auto  epsilon = params[wbhEPSILON];
        auto  gamma   = params[wbhGAMMA];
        auto  gamma_3 = 3.0/(gamma + 3.0);
        p[0] = gmx::power6(sigma);
        p[1] = gamma_3;
        // Prefactor of the dispersion
        p[2] = 2*epsilon/(1.0 - gamma_3);
        p[3] = gamma;
        p[4] = 1/sigma;
        return epsilon > 0 && gamma > 0 && sigma > 0;
    }
    //! \copydoc LJ12_6::compute
    static void compute(const double *p, SimdDouble rsq, SimdDouble *e, SimdDouble *fscal)
    {
        // Modified Buckingham: JCTC  Volume: 9  Page: 452  Year: 2012
        SimdDouble rinv     = simd::invsqrt(rsq);
        SimdDouble r        = rsq*rinv;
        SimdDouble r5       = rsq*rsq*r;
        SimdDouble inv6     = 1.0/(p[0] + r5*r);
        SimdDouble erep_exp = p[1]*simd::exp(p[3]*(1.0 - r*p[4]));
        SimdDouble vdisp    = -p[2]*p[0]*inv6;
        SimdDouble fdisp    = vdisp*6.0*r5*inv6;
        e[0]  += -vdisp*erep_exp;
        e[1]  += vdisp;
        *fscal = (-fdisp*erep_exp - vdisp*erep_exp*p[3]*p[4] + fdisp)*rinv;
    }
};

//! Buckingham
struct Buckingham
{
    //! \copydoc LJ12_6::prepare
    static bool prepare(const TopologyEntry &b, const std::vector<ActAtom> &, double, double *p)
    {
        p[0] = b.params()[bhA];
        p[1] = b.params()[bhB];
        p[2] = b.params()[bhC6];
        return p[0] > 0 && p[1] > 0 && p[2] > 0;
    }
    //! \copydoc LJ12_6::compute
    static void compute(const double *p, SimdDouble rsq, SimdDouble *e, SimdDouble *fscal)
    {
        SimdDouble rinv  = simd::invsqrt(rsq);
        SimdDouble rinv2 = rinv*rinv;
        SimdDouble erep  = p[0]*simd::exp(-p[1]*rsq*rinv);
        SimdDouble edisp = -p[2]*rinv2*rinv2*rinv2;
        e[0]  += erep;
        e[1]  += edisp;
        *fscal = (p[1]*erep + 6.0*edisp*rinv)*rinv;
    }
};

//! Tang-Toennies
struct TangToennies
{
    //! \copydoc LJ12_6::prepare
    static bool prepare(const TopologyEntry &b, const std::vector<ActAtom> &, double, double *p)
    {
        auto &params = b.params();
        p[0] = params[ttA];
        p[1] = params[ttB];
        p[2] = params[ttC6];
        p[3] = params[ttC8];
        p[4] = params[ttC10];
        return true;
    }
    //! \copydoc LJ12_6::compute
    static void compute(const double *p, SimdDouble rsq, SimdDouble *e, SimdDouble *fscal)
    {
        SimdDouble rinv = simd::invsqrt(rsq);
        SimdDouble r    = rsq*rinv;
        *fscal = 0.0;
        e[0]  += bornMayer(p[0], p[1], r, rinv, fscal);
        e[1]  += ttDispersion(p[1], p+2, r, rinv, fscal);
    }
};

//! Tang-Toennies with separate exponents
struct TT2b
{
    //! \copydoc LJ12_6::prepare
    static bool prepare(const TopologyEntry &b, const std::vector<ActAtom> &, double, double *p)
    {
        auto &params = b.params();
        p[0] = params[tt2bA];
        p[1] = params[tt2bBexch];
        p[2] = params[tt2bBdisp];
        p[3] = params[tt2bC6];
        p[4] = params[tt2bC8];
        p[5] = params[tt2bC10];
        return true;
    }
    //! \copydoc LJ12_6::compute
    static void compute(const double *p, SimdDouble rsq, SimdDouble *e, SimdDouble *fscal)
    {
        SimdDouble rinv = simd::invsqrt(rsq);
        SimdDouble r    = rsq*rinv;
        *fscal = 0.0;
        e[0]  += bornMayer(p[0], p[1], r, rinv, fscal);
        e[1]  += ttDispersion(p[2], p+3, r, rinv, fscal);
    }
};

//! Slater-ISA with Tang-Toennies dispersion
struct SlaterISATT
{
    //! \copydoc LJ12_6::prepare
    static bool prepare(const TopologyEntry &b, const std::vector<ActAtom> &, double, double *p)
    {
        auto &params = b.params();
        p[0] = params[SIttA];
        p[1] = params[SIttBexch];
        p[2] = params[SIttBdisp];
        p[3] = params[SIttC6];
        p[4] = params[SIttC8];
        p[5] = params[SIttC10];
        return true;
    }
    //! \copydoc LJ12_6::compute
    static void compute(const double *p, SimdDouble rsq, SimdDouble *e, SimdDouble *fscal)
    {
        SimdDouble rinv = simd::invsqrt(rsq);
        SimdDouble r    = rsq*rinv;
        *fscal = 0.0;
        e[0]  += slaterISA(p[0], p[1], r, rinv, fscal);
        e[1]  += ttDispersion(p[2], p+3, r, rinv, fscal);
    }
};

//! Coulomb between point or Gaussian charges
struct Coulomb
{
    //! \copydoc LJ12_6::prepare
    static bool prepare(const TopologyEntry &b, const std::vector<ActAtom> &atoms, double epsilonr, double *p)
    {
        auto &params   = b.params();
        auto &gindices = b.atomIndices();
        auto  izeta    = params[coulZETA];
        auto  jzeta    = params[coulZETA2];
        p[0] = ONE_4PI_EPS0*atoms[gindices[0]].charge()*atoms[gindices[1]].charge()/epsilonr;
        // Width of the combined Gaussian, as in Coulomb_GG, zero for point charges
        if (izeta == 0 || jzeta == 0)
        {
            p[1] = izeta + jzeta;
        }
        else
        {
            p[1] = izeta*jzeta/std::sqrt(izeta*izeta + jzeta*jzeta);
        }
        return true;
    }
    //! \copydoc LJ12_6::compute
    static void compute(const double *p, SimdDouble rsq, SimdDouble *e, SimdDouble *fscal)
    {
        auto       qq      = p[0];
        auto       zeta    = p[1];
        // Particles on top of each other only have an energy with Gaussian charges
        auto       overlap = !(rsq > 0.0);
        SimdDouble rinv    = simd::invsqrt(rsq);
        if (zeta == 0)
        {
            e[0]  += simd::select(overlap, 0.0, qq*rinv);
            *fscal = simd::select(overlap, 0.0, qq*rinv*rinv*rinv);
        }
        else
        {
            SimdDouble rz    = zeta*rsq*rinv;
            SimdDouble erfrz = simd::erf(rz);
            e[0]  += simd::select(overlap, qq*2*zeta*invsqrt_pi, qq*erfrz*rinv);
            *fscal = simd::select(overlap, 0.0,
                                  qq*(erfrz*rinv - 2*invsqrt_pi*zeta*simd::exp(-rz*rz))*rinv*rinv);
        }
    }
};

//! Polarization, a harmonic spring between shell and core
struct Polarization
{
    //! \copydoc LJ12_6::prepare
    static bool prepare(const TopologyEntry &b, const std::vector<ActAtom> &atoms, double, double *p)
    {
        auto &params = b.params();
        auto  alpha  = params[polALPHA];
        if (alpha == 0)
        {
            return false;
        }
        // Force constant from the shell charge
        auto q = atoms[b.atomIndices()[1]].charge();
        p[0]   = ONE_4PI_EPS0*q*q/alpha;
        p[1]   = params[polRHYPER];
        p[2]   = params[polFCHYPER];
        return true;
    }
    //! \copydoc LJ12_6::compute
    static void compute(const double *p, SimdDouble rsq, SimdDouble *e, SimdDouble *fscal)
    {
        auto ksh = p[0];
        auto fch = p[2];
        // A shell on top of its core has no energy
        auto overlap = !(rsq > 0.0);
        e[0]  += simd::select(overlap, 0.0, 0.5*ksh*rsq);
        *fscal = -ksh;
        if (fch > 0)
        {
            // Hyperpolarization beyond rhyper
            SimdDouble dr   = simd::sqrt(rsq);
            SimdDouble ddr  = simd::max(dr - p[1], 0.0);
            SimdDouble ddr2 = ddr*ddr;
            e[0]   += simd::select(overlap, 0.0, fch*ddr2*ddr2);
            *fscal -= simd::select(overlap, 0.0, 4.0*fch*ddr2*ddr/dr);
        }
    }
};

//! Harmonic bond
struct HarmonicBond
{
    //! \copydoc LJ12_6::prepare
    static bool prepare(const TopologyEntry &b, const std::vector<ActAtom> &, double, double *p)
    {
        p[0] = b.params()[bondLENGTH];
        p[1] = b.params()[bondKB];
        p[2] = b.params()[bondENERGY];
        return true;
    }
    //! \copydoc LJ12_6::compute
    static void compute(const double *p, SimdDouble rsq, SimdDouble *e, SimdDouble *fscal)
    {
        SimdDouble rinv = simd::invsqrt(rsq);
        SimdDouble dr   = rsq*rinv - p[0];
        e[0]  += 0.5*p[1]*dr*dr + p[2];
        *fscal = -p[1]*dr*rinv;
    }
};

//! Cubic bond
struct CubicBond
{
    //! \copydoc LJ12_6::prepare
    static bool prepare(const TopologyEntry &b, const std::vector<ActAtom> &, double, double *p)
    {
        auto &params = b.params();
        p[0] = params[cubicLENGTH];
        p[1] = params[cubicRMAX];
        p[2] = params[cubicKB];
        p[3] = params[cubicDE];
        // Do not allow the potential to go to -infinity
        auto rhi = (2*p[1] + p[0])/3;
        p[4] = rhi*rhi;
        return true;
    }
    //! \copydoc LJ12_6::compute
    static void compute(const double *p, SimdDouble rsq, SimdDouble *e, SimdDouble *fscal)
    {
        auto       b0   = p[0];
        auto       rmax = p[1];
        auto       kb   = p[2];
        SimdDouble r2   = simd::min(rsq, p[4]);
        SimdDouble rinv = simd::invsqrt(r2);
        SimdDouble r    = r2*rinv;
        e[0]  += kb*(r - b0)*(r - b0)*(rmax - r) - p[3];
        *fscal = rinv*kb*(b0 + 2*rmax - 3.0*r)*(b0 - r);
    }
};

//! Morse bond
struct MorseBond
{
    //! \copydoc LJ12_6::prepare
    static bool prepare(const TopologyEntry &b, const std::vector<ActAtom> &, double, double *p)
    {
        auto &params = b.params();
        p[0] = params[morseLENGTH];
        p[1] = params[morseBETA];
        p[2] = params[morseDE];
        p[3] = params[morseD0];
        return true;
    }
    //! \copydoc LJ12_6::compute
    static void compute(const double *p, SimdDouble rsq, SimdDouble *e, SimdDouble *fscal)
    {
        auto       beta    = p[1];
        auto       De      = p[2];
        SimdDouble dr      = simd::sqrt(rsq);
        SimdDouble expterm = simd::exp(-beta*(dr - p[0]));
        SimdDouble one_exp = 1.0 - expterm;
        e[0]  += De*(one_exp*one_exp - 1.0) + p[3];
        *fscal = 2*De*beta*expterm*(expterm - 1.0)/dr;
    }
};

//! Hua bond
struct HuaBond
{
    //! \copydoc LJ12_6::prepare
    static bool prepare(const TopologyEntry &b, const std::vector<ActAtom> &, double, double *p)
    {
        auto &params = b.params();
        p[0] = params[huaLENGTH];
        p[1] = params[huaB];
        p[2] = params[huaDE];
        p[3] = params[huaC];
        return true;
    }
    //! \copydoc LJ12_6::compute
    static void compute(const double *p, SimdDouble rsq, SimdDouble *e, SimdDouble *fscal)
    {
        auto       hb      = p[1];
        auto       De      = p[2];
        auto       c       = p[3];
        SimdDouble dr      = simd::sqrt(rsq);
        SimdDouble expterm = simd::exp(-hb*(dr - p[0]));
        SimdDouble exp2    = (1.0 - expterm)/(1.0 - c*expterm);
        SimdDouble denom   = c*expterm - 1.0;
        e[0]  += De*(exp2*exp2 - 1.0);
        *fscal = 2*De*hb*expterm*(expterm - 1.0)*(c - 1)/(dr*denom*denom*denom);
    }
};

/*! \brief Compute a pair potential for all conformations
 * \tparam       Kernel   The potential
 * \param[in]    entries  The interactions
 * \param[in]    atoms    The atoms
 * \param[in]    epsilonr Relative dielectric constant
 * \param[in]    stride   Number of conformations, a multiple of simdWidth
 * Without forces the compiler removes their computation from the kernel,
 * since the force is not used.
 * \tparam       Forces   Whether to compute forces
 * \param[in]    x        Coordinates
 * \param[inout] f        Forces, not used if Forces is false
 * \param[inout] energies Energies
 */
template <typename Kernel, bool Forces>
void computePairs(const TopologyEntryVector  &entries,
                  const std::vector<ActAtom> &atoms,
                  double                      epsilonr,
                  size_t                      stride,
                  const double               *x,
                  double                     *f,
                  double                     *energies)
{
    std::vector<PairParams> pairs;
    pairs.reserve(entries.size());
    for (const auto &b : entries)
    {
        const TopologyEntry &entry = *b.operator->();
        PairParams           pp;
        pp.i = DIM*entry.atomIndices()[0]*stride;
        pp.j = DIM*entry.atomIndices()[1]*stride;
        if (Kernel::prepare(entry, atoms, epsilonr, pp.p))
        {
            pairs.push_back(pp);
        }
    }
    for (const auto &pp : pairs)
    {
        for (size_t c = 0; c < stride; c += simd::simdWidth)
        {
            SimdDouble dx[DIM];
            SimdDouble rsq = 0.0;
            for (int m = 0; m < DIM; m++)
            {
                dx[m] = simd::load(x + pp.i + m*stride + c) - simd::load(x + pp.j + m*stride + c);
                rsq   = simd::fma(dx[m], dx[m], rsq);
            }
            SimdDouble e[explicitSimdEnergyTerms];
            for (int k = 0; k < explicitSimdEnergyTerms; k++)
            {
                e[k] = simd::load(energies + k*stride + c);
            }
            SimdDouble fscal;
            Kernel::compute(pp.p, rsq, e, &fscal);
            for (int k = 0; k < explicitSimdEnergyTerms; k++)
            {
                simd::store(energies + k*stride + c, e[k]);
            }
            if constexpr (Forces)
            {
                for (int m = 0; m < DIM; m++)
                {
                    double *fi = f + pp.i + m*stride + c;
                    double *fj = f + pp.j + m*stride + c;
                    simd::store(fi, simd::fma(fscal, dx[m], simd::load(fi)));
                    simd::store(fj, simd::load(fj) - fscal*dx[m]);
                }
            }
        }
    }
}

} // namespace

size_t explicitSimdWidth()
{
    return simd::simdWidth;
}

ExplicitSimd explicitSimdType()
{
#if ACT_SIMD_AVX_512
    return ExplicitSimd::Avx512;
#elif ACT_SIMD_AVX2
    return ExplicitSimd::Avx2;
#else
    return ExplicitSimd::ArmNeon;
#endif
}

bool computeExplicitSimd(Potential                   mypot,
                         const TopologyEntryVector  &entries,
                         const std::vector<ActAtom> &atoms,
                         double                      epsilonr,
                         size_t                      stride,
                         const double               *x,
                         double                     *f,
                         double                     *energies,
                         bool                        forces)
{
    auto pairs = [&](auto kernel)
    {
        if (forces)
        {
            computePairs<decltype(kernel), true>(entries, atoms, epsilonr, stride, x, f, energies);
        }
        else
        {
            computePairs<decltype(kernel), false>(entries, atoms, epsilonr, stride, x, f, energies);
        }
        return true;
    };
    switch (mypot)
    {
    case Potential::LJ12_6:               return pairs(LJ12_6{});
    case Potential::LJ12_6_4:             return pairs(LJ12_6_4{});
    case Potential::LJ8_6:                return pairs(LJ8_6{});
    case Potential::LJ14_7:               return pairs(LJ14_7{});
    case Potential::BORN_MAYER:           return pairs(BornMayer{});
    case Potential::SLATER_ISA:           return pairs(SlaterISA{});
    case Potential::MACDANIEL_SCHMIDT:    return pairs(DoubleExponential{});
    case Potential::WANG_BUCKINGHAM:      return pairs(WangBuckingham{});
    case Potential::BUCKINGHAM:           return pairs(Buckingham{});
    case Potential::TANG_TOENNIES:        return pairs(TangToennies{});
    case Potential::TT2b:                 return pairs(TT2b{});
    case Potential::SLATER_ISA_TT:        return pairs(SlaterISATT{});
    case Potential::COULOMB_POINT:
    case Potential::COULOMB_GAUSSIAN:     return pairs(Coulomb{});
    case Potential::POLARIZATION:         return pairs(Polarization{});
    case Potential::HARMONIC_BONDS:       return pairs(HarmonicBond{});
    case Potential::CUBIC_BONDS:          return pairs(CubicBond{});
    case Potential::MORSE_BONDS:          return pairs(MorseBond{});
    case Potential::HUA_BONDS:            return pairs(HuaBond{});
    default:
        return false;
    }
}

} // namespace alexandria

#else

namespace alexandria
{

size_t explicitSimdWidth()
{
    return 1;
}

ExplicitSimd explicitSimdType()
{
    return ExplicitSimd::None;
}

bool computeExplicitSimd(Potential, const TopologyEntryVector &, const std::vector<ActAtom> &,
                         double, size_t, const double *, double *, double *, bool)
{
    return false;
}

} // namespace alexandria

#endif

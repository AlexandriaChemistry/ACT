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
 * Pair potentials, harmonic angles, and Fourier and improper dihedrals
 * for many conformations at once with explicit SIMD.
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

#include <utility>

#include "pairunroll.h"

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
    //! The variant of the kernel, see pairVariants
    int    variant;
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
    /*! \copydoc LJ12_6::prepare
     * Each particle may have a second charge at the same position with its
     * own width, see ActAtom::chargeVS, such that there are up to four terms,
     * stored as pairs of the charge product and the combined width. Terms
     * without charge are stored with a zero charge product, except the first.
     */
    static bool prepare(const TopologyEntry &b, const std::vector<ActAtom> &atoms, double epsilonr, double *p)
    {
        auto &params   = b.params();
        auto &gindices = b.atomIndices();
        const auto &atomi = atoms[gindices[0]];
        const auto &atomj = atoms[gindices[1]];
        const double qi[2] = { atomi.charge(), atomi.chargeVS() };
        const double qj[2] = { atomj.charge(), atomj.chargeVS() };
        const double zi[2] = { params[coulZETA], atomi.zetaVS() };
        const double zj[2] = { params[coulZETA2], atomj.zetaVS() };
        for(int a = 0; a < 2; a++)
        {
            for(int c = 0; c < 2; c++)
            {
                double *pk = p + 2*(2*a + c);
                pk[0] = ONE_4PI_EPS0*qi[a]*qj[c]/epsilonr;
                // Width of the combined Gaussian, as in Coulomb_GG, zero for point charges
                if (zi[a] == 0 || zj[c] == 0)
                {
                    pk[1] = zi[a] + zj[c];
                }
                else
                {
                    pk[1] = zi[a]*zj[c]/std::sqrt(zi[a]*zi[a] + zj[c]*zj[c]);
                }
            }
        }
        return true;
    }
    /*! \brief The number of variants of compute, see pairVariants
     * Variant 0 is any pair, 1 a pair with one point charge each and 2 a
     * pair with one Gaussian charge each. Without second charges, the
     * common case, the variants 1 and 2 do not test the charges in the
     * loop over conformations.
     */
    static constexpr int nvariant = 3;
    /*! \brief The variant for an interaction
     * \param[in] p The parameters, see prepare
     * \return the variant
     */
    static int variant(const double *p)
    {
        if (p[2] != 0 || p[4] != 0 || p[6] != 0)
        {
            return 0;
        }
        return (p[1] == 0) ? 1 : 2;
    }
    /*! \copydoc LJ12_6::compute
     * \tparam V The variant, see nvariant. The variants 1 and 2 do the
     *           same operations as variant 0 does for them.
     */
    template <int V>
    static void compute(const double *p, SimdDouble rsq, SimdDouble *e, SimdDouble *fscal)
    {
        // Particles on top of each other only have an energy with Gaussian charges
        auto       overlap = !(rsq > 0.0);
        SimdDouble rinv    = simd::invsqrt(rsq);
        if constexpr (V == 1)
        {
            auto qq = p[0];
            e[0]  += simd::select(overlap, 0.0, qq*rinv);
            *fscal = simd::select(overlap, 0.0, qq*rinv*rinv*rinv);
            return;
        }
        else if constexpr (V == 2)
        {
            auto       qq   = p[0];
            auto       zeta = p[1];
            SimdDouble rz   = zeta*rsq*rinv;
            SimdDouble exprz2;
            SimdDouble erfrz = simd::erfAndExp(rz, &exprz2);
            e[0]  += simd::select(overlap, qq*2*zeta*invsqrt_pi, qq*erfrz*rinv);
            *fscal = simd::select(overlap, 0.0,
                                  qq*(erfrz*rinv - 2*invsqrt_pi*zeta*exprz2)*rinv*rinv);
            return;
        }
        for(int k = 0; k < 4; k++)
        {
            auto qq   = p[2*k];
            auto zeta = p[2*k+1];
            if (k > 0 && qq == 0)
            {
                continue;
            }
            SimdDouble ek, fk;
            if (zeta == 0)
            {
                ek = simd::select(overlap, 0.0, qq*rinv);
                fk = simd::select(overlap, 0.0, qq*rinv*rinv*rinv);
            }
            else
            {
                SimdDouble rz    = zeta*rsq*rinv;
                // The exponential for the force comes with the error function
                SimdDouble exprz2;
                SimdDouble erfrz = simd::erfAndExp(rz, &exprz2);
                ek = simd::select(overlap, qq*2*zeta*invsqrt_pi, qq*erfrz*rinv);
                fk = simd::select(overlap, 0.0,
                                  qq*(erfrz*rinv - 2*invsqrt_pi*zeta*exprz2)*rinv*rinv);
            }
            e[0]  += ek;
            *fscal = (k == 0) ? fk : *fscal + fk;
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

/*! \brief The number of variants of the compute function of a pair kernel
 * A kernel can have variants of compute for different kinds of
 * interactions, such that the loop over conformations does not test the
 * parameters, see Coulomb. It then has a constant nvariant, a function
 * variant giving the variant for the parameters and compute is a template
 * on the variant.
 * \tparam Kernel The potential
 * \return the number of variants, 1 if there are none
 */
template <typename Kernel>
static constexpr int pairVariants()
{
    if constexpr (requires { Kernel::nvariant; })
    {
        return Kernel::nvariant;
    }
    else
    {
        return 1;
    }
}

/*! \brief Compute one pair interaction for conformations c to c + simdWidth - 1
 * \tparam       Kernel   The potential
 * \tparam       Forces   Whether to compute forces
 * \tparam       V        The variant of the kernel, see pairVariants
 * \param[in]    pp       The interaction
 * \param[in]    stride   Number of conformations
 * \param[in]    c        The first conformation
 * \param[in]    x        Coordinates
 * \param[inout] f        Forces, not used if Forces is false
 * \param[inout] energies Energies
 */
template <typename Kernel, bool Forces, int V>
static inline void computePair(const PairParams &pp,
                               size_t            stride,
                               size_t            c,
                               const double     *x,
                               double           *f,
                               double           *energies)
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
    if constexpr (pairVariants<Kernel>() > 1)
    {
        Kernel::template compute<V>(pp.p, rsq, e, &fscal);
    }
    else
    {
        Kernel::compute(pp.p, rsq, e, &fscal);
    }
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

/*! \brief Compute one pair interaction for all conformations
 * The loop over conformations is unrolled pairUnroll times, such that the
 * instructions of that many SIMD vectors of conformations are available to
 * the execution units at once and the loop overhead is shared. Each
 * conformation is computed as without unrolling. The remaining SIMD
 * vectors, fewer than pairUnroll, are computed one at a time.
 * \tparam V The variant of the kernel, see pairVariants
 * The other parameters are those of computePair.
 */
#ifdef ACT_PAIR_UNROLL
//! The unroll factor of the pair loops, from the CMake option ACT_PAIR_UNROLL
constexpr size_t pairUnroll = ACT_PAIR_UNROLL;
#else
//! The unroll factor of the pair loops, the default of the SIMD type
constexpr size_t pairUnroll = simd::pairUnrollDefault;
#endif
static_assert(pairUnroll >= 1 && pairUnroll <= 16, "The unroll factor of the pair loops should be 1 to 16");

template <typename Kernel, bool Forces, int V>
static void computePairConformations(const PairParams &pp,
                                     size_t            stride,
                                     const double     *x,
                                     double           *f,
                                     double           *energies)
{
    constexpr size_t width = simd::simdWidth;
    size_t           c     = 0;
    for (; c + pairUnroll*width <= stride; c += pairUnroll*width)
    {
        // The vectors in increasing order
        [&]<size_t... u>(std::index_sequence<u...>)
        {
            (computePair<Kernel, Forces, V>(pp, stride, c + u*width, x, f, energies), ...);
        }(std::make_index_sequence<pairUnroll>{});
    }
    for (; c < stride; c += width)
    {
        computePair<Kernel, Forces, V>(pp, stride, c, x, f, energies);
    }
}

/*! \brief Compute one pair interaction with the variant of its kernel
 * \tparam V The first variant to test, see pairVariants
 * The other parameters are those of computePair.
 */
template <typename Kernel, bool Forces, int V = 0>
static void dispatchPairVariant(const PairParams &pp,
                                size_t            stride,
                                const double     *x,
                                double           *f,
                                double           *energies)
{
    if constexpr (V + 1 < pairVariants<Kernel>())
    {
        if (pp.variant != V)
        {
            dispatchPairVariant<Kernel, Forces, V + 1>(pp, stride, x, f, energies);
            return;
        }
    }
    computePairConformations<Kernel, Forces, V>(pp, stride, x, f, energies);
}

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
void computePairs(TopologyEntryRange          entries,
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
            if constexpr (pairVariants<Kernel>() > 1)
            {
                pp.variant = Kernel::variant(pp.p);
            }
            else
            {
                pp.variant = 0;
            }
            pairs.push_back(pp);
        }
    }
    // The interactions in the same order as without variants, such that
    // the forces and energies are summed in the same order
    for (const auto &pp : pairs)
    {
        dispatchPairVariant<Kernel, Forces>(pp, stride, x, f, energies);
    }
}

//! A harmonic angle with prepared parameters
struct AngleParams
{
    //! Offsets of the coordinates of the three atoms
    size_t i, j, k;
    //! The reference angle
    double theta0;
    //! The force constant
    double kt;
};

/*! \brief Compute harmonic angles for all conformations
 * As angles_simd in GROMACS, the cosine of the angle is computed from the
 * vectors and the angle with a SIMD arc cosine. The formulas are those of
 * angleOne in forces_simd.cpp.
 * \tparam       Forces   Whether to compute forces
 * \param[in]    entries  The interactions
 * \param[in]    stride   Number of conformations, a multiple of simdWidth
 * \param[in]    x        Coordinates
 * \param[inout] f        Forces, not used if Forces is false
 * \param[inout] energies Energies
 */
template <bool Forces>
void computeAngles(TopologyEntryRange         entries,
                   size_t                     stride,
                   const double              *x,
                   double                    *f,
                   double                    *energies)
{
    std::vector<AngleParams> angles;
    angles.reserve(entries.size());
    for (const auto &b : entries)
    {
        const TopologyEntry &entry   = *b.operator->();
        auto                &indices = entry.atomIndices();
        angles.push_back({ DIM*indices[0]*stride, DIM*indices[1]*stride, DIM*indices[2]*stride,
                           entry.params()[angleANGLE], entry.params()[angleKT] });
    }
    for (const auto &a : angles)
    {
        double halfk = 0.5*a.kt;
        for (size_t c = 0; c < stride; c += simd::simdWidth)
        {
            SimdDouble rij[DIM], rkj[DIM];
            SimdDouble ip = 0.0, ipa = 0.0, ipb = 0.0;
            for (int m = 0; m < DIM; m++)
            {
                SimdDouble xj = simd::load(x + a.j + m*stride + c);
                rij[m] = simd::load(x + a.i + m*stride + c) - xj;
                rkj[m] = simd::load(x + a.k + m*stride + c) - xj;
                ip     = simd::fma(rij[m], rkj[m], ip);
                ipa    = simd::fma(rij[m], rij[m], ipa);
                ipb    = simd::fma(rkj[m], rkj[m], ipb);
            }
            // As cos_angle, limited to [-1, 1]
            SimdDouble ipab  = ipa*ipb;
            SimdDouble costh = simd::select(ipab > 0.0, ip*simd::invsqrt(ipab), 1.0);
            costh            = simd::min(simd::max(costh, -1.0), 1.0);
            SimdDouble dth   = simd::acos(costh) - a.theta0;
            SimdDouble e     = simd::load(energies + c);
            simd::store(energies + c, simd::fma(halfk, dth*dth, e));
            if constexpr (Forces)
            {
                // There is no force for angles of 0 and 180 degrees
                SimdDouble costh2 = costh*costh;
                SimdDouble st     = simd::select(costh2 < 1.0,
                                                 (-a.kt*dth)*simd::invsqrt(1.0 - costh2), 0.0);
                SimdDouble sth    = st*costh;
                SimdDouble nrij_1 = simd::invsqrt(ipa);
                SimdDouble nrkj_1 = simd::invsqrt(ipb);
                SimdDouble cik    = st*nrij_1*nrkj_1;
                SimdDouble cii    = sth*nrij_1*nrij_1;
                SimdDouble ckk    = sth*nrkj_1*nrkj_1;
                for (int m = 0; m < DIM; m++)
                {
                    SimdDouble f_im = -(cik*rkj[m] - cii*rij[m]);
                    SimdDouble f_km = -(cik*rij[m] - ckk*rkj[m]);
                    SimdDouble f_jm = -f_im - f_km;
                    double    *fi   = f + a.i + m*stride + c;
                    double    *fj   = f + a.j + m*stride + c;
                    double    *fk   = f + a.k + m*stride + c;
                    simd::store(fi, simd::load(fi) + f_im);
                    simd::store(fj, simd::load(fj) + f_jm);
                    simd::store(fk, simd::load(fk) + f_km);
                }
            }
        }
    }
}

//! A dihedral with prepared parameters
struct DihedralParams
{
    //! Offsets of the coordinates of the four atoms
    size_t i, j, k, l;
    //! Number of parameters
    size_t nparam;
    //! The parameters
    double c[fdihNR];
};

//! Prepare the dihedrals, with at most fdihNR parameters
static std::vector<DihedralParams> prepareDihedrals(TopologyEntryRange         entries,
                                                    size_t                     stride)
{
    std::vector<DihedralParams> dihedrals;
    dihedrals.reserve(entries.size());
    for (const auto &b : entries)
    {
        const TopologyEntry &entry   = *b.operator->();
        auto                &indices = entry.atomIndices();
        auto                &params  = entry.params();
        DihedralParams       dp      = { DIM*indices[0]*stride, DIM*indices[1]*stride,
                                         DIM*indices[2]*stride, DIM*indices[3]*stride,
                                         std::min(params.size(), static_cast<size_t>(fdihNR)), {} };
        for (size_t p = 0; p < dp.nparam; p++)
        {
            dp.c[p] = params[p];
        }
        dihedrals.push_back(dp);
    }
    return dihedrals;
}

//! The vectors of a dihedral, as in dih_angle in forces_simd.cpp
struct DihedralVectors
{
    //! Vectors between the atoms
    SimdDouble r_ij[DIM], r_kj[DIM], r_kl[DIM];
    //! Normals of the planes, m = r_ij x r_kj and n = r_kj x r_kl
    SimdDouble m[DIM], n[DIM];
    //! Squared lengths of m, n and r_kj
    SimdDouble iprm, iprn, nrkj2;
    //! Inner products m.n and r_ij.n
    SimdDouble mn, ipr;
};

//! Compute the vectors of dihedral d for conformations c to c + simdWidth - 1
static inline void dihedralVectors(const DihedralParams &d,
                                   size_t                stride,
                                   size_t                c,
                                   const double         *x,
                                   DihedralVectors      *v)
{
    for (int m = 0; m < DIM; m++)
    {
        SimdDouble xj = simd::load(x + d.j + m*stride + c);
        SimdDouble xk = simd::load(x + d.k + m*stride + c);
        v->r_ij[m] = simd::load(x + d.i + m*stride + c) - xj;
        v->r_kj[m] = xk - xj;
        v->r_kl[m] = xk - simd::load(x + d.l + m*stride + c);
    }
    auto cprod = [](const SimdDouble *a, const SimdDouble *b, SimdDouble *axb)
    {
        axb[XX] = a[YY]*b[ZZ] - a[ZZ]*b[YY];
        axb[YY] = a[ZZ]*b[XX] - a[XX]*b[ZZ];
        axb[ZZ] = a[XX]*b[YY] - a[YY]*b[XX];
    };
    cprod(v->r_ij, v->r_kj, v->m);
    cprod(v->r_kj, v->r_kl, v->n);
    v->iprm  = 0.0;
    v->iprn  = 0.0;
    v->nrkj2 = 0.0;
    v->mn    = 0.0;
    v->ipr   = 0.0;
    for (int m = 0; m < DIM; m++)
    {
        v->iprm  = simd::fma(v->m[m], v->m[m], v->iprm);
        v->iprn  = simd::fma(v->n[m], v->n[m], v->iprn);
        v->nrkj2 = simd::fma(v->r_kj[m], v->r_kj[m], v->nrkj2);
        v->mn    = simd::fma(v->m[m], v->n[m], v->mn);
        v->ipr   = simd::fma(v->r_ij[m], v->n[m], v->ipr);
    }
}

/*! \brief Add the forces of a dihedral, as do_dih_fup_noshiftf in forces_simd.cpp
 * \param[in]    d      The dihedral
 * \param[in]    stride Number of conformations
 * \param[in]    c      The first conformation
 * \param[in]    v      The vectors of the dihedral
 * \param[in]    nrkj   Length of r_kj
 * \param[in]    s      1/(|m| |n|), finite
 * \param[in]    ddphi  Derivative of the energy with respect to the angle
 * \param[inout] f      The forces
 */
static inline void dihedralForces(const DihedralParams  &d,
                                  size_t                 stride,
                                  size_t                 c,
                                  const DihedralVectors &v,
                                  SimdDouble             nrkj,
                                  SimdDouble             s,
                                  SimdDouble             ddphi,
                                  double                *f)
{
    // There is no force if one of the planes is not defined
    SimdDouble toler  = v.nrkj2*GMX_REAL_EPS;
    simd::SimdBool ok = v.iprm > toler && v.iprn > toler;
    SimdDouble dd     = simd::select(ok, ddphi, 0.0)*nrkj;
    // 1/iprm and 1/iprn from s, such that there is no division
    SimdDouble s2     = s*s;
    SimdDouble a      = -dd*(s2*v.iprn);
    SimdDouble b      = dd*(s2*v.iprm);
    SimdDouble nrkj_2 = 1.0/simd::select(ok, v.nrkj2, 1.0);
    SimdDouble p      = 0.0, q = 0.0;
    for (int m = 0; m < DIM; m++)
    {
        p = simd::fma(v.r_ij[m], v.r_kj[m], p);
        q = simd::fma(v.r_kl[m], v.r_kj[m], q);
    }
    p *= nrkj_2;
    q *= nrkj_2;
    for (int m = 0; m < DIM; m++)
    {
        SimdDouble f_i  = a*v.m[m];
        SimdDouble f_l  = b*v.n[m];
        SimdDouble svec = p*f_i - q*f_l;
        double    *fi   = f + d.i + m*stride + c;
        double    *fj   = f + d.j + m*stride + c;
        double    *fk   = f + d.k + m*stride + c;
        double    *fl   = f + d.l + m*stride + c;
        simd::store(fi, simd::load(fi) + f_i);
        simd::store(fj, simd::load(fj) - (f_i - svec));
        simd::store(fk, simd::load(fk) - (f_l + svec));
        simd::store(fl, simd::load(fl) + f_l);
    }
}

/*! \brief 1/(|m| |n|) of a dihedral, 1 if one of the planes is not defined
 * \param[in] v The vectors of the dihedral
 */
static inline SimdDouble invNormProduct(const DihedralVectors &v)
{
    SimdDouble mn2 = v.iprm*v.iprn;
    return simd::invsqrt(simd::select(mn2 > 0.0, mn2, 1.0));
}

/*! \brief Compute Fourier dihedrals for all conformations
 * The energy is a polynomial in cos(phi), with phi = pi + psi and psi the
 * angle between the planes. Rather than computing psi with atan2 and then
 * its cosine and sine, as fourDihOne in forces_simd.cpp does, these are
 * computed from the vectors: cos(psi) = m.n/(|m| |n|) and, since
 * m x n = r_kj (r_ij.n), sin(psi) = |r_kj| (r_ij.n)/(|m| |n|). If one of
 * the planes is not defined psi is zero, as with atan2(0, 0).
 * \tparam       Forces   Whether to compute forces
 * \param[in]    entries  The interactions
 * \param[in]    stride   Number of conformations, a multiple of simdWidth
 * \param[in]    x        Coordinates
 * \param[inout] f        Forces, not used if Forces is false
 * \param[inout] energies Energies
 */
template <bool Forces>
void computeFourierDihedrals(TopologyEntryRange         entries,
                             size_t                     stride,
                             const double              *x,
                             double                    *f,
                             double                    *energies)
{
    for (const auto &d : prepareDihedrals(entries, stride))
    {
        for (size_t c = 0; c < stride; c += simd::simdWidth)
        {
            DihedralVectors v;
            dihedralVectors(d, stride, c, x, &v);
            SimdDouble s       = invNormProduct(v);
            SimdDouble defined = v.iprm*v.iprn;
            SimdDouble cos_phi = -simd::select(defined > 0.0, v.mn*s, 1.0);
            // As in fourDihOne, energy sum_j c_j cos^j(phi) and its derivative
            SimdDouble cosfac  = 1.0;
            SimdDouble ener    = d.c[0];
            SimdDouble dvdphi  = 0.0;
            for (size_t j = 1; j < d.nparam; j++)
            {
                dvdphi  = simd::fma(static_cast<double>(j)*d.c[j], cosfac, dvdphi);
                cosfac *= cos_phi;
                ener    = simd::fma(d.c[j], cosfac, ener);
            }
            simd::store(energies + c, simd::load(energies + c) + ener);
            if constexpr (Forces)
            {
                // -sin(phi) = sin(psi)
                SimdDouble nrkj = simd::sqrt(v.nrkj2);
                dvdphi *= nrkj*v.ipr*s;
                dihedralForces(d, stride, c, v, nrkj, s, dvdphi, f);
            }
        }
    }
}

/*! \brief Compute harmonic improper dihedrals for all conformations
 * The angle is computed with a SIMD arc tangent, as in GROMACS. With
 * m x n = r_kj (r_ij.n), atan2(|r_kj| (r_ij.n), m.n) is the angle with
 * the sign of r_ij.n that dih_angle in forces_simd.cpp computes.
 * \tparam       Forces   Whether to compute forces
 * \param[in]    entries  The interactions
 * \param[in]    stride   Number of conformations, a multiple of simdWidth
 * \param[in]    x        Coordinates
 * \param[inout] f        Forces, not used if Forces is false
 * \param[inout] energies Energies
 */
template <bool Forces>
void computeImproperDihedrals(TopologyEntryRange         entries,
                              size_t                     stride,
                              const double              *x,
                              double                    *f,
                              double                    *energies)
{
    for (const auto &d : prepareDihedrals(entries, stride))
    {
        double kA    = d.c[idihKPHI];
        double halfk = 0.5*kA;
        for (size_t c = 0; c < stride; c += simd::simdWidth)
        {
            DihedralVectors v;
            dihedralVectors(d, stride, c, x, &v);
            SimdDouble nrkj = simd::sqrt(v.nrkj2);
            SimdDouble phi  = simd::atan2(nrkj*v.ipr, v.mn);
            simd::store(energies + c, simd::fma(halfk*phi, phi, simd::load(energies + c)));
            if constexpr (Forces)
            {
                dihedralForces(d, stride, c, v, nrkj, invNormProduct(v), kA*phi, f);
            }
        }
    }
}

} // namespace

size_t explicitSimdWidth()
{
    return simd::simdWidth;
}

size_t explicitSimdPairUnroll()
{
    return pairUnroll;
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
                         TopologyEntryRange          entries,
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
    case Potential::HARMONIC_ANGLES:
        if (forces)
        {
            computeAngles<true>(entries, stride, x, f, energies);
        }
        else
        {
            computeAngles<false>(entries, stride, x, f, energies);
        }
        return true;
    case Potential::FOURIER_DIHEDRALS:
        if (forces)
        {
            computeFourierDihedrals<true>(entries, stride, x, f, energies);
        }
        else
        {
            computeFourierDihedrals<false>(entries, stride, x, f, energies);
        }
        return true;
    case Potential::HARMONIC_DIHEDRALS:
        if (forces)
        {
            computeImproperDihedrals<true>(entries, stride, x, f, energies);
        }
        else
        {
            computeImproperDihedrals<false>(entries, stride, x, f, energies);
        }
        return true;
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

size_t explicitSimdPairUnroll()
{
    return 1;
}

ExplicitSimd explicitSimdType()
{
    return ExplicitSimd::None;
}

bool computeExplicitSimd(Potential, TopologyEntryRange, const std::vector<ActAtom> &,
                         double, size_t, const double *, double *, double *, bool)
{
    return false;
}

} // namespace alexandria

#endif

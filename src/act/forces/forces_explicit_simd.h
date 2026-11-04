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
 * Potentials computed for many conformations at once with explicit SIMD.
 *
 * The SIMD lanes hold different conformations. The implementation depends
 * on the processor, at present there are ones for ARM NEON, AVX2 and AVX-512,
 * see simd_math.h. All functions are available in double and single
 * precision, in every build, such that code can use double precision
 * where it needs it. Where no
 * implementation is available, the functions here report that, and the
 * code in forces_simd.cpp is used instead. Whether the explicit SIMD can
 * and should be used on the processor at hand is decided by
 * useExplicitSimd() in simd_information.h.
 * \author David van der Spoel <david.vanderspoel@icm.uu.se>
 */
#ifndef ACT_FORCES_EXPLICIT_SIMD_H
#define ACT_FORCES_EXPLICIT_SIMD_H

#include <array>
#include <cstddef>
#include <vector>

#include "act/basics/coordinate.h"
#include "act/forcefield/potential.h"
#include "act/molprop/topologyentry.h"
#include "act/topology/actatom.h"

namespace alexandria
{

//! Number of energy terms a potential can have
constexpr int explicitSimdEnergyTerms = 3;

//! The instruction sets that the explicit SIMD can be compiled for
enum class ExplicitSimd
{
    None,
    ArmNeon,
    Avx2,
    Avx512
};

/*! \brief Terms of the Tang-Toennies damping functions
 * The damping functions for C6, C8 and C10 use the terms (br)^k/k! up to
 * ttTerms. For br < ttTailMax, they are computed from the terms beyond,
 * up to ttTailTerms, which then have a relative error below 1e-17. For
 * br >= ttTailMax, the relative error of the difference 1 - exp(-br)
 * sum_k (br)^k/k! is below 1e-14 in double precision.
 */
constexpr int    ttTerms     = 10;
//! \copydoc ttTerms
constexpr int    ttTailTerms = 38;
//! \copydoc ttTerms
constexpr double ttTailMax   = 6;

/*! \brief Coefficients of the terms beyond ttTerms of the damping functions
 * The sum of (br)^k/k! for k = ttTerms+1 ... ttTailTerms is
 * (br)^(ttTerms+1)/(ttTerms+1)! times the polynomial in br with these
 * coefficients, (ttTerms+1)!/k!.
 */
constexpr std::array<double, ttTailTerms-ttTerms> ttTailCoefficients = []()
{
    std::array<double, ttTailTerms-ttTerms> c{};
    double                                  f = 1;
    for (int j = 0; j < ttTailTerms-ttTerms; j++)
    {
        if (j > 0)
        {
            f /= ttTerms+1+j;
        }
        c[j] = f;
    }
    return c;
}();

/*! \brief The instruction set the explicit SIMD was compiled for
 * Since the code is compiled for this instruction set, this function is
 * safe to call on any processor, but computeExplicitSimd is not.
 * \return the instruction set
 */
ExplicitSimd explicitSimdType();

/*! \brief The number of conformations computed at once
 * The arrays passed to computeExplicitSimd have to hold a multiple of
 * this number of conformations.
 * \tparam T The floating point type of the arrays, float or double
 * \return the number of conformations computed at once
 */
template <typename T>
size_t explicitSimdWidth();

/*! \brief How often the loops over conformations of pair interactions are unrolled
 * \return the unroll factor, the default of the instruction set or that of
 *         the CMake option ACT_PAIR_UNROLL
 */
size_t explicitSimdPairUnroll();

/*! \brief Compute a potential for all conformations with explicit SIMD
 * Coordinate m of atom a in slot c is stored at x[(DIM*a+m)*stride + c],
 * in units of coordinateUnit, see coordinate.h,
 * the forces likewise and energy term k at energies[k*stride + c]. Forces
 * and energies are added to. Slots 0 to nslot - 1 are computed.
 * \param[in]    mypot    The potential
 * \param[in]    entries  The interactions
 * \param[in]    atoms    The atoms
 * \param[in]    epsilonr Relative dielectric constant
 * \param[in]    stride   Distance between the rows of the arrays
 * \param[in]    nslot    Number of slots to compute, a multiple of explicitSimdWidth<double>(),
 *                        at most stride
 * \param[in]    x        Coordinates
 * \param[inout] f        Forces, not used if forces is false
 * \param[inout] energies Energies
 * \param[in]    forces   Whether to compute forces
 * Only call this if useExplicitSimd() returns true, since the processor
 * may not support the instructions otherwise.
 * \return false if this potential is not available with explicit SIMD
 */
bool computeExplicitSimd(Potential                   mypot,
                         TopologyEntryRange          entries,
                         const std::vector<ActAtom> &atoms,
                         double                      epsilonr,
                         size_t                      stride,
                         size_t                      nslot,
                         const Coordinate           *x,
                         double                     *f,
                         double                     *energies,
                         bool                        forces);

/*! \brief Compute a potential for all conformations with explicit SIMD in single precision
 * As the double precision version, with a number of slots that is a
 * multiple of explicitSimdWidth<float>(). The parameters of the interactions are
 * prepared in double precision and rounded when they are used. The energy
 * of each interaction is computed in single precision and added to the
 * energies in double precision, such that sums of many interactions, and
 * differences of such sums, keep the accuracy of each interaction.
 */
bool computeExplicitSimd(Potential                   mypot,
                         TopologyEntryRange          entries,
                         const std::vector<ActAtom> &atoms,
                         double                      epsilonr,
                         size_t                      stride,
                         size_t                      nslot,
                         const Coordinate           *x,
                         float                      *f,
                         double                     *energies,
                         bool                        forces);

/*! \brief The Coulomb energies of interactions for charges of one with explicit SIMD
 * The interaction between two particles, each with a charge and possibly a
 * second one at the same position, has four terms, see ActAtom::chargeVS.
 * Term k = 2*a + b is the energy between charge a of the first particle
 * and charge b of the second one, both one, without ONE_4PI_EPS0 and the
 * dielectric constant, as COULOMB_POINT and COULOMB_GAUSSIAN compute it.
 * \param[in]  entries The interactions
 * \param[in]  atoms   The atoms, for the widths of the second charges
 * \param[in]  stride  Distance between the rows of the arrays
 * \param[in]  nslot   Number of slots to compute, a multiple of explicitSimdWidth<double>(),
 *                     at most stride
 * \param[in]  x       Coordinates, as in computeExplicitSimd
 * \param[in]  terms   For each interaction, bit k tells whether to compute term k
 * \param[out] kernel  Term k of interaction e of conformation c at (4*e + k)*stride + c,
 *                     the terms that are not computed are not changed
 * Only call this if useExplicitSimd() returns true.
 * \return false if there is no explicit SIMD
 */
bool coulombKernelsExplicitSimd(TopologyEntryRange          entries,
                                const std::vector<ActAtom> &atoms,
                                size_t                      stride,
                                size_t                      nslot,
                                const Coordinate           *x,
                                const unsigned char        *terms,
                                double                     *kernel);

/*! \brief The Coulomb energies for charges of one with explicit SIMD in single precision
 * As the double precision version, with a number of slots that is a
 * multiple of explicitSimdWidth<float>().
 */
bool coulombKernelsExplicitSimd(TopologyEntryRange          entries,
                                const std::vector<ActAtom> &atoms,
                                size_t                      stride,
                                size_t                      nslot,
                                const Coordinate           *x,
                                const unsigned char        *terms,
                                float                      *kernel);

/*! \brief Add the terms of one dimer to the sums of the second virial with explicit SIMD
 * The Boltzmann factor at temperature t is g = exp(-energy*beta[t]), zero
 * for overlapping atoms or if energy*beta[t] is 700 or more. Row k of the
 * sums is s[k*stride ...]: row 0 gets weight*(g - 1), rows 1 to 8 get
 * g*wx[k-1] and row 9 gets weight. Temperatures with a beta of zero get no
 * terms.
 * \param[in]    nT       Number of temperatures, a multiple of explicitSimdWidth<double>()
 * \param[in]    beta     1/(kB T) of each temperature
 * \param[in]    energy   The energy of the dimer
 * \param[in]    overlap  Whether the atoms overlap, then g is zero
 * \param[in]    weight   The weight of the dimer
 * \param[in]    wx       The weight times the squared forces and torque
 *                         components, all finite
 * \param[inout] s        The sums
 * \param[in]    stride   Distance between the rows of the sums
 * Only call this if useExplicitSimd() returns true.
 * \return false if there is no explicit SIMD
 */
bool addB2TermsExplicitSimd(size_t        nT,
                            const double *beta,
                            double        energy,
                            bool          overlap,
                            double        weight,
                            const double  wx[8],
                            double       *s,
                            size_t        stride);

} // namespace alexandria

#endif

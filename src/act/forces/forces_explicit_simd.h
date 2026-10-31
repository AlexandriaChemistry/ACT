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
 * see simd_math.h. Where no
 * implementation is available, the functions here report that, and the
 * code in forces_simd.cpp is used instead. Whether the explicit SIMD can
 * and should be used on the processor at hand is decided by
 * useExplicitSimd() in simd_information.h.
 * \author David van der Spoel <david.vanderspoel@icm.uu.se>
 */
#ifndef ACT_FORCES_EXPLICIT_SIMD_H
#define ACT_FORCES_EXPLICIT_SIMD_H

#include <cstddef>
#include <vector>

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

/*! \brief The instruction set the explicit SIMD was compiled for
 * Since the code is compiled for this instruction set, this function is
 * safe to call on any processor, but computeExplicitSimd is not.
 * \return the instruction set
 */
ExplicitSimd explicitSimdType();

/*! \brief The number of conformations computed at once
 * The arrays passed to computeExplicitSimd have to hold a multiple of
 * this number of conformations.
 * \return the number of conformations computed at once
 */
size_t explicitSimdWidth();

/*! \brief How often the loops over conformations of pair interactions are unrolled
 * \return the unroll factor, the default of the instruction set or that of
 *         the CMake option ACT_PAIR_UNROLL
 */
size_t explicitSimdPairUnroll();

/*! \brief Compute a potential for all conformations with explicit SIMD
 * Coordinate m of atom a in conformation c is stored at
 * x[(DIM*a+m)*stride + c], the forces likewise and energy term k at
 * energies[k*stride + c]. Forces and energies are added to.
 * \param[in]    mypot    The potential
 * \param[in]    entries  The interactions
 * \param[in]    atoms    The atoms
 * \param[in]    epsilonr Relative dielectric constant
 * \param[in]    stride   Number of conformations, a multiple of explicitSimdWidth()
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
                         const double               *x,
                         double                     *f,
                         double                     *energies,
                         bool                        forces);

} // namespace alexandria

#endif

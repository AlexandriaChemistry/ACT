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
 * on the processor, at present there are ones for ARM NEON and AVX-512,
 * see simd_math.h. Where no
 * implementation is available, the functions here report that, and the
 * code in forces_simd.cpp is used instead. That code is also used when
 * the environment variable ACT_NO_EXPLICIT_SIMD is set, which allows
 * comparing the performance and results of the two.
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

/*! \brief The number of conformations computed at once
 * The arrays passed to computeExplicitSimd have to hold a multiple of
 * this number of conformations.
 */
size_t explicitSimdWidth();

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
 * \param[inout] f        Forces
 * \param[inout] energies Energies
 * \return false if this potential is not available with explicit SIMD,
 *         or if the environment variable ACT_NO_EXPLICIT_SIMD is set
 */
bool computeExplicitSimd(Potential                   mypot,
                         const TopologyEntryVector  &entries,
                         const std::vector<ActAtom> &atoms,
                         double                      epsilonr,
                         size_t                      stride,
                         const double               *x,
                         double                     *f,
                         double                     *energies);

} // namespace alexandria

#endif

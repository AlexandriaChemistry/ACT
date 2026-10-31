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
 * Runtime detection of the SIMD instructions of the processor.
 *
 * Decides whether the explicit SIMD in forces_explicit_simd.cpp is used
 * and reports which SIMD instructions the processor has and which ones
 * the code was compiled for. This file is compiled with the normal
 * compiler flags, such that it runs on any processor.
 * \author David van der Spoel <david.vanderspoel@icm.uu.se>
 */
#ifndef ACT_FORCES_SIMD_INFORMATION_H
#define ACT_FORCES_SIMD_INFORMATION_H

#include "act/utility/timings.h"

namespace gmx
{
class TextWriter;
}

namespace alexandria
{

/*! \brief Whether to compute potentials with explicit SIMD
 * \return true if explicit SIMD was compiled in, the processor supports
 *         it and the environment variable ACT_NO_EXPLICIT_SIMD is not
 *         set. Setting that variable is useful for comparing the
 *         performance and results with and without explicit SIMD.
 */
bool useExplicitSimd();

/*! \brief Turn explicit SIMD on or off, for tests comparing the routines
 * Explicit SIMD is only turned on if useExplicitSimd would return true
 * without ACT_NO_EXPLICIT_SIMD, that is if it was compiled in and the
 * processor supports it.
 * \param[in] use Whether to use explicit SIMD
 * \return whether explicit SIMD was used before
 */
bool setUseExplicitSimd(bool use);

/*! \brief Whether explicit SIMD can be used on this processor
 * \return true if it was compiled in and the processor supports it,
 *         regardless of ACT_NO_EXPLICIT_SIMD
 */
bool explicitSimdAvailable();

/*! \brief The theoretical double precision floating point performance of one core
 * It is the maximum clock frequency times the floating point operations
 * per cycle of the widest SIMD instructions of the processor with fused
 * multiply-add, assuming the usual number of units for the processor.
 * The environment variable ACT_PEAK_GFLOPS, in GFlop/s, overrides it.
 * \return the peak performance, with zero GFlop/s if it is not known
 */
PeakPerformance theoreticalPeakPerCore();

/*! \brief Write the SIMD support of the processor and the code
 * \param[in] writer Where to write to, one line per item
 */
void writeSimdInformation(gmx::TextWriter *writer);

} // namespace alexandria

#endif

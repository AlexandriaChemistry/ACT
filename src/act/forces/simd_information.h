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

/*! \brief Write the SIMD support of the processor and the code
 * \param[in] writer Where to write to, one line per item
 */
void writeSimdInformation(gmx::TextWriter *writer);

} // namespace alexandria

#endif

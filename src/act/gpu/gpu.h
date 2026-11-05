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
 * Interface to the GPU code, without SYCL types, such that code compiled
 * with any C++ compiler can call it, see cmake/actManageGpu.cmake.
 * \author David van der Spoel <david.vanderspoel@icm.uu.se>
 */
#ifndef ACT_GPU_GPU_H
#define ACT_GPU_GPU_H

#include <string>
#include <vector>

namespace alexandria
{

namespace gpu
{

/*! \brief The programming model of the GPU code
 * \return for instance "SYCL (AdaptiveCpp 25.10.0)", or an empty string if
 *         ACT was built without GPU support
 */
std::string implementation();

/*! \brief The GPUs that the programming model finds
 * The first call initializes the runtime of the programming model.
 * AdaptiveCpp warns about each backend that it cannot load, for instance
 * that of AMD GPUs on a computer with NVIDIA GPUs only, unless the
 * environment variable ACPP_DEBUG_LEVEL is 1 or 0, as ACTRC sets it.
 * \return the names of the GPUs, empty without GPU support or without GPUs
 */
std::vector<std::string> devices();

} // namespace gpu

} // namespace alexandria

#endif

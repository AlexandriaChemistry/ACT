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
 * The GPUs that SYCL finds and the one that is used. This file is compiled with acpp, see
 * cmake/actManageGpu.cmake, and includes no ACT headers but gpu.h.
 * \author David van der Spoel <david.vanderspoel@icm.uu.se>
 */
#include "gpu.h"

#include <cstdlib>

#include <sycl/sycl.hpp>

#include "sycl/context.h"

namespace alexandria
{

namespace gpu
{

std::string implementation()
{
    return std::string("SYCL (AdaptiveCpp ") + ACT_ACPP_VERSION + ")";
}

std::vector<std::string> devices()
{
    std::vector<std::string> names;
    try
    {
        for (const auto &d : sycl::device::get_devices(sycl::info::device_type::gpu))
        {
            names.push_back(d.get_info<sycl::info::device::name>());
        }
    }
    catch (const sycl::exception &)
    {
        // No usable GPU runtime, as without GPUs
        names.clear();
    }
    return names;
}

sycl::queue *gpuQueue()
{
    // Never deleted, since the order of destruction at the end of the
    // program is unknown relative to that of the SYCL runtime
    static sycl::queue *queue = []()
    {
        sycl::queue *q = nullptr;
        try
        {
            auto gpus = sycl::device::get_devices(sycl::info::device_type::gpu);
            if (!gpus.empty())
            {
                q = new sycl::queue(gpus[0], sycl::property::queue::in_order());
            }
        }
        catch (const sycl::exception &)
        {
            q = nullptr;
        }
        return q;
    }();
    return queue;
}

//! Whether the GPU is used, see useGpu and setUseGpu
static bool &gpuInUse()
{
    static bool use = (gpuQueue() != nullptr && std::getenv("ACT_NO_GPU") == nullptr);
    return use;
}

bool useGpu()
{
    return gpuInUse();
}

bool setUseGpu(bool use)
{
    bool before = gpuInUse();
    gpuInUse()  = use && gpuQueue() != nullptr;
    return before;
}

} // namespace gpu

} // namespace alexandria

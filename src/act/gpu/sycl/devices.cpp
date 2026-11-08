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

#include <cstdio>
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

/*! \brief The single precision lanes of a compute unit of a GPU
 * NVIDIA streaming multiprocessors have 128 lanes, those of P100, V100 and
 * A100 64. AMD compute units have 64 lanes, without packed instructions.
 * \param[in] vendor The vendor of the GPU
 * \param[in] name   The name of the GPU
 * \return the number of lanes, 0 if not known
 */
static int lanesPerComputeUnit(const std::string &vendor, const std::string &name)
{
    auto has = [](const std::string &s, const char *part) { return s.find(part) != std::string::npos; };
    if (has(vendor, "NVIDIA") || has(name, "NVIDIA"))
    {
        for (const char *model : { "P100", "V100", "A100", "A30" })
        {
            if (has(name, model))
            {
                return 64;
            }
        }
        return 128;
    }
    if (has(vendor, "AMD") || has(name, "AMD") || has(name, "gfx"))
    {
        return 64;
    }
    return 0;
}

std::vector<DevicePeak> devicePeaks()
{
    std::vector<DevicePeak> peaks;
    try
    {
        for (const auto &d : sycl::device::get_devices(sycl::info::device_type::gpu))
        {
            DevicePeak  peak;
            std::string name   = d.get_info<sycl::info::device::name>();
            std::string vendor = d.get_info<sycl::info::device::vendor>();
            unsigned    units  = d.get_info<sycl::info::device::max_compute_units>();
            double      ghz    = d.get_info<sycl::info::device::max_clock_frequency>()*1e-3;
            int         lanes  = lanesPerComputeUnit(vendor, name);
            if (units > 0 && ghz > 0 && lanes > 0)
            {
                peak.gflops = 2.0*lanes*units*ghz;
                peak.basis  = std::to_string(units) + " compute units x " + std::to_string(lanes) +
                              " lanes x 2 flop (FMA) x ";
                char buf[32];
                std::snprintf(buf, sizeof(buf), "%.2f GHz", ghz);
                peak.basis += buf;
            }
            else
            {
                peak.basis = "unknown lanes or frequency";
            }
            peaks.push_back(peak);
        }
    }
    catch (const sycl::exception &)
    {
        peaks.clear();
    }
    return peaks;
}

/*! \brief The number of this process among those on the same node
 * From the environment variables that Slurm and the MPI libraries set,
 * zero if there is none.
 * \return the local rank
 */
static int localRank()
{
    for (const char *name : { "SLURM_LOCALID", "OMPI_COMM_WORLD_LOCAL_RANK", "MPI_LOCALRANKID",
                              "PMI_LOCAL_RANK", "MV2_COMM_WORLD_LOCAL_RANK" })
    {
        const char *value = std::getenv(name);
        if (value && *value)
        {
            return std::atoi(value);
        }
    }
    return 0;
}

int gpuIndex()
{
    static const int index = []()
    {
        try
        {
            auto gpus = sycl::device::get_devices(sycl::info::device_type::gpu);
            return gpus.empty() ? -1 : static_cast<int>(localRank() % gpus.size());
        }
        catch (const sycl::exception &)
        {
            return -1;
        }
    }();
    return index;
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
                // The processes on a node are spread over its GPUs
                // In order, and without a backend event for each operation,
                // since ACT waits for the queue and never for single events.
                // This shortens the gaps between the many short kernels.
                q = new sycl::queue(gpus[gpuIndex()], sycl::property_list{ sycl::property::queue::in_order(),
                                                                           sycl::property::queue::AdaptiveCpp_coarse_grained_events() });
                // A GPU may take only a limited number of processes, as with
                // the NVIDIA multi-process service, then a kernel fails, and
                // this process computes without it
                int *d     = sycl::malloc_device<int>(1, *q);
                int  value = 0;
                if (d)
                {
                    q->parallel_for(sycl::range<1>(1), [=](sycl::id<1>) { *d = 1993; });
                    q->memcpy(&value, d, sizeof(int)).wait_and_throw();
                    sycl::free(d, *q);
                }
                if (value != 1993)
                {
                    q = nullptr;
                }
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
    // Without starting the GPU when it is turned off, since a GPU may
    // only take a limited number of processes
    static bool use = (std::getenv("ACT_NO_GPU") == nullptr && gpuQueue() != nullptr);
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

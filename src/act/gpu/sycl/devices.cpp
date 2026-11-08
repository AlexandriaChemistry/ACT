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

#include <algorithm>
#include <cstdio>
#include <cstdlib>
#include <filesystem>
#include <fstream>

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

/*! \brief The directory of the AdaptiveCpp runtime library this process loaded
 * \return the directory, empty if not known, for instance on macOS
 */
static std::string acppLibraryDirectory()
{
    std::ifstream maps("/proc/self/maps");
    std::string   line;
    while (std::getline(maps, line))
    {
        auto pos = line.find('/');
        if (pos != std::string::npos && line.find("libacpp-rt") != std::string::npos)
        {
            return std::filesystem::path(line.substr(pos)).parent_path().string();
        }
    }
    return "";
}

/*! \brief The backends of AdaptiveCpp, the plugins next to its runtime library
 * \param[in] libdir The directory of the runtime library
 * \return the names of the backends, such as cuda, hip, ocl and omp
 */
static std::vector<std::string> acppBackends(const std::string &libdir)
{
    std::vector<std::string> names;
    std::error_code          ec;
    const std::string        prefix = "librt-backend-";
    for (const auto &entry : std::filesystem::directory_iterator(libdir + "/hipSYCL", ec))
    {
        std::string file = entry.path().stem().string();
        if (file.compare(0, prefix.size(), prefix) == 0)
        {
            names.push_back(file.substr(prefix.size()));
        }
    }
    std::sort(names.begin(), names.end());
    return names;
}

std::vector<std::string> noGpuReasons()
{
    std::vector<std::string> reasons;
    auto exists = [](const char *file)
    {
        std::error_code ec;
        return std::filesystem::exists(file, ec);
    };
    bool nvidia = exists("/dev/nvidiactl") || exists("/dev/nvidia0");
    bool amd    = exists("/dev/kfd");
    // A queue system that gives a job no GPU may hide them like this
    for (const char *name : { "CUDA_VISIBLE_DEVICES", "ROCR_VISIBLE_DEVICES", "HIP_VISIBLE_DEVICES" })
    {
        const char *value = std::getenv(name);
        if (value && (std::string(value).empty() || std::string(value) == "NoDevFiles" ||
                      std::string(value) == "-1"))
        {
            reasons.push_back(std::string(name) + " is \"" + value + "\", the queue system gave this job no GPU,");
            reasons.push_back("request one, for instance with sbatch --gres=gpu:1");
        }
    }
    std::string libdir   = acppLibraryDirectory();
    auto        backends = libdir.empty() ? std::vector<std::string>() : acppBackends(libdir);
    auto        has      = [&backends](const char *name)
    {
        return std::find(backends.begin(), backends.end(), name) != backends.end();
    };
    std::string list;
    for (const auto &b : backends)
    {
        list += (list.empty() ? "" : " ") + b;
    }
    if (!backends.empty())
    {
        if (nvidia && !has("cuda"))
        {
            reasons.push_back("this node has NVIDIA GPUs, but AdaptiveCpp has no CUDA backend, build");
            reasons.push_back("ACT with an AdaptiveCpp that has one, for instance from conda-forge");
            reasons.push_back("\"adaptivecpp=*=cuda*\", given with -DAdaptiveCpp_DIR");
        }
        if (amd && !has("hip"))
        {
            reasons.push_back("this node has AMD GPUs, but AdaptiveCpp has no HIP backend, build ACT");
            reasons.push_back("with an AdaptiveCpp that has one, given with -DAdaptiveCpp_DIR");
        }
    }
    if (!nvidia && !amd)
    {
        reasons.push_back("this node has no NVIDIA or AMD GPU device files (/dev/nvidiactl, /dev/kfd),");
        reasons.push_back("run on a node with GPUs");
    }
    if (!libdir.empty())
    {
        reasons.push_back("AdaptiveCpp in " + libdir + " has the backends: " + (list.empty() ? "none" : list));
    }
    reasons.push_back("acpp-info -l lists the devices that AdaptiveCpp finds, and the environment");
    reasons.push_back("variable ACPP_DEBUG_LEVEL=2 shows its messages");
    return reasons;
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

//! Whether turnOffGpu was called
static bool &gpuOff()
{
    static bool off = false;
    return off;
}

//! Whether the GPU is used, see useGpu and setUseGpu
static bool &gpuInUse()
{
    // Without starting the GPU when it is turned off, since a GPU may
    // only take a limited number of processes
    static bool use = (!gpuOff() && gpuQueue() != nullptr);
    return use;
}

void turnOffGpu()
{
    gpuOff()   = true;
    gpuInUse() = false;
}

bool gpuTurnedOff()
{
    return gpuOff();
}

bool useGpu()
{
    return gpuInUse();
}

bool setUseGpu(bool use)
{
    bool before = gpuInUse();
    gpuInUse()  = use && !gpuOff() && gpuQueue() != nullptr;
    return before;
}

} // namespace gpu

} // namespace alexandria

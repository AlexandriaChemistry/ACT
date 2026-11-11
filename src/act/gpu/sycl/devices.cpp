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
    NoGpuFacts facts;
    auto       exists = [](const char *file)
    {
        std::error_code ec;
        return std::filesystem::exists(file, ec);
    };
    facts.nvidiaDevices = exists("/dev/nvidiactl") || exists("/dev/nvidia0");
    facts.amdDevices    = exists("/dev/kfd");
    for (const char *name : { "CUDA_VISIBLE_DEVICES", "ROCR_VISIBLE_DEVICES", "HIP_VISIBLE_DEVICES" })
    {
        const char *value = std::getenv(name);
        if (value)
        {
            facts.visibleDevices.push_back({ name, value });
        }
    }
    facts.runtimeDirectory = acppLibraryDirectory();
    if (!facts.runtimeDirectory.empty())
    {
        facts.backends = acppBackends(facts.runtimeDirectory);
    }
    return explainNoGpu(facts);
}

std::vector<DevicePeak> devicePeaks()
{
    std::vector<DevicePeak> peaks;
    try
    {
        for (const auto &d : sycl::device::get_devices(sycl::info::device_type::gpu))
        {
            std::string name   = d.get_info<sycl::info::device::name>();
            std::string vendor = d.get_info<sycl::info::device::vendor>();
            unsigned    units  = d.get_info<sycl::info::device::max_compute_units>();
            double      ghz    = d.get_info<sycl::info::device::max_clock_frequency>()*1e-3;
            peaks.push_back(estimatePeak(vendor, name, units, ghz));
        }
    }
    catch (const sycl::exception &)
    {
        peaks.clear();
    }
    return peaks;
}

/*! \brief Let AdaptiveCpp report only errors, unless ACPP_DEBUG_LEVEL is set
 * Otherwise it warns about each GPU backend that it cannot load, for instance
 * that of AMD GPUs on a computer with NVIDIA GPUs only. A static object in the
 * headers of AdaptiveCpp reads its settings when the static objects of this
 * library are made, before main, so this runs before them, as a constructor
 * with a priority.
 */
__attribute__((constructor(101))) static void quietAdaptiveCpp()
{
    setenv("ACPP_DEBUG_LEVEL", "1", 0);
}

/*! \brief The number of this process among those on the same node
 * From the environment variables that Slurm and the MPI libraries set,
 * zero if there is none.
 * \return the local rank
 */
static int localRank()
{
    // Those of the MPI libraries first: mpirun in a Slurm job inherits the
    // SLURM_LOCALID of the job script, which is then zero in all processes
    for (const char *name : { "OMPI_COMM_WORLD_LOCAL_RANK", "MPI_LOCALRANKID", "PMI_LOCAL_RANK",
                              "MV2_COMM_WORLD_LOCAL_RANK", "SLURM_LOCALID" })
    {
        const char *value = std::getenv(name);
        if (value && *value)
        {
            return std::atoi(value);
        }
    }
    return 0;
}

//! The GPU of the node chosen by selectGpuOfProcess
static std::string &chosenGpu()
{
    static std::string chosen;
    return chosen;
}

void selectGpuOfProcess()
{
    int  rank   = localRank();
    bool listed = false;
    for (const char *name : { "CUDA_VISIBLE_DEVICES", "ROCR_VISIBLE_DEVICES", "HIP_VISIBLE_DEVICES" })
    {
        const char *value = std::getenv(name);
        if (!value || !*value)
        {
            continue;
        }
        listed = true;
        std::string              all(value);
        std::vector<std::string> ids;
        size_t                   start = 0;
        while (start <= all.size())
        {
            size_t end = all.find(',', start);
            end        = end == std::string::npos ? all.size() : end;
            if (end > start)
            {
                ids.push_back(all.substr(start, end - start));
            }
            start = end + 1;
        }
        if (ids.size() > 1)
        {
            const std::string &id = ids[rank % ids.size()];
            chosenGpu() = id + " of " + all;
            setenv(name, id.c_str(), 1);
        }
    }
    if (!listed)
    {
        // The NVIDIA GPUs of the node
        int n = 0;
        while (std::filesystem::exists("/dev/nvidia" + std::to_string(n)))
        {
            n++;
        }
        if (n > 1)
        {
            std::string id = std::to_string(rank % n);
            chosenGpu()    = id + " of " + std::to_string(n);
            setenv("CUDA_VISIBLE_DEVICES", id.c_str(), 1);
        }
    }
}

std::string gpuOfProcess()
{
    return chosenGpu();
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

//! The reason given to turnOffGpu
static std::string &offReason()
{
    static std::string reason;
    return reason;
}

void turnOffGpu(const std::string &reason)
{
    gpuOff()     = true;
    offReason()  = reason;
    gpuInUse()   = false;
}

bool gpuTurnedOff()
{
    return gpuOff();
}

std::string gpuOffReason()
{
    return offReason();
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

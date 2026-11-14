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
 * The parts of the GPU information of gpu.h that need no GPU runtime, such
 * that they can be tested in every build: the peak estimate and the
 * explanation of why no GPU was found. Compiled into the GPU library with
 * GPU support, and into the ACT library without.
 * \author David van der Spoel <david.vanderspoel@icm.uu.se>
 */
#include "gpu.h"

#include <algorithm>
#include <cctype>
#include <cstdio>
#include <cstdlib>
#include <filesystem>
#include <fstream>
#include <system_error>

namespace alexandria
{

namespace gpu
{

/*! \brief The single precision lanes of a compute unit of a GPU, see estimatePeak
 * \param[in] vendor The vendor of the GPU
 * \param[in] name   The name of the GPU
 * \return the number of lanes, 0 if not known
 */
static int lanesPerComputeUnit(const std::string &vendor, const std::string &name)
{
    const auto has = [](const std::string &s, const char *part) { return s.find(part) != std::string::npos; };
    // A model name by itself, such that the RTX A3000 is not an A30
    const auto model = [](const std::string &s, const std::string &part)
    {
        for (auto pos = s.find(part); pos != std::string::npos; pos = s.find(part, pos + 1))
        {
            auto end = pos + part.size();
            if ((pos == 0 || !std::isalnum(static_cast<unsigned char>(s[pos-1]))) &&
                (end == s.size() || !std::isalnum(static_cast<unsigned char>(s[end]))))
            {
                return true;
            }
        }
        return false;
    };
    if (has(vendor, "NVIDIA") || has(name, "NVIDIA"))
    {
        for (const char *m : { "P100", "V100", "A100", "A30" })
        {
            if (model(name, m))
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

DevicePeak estimatePeak(const std::string &vendor, const std::string &name,
                        unsigned units, double ghz)
{
    DevicePeak peak;
    int        lanes = lanesPerComputeUnit(vendor, name);
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
    return peak;
}

/*! \brief The first number in an environment variable
 * \param[in] name The variable
 * \return the number, 0 if the variable is not set or does not start with one,
 *         such as SLURM_STEP_TASKS_PER_NODE, which can be "32(x2)"
 */
static int numberInEnvironment(const char *name)
{
    const char *value = std::getenv(name);
    return value && std::isdigit(static_cast<unsigned char>(*value)) ? std::atoi(value) : 0;
}

OutOfMemoryFacts outOfMemoryFacts(const std::string &device)
{
    OutOfMemoryFacts facts;
    // Those of the MPI libraries first, as in localRank
    for (const char *name : { "OMPI_COMM_WORLD_LOCAL_SIZE", "MPI_LOCALNRANKS", "PMI_LOCAL_SIZE",
                              "MV2_COMM_WORLD_LOCAL_SIZE", "SLURM_STEP_TASKS_PER_NODE" })
    {
        facts.nodeProcesses = numberInEnvironment(name);
        if (facts.nodeProcesses > 0)
        {
            break;
        }
    }
    // selectGpuOfProcess tells "1 of 0,1,2,3" or "1 of 4" if the node has several
    const auto chosen = gpuOfProcess();
    auto of     = chosen.find(" of ");
    if (of != std::string::npos)
    {
        auto list = chosen.substr(of + 4);
        facts.nodeGpus = list.find(',') != std::string::npos
            ? static_cast<int>(std::count(list.begin(), list.end(), ',')) + 1
            : std::atoi(list.c_str());
    }
    else
    {
        std::error_code error;
        while (std::filesystem::exists("/dev/nvidia" + std::to_string(facts.nodeGpus), error))
        {
            facts.nodeGpus++;
        }
    }
    facts.nvidia = device.find("NVIDIA") != std::string::npos;
    if (facts.nvidia)
    {
        // The control daemon writes its process id there while it runs; the
        // pipes stay after it stopped
        const char   *pipes = std::getenv("CUDA_MPS_PIPE_DIRECTORY");
        std::ifstream pidFile(std::string(pipes && *pipes ? pipes : "/tmp/nvidia-mps") +
                              "/nvidia-cuda-mps-control.pid");
        long          pid = 0;
        if (pidFile >> pid && pid > 0)
        {
            std::error_code error;
            facts.mpsRuns = std::filesystem::exists("/proc/" + std::to_string(pid), error);
        }
    }
    return facts;
}

std::string outOfMemoryMessage(const std::string &device, double totalBytes,
                               size_t allocatedBytes, size_t requestedBytes,
                               const OutOfMemoryFacts &facts)
{
    char buf[512];
    std::snprintf(buf, sizeof(buf),
                  "The GPU (%s, %.1f GB) has no memory left: this process had %.1f MB on it "
                  "and asked for %.1f MB more.",
                  device.c_str(), totalBytes/1e9, allocatedBytes/1e6, requestedBytes/1e6);
    std::string message(buf);
    if (facts.nodeProcesses > 0)
    {
        int gpus = std::max(1, facts.nodeGpus);
        int per  = (facts.nodeProcesses + gpus - 1)/gpus;
        if (gpus == 1)
        {
            std::snprintf(buf, sizeof(buf), " All %d processes on this node use this GPU and share "
                          "its memory.", facts.nodeProcesses);
        }
        else
        {
            std::snprintf(buf, sizeof(buf), " The %d processes on this node share its %d GPUs, "
                          "up to %d processes per GPU.", facts.nodeProcesses, gpus, per);
        }
        message += buf;
    }
    else
    {
        message += " All processes that use this GPU share its memory.";
    }
    if (facts.nvidia && !facts.mpsRuns)
    {
        message += "\n\nThe NVIDIA multi-process service (MPS) does not run, so each process "
            "makes a context of its own on the GPU, which takes several hundred MB. "
            "Start MPS in the job script before mpirun or srun (see the section GPU support "
            "with SYCL of the installation chapter of the manual):\n"
            "  export CUDA_MPS_PIPE_DIRECTORY=/tmp/mps-$SLURM_JOB_ID\n"
            "  export CUDA_MPS_LOG_DIRECTORY=/tmp/mpslog-$SLURM_JOB_ID\n"
            "  mkdir -p $CUDA_MPS_PIPE_DIRECTORY $CUDA_MPS_LOG_DIRECTORY\n"
            "  nvidia-cuda-mps-control -d\n"
            "and stop it at the end with: echo quit | nvidia-cuda-mps-control";
    }
    message += "\n\nTo use less GPU memory:\n"
        "- run fewer processes per GPU: a smaller -n of mpirun or srun, or more GPUs "
        "(sbatch --gres=gpu:N, or more nodes);\n"
        "- use a smaller train_ff -lockstep, down to 1: a process keeps a copy of its "
        "molecules on the GPU for each individual that it computes at once;\n"
        "- in train_ff, use -nogpu_master, such that the master process does not use the GPU;\n"
        "- or train with -nogpu_resident, which keeps no molecules on the GPU.";
    return message;
}

OutOfMemory::OutOfMemory(const std::string &device, double totalBytes,
                         size_t allocatedBytes, size_t requestedBytes)
    : message_(outOfMemoryMessage(device, totalBytes, allocatedBytes, requestedBytes,
                                  outOfMemoryFacts(device)))
{
}

std::vector<std::string> explainNoGpu(const NoGpuFacts &facts)
{
    std::vector<std::string> reasons;
    // A queue system that gives a job no GPU may hide them like this
    for (const auto &vd : facts.visibleDevices)
    {
        if (vd.second.empty() || vd.second == "NoDevFiles" || vd.second == "-1")
        {
            reasons.push_back(vd.first + " is \"" + vd.second + "\", the queue system gave this job no GPU,");
            reasons.emplace_back("request one, for instance with sbatch --gres=gpu:1");
        }
    }
    const auto &backends = facts.backends;
    const auto  has      = [&backends](const char *name)
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
        if (facts.nvidiaDevices && !has("cuda"))
        {
            reasons.emplace_back("this node has NVIDIA GPUs, but AdaptiveCpp has no CUDA backend, build");
            reasons.emplace_back("ACT with an AdaptiveCpp that has one, for instance from conda-forge");
            reasons.emplace_back("\"adaptivecpp=*=cuda*\", given with -DAdaptiveCpp_DIR");
        }
        if (facts.amdDevices && !has("hip"))
        {
            reasons.emplace_back("this node has AMD GPUs, but AdaptiveCpp has no HIP backend, build ACT");
            reasons.emplace_back("with an AdaptiveCpp that has one, given with -DAdaptiveCpp_DIR");
        }
    }
    if (!facts.nvidiaDevices && !facts.amdDevices)
    {
        reasons.emplace_back("this node has no NVIDIA or AMD GPU device files (/dev/nvidiactl, /dev/kfd),");
        reasons.emplace_back("run on a node with GPUs");
    }
    if (!facts.runtimeDirectory.empty())
    {
        reasons.push_back("AdaptiveCpp in " + facts.runtimeDirectory + " has the backends: " +
                          (list.empty() ? "none" : list));
    }
    reasons.emplace_back("acpp-info -l lists the devices that AdaptiveCpp finds, and the environment");
    reasons.emplace_back("variable ACPP_DEBUG_LEVEL=2 shows its messages");
    return reasons;
}

} // namespace gpu

} // namespace alexandria

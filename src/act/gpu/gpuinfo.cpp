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
    auto has = [](const std::string &s, const char *part) { return s.find(part) != std::string::npos; };
    // A model name by itself, such that the RTX A3000 is not an A30
    auto model = [](const std::string &s, const std::string &part)
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

std::vector<std::string> explainNoGpu(const NoGpuFacts &facts)
{
    std::vector<std::string> reasons;
    // A queue system that gives a job no GPU may hide them like this
    for (const auto &vd : facts.visibleDevices)
    {
        if (vd.second.empty() || vd.second == "NoDevFiles" || vd.second == "-1")
        {
            reasons.push_back(vd.first + " is \"" + vd.second + "\", the queue system gave this job no GPU,");
            reasons.push_back("request one, for instance with sbatch --gres=gpu:1");
        }
    }
    const auto &backends = facts.backends;
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
        if (facts.nvidiaDevices && !has("cuda"))
        {
            reasons.push_back("this node has NVIDIA GPUs, but AdaptiveCpp has no CUDA backend, build");
            reasons.push_back("ACT with an AdaptiveCpp that has one, for instance from conda-forge");
            reasons.push_back("\"adaptivecpp=*=cuda*\", given with -DAdaptiveCpp_DIR");
        }
        if (facts.amdDevices && !has("hip"))
        {
            reasons.push_back("this node has AMD GPUs, but AdaptiveCpp has no HIP backend, build ACT");
            reasons.push_back("with an AdaptiveCpp that has one, given with -DAdaptiveCpp_DIR");
        }
    }
    if (!facts.nvidiaDevices && !facts.amdDevices)
    {
        reasons.push_back("this node has no NVIDIA or AMD GPU device files (/dev/nvidiactl, /dev/kfd),");
        reasons.push_back("run on a node with GPUs");
    }
    if (!facts.runtimeDirectory.empty())
    {
        reasons.push_back("AdaptiveCpp in " + facts.runtimeDirectory + " has the backends: " +
                          (list.empty() ? "none" : list));
    }
    reasons.push_back("acpp-info -l lists the devices that AdaptiveCpp finds, and the environment");
    reasons.push_back("variable ACPP_DEBUG_LEVEL=2 shows its messages");
    return reasons;
}

} // namespace gpu

} // namespace alexandria

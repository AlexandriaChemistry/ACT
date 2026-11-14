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
 * Tests of the GPU information of gpu.h: the peak estimate, the explanation
 * of why no GPU was found, and that what the log file shows is consistent.
 * \author David van der Spoel <david.vanderspoel@icm.uu.se>
 */
#include "../gpu.h"

#include <string>
#include <vector>

#include <gtest/gtest.h>

namespace alexandria
{

namespace
{

//! \return whether one of the lines contains the text
bool contains(const std::vector<std::string> &lines, const std::string &text)
{
    for (const auto &line : lines)
    {
        if (line.find(text) != std::string::npos)
        {
            return true;
        }
    }
    return false;
}

// Test: An AMD GPU has 64 lanes per compute unit, as the log shows for the MI250X
TEST(GpuPeakTest, AmdInstinct)
{
    const auto peak = gpu::estimatePeak("AMD", "AMD Instinct MI250X", 110, 1.7);
    EXPECT_NEAR(2*64*110*1.7, peak.gflops, 1e-9);
    EXPECT_EQ("110 compute units x 64 lanes x 2 flop (FMA) x 1.70 GHz", peak.basis);
}

// Test: The architecture name of an AMD GPU is enough, the vendor may be missing
TEST(GpuPeakTest, AmdArchitectureName)
{
    EXPECT_NEAR(2*64*104*1.7, gpu::estimatePeak("", "gfx90a:sramecc+:xnack-", 104, 1.7).gflops, 1e-9);
}

// Test: NVIDIA GPUs have 128 lanes per streaming multiprocessor
TEST(GpuPeakTest, NvidiaWith128Lanes)
{
    EXPECT_NEAR(2*128*46*1.65, gpu::estimatePeak("NVIDIA Corporation", "NVIDIA RTX A4500", 46, 1.65).gflops, 1e-9);
    EXPECT_NEAR(2*128*132*1.98, gpu::estimatePeak("NVIDIA", "NVIDIA GH200 480GB", 132, 1.98).gflops, 1e-9);
}

// Test: P100, V100, A100 and A30 have 64 lanes per streaming multiprocessor
TEST(GpuPeakTest, NvidiaWith64Lanes)
{
    for (const char *name : { "Tesla P100-PCIE-16GB", "Tesla V100-SXM2-32GB",
                              "NVIDIA A100-SXM4-40GB", "NVIDIA A100 80GB PCIe", "NVIDIA A30" })
    {
        SCOPED_TRACE(name);
        EXPECT_NEAR(2*64*108*1.41, gpu::estimatePeak("NVIDIA", name, 108, 1.41).gflops, 1e-9);
    }
}

// Test: The names of other GPUs that contain those of the 64-lane models have 128 lanes
TEST(GpuPeakTest, NvidiaModelNamesNotParts)
{
    for (const char *name : { "NVIDIA RTX A3000 Laptop GPU", "NVIDIA RTX A1000", "NVIDIA RTX A300X" })
    {
        SCOPED_TRACE(name);
        EXPECT_NEAR(2*128*16*1.5, gpu::estimatePeak("NVIDIA", name, 16, 1.5).gflops, 1e-9);
    }
}

// Test: Other vendors, and GPUs without compute units or frequency, have no estimate
TEST(GpuPeakTest, Unknown)
{
    for (const auto &peak : { gpu::estimatePeak("Intel(R) Corporation", "Intel(R) Data Center GPU Max 1550", 128, 1.6),
                              gpu::estimatePeak("AMD", "AMD Instinct MI250X", 0, 1.7),
                              gpu::estimatePeak("NVIDIA", "NVIDIA RTX A4500", 46, 0) })
    {
        EXPECT_EQ(0, peak.gflops);
        EXPECT_EQ("unknown lanes or frequency", peak.basis);
    }
}

// Test: A node without GPU device files gets that reason, and always how to look further
TEST(NoGpuTest, NodeWithoutGpus)
{
    const auto reasons = gpu::explainNoGpu(gpu::NoGpuFacts());
    EXPECT_TRUE(contains(reasons, "no NVIDIA or AMD GPU device files"));
    EXPECT_TRUE(contains(reasons, "run on a node with GPUs"));
    EXPECT_TRUE(contains(reasons, "acpp-info -l"));
    EXPECT_TRUE(contains(reasons, "ACPP_DEBUG_LEVEL=2"));
    EXPECT_FALSE(contains(reasons, "queue system"));
    EXPECT_FALSE(contains(reasons, "backend"));
}

// Test: A job that the queue system gave no GPU, as the values of the variables show
TEST(NoGpuTest, JobWithoutGpu)
{
    for (const char *value : { "", "NoDevFiles", "-1" })
    {
        SCOPED_TRACE(value);
        gpu::NoGpuFacts facts;
        facts.nvidiaDevices = true;
        facts.visibleDevices.emplace_back("CUDA_VISIBLE_DEVICES", value);
        const auto reasons = gpu::explainNoGpu(facts);
        EXPECT_TRUE(contains(reasons, std::string("CUDA_VISIBLE_DEVICES is \"") + value +
                                      "\", the queue system gave this job no GPU"));
        EXPECT_TRUE(contains(reasons, "sbatch --gres=gpu:1"));
        EXPECT_FALSE(contains(reasons, "device files"));
    }
}

// Test: A list of GPUs in the variables is no reason
TEST(NoGpuTest, JobWithGpus)
{
    gpu::NoGpuFacts facts;
    facts.amdDevices = true;
    facts.visibleDevices.emplace_back("ROCR_VISIBLE_DEVICES", "0,1");
    EXPECT_FALSE(contains(gpu::explainNoGpu(facts), "queue system"));
}

// Test: NVIDIA GPUs but no CUDA backend, with the backends that there are
TEST(NoGpuTest, NoCudaBackend)
{
    gpu::NoGpuFacts facts;
    facts.nvidiaDevices    = true;
    facts.runtimeDirectory = "/opt/acpp/lib";
    facts.backends         = { "hip", "omp" };
    const auto reasons = gpu::explainNoGpu(facts);
    EXPECT_TRUE(contains(reasons, "NVIDIA GPUs, but AdaptiveCpp has no CUDA backend"));
    EXPECT_TRUE(contains(reasons, "adaptivecpp=*=cuda*"));
    EXPECT_TRUE(contains(reasons, "AdaptiveCpp in /opt/acpp/lib has the backends: hip omp"));
    EXPECT_FALSE(contains(reasons, "HIP backend"));
    EXPECT_FALSE(contains(reasons, "device files"));
}

// Test: AMD GPUs but no HIP backend
TEST(NoGpuTest, NoHipBackend)
{
    gpu::NoGpuFacts facts;
    facts.amdDevices = true;
    facts.backends   = { "cuda", "omp" };
    const auto reasons = gpu::explainNoGpu(facts);
    EXPECT_TRUE(contains(reasons, "AMD GPUs, but AdaptiveCpp has no HIP backend"));
    EXPECT_FALSE(contains(reasons, "CUDA backend"));
}

// Test: The backend for the GPUs of the node is there, so no backend is blamed
TEST(NoGpuTest, MatchingBackend)
{
    gpu::NoGpuFacts facts;
    facts.nvidiaDevices = true;
    facts.amdDevices    = true;
    facts.backends      = { "cuda", "hip", "omp" };
    const auto reasons = gpu::explainNoGpu(facts);
    EXPECT_FALSE(contains(reasons, "no CUDA backend"));
    EXPECT_FALSE(contains(reasons, "no HIP backend"));
    EXPECT_TRUE(contains(reasons, "acpp-info -l"));
}

// Test: A runtime without backends says so
TEST(NoGpuTest, NoBackends)
{
    gpu::NoGpuFacts facts;
    facts.nvidiaDevices    = true;
    facts.runtimeDirectory = "/opt/acpp/lib";
    const auto reasons = gpu::explainNoGpu(facts);
    EXPECT_TRUE(contains(reasons, "AdaptiveCpp in /opt/acpp/lib has the backends: none"));
}

// Test: The message of running out of GPU memory tells the sizes and what to
// do, and the program stops with it as with other failed allocations
TEST(GpuOutOfMemoryTest, Message)
{
    try
    {
        throw gpu::OutOfMemory("NVIDIA RTX A4500", 20.0e9, 350000000, 1200000000);
    }
    catch (const std::bad_alloc &ex)
    {
        const std::string what = ex.what();
        EXPECT_NE(std::string::npos, what.find("NVIDIA RTX A4500, 20.0 GB")) << what;
        EXPECT_NE(std::string::npos, what.find("had 350.0 MB")) << what;
        EXPECT_NE(std::string::npos, what.find("1200.0 MB more")) << what;
        EXPECT_NE(std::string::npos, what.find("-lockstep")) << what;
        EXPECT_NE(std::string::npos, what.find("-nogpu_resident")) << what;
        EXPECT_NE(std::string::npos, what.find("-nogpu_master")) << what;
        return;
    }
    FAIL() << "Not caught as std::bad_alloc";
}

// Test: The message tells how many processes share each GPU, and how to start
// the NVIDIA multi-process service only if it does not run on an NVIDIA GPU
TEST(GpuOutOfMemoryTest, Facts)
{
    gpu::OutOfMemoryFacts facts;
    facts.nodeProcesses = 33;
    facts.nodeGpus      = 1;
    facts.nvidia        = true;
    const auto noMps = gpu::outOfMemoryMessage("NVIDIA RTX A4500", 21.2e9, 63600000, 122900000, facts);
    EXPECT_NE(std::string::npos, noMps.find("All 33 processes on this node use this GPU")) << noMps;
    EXPECT_NE(std::string::npos, noMps.find("(MPS) does not run")) << noMps;
    EXPECT_NE(std::string::npos, noMps.find("nvidia-cuda-mps-control -d")) << noMps;
    EXPECT_NE(std::string::npos, noMps.find("a smaller -n of mpirun or srun")) << noMps;

    facts.mpsRuns = true;
    const auto mps = gpu::outOfMemoryMessage("NVIDIA RTX A4500", 21.2e9, 63600000, 122900000, facts);
    EXPECT_EQ(std::string::npos, mps.find("does not run")) << mps;

    facts.nodeProcesses = 64;
    facts.nodeGpus      = 8;
    facts.nvidia        = false;
    facts.mpsRuns       = false;
    const auto amd = gpu::outOfMemoryMessage("AMD Instinct MI250X", 68.7e9, 1, 1, facts);
    EXPECT_NE(std::string::npos, amd.find("The 64 processes on this node share its 8 GPUs, up to 8 processes per GPU")) << amd;
    EXPECT_EQ(std::string::npos, amd.find("MPS")) << amd;

    facts.nodeProcesses = 0;
    const auto unknown = gpu::outOfMemoryMessage("AMD Instinct MI250X", 68.7e9, 1, 1, facts);
    EXPECT_NE(std::string::npos, unknown.find("All processes that use this GPU share")) << unknown;
}

// Test: What the log file shows about the GPUs of this computer is consistent:
// a peak estimate for each GPU, reasons only without GPUs and only with GPU
// support, and no GPU without GPU support
TEST(GpuInformationTest, Consistent)
{
    const auto names = gpu::devices();
    const auto peaks = gpu::devicePeaks();
    EXPECT_EQ(names.size(), peaks.size());
    if (gpu::implementation().empty())
    {
        EXPECT_TRUE(names.empty());
        EXPECT_TRUE(gpu::noGpuReasons().empty());
        EXPECT_FALSE(gpu::useGpu());
        EXPECT_EQ(-1, gpu::gpuIndex());
    }
    else if (names.empty())
    {
        EXPECT_TRUE(contains(gpu::noGpuReasons(), "acpp-info -l"));
        EXPECT_FALSE(gpu::useGpu());
    }
}

} // namespace

} // namespace alexandria

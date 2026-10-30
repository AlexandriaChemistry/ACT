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
 * \author David van der Spoel <david.vanderspoel@icm.uu.se>
 */
#include "actpre.h"

#include "simd_information.h"

#include <algorithm>
#include <cstdint>
#include <cstdlib>
#include <cstring>
#include <fstream>
#include <string>
#include <utility>
#include <vector>

#include "act/forces/forces_explicit_simd.h"
#include "gromacs/utility/stringutil.h"
#include "gromacs/utility/textwriter.h"

#if (defined(__x86_64__) || defined(__i386__)) && defined(__GNUC__)
#    include <cpuid.h>
#    define ACT_X86 1
#endif
#if defined(__aarch64__) && defined(__linux__)
#    include <sys/auxv.h>
#endif
#if defined(__APPLE__)
#    include <sys/sysctl.h>
#endif

namespace alexandria
{

namespace
{

//! \return The name of the processor, or an empty string if unknown
std::string processorName()
{
    char name[256] = { 0 };
#if ACT_X86
    unsigned int regs[12];
    if (__get_cpuid(0x80000000, &regs[0], &regs[1], &regs[2], &regs[3]) && regs[0] >= 0x80000004)
    {
        for (unsigned int i = 0; i < 3; i++)
        {
            __get_cpuid(0x80000002 + i, &regs[4*i], &regs[4*i+1], &regs[4*i+2], &regs[4*i+3]);
        }
        std::memcpy(name, regs, sizeof(regs));
    }
#elif defined(__APPLE__)
    size_t len = sizeof(name) - 1;
    if (sysctlbyname("machdep.cpu.brand_string", name, &len, nullptr, 0) != 0)
    {
        name[0] = '\0';
    }
#endif
    return gmx::stripString(name);
}

//! \return The SIMD instructions of the processor, detected at runtime
std::string processorSimd()
{
    std::string simd;
#if ACT_X86
    __builtin_cpu_init();
    // The names have to be literals, therefore no loop
    if (__builtin_cpu_supports("sse2"))    { simd += " SSE2"; }
    if (__builtin_cpu_supports("sse4.1"))  { simd += " SSE4.1"; }
    if (__builtin_cpu_supports("avx"))     { simd += " AVX"; }
    if (__builtin_cpu_supports("avx2"))    { simd += " AVX2"; }
    if (__builtin_cpu_supports("fma"))     { simd += " FMA"; }
    if (__builtin_cpu_supports("avx512f")) { simd += " AVX-512"; }
#elif defined(__aarch64__)
    // NEON is part of every aarch64 processor
    simd += " ARM-NEON";
#    if defined(__linux__) && defined(HWCAP_SVE)
    if (getauxval(AT_HWCAP) & HWCAP_SVE) { simd += " SVE"; }
#    endif
#endif
    return simd.empty() ? "none" : simd.substr(1);
}

#if ACT_X86
//! \return The vendor of an x86 processor, e.g. GenuineIntel or AuthenticAMD
std::string x86Vendor()
{
    unsigned int eax, ebx, ecx, edx;
    char         vendor[13] = { 0 };
    if (__get_cpuid(0, &eax, &ebx, &ecx, &edx))
    {
        std::memcpy(vendor, &ebx, 4);
        std::memcpy(vendor + 4, &edx, 4);
        std::memcpy(vendor + 8, &ecx, 4);
    }
    return vendor;
}

//! \return The family of an x86 processor, including the extended family
unsigned int x86Family()
{
    unsigned int eax = 0, ebx, ecx, edx;
    if (!__get_cpuid(1, &eax, &ebx, &ecx, &edx))
    {
        return 0;
    }
    unsigned int family = (eax >> 8) & 0xF;
    return family == 0xF ? family + ((eax >> 20) & 0xFF) : family;
}
#endif

/*! \brief The maximum clock frequency
 * \param[in]  name   Name of the processor
 * \param[out] source How it was determined
 * \return the frequency in GHz, or zero if not known
 */
double maxFrequency(const std::string &name, std::string *source)
{
#if defined(__linux__)
    {
        std::ifstream fmax("/sys/devices/system/cpu/cpu0/cpufreq/cpuinfo_max_freq");
        double        khz = 0;
        if (fmax >> khz && khz > 0)
        {
            *source = "maximum frequency from cpufreq";
            return 1e-6*khz;
        }
    }
#endif
#if defined(__APPLE__)
    uint64_t hz  = 0;
    size_t   len = sizeof(hz);
    if (sysctlbyname("hw.cpufrequency_max", &hz, &len, nullptr, 0) == 0 && hz > 0)
    {
        *source = "maximum frequency from sysctl";
        return 1e-9*hz;
    }
    // Apple silicon does not report its frequency, the nominal maximum of
    // the performance cores, the more specific names first
    const std::vector<std::pair<std::string, double>> apple = {
        { "Apple M4 Max", 4.51 }, { "Apple M4 Pro", 4.51 }, { "Apple M4", 4.4 },
        { "Apple M3", 4.05 },
        { "Apple M2 Ultra", 3.504 }, { "Apple M2 Max", 3.504 }, { "Apple M2 Pro", 3.504 }, { "Apple M2", 3.49 },
        { "Apple M1 Ultra", 3.228 }, { "Apple M1 Max", 3.228 }, { "Apple M1 Pro", 3.228 }, { "Apple M1", 3.204 }
    };
    for (const auto &chip : apple)
    {
        if (name.compare(0, chip.first.size(), chip.first) == 0)
        {
            *source = "nominal maximum of the performance cores";
            return chip.second;
        }
    }
#endif
    // Intel processor names end with the base frequency, e.g. @ 2.40GHz
    auto at = name.rfind('@');
    if (at != std::string::npos)
    {
        const char *begin = name.c_str() + at + 1;
        char       *end   = nullptr;
        double      ghz   = std::strtod(begin, &end);
        if (end != begin && ghz > 0 && std::strncmp(end, "GHz", 3) == 0)
        {
            *source = "base frequency from the processor name";
            return ghz;
        }
    }
#if defined(__linux__)
    // The current frequency, the highest of all cores
    std::ifstream cpuinfo("/proc/cpuinfo");
    std::string   line;
    double        mhz = 0;
    while (std::getline(cpuinfo, line))
    {
        if (line.compare(0, 7, "cpu MHz") == 0)
        {
            auto colon = line.find(':');
            if (colon != std::string::npos)
            {
                mhz = std::max(mhz, std::strtod(line.c_str() + colon + 1, nullptr));
            }
        }
    }
    if (mhz > 0)
    {
        *source = "highest current frequency in /proc/cpuinfo";
        return 1e-3*mhz;
    }
#endif
    return 0;
}

/*! \brief Double precision floating point operations per cycle of one core
 * Fused multiply-adds count as two operations, assuming two FMA units
 * on x86 processors, four on Apple silicon and two on other ARM processors.
 * \param[in]  name Name of the processor
 * \param[out] simd The instructions and units it is based on
 * \return the number of operations, or zero if not known
 */
int flopsPerCycle(const std::string &name, std::string *simd)
{
#if ACT_X86
    (void) name;
    __builtin_cpu_init();
    if (__builtin_cpu_supports("avx512f"))
    {
        // Zen 4 executes 512-bit instructions in two halves
        if (x86Vendor() == "AuthenticAMD" && x86Family() < 0x1A)
        {
            *simd = "AVX-512 executed as two 256-bit halves, two FMA units";
            return 16;
        }
        *simd = "AVX-512, two 512-bit FMA units";
        return 32;
    }
    if (__builtin_cpu_supports("avx2") && __builtin_cpu_supports("fma"))
    {
        *simd = "AVX2, two 256-bit FMA units";
        return 16;
    }
    if (__builtin_cpu_supports("avx"))
    {
        *simd = "AVX, 256-bit addition and multiplication";
        return 8;
    }
    *simd = "SSE2, 128-bit addition and multiplication";
    return 4;
#elif defined(__aarch64__)
    if (name.compare(0, 5, "Apple") == 0)
    {
        *simd = "ARM-NEON, four 128-bit FMA units";
        return 16;
    }
    *simd = "ARM-NEON, two 128-bit FMA units";
    return 8;
#else
    (void) name;
    *simd = "unknown instructions";
    return 0;
#endif
}

//! \return Whether the processor can run code for this instruction set
bool processorSupports(ExplicitSimd type)
{
    switch (type)
    {
    case ExplicitSimd::ArmNeon:
#if defined(__aarch64__)
        return true;
#else
        return false;
#endif
    case ExplicitSimd::Avx2:
#if ACT_X86
        __builtin_cpu_init();
        return __builtin_cpu_supports("avx2") && __builtin_cpu_supports("fma");
#else
        return false;
#endif
    case ExplicitSimd::Avx512:
#if ACT_X86
        // This also checks that the operating system saves the AVX-512 registers
        __builtin_cpu_init();
        return __builtin_cpu_supports("avx512f");
#else
        return false;
#endif
    default:
        return false;
    }
}

//! \return The name of an instruction set
const char *simdName(ExplicitSimd type)
{
    switch (type)
    {
    case ExplicitSimd::ArmNeon: return "ARM-NEON";
    case ExplicitSimd::Avx2:    return "AVX2";
    case ExplicitSimd::Avx512:  return "AVX-512";
    default:                    return "none";
    }
}

//! \return The SIMD instructions the compiler may use for the rest of the code
const char *compilerSimd()
{
#if defined(__AVX512F__)
    return "AVX-512";
#elif defined(__AVX2__)
    return "AVX2";
#elif defined(__AVX__)
    return "AVX";
#elif defined(__SSE4_1__)
    return "SSE4.1";
#elif defined(__SSE2__)
    return "SSE2";
#elif defined(__ARM_NEON)
    return "ARM-NEON";
#else
    return "none";
#endif
}

//! \return Whether the environment variable turns off the explicit SIMD
bool explicitSimdTurnedOff()
{
    return std::getenv("ACT_NO_EXPLICIT_SIMD") != nullptr;
}

} // namespace

bool useExplicitSimd()
{
    static const bool use = (processorSupports(explicitSimdType()) && !explicitSimdTurnedOff());
    return use;
}

PeakPerformance theoreticalPeakPerCore()
{
    PeakPerformance peak;
    const char     *env = std::getenv("ACT_PEAK_GFLOPS");
    if (env)
    {
        char   *end    = nullptr;
        double  gflops = std::strtod(env, &end);
        if (end != env && gflops > 0)
        {
            peak.gflops = gflops;
            peak.basis  = "given by ACT_PEAK_GFLOPS";
            return peak;
        }
    }
    auto        name = processorName();
    std::string simd, source;
    int         fpc  = flopsPerCycle(name, &simd);
    double      ghz  = maxFrequency(name, &source);
    std::string cpu  = name.empty() ? "unknown processor" : name;
    if (fpc > 0 && ghz > 0)
    {
        peak.gflops = fpc*ghz;
        peak.basis  = gmx::formatString("%s, %d flop/cycle (%s) at %.2f GHz (%s)",
                                        cpu.c_str(), fpc, simd.c_str(), ghz, source.c_str());
    }
    else
    {
        peak.basis  = gmx::formatString("%s, %s", cpu.c_str(),
                                        ghz > 0 ? "unknown operations per cycle" : "unknown frequency");
    }
    return peak;
}

void writeSimdInformation(gmx::TextWriter *writer)
{
    std::string name = processorName();
    if (!name.empty())
    {
        writer->writeLine(gmx::formatString("Processor:          %s", name.c_str()));
    }
    writer->writeLine(gmx::formatString("Processor SIMD:     %s (detected at runtime)", processorSimd().c_str()));
    writer->writeLine(gmx::formatString("Compiler SIMD:      %s (loops vectorized by the compiler)", compilerSimd()));
    ExplicitSimd type = explicitSimdType();
    std::string  status;
    if (type == ExplicitSimd::None)
    {
        status = "not compiled in, potentials vectorized by the compiler";
        if (processorSupports(ExplicitSimd::ArmNeon))
        {
            status += ", configure with -DACT_SIMD=ARM_NEON to use ARM-NEON";
        }
        else if (processorSupports(ExplicitSimd::Avx512))
        {
            status += ", configure with -DACT_SIMD=AVX_512 to use AVX-512";
        }
        else if (processorSupports(ExplicitSimd::Avx2))
        {
            status += ", configure with -DACT_SIMD=AVX2 to use AVX2";
        }
    }
    else if (!processorSupports(type))
    {
        status = gmx::formatString("%s compiled in, but not supported by this processor, not used", simdName(type));
    }
    else if (explicitSimdTurnedOff())
    {
        status = gmx::formatString("%s, turned off by ACT_NO_EXPLICIT_SIMD", simdName(type));
    }
    else
    {
        status = gmx::formatString("%s, %zu conformations at once, in use", simdName(type), explicitSimdWidth());
    }
    writer->writeLine(gmx::formatString("Explicit SIMD:      %s", status.c_str()));
    auto peak = theoreticalPeakPerCore();
    if (peak.gflops > 0)
    {
        writer->writeLine(gmx::formatString("Peak per core:      %.2f GFlop/s double precision, %s", peak.gflops, peak.basis.c_str()));
    }
    else
    {
        writer->writeLine(gmx::formatString("Peak per core:      not known, %s", peak.basis.c_str()));
    }
}

} // namespace alexandria

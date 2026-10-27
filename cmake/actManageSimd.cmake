#
# This source file is part of the Alexandria Chemistry Toolkit.
#
# Copyright (C) 2026
#
# Developers:
#             Mohammad Mehdi Ghahremanpour,
#             Julian Marrades,
#             Marie-Madeleine Walz,
#             Paul J. van Maaren,
#             David van der Spoel (Project leader)
#
# This program is free software; you can redistribute it and/or
# modify it under the terms of the GNU General Public License
# as published by the Free Software Foundation; either version 2
# of the License, or (at your option) any later version.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with this program; if not, write to the Free Software
# Foundation, Inc., 51 Franklin Street, Fifth Floor,
# Boston, MA  02110-1301, USA.
#

# Explicit SIMD for the force routines
#
# The force routines have explicit SIMD code for ARM NEON, AVX2 with FMA and
# AVX-512. The option ACT_SIMD selects one of them. The default AUTO selects
# the best one for the processor CMake runs on, NONE turns explicit SIMD off.
# Only forces_explicit_simd.cpp is compiled with the flags for the SIMD, and
# it is only used if the processor running the code supports them, otherwise
# the code vectorized by the compiler is used. Therefore the programs also
# run on processors with less SIMD than the one they were compiled on.
#
# Sets
#   ACT_SIMD_FLAGS      Compiler flags for the files with explicit SIMD
#   ACT_SIMD_ARM_NEON   1 if ARM NEON is used, otherwise 0, for config.h
#   ACT_SIMD_AVX2       Likewise for AVX2
#   ACT_SIMD_AVX_512    Likewise for AVX-512

include(CheckCXXCompilerFlag)

set(ACT_SIMD "AUTO" CACHE STRING "Explicit SIMD for the force routines: AUTO (the best for this processor), NONE, ARM_NEON, AVX2 (with FMA, for instance AMD Zen 2/3) or AVX_512")
set_property(CACHE ACT_SIMD PROPERTY STRINGS AUTO NONE ARM_NEON AVX2 AVX_512)

# Compiler flags for each choice
set(_act_simd_flags_NONE     "")
set(_act_simd_flags_ARM_NEON "")
set(_act_simd_flags_AVX2     "-mavx2;-mfma")
set(_act_simd_flags_AVX_512  "-mavx512f")

# Detect the best SIMD of the processor CMake runs on by running a small
# program. The result is stored in the cache, such that this is done once.
function(act_detect_simd result)
    if(NOT DEFINED ACT_SIMD_DETECTED)
        set(_detected "NONE")
        if(CMAKE_CROSSCOMPILING)
            message(STATUS "Cannot detect the SIMD of the processor when cross-compiling, select it with ACT_SIMD")
        else()
            set(_dir ${CMAKE_BINARY_DIR}${CMAKE_FILES_DIRECTORY}/ActDetectSimd)
            file(WRITE ${_dir}/detectsimd.cpp "
#include <cstdio>
int main()
{
#if (defined(__x86_64__) || defined(__i386__)) && defined(__GNUC__)
    __builtin_cpu_init();
    if (__builtin_cpu_supports(\"avx512f\"))
    {
        std::printf(\"AVX_512\");
    }
    else if (__builtin_cpu_supports(\"avx2\") && __builtin_cpu_supports(\"fma\"))
    {
        std::printf(\"AVX2\");
    }
    else
    {
        std::printf(\"NONE\");
    }
#elif defined(__aarch64__) && defined(__ARM_NEON)
    std::printf(\"ARM_NEON\");
#else
    std::printf(\"NONE\");
#endif
    return 0;
}
")
            try_run(_run_result _compile_result ${_dir} ${_dir}/detectsimd.cpp
                    RUN_OUTPUT_VARIABLE _output)
            if(_compile_result AND _run_result EQUAL 0)
                string(STRIP "${_output}" _detected)
            else()
                message(STATUS "Could not detect the SIMD of the processor, select it with ACT_SIMD")
            endif()
        endif()
        set(ACT_SIMD_DETECTED "${_detected}" CACHE INTERNAL "The best explicit SIMD for the processor CMake runs on")
        message(STATUS "Detected the best SIMD for this processor: ${_detected}")
    endif()
    set(${result} "${ACT_SIMD_DETECTED}" PARENT_SCOPE)
endfunction()

# Check whether the compiler can produce code for a SIMD choice
function(act_compiler_supports_simd simd result)
    set(_flags "${_act_simd_flags_${simd}}")
    if(simd STREQUAL "ARM_NEON")
        # NEON is part of every aarch64 processor, no flags are needed
        if(CMAKE_SYSTEM_PROCESSOR MATCHES "^(aarch64|arm64|ARM64)$")
            set(${result} TRUE PARENT_SCOPE)
        else()
            set(${result} FALSE PARENT_SCOPE)
        endif()
    elseif(_flags)
        string(REPLACE ";" " " _flags_string "${_flags}")
        string(MAKE_C_IDENTIFIER "CXX_FLAG_ACT_SIMD_${simd}" _var)
        check_cxx_compiler_flag("${_flags_string}" ${_var})
        set(${result} ${${_var}} PARENT_SCOPE)
    else()
        set(${result} TRUE PARENT_SCOPE)
    endif()
endfunction()

if(ACT_SIMD STREQUAL "AUTO")
    act_detect_simd(_act_detected)
    # Use less SIMD if the compiler cannot produce code for the best one
    if(_act_detected STREQUAL "AVX_512")
        set(_act_choices AVX_512 AVX2)
    else()
        set(_act_choices ${_act_detected})
    endif()
    set(_act_simd "NONE")
    foreach(_simd ${_act_choices})
        act_compiler_supports_simd(${_simd} _supported)
        if(_supported)
            set(_act_simd ${_simd})
            break()
        endif()
        message(STATUS "The compiler cannot produce ${_simd} code")
    endforeach()
elseif(ACT_SIMD MATCHES "^(NONE|ARM_NEON|AVX2|AVX_512)$")
    set(_act_simd ${ACT_SIMD})
    act_compiler_supports_simd(${_act_simd} _supported)
    if(NOT _supported)
        message(FATAL_ERROR "ACT_SIMD is ${ACT_SIMD}, but the compiler cannot produce code for it for ${CMAKE_SYSTEM_PROCESSOR}")
    endif()
else()
    message(FATAL_ERROR "ACT_SIMD should be AUTO, NONE, ARM_NEON, AVX2 or AVX_512, not ${ACT_SIMD}")
endif()

set(ACT_SIMD_FLAGS "${_act_simd_flags_${_act_simd}}")
foreach(_simd ARM_NEON AVX2 AVX_512)
    if(_act_simd STREQUAL _simd)
        set(ACT_SIMD_${_simd} 1)
    else()
        set(ACT_SIMD_${_simd} 0)
    endif()
endforeach()
if(_act_simd STREQUAL "NONE")
    message(STATUS "No explicit SIMD for the force routines")
else()
    message(STATUS "Using ${_act_simd} for the explicit SIMD force routines")
endif()

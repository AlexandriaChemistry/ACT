#
# This file is part of the GROMACS molecular simulation package.
#
# Copyright (c) 2012,2013,2014,2015,2016,2017, by the GROMACS development team, led by
# Mark Abraham, David van der Spoel, Berk Hess, and Erik Lindahl,
# and including many others, as listed in the AUTHORS file in the
# top-level source directory and at http://www.gromacs.org.
#
# GROMACS is free software; you can redistribute it and/or
# modify it under the terms of the GNU Lesser General Public License
# as published by the Free Software Foundation; either version 2.1
# of the License, or (at your option) any later version.
#
# GROMACS is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
# Lesser General Public License for more details.
#
# You should have received a copy of the GNU Lesser General Public
# License along with GROMACS; if not, see
# http://www.gnu.org/licenses, or write to the Free Software Foundation,
# Inc., 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301  USA.
#
# If you want to redistribute modifications to GROMACS, please
# consider that scientific software is very special. Version
# control is crucial - bugs must be traceable. We will be happy to
# consider code for inclusion in the official distribution, but
# derived work must not be called official GROMACS. Details are found
# in the README & COPYING files - if they are missing, get the
# official version at http://www.gromacs.org.
#
# To help us fund GROMACS development, we humbly ask that you cite
# the research papers on the package. Check out http://www.gromacs.org.

# Manage the OpenMP setup. This wrapper file checks for some known-bad compiler
# versions before trying to detect OpenMP with the standard find-package-module,
# and then does some additional tests for flags afterwards.

if(GMX_OPENMP)
    # We should do OpenMP detection if we get here
    # OpenMP check must come before other CFLAGS!
    find_package(OpenMP)
    if(OPENMP_FOUND)
        # CMake on Windows doesn't support linker flags passed to target_link_libraries
        # (i.e. it treats /openmp as \openmp library file). Also, no OpenMP linker flags are needed.
        if(NOT (WIN32 AND NOT MINGW))
            if(CMAKE_COMPILER_IS_GNUCC AND GMX_PREFER_STATIC_OPENMP AND NOT APPLE)
                set(OpenMP_LINKER_FLAGS "-Wl,-static -lgomp -lrt -Wl,-Bdynamic -lpthread")
                set(OpenMP_SHARED_LINKER_FLAGS "")
            else()
                # Only set a linker flag if the user didn't set them manually
                if(NOT DEFINED OpenMP_LINKER_FLAGS)
                    set(OpenMP_LINKER_FLAGS "${OpenMP_C_FLAGS}")
                endif()
                if(NOT DEFINED OpenMP_SHARED_LINKER_FLAGS)
                    set(OpenMP_SHARED_LINKER_FLAGS "${OpenMP_C_FLAGS}")
                endif()
            endif()
        endif()
        if(MINGW)
            #GCC Bug 48659
            set(OpenMP_C_FLAGS "${OpenMP_C_FLAGS} -mstackrealign")
        endif()
    else()
        message(WARNING
                "The compiler you are using does not support OpenMP parallelism. This might hurt your performance a lot, in particular with GPUs. Try using a more recent version, or a different compiler. For now, we are proceeding by turning off OpenMP.")
        set(GMX_OPENMP OFF CACHE STRING "Whether GROMACS will use OpenMP parallelism." FORCE)
    endif()
endif()
# Compiler vectorization of loops marked with "#pragma omp simd".
# This only enables the OpenMP SIMD directives, not OpenMP multithreading,
# so it can be used independently of GMX_OPENMP. The flags are collected in
# GMX_SIMD_VECTORIZATION_C_FLAGS, GMX_SIMD_VECTORIZATION_CXX_FLAGS and
# GMX_SIMD_VECTORIZATION_LINKER_FLAGS and added to the compiler flags in the
# main CMakeLists.txt.
option(GMX_OPENMP_SIMD "Enable OpenMP SIMD directives for compiler vectorization of loops" ON)
option(GMX_VECTOR_MATH "Use a vector math library for math functions (exp, erf, ...) in vectorized loops" ON)
mark_as_advanced(GMX_OPENMP_SIMD GMX_VECTOR_MATH)
set(GMX_SIMD_VECTORIZATION_C_FLAGS "")
set(GMX_SIMD_VECTORIZATION_CXX_FLAGS "")
set(GMX_SIMD_VECTORIZATION_LINKER_FLAGS "")
if(GMX_OPENMP_SIMD)
    include(CheckCXXCompilerFlag)
    include(CheckCXXSourceCompiles)
    # Flag for the OpenMP SIMD directives
    set(_omp_simd_flag "")
    foreach(_flag -fopenmp-simd -qopenmp-simd)
        string(MAKE_C_IDENTIFIER "CXX_FLAG${_flag}" _var)
        check_cxx_compiler_flag(${_flag} ${_var})
        if(${_var})
            set(_omp_simd_flag ${_flag})
            break()
        endif()
    endforeach()
    if(_omp_simd_flag)
        set(GMX_SIMD_VECTORIZATION_C_FLAGS   "${_omp_simd_flag}")
        set(GMX_SIMD_VECTORIZATION_CXX_FLAGS "${_omp_simd_flag}")
        # A loop marked for vectorization may still not be vectorized, for
        # instance when it calls math functions without a vector math library.
        # The compiler then warns, which should not stop the build with -Werror.
        check_cxx_compiler_flag(-Wno-error=pass-failed CXX_FLAG_WNO_ERROR_PASS_FAILED)
        if(CXX_FLAG_WNO_ERROR_PASS_FAILED)
            set(GMX_SIMD_VECTORIZATION_CXX_FLAGS "${GMX_SIMD_VECTORIZATION_CXX_FLAGS} -Wno-error=pass-failed")
        endif()
        if(GMX_VECTOR_MATH)
            # Vector math library, the first one that compiles and links a
            # vectorized loop calling exp and erf is used.
            if(APPLE)
                set(_veclibs "Darwin_libsystem_m")
            else()
                set(_veclibs "libmvec;SVML")
            endif()
            set(_vector_math_found FALSE)
            foreach(_veclib ${_veclibs})
                set(_veclib_libs "")
                if(_veclib STREQUAL "libmvec")
                    set(_veclib_libs "-lmvec")
                endif()
                set(CMAKE_REQUIRED_FLAGS "${_omp_simd_flag} -fveclib=${_veclib}")
                set(CMAKE_REQUIRED_LIBRARIES "${_veclib_libs}")
                string(MAKE_C_IDENTIFIER "CXX_VECLIB_${_veclib}" _var)
                check_cxx_source_compiles("
#include <cmath>
int main(int argc, char *argv[])
{
    double x[64], y[64];
    for (int i = 0; i < 64; i++) { x[i] = argc*0.01*i; }
#pragma omp simd
    for (int i = 0; i < 64; i++) { y[i] = std::exp(-x[i]) + std::erf(x[i]); }
    return y[argc] > 0 ? 0 : 1;
}" ${_var})
                unset(CMAKE_REQUIRED_FLAGS)
                unset(CMAKE_REQUIRED_LIBRARIES)
                if(${_var})
                    set(GMX_SIMD_VECTORIZATION_C_FLAGS   "${GMX_SIMD_VECTORIZATION_C_FLAGS} -fveclib=${_veclib}")
                    set(GMX_SIMD_VECTORIZATION_CXX_FLAGS "${GMX_SIMD_VECTORIZATION_CXX_FLAGS} -fveclib=${_veclib}")
                    set(GMX_SIMD_VECTORIZATION_LINKER_FLAGS "${_veclib_libs}")
                    set(_vector_math_found TRUE)
                    message(STATUS "Using vector math library ${_veclib} for vectorized loops")
                    break()
                endif()
            endforeach()
            if(NOT _vector_math_found)
                message(STATUS "No vector math library found, loops calling math functions will not be vectorized")
            endif()
        endif()
        message(STATUS "Flags for compiler vectorization: ${GMX_SIMD_VECTORIZATION_CXX_FLAGS}")
    else()
        message(STATUS "The compiler does not support OpenMP SIMD directives, loops marked with #pragma omp simd may not be vectorized")
    endif()
endif()

gmx_dependent_cache_variable(GMX_OPENMP_MAX_THREADS
    "Maximum number of OpenMP Threads supported. Has to be 32 or a multiple of 64."
    STRING 64 GMX_OPENMP)
mark_as_advanced(GMX_OPENMP_MAX_THREADS)
math(EXPR MAX_THREAD_MOD "${GMX_OPENMP_MAX_THREADS} % 64")
if (NOT GMX_OPENMP_MAX_THREADS EQUAL 32 AND NOT ${MAX_THREAD_MOD} EQUAL 0)
    message(FATAL_ERROR "Only 32 or multiples of 64 supported for GMX_OPENMP_MAX_THREADS.")
endif()

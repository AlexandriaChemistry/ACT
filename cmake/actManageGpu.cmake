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

# GPU acceleration of the force routines
#
# The option ACT_GPU selects the programming model for GPUs: NONE (the
# default) or SYCL. SYCL needs AdaptiveCpp (https://adaptivecpp.github.io),
# for instance from conda-forge, found with CMAKE_PREFIX_PATH or
# AdaptiveCpp_DIR. Its generic target compiles the kernels for the GPU at
# run time, such that one build runs on NVIDIA and AMD GPUs and, without a
# GPU, on the processor. Only the library actGpu in src/act/gpu/sycl is
# compiled with acpp, the rest of ACT with the C++ compiler, and the
# interface between them has no SYCL types. The libraries of AdaptiveCpp
# need the C++ library of GCC 13 or newer at run time.
#
# Sets
#   ACT_GPU_SYCL        1 if SYCL is used, otherwise 0, for config.h
#   ACT_GPU_DESCRIPTION The GPU support for the log file, for config.h

set(ACT_GPU "NONE" CACHE STRING "GPU acceleration of the force routines: NONE or SYCL (with AdaptiveCpp)")
set_property(CACHE ACT_GPU PROPERTY STRINGS NONE SYCL)

set(ACT_GPU_SYCL 0)
set(ACT_GPU_DESCRIPTION "disabled")
if(ACT_GPU STREQUAL "SYCL")
    find_package(AdaptiveCpp CONFIG REQUIRED)
    set(ACT_GPU_SYCL 1)
    # The backends of AdaptiveCpp are plugins next to its runtime library,
    # for instance conda-forge has variants without the CUDA backend, which
    # cannot use NVIDIA GPUs
    get_filename_component(_acpp_libdir "${AdaptiveCpp_DIR}/../.." ABSOLUTE)
    file(GLOB _acpp_backend_libs "${_acpp_libdir}/hipSYCL/librt-backend-*")
    set(_acpp_backends "")
    foreach(_lib ${_acpp_backend_libs})
        get_filename_component(_name "${_lib}" NAME_WE)
        string(REPLACE "librt-backend-" "" _name "${_name}")
        list(APPEND _acpp_backends ${_name})
    endforeach()
    list(SORT _acpp_backends)
    list(JOIN _acpp_backends " " _acpp_backend_list)
    if(_acpp_backends)
        set(ACT_GPU_DESCRIPTION "SYCL (AdaptiveCpp ${AdaptiveCpp_VERSION}, backends ${_acpp_backend_list})")
    else()
        set(ACT_GPU_DESCRIPTION "SYCL (AdaptiveCpp ${AdaptiveCpp_VERSION})")
    endif()
    message(STATUS "Using SYCL with AdaptiveCpp ${AdaptiveCpp_VERSION} (${ACPP_COMPILER}) for the GPU force routines, backends: ${_acpp_backend_list}")
    if(_acpp_backends AND NOT "cuda" IN_LIST _acpp_backends AND NOT "hip" IN_LIST _acpp_backends)
        message(WARNING "AdaptiveCpp in ${_acpp_libdir} has neither the CUDA nor the HIP backend, "
                        "such that ACT cannot use NVIDIA or AMD GPUs. Install a variant with a GPU "
                        "backend, for instance from conda-forge \"adaptivecpp=*=cuda*\" for NVIDIA GPUs, "
                        "and give its AdaptiveCpp_DIR.")
    elseif(_acpp_backends AND NOT "cuda" IN_LIST _acpp_backends)
        message(WARNING "AdaptiveCpp in ${_acpp_libdir} has no CUDA backend, such that ACT cannot use "
                        "NVIDIA GPUs, only AMD GPUs with its HIP backend. For NVIDIA GPUs, install a "
                        "variant with it, for instance from conda-forge \"adaptivecpp=*=cuda*\", and "
                        "give its AdaptiveCpp_DIR.")
    endif()
elseif(NOT ACT_GPU STREQUAL "NONE")
    message(FATAL_ERROR "ACT_GPU should be NONE or SYCL, not ${ACT_GPU}")
endif()

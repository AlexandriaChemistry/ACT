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
    set(ACT_GPU_DESCRIPTION "SYCL (AdaptiveCpp ${AdaptiveCpp_VERSION})")
    message(STATUS "Using SYCL with AdaptiveCpp ${AdaptiveCpp_VERSION} (${ACPP_COMPILER}) for the GPU force routines")
elseif(NOT ACT_GPU STREQUAL "NONE")
    message(FATAL_ERROR "ACT_GPU should be NONE or SYCL, not ${ACT_GPU}")
endif()

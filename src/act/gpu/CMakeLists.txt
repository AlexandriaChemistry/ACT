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

# The interface in gpu.h is implemented with SYCL in the library actGpu,
# which is compiled with acpp, see cmake/actManageGpu.cmake, or without a
# GPU in nogpu.cpp.
if(ACT_GPU_SYCL)
    # acpp is clang, which does not know some of the flags for GCC that
    # cmake/gmxCFlags.cmake adds, therefore the targets of this directory
    # only have the flags of the cache, those of the build type included
    set(CMAKE_CXX_FLAGS "$CACHE{CMAKE_CXX_FLAGS} -Wall -Wextra")
    string(TOUPPER "${CMAKE_BUILD_TYPE}" _act_build_type)
    if(_act_build_type)
        set(CMAKE_CXX_FLAGS_${_act_build_type} "$CACHE{CMAKE_CXX_FLAGS_${_act_build_type}}")
    endif()
    set(ACT_GPU_SYCL_SOURCES sycl/coulomb.cpp sycl/devices.cpp)
    add_library(actGpu SHARED ${ACT_GPU_SYCL_SOURCES})
    target_include_directories(actGpu PRIVATE ${CMAKE_CURRENT_SOURCE_DIR})
    target_compile_definitions(actGpu PRIVATE ACT_ACPP_VERSION="${AdaptiveCpp_VERSION}")
    add_sycl_to_target(TARGET actGpu SOURCES ${ACT_GPU_SYCL_SOURCES})
    target_link_libraries(actLib PUBLIC actGpu)
    install(TARGETS actGpu LIBRARY DESTINATION ${CMAKE_INSTALL_LIBDIR})
else()
    target_sources(actLib PRIVATE nogpu.cpp)
endif()
target_sources(actLib PRIVATE gpu.h)

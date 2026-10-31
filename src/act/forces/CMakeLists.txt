#
#
#  This source file is part of the Alexandria program.
#
#  Copyright (C) 2022
#
#  Developers:
#             Mohammad Mehdi Ghahremanpour,
#             Paul J. van Maaren,
#             David van der Spoel (Project leader)
#
#  This program is free software; you can redistribute it and/or
#  modify it under the terms of the GNU General Public License
#  as published by the Free Software Foundation; either version 2
#  of the License, or (at your option) any later version.
#
#  This program is distributed in the hope that it will be useful,
#  but WITHOUT ANY WARRANTY; without even the implied warranty of
#  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#  GNU General Public License for more details.
#
#  You should have received a copy of the GNU General Public License
#  along with this program; if not, write to the Free Software
#  Foundation, Inc., 51 Franklin Street, Fifth Floor,
#  Boston, MA  02110-1301, USA.
#

file(GLOB UTILITY_SOURCES *.cpp *.h)
target_sources(actLib PRIVATE ${UTILITY_SOURCES})
if(ACT_SIMD_FLAGS)
    set_source_files_properties(forces_explicit_simd.cpp TARGET_DIRECTORY actLib
                                PROPERTIES COMPILE_OPTIONS "${ACT_SIMD_FLAGS}")
endif()
# For benchmarking. The value is in a generated header that only the file
# with the explicit SIMD includes, such that only that file is compiled again
# when it changes.
set(ACT_PAIR_UNROLL "" CACHE STRING
    "Unroll factor of the loops over conformations of pair interactions with explicit SIMD, 1 to 16, empty for the default of the instruction set")
if(ACT_PAIR_UNROLL AND NOT ACT_PAIR_UNROLL MATCHES "^([1-9]|1[0-6])$")
    message(FATAL_ERROR "ACT_PAIR_UNROLL should be 1 to 16, not ${ACT_PAIR_UNROLL}")
endif()
configure_file(pairunroll.h.cmakein ${CMAKE_BINARY_DIR}/src/pairunroll.h)
	
if (BUILD_TESTING)
    add_subdirectory(tests)
endif()


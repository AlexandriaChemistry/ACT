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

# The suffix of the names of the programs and libraries. A double precision
# build gets the suffix _d, such that it can be installed next to a single
# precision build, in the same place: alexandria and alexandria_d, and
# libactgromacs and libactgromacs_d. With ACT_DEFAULT_SUFFIX=OFF the suffixes
# are ACT_BINARY_SUFFIX and ACT_LIBS_SUFFIX, which may then be set freely.

option(ACT_DEFAULT_SUFFIX "Add _d to the names of the programs and libraries of a double precision build" ON)
mark_as_advanced(ACT_DEFAULT_SUFFIX)

if(ACT_DEFAULT_SUFFIX)
    if(ACT_DOUBLE)
        set(_act_suffix "_d")
    else()
        set(_act_suffix "")
    endif()
    # Follow ACT_DOUBLE, also when it changes in an existing build directory
    set(ACT_BINARY_SUFFIX "${_act_suffix}" CACHE STRING "Suffix of the names of the programs" FORCE)
    set(ACT_LIBS_SUFFIX   "${_act_suffix}" CACHE STRING "Suffix of the names of the libraries" FORCE)
    unset(_act_suffix)
else()
    set(ACT_BINARY_SUFFIX "" CACHE STRING "Suffix of the names of the programs")
    set(ACT_LIBS_SUFFIX   "" CACHE STRING "Suffix of the names of the libraries")
endif()
mark_as_advanced(ACT_BINARY_SUFFIX ACT_LIBS_SUFFIX)

if(ACT_BINARY_SUFFIX OR ACT_LIBS_SUFFIX)
    message(STATUS "The programs get the suffix '${ACT_BINARY_SUFFIX}' and the libraries '${ACT_LIBS_SUFFIX}'")
endif()

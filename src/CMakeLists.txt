#
# This file is part of the GROMACS molecular simulation package.
#
# Copyright (c) 2009,2010,2011,2012,2013,2014,2015,2016,2017,2018, by the GROMACS development team, led by
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

######################################
# Output compiler and CFLAGS used
######################################
include(GetCompilerInfo.cmake)
get_compiler_info(C BUILD_C_COMPILER BUILD_CFLAGS)
get_compiler_info(CXX BUILD_CXX_COMPILER BUILD_CXXFLAGS)
if(GMX_USE_CUDA)
    if(NOT GMX_CLANG_CUDA)
        GMX_SET_CUDA_NVCC_FLAGS()
    endif()

    get_cuda_compiler_info(CUDA_COMPILER_INFO CUDA_COMPILER_FLAGS)
endif()

string(TOUPPER "${CMAKE_BUILD_TYPE}" CMAKE_BUILD_TYPE_UPPER)
configure_file(config.h.cmakein config.h)
configure_file(actpre-config.h.cmakein actpre-config.h)
configure_file(buildinfo.h.cmakein buildinfo.h ESCAPE_QUOTES)

####
set(IGNORED_CLANG_ALL_WARNINGS
    "-Werror"
    "-Wno-c++98-compat -Wno-c++98-compat-pedantic" #No intention of C++98 compability
    "-Wno-source-uses-openmp" #Don't warn for no-omp build
    "-Wno-documentation-unknown-command" #Custom commands are used
    "-Wno-covered-switch-default" #GCC gives maybe-uninitialized without default label and checks for illegal enum values.
    "-Wno-switch-enum" # default statement for enum is OK

    #Following ones are undecided/TODO
    "-Wno-disabled-macro-expansion"
    "-Wno-deprecated-declarations"
    "-Wno-deprecated-copy-with-dtor"
    "-Wno-cast-align"
    "-Wno-reserved-id-macro"
    "-Wno-global-constructors"
    "-Wno-exit-time-destructors"
    "-Wno-unused-macros"
    "-Wno-weak-vtables"
    "-Wno-conditional-uninitialized"
    "-Wno-format-nonliteral"
    "-Wno-shadow"
    "-Wno-cast-qual"
    "-Wno-documentation"
    "-Wno-used-but-marked-unused"
    "-Wno-padded"
    "-Wno-float-equal"
    "-Wno-old-style-cast"
    "-Wno-conversion"
    "-Wno-double-promotion")
string(REPLACE " " ";" IGNORED_CLANG_ALL_WARNINGS "${IGNORED_CLANG_ALL_WARNINGS}")

if (GMX_CLANG_TIDY)
   set(CLANG_TIDY "clang-tidy" CACHE STRING "Name of clang-tidy executable")
   find_program(CLANG_TIDY_EXE NAMES "${CLANG_TIDY}"
       DOC "Path to clang-tidy executable")
   if(NOT CLANG_TIDY_EXE)
       message(FATAL_ERROR "clang-tidy not found.")
   endif()
   if(CMAKE_VERSION VERSION_LESS "3.6")
       message(FATAL_ERROR "clang-tidy support requires cmake 3.6.")
   endif()
endif()
#####

if (BUILD_TESTING)
    if(NOT GMX_DEVELOPER_BUILD)
        set(UNITTEST_TARGET_OPTIONS EXCLUDE_FROM_ALL)
    endif()
    if (GMX_BUILD_UNITTESTS)
        add_subdirectory(external/googletest)
    endif()
    include(testutils/TestMacros.cmake)
    add_subdirectory(testutils)
endif()
add_subdirectory(gromacs)

# The ACT library file
add_library(actLib "")
# And the source code
add_subdirectory(act)

# External code that is included in the ACT library
add_subdirectory(external/quasirandom_sequences)
add_subdirectory(external/stlbfgs)



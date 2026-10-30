#
# This file is part of the GROMACS molecular simulation package.
#
# Copyright (c) 2010,2011,2012,2013,2014,2015,2016,2017,2018, by the GROMACS development team, led by
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

set(LIBGROMACS_SOURCES)

set_property(GLOBAL PROPERTY GMX_LIBGROMACS_SOURCES)
set_property(GLOBAL PROPERTY GMX_LIBGROMACS_GPU_IMPL_SOURCES)
set_property(GLOBAL PROPERTY GMX_AVX_512_SOURCE)

function (_gmx_add_files_to_property PROPERTY)
    foreach (_file ${ARGN})
        if (IS_ABSOLUTE "${_file}")
            set_property(GLOBAL APPEND PROPERTY ${PROPERTY} ${_file})
        else()
            set_property(GLOBAL APPEND PROPERTY ${PROPERTY}
                         ${CMAKE_CURRENT_LIST_DIR}/${_file})
        endif()
    endforeach()
endfunction ()

function (gmx_add_libgromacs_sources)
    _gmx_add_files_to_property(GMX_LIBGROMACS_SOURCES ${ARGN})
endfunction ()

add_subdirectory(gmxlib)
add_subdirectory(commandline)
add_subdirectory(math)
add_subdirectory(mdtypes)
add_subdirectory(onlinehelp)
add_subdirectory(options)
add_subdirectory(pbcutil)
add_subdirectory(random)
add_subdirectory(topology)
add_subdirectory(utility)
add_subdirectory(fileio)
add_subdirectory(compat)
add_subdirectory(gmxpreprocess)

get_property(PROPERTY_SOURCES GLOBAL PROPERTY GMX_LIBGROMACS_SOURCES)
list(APPEND LIBGROMACS_SOURCES ${GMXLIB_SOURCES} ${MDLIB_SOURCES} ${PROPERTY_SOURCES})

configure_file(version.h.cmakein version.h)

# This code is here instead of utility/CMakeLists.txt, because CMake
# custom commands and source file properties can only be set in the directory
# that contains the target that uses them.
# TODO: Generate a header instead that can be included from baseversion.c.
# That probably simplifies things somewhat.
set(GENERATED_VERSION_FILE utility/baseversion-gen.cpp)
gmx_configure_version_file(
    utility/baseversion-gen.cpp.cmakein ${GENERATED_VERSION_FILE}
    REMOTE_HASH
    EXTRA_VARS
        ACT_SOURCE_DOI
    )
list(APPEND LIBGROMACS_SOURCES ${GENERATED_VERSION_FILE}
)

add_library(libactgromacs ${LIBGROMACS_SOURCES})

# MPI is public, since the headers include mpi.h
target_link_libraries(libactgromacs
                      PRIVATE
                      ${GMX_EXTRA_LIBRARIES}
                      ${GMX_SHARED_LINKER_FLAGS}
                      ${GMX_STDLIB_LIBRARIES}
                      PUBLIC
                      ${GMX_COMMON_LIBRARIES}
                      ${GMX_PUBLIC_LIBRARIES}
                      )
set_target_properties(libactgromacs PROPERTIES
                      OUTPUT_NAME "actgromacs${GMX_LIBS_SUFFIX}"
                      SOVERSION ${LIBRARY_SOVERSION_MAJOR}
                      VERSION ${LIBRARY_VERSION}
                      COMPILE_FLAGS "${OpenMP_C_FLAGS}")

if (GMX_CLANG_TIDY)
   set_target_properties(libactgromacs PROPERTIES CXX_CLANG_TIDY
       "${CLANG_TIDY_EXE};-warnings-as-errors=*")
endif()

install(TARGETS libactgromacs
        EXPORT libactgromacs
        LIBRARY DESTINATION ${CMAKE_INSTALL_LIBDIR}
        RUNTIME DESTINATION ${CMAKE_INSTALL_BINDIR}
        ARCHIVE DESTINATION ${CMAKE_INSTALL_LIBDIR}
        COMPONENT libraries)

include(InstallLibInfo.cmake)


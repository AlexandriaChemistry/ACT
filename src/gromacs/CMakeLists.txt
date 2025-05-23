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
add_subdirectory(mdlib)
add_subdirectory(applied-forces)
add_subdirectory(listed-forces)
add_subdirectory(commandline)
add_subdirectory(gpu_utils)
add_subdirectory(math)
add_subdirectory(mdrunutility)
add_subdirectory(mdtypes)
add_subdirectory(onlinehelp)
add_subdirectory(options)
add_subdirectory(pbcutil)
add_subdirectory(random)
add_subdirectory(tables)
add_subdirectory(timing)
add_subdirectory(topology)
add_subdirectory(utility)
add_subdirectory(fileio)
add_subdirectory(compat)
add_subdirectory(gmxpreprocess)
add_subdirectory(statistics)
#add_subdirectory(selection)

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

# Recent versions of gcc and clang give warnings on scanner.cpp, which
# is a generated source file. These are awkward to suppress inline, so
# we do it in the compilation command (after testing that the compiler
# supports the suppressions). Same issue exists for nonbonded kernels
# so we supress them for all generated files.
include(CheckCXXCompilerFlag)
check_cxx_compiler_flag("-Wno-unused -Wno-unused-parameter" HAS_NO_UNUSED)
check_cxx_compiler_flag(-Wno-missing-declarations HAS_NO_MISSING_DECL)
check_cxx_compiler_flag(-Wno-missing-prototypes HAS_NO_MISSING_PROTO)
check_cxx_compiler_flag(/wd4101 HAS_NO_MSVC_UNUSED)
if (NOT MSVC)
    check_cxx_compiler_flag(-wd1419 HAS_DECL_IN_SOURCE)
endif()
# TODO The group scheme kernels don't use proper function prototype
# declarations, and clang warns about such use, which we suppress
# rather than fix. We would prefer to use no suppressions. However
# other compilers do not support such a warning suppression for C++
# source files, and issue warnings about that. Remove the use of
# -Wno-missing-prototypes here and above when the group scheme is
# removed.

if(SIMD_AVX_512_CXX_SUPPORTED AND NOT ("${GMX_SIMD_ACTIVE}" STREQUAL "AVX_512_KNL"))
    # Since we might be overriding -march=core-avx2, add a flag so we don't warn for this specific file.
    # On KNL this can cause illegal instruction because the compiler might use non KNL AVX instructions
    # with the SIMD_AVX_512_CXX_FLAGS flags.
    set_source_files_properties(hardware/identifyavx512fmaunits.cpp PROPERTIES COMPILE_FLAGS "${SIMD_AVX_512_CXX_FLAGS} ${CXX_NO_UNUSED_OPTION_WARNING_FLAGS}")
endif()

target_link_libraries(libactgromacs
                      PRIVATE
                      ${EXTRAE_LIBRARIES}
                      ${GMX_EXTRA_LIBRARIES}
                      ${GMX_COMMON_LIBRARIES}
                      ${FFT_LIBRARIES} ${LINEAR_ALGEBRA_LIBRARIES}
                      ${THREAD_LIB} ${GMX_SHARED_LINKER_FLAGS}
                      ${OpenCL_LIBRARIES}
                      ${GMX_STDLIB_LIBRARIES}
                      PUBLIC
                      ${GMX_PUBLIC_LIBRARIES}
                      )
set_target_properties(libactgromacs PROPERTIES
                      OUTPUT_NAME "actgromacs${GMX_LIBS_SUFFIX}"
                      SOVERSION ${LIBRARY_SOVERSION_MAJOR}
                      VERSION ${LIBRARY_VERSION}
                      COMPILE_FLAGS "${OpenMP_C_FLAGS}")

if (CMAKE_CXX_COMPILER_ID MATCHES "Clang" AND CMAKE_CXX_COMPILER_VERSION MATCHES "^6\.0")
   target_compile_options(libactgromacs PRIVATE $<$<COMPILE_LANGUAGE:CXX>:-Weverything ${IGNORED_CLANG_ALL_WARNINGS}>)
endif()
if (CMAKE_CXX_COMPILER_ID STREQUAL "MSVC")
   target_compile_options(libactgromacs PRIVATE $<$<COMPILE_LANGUAGE:CXX>:/analyze /analyze:stacksize 70000
     #Control flow warnings are disabled because the commond line output is insufficient. There is no tool
     #to convert the xml report to e.g. HTML and even in Visual Studio the viewer doesn't work with cmake support.
     /wd6001  #unitialized memory
     /wd6011  #derefencing NULL
     /wd6053  #prior call not zero-terminate
     /wd6054  #might not be zero-terminated
     /wd6385  #reading invalid data
     /wd6386  #buffer overrun
     /wd6387  #could be '0'
     /wd28199 #uninitialized memory
     # For compile time constant (e.g. templates) the following warnings have flase postives
     /wd6239  #(<non-zero> && <expr>)
     /wd6240  #(<expr> && <non-zero>)
     /wd6294  #Ill-defined for-loop
     /wd6326  #comparison of constant with other constant
     /wd28020 #expression involving paramter is not true
     # Misc
     /wd6330  #incorrect type to function (warns for char (instead of unsigned) for isspace/isalpha/isdigit/..))
     /wd6993  #OpenMP ignored
     #TODO
     /wd6031  #return value ignored (important - mostly warnigns about sscanf)
     /wd6244  #hides declaration (known issue - we ingore similar warnings for other compilers)
     /wd6246  #hides declaration
     >
   )
endif()

if (GMX_CLANG_TIDY)
   set_target_properties(libactgromacs PROPERTIES CXX_CLANG_TIDY
       "${CLANG_TIDY_EXE};-warnings-as-errors=*")
endif()

# Only install the library in mdrun-only mode if it is actually necessary
# for the binary
if (NOT GMX_BUILD_MDRUN_ONLY OR BUILD_SHARED_LIBS)
    install(TARGETS libactgromacs
            EXPORT libactgromacs
            LIBRARY DESTINATION ${CMAKE_INSTALL_LIBDIR}
            RUNTIME DESTINATION ${CMAKE_INSTALL_BINDIR}
            ARCHIVE DESTINATION ${CMAKE_INSTALL_LIBDIR}
            COMPONENT libraries)
endif()

if (NOT GMX_BUILD_MDRUN_ONLY)
    include(InstallLibInfo.cmake)
endif()


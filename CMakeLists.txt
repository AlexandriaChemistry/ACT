#
# Main file for the Alexandria Chemistry Toolkit which is based
# on a fork off of GROMACS.
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

# 3.20 for cmake_path and 3.24 for CMAKE_COMPILE_WARNING_AS_ERROR. With
# newer versions the policies up to 3.31 are used, such that, e.g.,
# duplicate libraries are removed from the link lines (CMP0156).
cmake_minimum_required(VERSION 3.24...3.31)

# CMake modules/macros are in a subdirectory to keep this file cleaner
# This needs to be set before project() in order to pick up toolchain files
list(APPEND CMAKE_MODULE_PATH ${CMAKE_CURRENT_SOURCE_DIR}/cmake ${CMAKE_CURRENT_SOURCE_DIR}/cmake/Platform)

project(ACT)

include(gmxManageCcache)

set(CMAKE_LIBRARY_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/lib)
set(CMAKE_ARCHIVE_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/lib)
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/bin)

find_package (Python3 COMPONENTS Interpreter Development)
if (NOT CMAKE_INSTALL_PYTHON_LIBDIR)
  set(CMAKE_INSTALL_PYTHON_LIBDIR
    "lib/python${Python3_VERSION_MAJOR}.${Python3_VERSION_MINOR}/site-packages")
endif ()

# Set up common version variables, as well as general information about
# the build tree (whether the build is from a source package or from a git
# repository).  Also declares a few functions that will be used for generating
# version info files later.
include(gmxBuildTreeInfo)
include(actVersionInfo)

if(CMAKE_INSTALL_PREFIX_INITIALIZED_TO_DEFAULT AND UNIX)
    set(CMAKE_INSTALL_PREFIX "/usr/local/ACT" CACHE STRING "Installation prefix (installation will need write permissions here)" FORCE)
endif()
if("${CMAKE_INSTALL_PREFIX}" STREQUAL "${CMAKE_BINARY_DIR}")
    message(FATAL_ERROR "ACT cannot be installed into the build tree, choose a different location for CMAKE_INSTALL_PREFIX")
endif()

include(gmxBuildTypeReference)
include(gmxBuildTypeProfile)
include(gmxBuildTypeTSAN)
include(gmxBuildTypeASAN)
include(gmxBuildTypeMSAN)
include(gmxBuildTypeReleaseWithAssert)

if(NOT CMAKE_BUILD_TYPE)
    set(CMAKE_BUILD_TYPE "Release" CACHE STRING "Choose the type of build, options are: Debug Release RelWithDebInfo MinSizeRel Reference RelWithAssert Profile TSAN ASAN MSAN." FORCE)
    # Set the possible values of build type for cmake-gui
    set_property(CACHE CMAKE_BUILD_TYPE PROPERTY STRINGS "Debug" "Release"
        "MinSizeRel" "RelWithDebInfo" "Reference" "RelWithAssert" "Profile" "TSAN" "ASAN" "MSAN")
endif()
if(CMAKE_CONFIGURATION_TYPES)
    # Add appropriate ACT-specific build types for the Visual
    # Studio generator (Debug, Release, MinSizeRel and RelWithDebInfo
    # are already present by default).
    list(APPEND CMAKE_CONFIGURATION_TYPES "RelWithAssert" "Reference")
    list(REMOVE_DUPLICATES CMAKE_CONFIGURATION_TYPES)
    set(CMAKE_CONFIGURATION_TYPES "${CMAKE_CONFIGURATION_TYPES}" CACHE STRING
        "List of configuration types"
        FORCE)
endif()
set(build_types_with_explicit_flags RELEASE DEBUG RELWITHDEBINFO RELWITHASSERT MINSIZEREL PROFILE)

set_property(GLOBAL PROPERTY FIND_LIBRARY_USE_LIB64_PATHS ON)

include(gmxCTestUtilities)
gmx_ctest_init()

include(gmxCPackUtilities)
gmx_cpack_init()

# Variables that accumulate stuff influencing the installed headers
set(INSTALLED_HEADER_INCLUDE_DIRS "")
set(INSTALLED_HEADER_DEFINITIONS "")

########################################################################
# Global non-cache variables for implementing the build system
########################################################################

# These variables collect libraries that ACT requires for
# linking. They should be appended to with list(APPEND ${name}
# new-library) calls. They are:
#  - Libraries that are required for libactgromacs (only)
set(GMX_EXTRA_LIBRARIES "")
#  - Libraries that are required for all code in the repository
set(GMX_COMMON_LIBRARIES "")

########################################################################
# Check and warn if cache generated on a different host is being reused
########################################################################
if(CMAKE_HOST_UNIX)
    execute_process(COMMAND hostname
                    OUTPUT_VARIABLE TMP_HOSTNAME
                    OUTPUT_STRIP_TRAILING_WHITESPACE)
    if(GMX_BUILD_HOSTNAME AND NOT "${GMX_BUILD_HOSTNAME}" STREQUAL "${TMP_HOSTNAME}")
        message(WARNING "
            The CMake cache, probably generated on a different host (${GMX_BUILD_HOSTNAME}),
            is being reused! This could lead to inconsistencies; therefore, it is
            recommended to regenerate the cache!")
    endif()
    set(GMX_BUILD_HOSTNAME "${TMP_HOSTNAME}" CACHE INTERNAL
            "Hostname of the machine where the cache was generated.")
endif()

########################################################################
# Detect architecture before setting options so we can alter defaults
########################################################################
# Detect the architecture the compiler is targetting, detect
# SIMD instructions possibilities on that hardware, suggest SIMD instruction set
# to use if none is specified, and populate the cache option for CPU
# SIMD.
include(gmxDetectTargetArchitecture)
gmx_detect_target_architecture()

# Permit the user to specify a particular standard library, e.g. compiling
# with "-stdlib=libc++" and linking with "-lc++abi -lc++" to get clang's libcxx.
set(GMX_STDLIB_CXX_FLAGS "" CACHE STRING "Compiler flag for a C++ standard library flavour")
set(GMX_STDLIB_LIBRARIES "" CACHE STRING "Linker libraries for a particular C++ standard library")
mark_as_advanced(GMX_STDLIB_CXX_FLAGS)
mark_as_advanced(GMX_STDLIB_LIBRARIES)

########################################################################
# Detect CXX20 support and flags
########################################################################
include(CheckCCompilerFlag)
include(CheckCXXCompilerFlag)

# All C++ code is compiled as C++20, which also applies to the
# configuration tests, with the expected standard library.
set(CMAKE_CXX_STANDARD 20)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${GMX_STDLIB_CXX_FLAGS}")

########################################################################
# User input options                                                   #
########################################################################
include(gmxOptionUtilities)
set(CMAKE_COMPILE_WARNING_AS_ERROR ON)
set(CMAKE_PREFIX_PATH "" CACHE STRING "Extra locations to search for external libraries and tools (give directory without lib, bin, or include)")

gmx_set_boolean(GMX_DOUBLE_DEFAULT ON)
option(GMX_DOUBLE "Use double precision (much slower, use only if you really need it)" ${GMX_DOUBLE_DEFAULT})

option(GMX_MPI    "Build a parallel (message-passing) version of ACT" ON)
gmx_dependent_option(
    GMX_MPI_IN_PLACE
    "Enable MPI_IN_PLACE for MPIs that have it defined"
    ON
    GMX_MPI)
mark_as_advanced(GMX_MPI_IN_PLACE)

# There is no GPU support, but the code tests for it
set(GMX_GPU_ACCELERATION_FRAMEWORK "GMX_GPU_NONE")

include(actManageSimd)

option(GMX_BROKEN_CALLOC "Work around broken calloc()" OFF)
mark_as_advanced(GMX_BROKEN_CALLOC)

option(GMX_OPENMP "Enable OpenMP-based multithreading" OFF)

option(GMX_SKIP_DEFAULT_CFLAGS "Don't automatically add suggested/required Compiler flags." OFF)
mark_as_advanced(GMX_SKIP_DEFAULT_CFLAGS)

gmx_set_boolean(GMX_COMPILER_WARNINGS_DEFAULT "NOT SOURCE_IS_SOURCE_DISTRIBUTION")
option(GMX_COMPILER_WARNINGS
    "Enable a default set of compiler warnings"
    ${GMX_COMPILER_WARNINGS_DEFAULT})
mark_as_advanced(GMX_COMPILER_WARNINGS)

option(GMX_BUILD_SHARED_EXE
    "Build exectuables as shared binaries. If not set, this disables rpath and dynamic linker flags in an attempt to build a static binary, but this may require setting up the toolchain properly and making appropriate libraries available."
    ON)
mark_as_advanced(GMX_BUILD_SHARED_EXE)

######################################################################
# Detect OpenMP support
######################################################################
# The OpenMP detection _must_ come before tests for other CFLAGS.
include(gmxManageOpenMP)



######################################################################
# Compiler tests
# These need to be done early (before further tests).
#####################################################################

include(gmxCFlags)
gmx_c_flags()

# This variable should be used for additional compiler flags which are not
# generated in gmxCFlags nor are SIMD or MPI related.
#
# TODO These variables should be consolidated into
# EXTRA_COMPILER_FLAGS so that we we don't perpetrate bugs where
# things that work in C compilation (e.g. merging from old branches)
# also work for C++ compilation.
set(EXTRA_C_FLAGS "")
set(EXTRA_CXX_FLAGS "")

# Run through a number of tests for buggy compilers and other issues
include(gmxTestCompilerProblems)
gmx_test_compiler_problems()

# Implement double-precision option. This is complicated because we
# need installed headers to use the precision mode of the build that
# produced the library, but cannot use config.h in that case. We also
# want such variables to always have a definition, because #if is more
# robust than #ifdef. So, we put this value on the compiler command
# line in all cases.

if(GMX_DOUBLE)
    set(GMX_DOUBLE_VALUE 1)
else()
    set(GMX_DOUBLE_VALUE 0)
endif()
add_compile_definitions(GMX_DOUBLE=${GMX_DOUBLE_VALUE})
list(APPEND INSTALLED_HEADER_DEFINITIONS "-DGMX_DOUBLE=${GMX_DOUBLE_VALUE}")

if(WIN32)
    list(APPEND GMX_EXTRA_LIBRARIES "wsock32")
    add_compile_definitions(GMX_HAVE_WINSOCK)
endif()



########################################################################
# Basic system tests (standard libraries, headers, functions, types)   #
########################################################################
include(CheckIncludeFiles)
include(CheckIncludeFileCXX)
check_include_files(unistd.h     HAVE_UNISTD_H)
check_include_files(pwd.h        HAVE_PWD_H)
check_include_files(dirent.h     HAVE_DIRENT_H)
check_include_files(sys/time.h   HAVE_SYS_TIME_H)
check_include_files(io.h         HAVE_IO_H)
check_include_files(sched.h      HAVE_SCHED_H)

check_include_files(regex.h      HAVE_POSIX_REGEX)
# TODO: It could be nice to inform the user if no regex support is found,
# as selections won't be fully functional.

include(CheckCXXSymbolExists)
check_cxx_symbol_exists(sysconf           unistd.h     HAVE_SYSCONF)
check_cxx_symbol_exists(nice              unistd.h     HAVE_NICE)
check_cxx_symbol_exists(fsync             unistd.h     HAVE_FSYNC)
check_cxx_symbol_exists(_fileno           stdio.h      HAVE__FILENO)
check_cxx_symbol_exists(fileno            stdio.h      HAVE_FILENO)
check_cxx_symbol_exists(_commit           io.h         HAVE__COMMIT)

# We cannot check for the __builtins as symbols, but check if code compiles
check_cxx_source_compiles("int main(){ return __builtin_clz(1);}"   HAVE_BUILTIN_CLZ)
check_cxx_source_compiles("int main(){ return __builtin_clzll(1);}" HAVE_BUILTIN_CLZLL)

include(CheckLibraryExists)
find_library(HAVE_LIBM m)
mark_as_advanced(HAVE_LIBM)
check_library_exists(m feenableexcept "" HAVE_FEENABLEEXCEPT)
check_library_exists(m fedisableexcept "" HAVE_FEDISABLEEXCEPT)

include(TestSchedAffinity)
test_sched_affinity(HAVE_SCHED_AFFINITY)

# Aligned memory allocation. We need to check for both mm_malloc(),
# posix_memalign(), memalign(), and on windows also _aligned_malloc()
include(gmxTestMMMalloc)
gmx_test_mm_malloc(HAVE__MM_MALLOC)
check_cxx_symbol_exists(posix_memalign    stdlib.h     HAVE_POSIX_MEMALIGN)
check_cxx_symbol_exists(memalign          stdlib.h     HAVE_MEMALIGN)

include(TestBigEndian)
test_big_endian(GMX_INTEGER_BIG_ENDIAN)

gmx_set_boolean(GMX_USE_NICE "HAVE_UNISTD_H AND HAVE_NICE")

########################################################################
#Process MPI settings
########################################################################
include(gmxManageMPI)

########################################################################
#Process shared/static library settings
########################################################################
include(gmxManageSharedLibraries)


########################################################################
# Find external packages                                               #
########################################################################

# The libraries are linked as imported targets, which also provide the
# include directories.
find_package(LibXml2 REQUIRED)
list(APPEND GMX_EXTRA_LIBRARIES LibXml2::LibXml2)

find_package(Eigen3 5 REQUIRED NO_MODULE)
list(APPEND GMX_EXTRA_LIBRARIES Eigen3::Eigen)

# When compiling alexandria some extra libraries become optional
# CLN, sqlite3
option(ACT_CLN "Use CLN library to compute Slater integrals" OFF)
if(ACT_CLN)
    find_package(cln REQUIRED)
    if(CLN_FOUND)
      include_directories(${CLN_INCLUDE_DIR})
      set(HAVE_LIBCLN 1)
      list(APPEND	GMX_EXTRA_LIBRARIES ${CLN_LIBRARIES})
    else(CLN_FOUND)
      message(FATAL_ERROR "CLN library not found. Please install the library or use -DACT_CLN=OFF.")
    endif(CLN_FOUND)
else()
    set(HAVE_LIBCLN 0)
endif()

find_package(SQLite3)
if(SQLite3_FOUND)
    set(HAVE_LIBSQLITE3 1)
    list(APPEND GMX_EXTRA_LIBRARIES SQLite::SQLite3)
endif()

# nlohmann/json for JSON reading/writing (header-only library)
find_package(nlohmann_json REQUIRED)
list(APPEND GMX_EXTRA_LIBRARIES nlohmann_json::nlohmann_json)

# We need Boost since RDKit needs it
find_package(Boost CONFIG REQUIRED)
list(APPEND GMX_EXTRA_LIBRARIES Boost::headers)

# These are the necessary RDKit libraries
set(rdkitlibs RDGeneral SmilesParse FileParsers GraphMol Inchi RDInchiLib SubstructMatch SubstructLibrary GeneralizedSubstruct  FilterCatalog DetermineBonds Descriptors )
find_package(RDKit REQUIRED COMPONENTS ${rdkitlibs} )
foreach(lib IN ITEMS ${rdkitlibs})
    list(APPEND GMX_EXTRA_LIBRARIES RDKit::${lib})
endforeach()

if(CYGWIN)
    set(GMX_CYGWIN 1)
endif()

if(WIN32)
    set(GMX_NATIVE_WINDOWS 1)
    # This makes windows.h not declare min/max as macros that would break
    # C++ code using std::min/std::max.
    add_compile_definitions(NOMINMAX)
endif()

option(GMX_BUILD_UNITTESTS "Build unit tests with BUILD_TESTING" ON)
mark_as_advanced(GMX_BUILD_UNITTESTS)
gmx_add_cache_dependency(GMX_BUILD_UNITTESTS BOOL BUILD_TESTING OFF)

########################################################################
# Our own GROMACS tests
########################################################################

add_compile_definitions(HAVE_CONFIG_H)
include_directories(BEFORE ${CMAKE_SOURCE_DIR}/src)
# TODO required at high level because both libactgromacs and progs/mdrun
# require it, both for thread-MPI and its atomics and mutexes.
include_directories(BEFORE SYSTEM ${CMAKE_SOURCE_DIR}/src/external/thread_mpi/include)
# Required for config.h, maybe should only be set in src/CMakeLists.txt
include_directories(BEFORE ${CMAKE_BINARY_DIR}/src)

include(gmxTestInlineASM)
gmx_test_inline_asm_gcc_x86(GMX_X86_GCC_INLINE_ASM)

include(gmxSetBuildInformation)
gmx_set_build_information()

gmx_option_multichoice(
    GMX_USE_RDTSCP
    "Use low-latency RDTSCP instruction for CPU-based timers for mdrun execution; might need to be off when compiling for heterogeneous environments)"
    "AUTO"
    OFF ON AUTO DETECT)
mark_as_advanced(GMX_USE_RDTSCP)

macro(gmx_check_rdtscp)
    if (CPU_DETECTION_FEATURES MATCHES "rdtscp")
        set(HAVE_RDTSCP 1)
        set(RDTSCP_DETECTION_MESSAGE " - detected on the build host")
    else()
        set(RDTSCP_DETECTION_MESSAGE " - not detected on the build host")
    endif()
endmacro()

set(HAVE_RDTSCP 0)
if (GMX_USE_RDTSCP STREQUAL "ON")
    set(HAVE_RDTSCP 1)
elseif(GMX_USE_RDTSCP STREQUAL "DETECT")
    gmx_check_rdtscp()
endif()
gmx_check_if_changed(HAVE_RDTSCP_CHANGED HAVE_RDTSCP)
if (HAVE_RDTSCP_CHANGED)
    if (HAVE_RDTSCP)
        message(STATUS "Enabling RDTSCP support${RDTSCP_DETECTION_MESSAGE}")
    else()
        message(STATUS "Disabling RDTSCP support${RDTSCP_DETECTION_MESSAGE}")
    endif()
endif()

include(gmxTestPipes)
gmx_test_pipes(HAVE_PIPES)

include(gmxTestXDR)
gmx_test_xdr(GMX_SYSTEM_XDR)
# Darwin has system XDR, but it uses a three-argument flavour of
# xdrproc_t that it guarantees will still work if you pass the normal
# two-argument xdr filters, but gcc 8 warns about the cast necessary
# to do that, so it's simpler to just use our own XDR library.
#
# TODO It would be better to craft a cmake test which fails if such
# XDR operations cause warnings, and succeeds otherwise, because it is
# generally preferable to use system libraries where possible.
if(NOT GMX_SYSTEM_XDR OR CMAKE_SYSTEM_NAME STREQUAL "Darwin")
    set(GMX_INTERNAL_XDR 1)
endif()


include(gmxManageCycleCounters)
gmx_manage_cycle_counters()

# Math and thread libraries must often come after all others when linking...
if (HAVE_LIBM)
    list(APPEND GMX_PUBLIC_LIBRARIES m)
endif()

# # # # # # # # # # NO MORE TESTS AFTER THIS LINE! # # # # # # # # # # #
# these are set after everything else
if (NOT GMX_SKIP_DEFAULT_CFLAGS)
    set(CMAKE_C_FLAGS "${GMX_SIMD_VECTORIZATION_C_FLAGS} ${MPI_COMPILE_FLAGS} ${EXTRA_C_FLAGS} ${CMAKE_C_FLAGS}")
    set(CMAKE_CXX_FLAGS "${GMX_SIMD_VECTORIZATION_CXX_FLAGS} ${MPI_COMPILE_FLAGS} ${EXTRA_CXX_FLAGS} ${CMAKE_CXX_FLAGS}")
    set(CMAKE_EXE_LINKER_FLAGS "${MPI_LINKER_FLAGS} ${GMX_SIMD_VECTORIZATION_LINKER_FLAGS} ${CMAKE_EXE_LINKER_FLAGS}")
    set(CMAKE_SHARED_LINKER_FLAGS "${MPI_LINKER_FLAGS} ${GMX_SIMD_VECTORIZATION_LINKER_FLAGS} ${CMAKE_SHARED_LINKER_FLAGS}")
else()
    message("Recommended flags which are not added because GMX_SKIP_DEFAULT_CFLAGS=yes:")
    message("CMAKE_C_FLAGS: ${GMX_SIMD_VECTORIZATION_C_FLAGS} ${MPI_COMPILE_FLAGS} ${EXTRA_C_FLAGS} ${GMXC_CFLAGS}")
    message("CMAKE_C_FLAGS_RELEASE: ${GMXC_CFLAGS_RELEASE}")
    message("CMAKE_C_FLAGS_DEBUG: ${GMXC_CFLAGS_DEBUG}")
    message("CMAKE_CXX_FLAGS: ${GMX_SIMD_VECTORIZATION_CXX_FLAGS} ${MPI_COMPILE_FLAGS} ${EXTRA_CXX_FLAGS} ${GMXC_CXXFLAGS}")
    message("CMAKE_CXX_FLAGS_RELEASE: ${GMXC_CXXFLAGS_RELEASE}")
    message("CMAKE_CXX_FLAGS_DEBUG: ${GMXC_CXXFLAGS_DEBUG}")
    message("CMAKE_EXE_LINKER_FLAGS: ${MPI_LINKER_FLAGS} ${GMX_SIMD_VECTORIZATION_LINKER_FLAGS}")
    message("CMAKE_SHARED_LINKER_FLAGS: ${MPI_LINKER_FLAGS} ${GMX_SIMD_VECTORIZATION_LINKER_FLAGS}")
endif()

if(NOT GMX_OPENMP)
    #Unset all OpenMP flags in case OpenMP was disabled either by the user
    #or because it was only partially detected (e.g. only for C but not C++ compiler)
    unset(OpenMP_C_FLAGS CACHE)
    unset(OpenMP_CXX_FLAGS CACHE)
else()
    set(GMX_EXE_LINKER_FLAGS ${GMX_EXE_LINKER_FLAGS} ${OpenMP_LINKER_FLAGS})
    set(GMX_SHARED_LINKER_FLAGS ${GMX_SHARED_LINKER_FLAGS} ${OpenMP_SHARED_LINKER_FLAGS})
endif()


########################################################################
# Specify install locations
########################################################################
# Use GNUInstallDirs to set paths on multiarch systems.
include(GNUInstallDirs)

set(ACT_INSTALL_DATASUBDIR "act" CACHE STRING "Subdirectory for ACT data under CMAKE_INSTALL_DATADIR")
mark_as_advanced(ACT_INSTALL_DATASUBDIR)

# Internal convenience so we do not have to join two path segments in the code
set(ACT_INSTALL_ACTDATADIR ${CMAKE_INSTALL_DATADIR}/${ACT_INSTALL_DATASUBDIR})

# If the nesting level wrt. the installation root is changed,
# gromacs-config.cmake.cmakein needs to be adapted.
set(GMX_INSTALL_CMAKEDIR  ${CMAKE_INSTALL_DATAROOTDIR}/cmake)

# TODO: Make GMXRC adapt if this is changed
set(GMX_INSTALL_PKGCONFIGDIR ${CMAKE_INSTALL_LIBDIR}/pkgconfig)

list(APPEND INSTALLED_HEADER_INCLUDE_DIRS ${CMAKE_INSTALL_INCLUDEDIR})

################################################################
# Shared library load path settings
################################################################
if(NOT GMX_BUILD_SHARED_EXE)
    # No rpath
    set(CMAKE_SKIP_RPATH TRUE)
    set(CMAKE_EXE_LINK_DYNAMIC_C_FLAGS) # remove -Wl,-Bdynamic
    set(CMAKE_EXE_LINK_DYNAMIC_CXX_FLAGS)
else()
    # The build folder always has bin/ and lib/; if we are also going to
    # install to lib/, then the installation RPATH works also in the build
    # tree.  This makes installation slightly faster (no need to rewrite the
    # RPATHs), and makes the binaries in the build tree relocatable.
    if(CMAKE_INSTALL_LIBDIR STREQUAL "lib")
        set(CMAKE_BUILD_WITH_INSTALL_RPATH TRUE)
        set(CMAKE_BUILD_WITH_INSTALL_NAME_DIR TRUE)
    endif()
    # Set the RPATH as relative to the executable location to make the
    # binaries relocatable.
    if(CMAKE_SYSTEM_NAME STREQUAL "Darwin") #Assume OS X >=10.5
        set(CMAKE_INSTALL_RPATH "@loader_path/../${CMAKE_INSTALL_LIBDIR}")
        set(CMAKE_INSTALL_NAME_DIR ${CMAKE_INSTALL_RPATH})
    else()
        set(CMAKE_INSTALL_RPATH "\$ORIGIN/../${CMAKE_INSTALL_LIBDIR}")
    endif()
    set(CMAKE_INSTALL_RPATH_USE_LINK_PATH TRUE)
    set(CMAKE_MACOSX_RPATH 1)
endif()

#COPYING file: Only necessary for binary distributions.
#Simpler to always install.
install(FILES COPYING DESTINATION ${ACT_INSTALL_ACTDATADIR} COMPONENT data)

if (BUILD_TESTING)
    include(tests/CheckTarget.cmake)
endif()

find_package(ImageMagick QUIET COMPONENTS convert)
include(gmxTestImageMagick)
GMX_TEST_IMAGEMAGICK(IMAGE_CONVERT_POSSIBLE)
add_subdirectory(docs)
add_subdirectory(share)
add_subdirectory(scripts)
add_subdirectory(src)

gmx_cpack_write_config()

#######################
## uninstall target
#######################
configure_file("${CMAKE_CURRENT_SOURCE_DIR}/cmake/cmake_uninstall.cmake.in"
               "${CMAKE_CURRENT_BINARY_DIR}/cmake/cmake_uninstall.cmake"
               @ONLY)
add_custom_target(uninstall
                  "${CMAKE_COMMAND}" -P
                  "${CMAKE_CURRENT_BINARY_DIR}/cmake/cmake_uninstall.cmake")
set_directory_properties(PROPERTIES
            ADDITIONAL_CLEAN_FILES "install_manifest.txt")

#
# This source file is part of the Alexandria Chemistry Toolkit.
#
# This program is free software; you can redistribute it and/or
# modify it under the terms of the GNU General Public License
# as published by the Free Software Foundation; either version 2
# of the License, or (at your option) any later version.
#

# The Python interface to ACT, the package act. It runs the alexandria
# program, so it needs nothing to be compiled. The tests use the standard
# unittest module, such that they run with any Python 3.9 or newer without
# other packages; pytest runs them as well.
#
# The package is installed in site-packages/act of the Python that CMake
# found. ACT-utilities, which had its modules there until version 1.9b, has
# them in site-packages/act_utilities now. A build for a package manager that
# puts the Python code in a package of its own can switch the installation
# off.
option(ACT_INSTALL_PYTHON_PACKAGE "Install the Python package act, the interface to the alexandria program" ON)

if(ACT_INSTALL_PYTHON_PACKAGE)
    set(ACT_PYTHON_PACKAGE_DIR ${CMAKE_INSTALL_PYTHON_LIBDIR}/act)
    file(GLOB ACT_PYTHON_PACKAGE_FILES CONFIGURE_DEPENDS ${CMAKE_CURRENT_SOURCE_DIR}/act/*.py)
    install(FILES ${ACT_PYTHON_PACKAGE_FILES}
            DESTINATION ${ACT_PYTHON_PACKAGE_DIR} COMPONENT python)
    # Where the program is, relative to the package, such that the package
    # runs the alexandria it was installed with, and the installation can
    # be moved as a whole.
    file(RELATIVE_PATH ACT_PYTHON_TO_BINDIR
         ${CMAKE_INSTALL_PREFIX}/${ACT_PYTHON_PACKAGE_DIR}
         ${CMAKE_INSTALL_PREFIX}/${CMAKE_INSTALL_BINDIR})
    file(GENERATE OUTPUT ${CMAKE_CURRENT_BINARY_DIR}/_installed.py
         CONTENT "\"\"\"!
@file
@brief Written by the installation of ACT: where the program is.
\"\"\"
## The alexandria program, relative to the directory of this package
EXECUTABLE = \"${ACT_PYTHON_TO_BINDIR}/$<TARGET_FILE_NAME:alexandria>\"
")
    install(FILES ${CMAKE_CURRENT_BINARY_DIR}/_installed.py
            DESTINATION ${ACT_PYTHON_PACKAGE_DIR} COMPONENT python)
endif()

# The classes of ACT in Python: the module act._core, see src/core.cpp.
# It is built next to a copy of the package in the build tree, such that
# the package can be used and tested from there with
# PYTHONPATH=<build>/python.
if(ACT_PYTHON)
    if(NOT Python3_Development_FOUND)
        message(FATAL_ERROR "ACT_PYTHON needs the headers of Python ${Python3_VERSION}, "
                            "which were not found for ${Python3_EXECUTABLE}")
    endif()
    # The same Python as the rest of ACT uses
    set(PYBIND11_FINDPYTHON ON)
    set(Python_EXECUTABLE ${Python3_EXECUTABLE})
    find_package(pybind11 2.10 CONFIG REQUIRED)
    set(ACT_PYTHON_BUILD_PACKAGE ${CMAKE_BINARY_DIR}/python/act)
    pybind11_add_module(act_python_core MODULE src/core.cpp)
    set_target_properties(act_python_core PROPERTIES
                          OUTPUT_NAME _core
                          LIBRARY_OUTPUT_DIRECTORY ${ACT_PYTHON_BUILD_PACKAGE})
    # dladdr, to find where the module is, is in a library of its own on
    # systems with an older C library
    target_link_libraries(act_python_core PRIVATE actLib ${CMAKE_DL_LIBS})
    # From the package where the module is installed to the installation
    # prefix, where the data files of ACT are
    file(RELATIVE_PATH ACT_PYTHON_TO_PREFIX
         ${CMAKE_INSTALL_PREFIX}/${CMAKE_INSTALL_PYTHON_LIBDIR}/act ${CMAKE_INSTALL_PREFIX})
    target_compile_definitions(act_python_core PRIVATE ACT_PYTHON_TO_PREFIX="${ACT_PYTHON_TO_PREFIX}")
    # With shared libraries the module has to find those of ACT from where
    # it is: python/act in the build tree, the package directory when
    # installed. The programs find them from bin, see CMAKE_INSTALL_RPATH.
    if(BUILD_SHARED_LIBS)
        if(APPLE)
            set(ACT_PYTHON_ORIGIN "@loader_path")
        else()
            set(ACT_PYTHON_ORIGIN "\$ORIGIN")
        endif()
        file(RELATIVE_PATH ACT_PYTHON_TO_LIBDIR
             ${CMAKE_INSTALL_PREFIX}/${CMAKE_INSTALL_PYTHON_LIBDIR}/act
             ${CMAKE_INSTALL_PREFIX}/${CMAKE_INSTALL_LIBDIR})
        set_property(TARGET act_python_core APPEND PROPERTY INSTALL_RPATH
                     "${ACT_PYTHON_ORIGIN}/../../${CMAKE_INSTALL_LIBDIR}"
                     "${ACT_PYTHON_ORIGIN}/${ACT_PYTHON_TO_LIBDIR}")
        set_property(TARGET act_python_core APPEND PROPERTY BUILD_RPATH
                     "${ACT_PYTHON_ORIGIN}/../../${CMAKE_INSTALL_LIBDIR}")
        if(APPLE)
            # The libraries of ACT are named after where the programs find
            # them, see CMAKE_INSTALL_NAME_DIR; the module looks them up
            # in its own run path instead
            foreach(library actLib libactgromacs)
                add_custom_command(TARGET act_python_core POST_BUILD
                    COMMAND install_name_tool -change
                            "${CMAKE_INSTALL_NAME_DIR}/$<TARGET_SONAME_FILE_NAME:${library}>"
                            "@rpath/$<TARGET_SONAME_FILE_NAME:${library}>"
                            $<TARGET_FILE:act_python_core>
                    VERBATIM)
            endforeach()
        endif()
    endif()
    # The package itself in the build tree, next to the module
    file(GLOB ACT_PYTHON_SOURCE_FILES CONFIGURE_DEPENDS ${CMAKE_CURRENT_SOURCE_DIR}/act/*.py)
    set(ACT_PYTHON_BUILD_FILES)
    foreach(source ${ACT_PYTHON_SOURCE_FILES})
        get_filename_component(name ${source} NAME)
        add_custom_command(OUTPUT ${ACT_PYTHON_BUILD_PACKAGE}/${name}
                           COMMAND ${CMAKE_COMMAND} -E copy_if_different ${source}
                                   ${ACT_PYTHON_BUILD_PACKAGE}/${name}
                           DEPENDS ${source}
                           COMMENT "Copying ${name} of the Python package")
        list(APPEND ACT_PYTHON_BUILD_FILES ${ACT_PYTHON_BUILD_PACKAGE}/${name})
    endforeach()
    add_custom_target(act_python_package ALL DEPENDS ${ACT_PYTHON_BUILD_FILES})
    add_dependencies(act_python_core act_python_package)
    if(ACT_INSTALL_PYTHON_PACKAGE)
        install(TARGETS act_python_core
                LIBRARY DESTINATION ${ACT_PYTHON_PACKAGE_DIR} COMPONENT python)
    endif()
endif()

if(BUILD_TESTING AND Python3_Interpreter_FOUND)
    add_test(NAME PythonInterfaceTests
             COMMAND ${Python3_EXECUTABLE} -m unittest discover
                     -s ${CMAKE_CURRENT_SOURCE_DIR}/tests
                     -t ${CMAKE_CURRENT_SOURCE_DIR}
             WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR})
    # The tests with the real program find it and the test data through
    # these variables, and are skipped without them.
    set_tests_properties(PythonInterfaceTests PROPERTIES
        ENVIRONMENT "ACT_ALEXANDRIA=$<TARGET_FILE:alexandria>;ACT_TESTDATA=${CMAKE_SOURCE_DIR}/src/act/testdata"
        LABELS "PythonTest"
        TIMEOUT 300)
    if(ACT_PYTHON)
        # The classes, with the package of the build tree that has the module
        add_test(NAME PythonClassTests
                 COMMAND ${Python3_EXECUTABLE}
                         ${CMAKE_CURRENT_SOURCE_DIR}/tests/test_core.py -v
                 WORKING_DIRECTORY ${CMAKE_BINARY_DIR}/python)
        set_tests_properties(PythonClassTests PROPERTIES
            ENVIRONMENT "ACT_PYTHON_PATH=${CMAKE_BINARY_DIR}/python;ACT_ALEXANDRIA=$<TARGET_FILE:alexandria>;ACT_TESTDATA=${CMAKE_SOURCE_DIR}/src/act/testdata;ACT_REQUIRE_CLASSES=1"
            LABELS "PythonTest"
            TIMEOUT 300)
    endif()
    if(ACT_INSTALL_PYTHON_PACKAGE)
        # Installs the package in a directory of the build tree and imports
        # it from there. With the module of the classes, if that has the
        # libraries of ACT in it: shared libraries are not in that directory.
        set(ACT_PYTHON_INSTALL_TEST_CLASSES OFF)
        if(ACT_PYTHON AND NOT BUILD_SHARED_LIBS)
            set(ACT_PYTHON_INSTALL_TEST_CLASSES ON)
        endif()
        add_test(NAME PythonInstallTest
                 COMMAND ${CMAKE_COMMAND}
                         -DACT_BINARY_DIR=${CMAKE_BINARY_DIR}
                         -DACT_PYTHON=${Python3_EXECUTABLE}
                         -DACT_SITE_DIR=${CMAKE_INSTALL_PYTHON_LIBDIR}
                         -DACT_BIN_DIR=${CMAKE_INSTALL_BINDIR}
                         -DACT_PROGRAM=$<TARGET_FILE_NAME:alexandria>
                         -DACT_CLASSES=${ACT_PYTHON_INSTALL_TEST_CLASSES}
                         -DACT_CHECK=${CMAKE_CURRENT_SOURCE_DIR}/tests/check_installed.py
                         -P ${CMAKE_CURRENT_SOURCE_DIR}/tests/install_test.cmake)
        set_tests_properties(PythonInstallTest PROPERTIES
            LABELS "PythonTest"
            TIMEOUT 120)
    endif()
endif()

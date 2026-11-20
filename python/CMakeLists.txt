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
    if(ACT_INSTALL_PYTHON_PACKAGE)
        # Installs the package in a directory of the build tree and imports
        # it from there
        add_test(NAME PythonInstallTest
                 COMMAND ${CMAKE_COMMAND}
                         -DACT_BINARY_DIR=${CMAKE_BINARY_DIR}
                         -DACT_PYTHON=${Python3_EXECUTABLE}
                         -DACT_SITE_DIR=${CMAKE_INSTALL_PYTHON_LIBDIR}
                         -DACT_BIN_DIR=${CMAKE_INSTALL_BINDIR}
                         -DACT_PROGRAM=$<TARGET_FILE_NAME:alexandria>
                         -DACT_CHECK=${CMAKE_CURRENT_SOURCE_DIR}/tests/check_installed.py
                         -P ${CMAKE_CURRENT_SOURCE_DIR}/tests/install_test.cmake)
        set_tests_properties(PythonInstallTest PROPERTIES
            LABELS "PythonTest"
            TIMEOUT 120)
    endif()
endif()

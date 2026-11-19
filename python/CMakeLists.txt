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
# The package is not installed yet: the directory site-packages/act is where
# ACT-utilities puts its modules, and it adds that directory to the path of
# Python, which would make the modules of this package importable under
# their bare names (commands, options, ...) as well. See
# docs/python-api-gui-plan.md.

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
endif()

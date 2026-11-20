#
# This source file is part of the Alexandria Chemistry Toolkit.
#
# This program is free software; you can redistribute it and/or
# modify it under the terms of the GNU General Public License
# as published by the Free Software Foundation; either version 2
# of the License, or (at your option) any later version.
#
# Test of the installation of the Python package act: installs the package,
# and nothing else of ACT, in a new directory in the build tree and lets
# check_installed.py import it from there. Run by ctest, see
# python/CMakeLists.txt for the variables.
#
set(_prefix "${ACT_BINARY_DIR}/python/install-test")
file(REMOVE_RECURSE "${_prefix}")
execute_process(
    COMMAND ${CMAKE_COMMAND} --install ${ACT_BINARY_DIR} --component python --prefix ${_prefix}
    RESULT_VARIABLE _result
    OUTPUT_VARIABLE _output
    ERROR_VARIABLE  _output)
if(NOT _result EQUAL 0)
    message(FATAL_ERROR "Installing the Python package failed:\n${_output}")
endif()
# -E -s: without PYTHONPATH and the packages of the user, which may have
# another act
execute_process(
    COMMAND ${ACT_PYTHON} -E -s ${ACT_CHECK}
            ${_prefix}/${ACT_SITE_DIR} ${_prefix}/${ACT_BIN_DIR}/${ACT_PROGRAM}
    WORKING_DIRECTORY ${_prefix}
    RESULT_VARIABLE _result)
if(NOT _result EQUAL 0)
    message(FATAL_ERROR "The installed Python package does not work")
endif()

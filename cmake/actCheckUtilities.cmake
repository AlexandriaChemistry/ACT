#
# This source file is part of the Alexandria Chemistry Toolkit.
#
# This program is free software; you can redistribute it and/or
# modify it under the terms of the GNU General Public License
# as published by the Free Software Foundation; either version 2
# of the License, or (at your option) any later version.
#

# ACT requires the Python modules and scripts of ACT-utilities,
# https://github.com/AlexandriaChemistry/ACT-utilities, which were in
# src/act/python before. Check that the Python found by CMake can import
# them, either from its own site-packages, for instance in a conda
# environment, or from the installation prefix of ACT, where ACTRC adds them
# to PYTHONPATH.

if(NOT Python3_Interpreter_FOUND)
    message(FATAL_ERROR "ACT needs Python 3 with the modules of ACT-utilities, but no Python 3 was found")
endif()

set(_act_utilities_path "${CMAKE_INSTALL_PREFIX}/${CMAKE_INSTALL_PYTHON_LIBDIR}/act")
if(DEFINED ENV{PYTHONPATH} AND NOT "$ENV{PYTHONPATH}" STREQUAL "")
    string(APPEND _act_utilities_path ":$ENV{PYTHONPATH}")
endif()
execute_process(
    COMMAND ${CMAKE_COMMAND} -E env "PYTHONPATH=${_act_utilities_path}"
            ${Python3_EXECUTABLE} -c "import actutils, molprops; print(actutils.__file__)"
    RESULT_VARIABLE _act_utilities_result
    OUTPUT_VARIABLE _act_utilities_file
    ERROR_QUIET
    OUTPUT_STRIP_TRAILING_WHITESPACE)
if(NOT _act_utilities_result EQUAL 0)
    message(FATAL_ERROR
        "ACT requires the Python modules and scripts of ACT-utilities, which "
        "${Python3_EXECUTABLE} cannot import. Install them first, either with conda:\n"
        "  conda install -c conda-forge act-utilities\n"
        "or from source, into the same prefix as ACT:\n"
        "  git clone https://github.com/AlexandriaChemistry/ACT-utilities.git\n"
        "  cmake -S ACT-utilities -B ACT-utilities/build -DCMAKE_INSTALL_PREFIX=${CMAKE_INSTALL_PREFIX}\n"
        "  cmake --install ACT-utilities/build\n")
endif()
get_filename_component(_act_utilities_dir "${_act_utilities_file}" DIRECTORY)
message(STATUS "Found the Python modules of ACT-utilities in ${_act_utilities_dir}")

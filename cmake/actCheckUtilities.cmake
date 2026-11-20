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
#
# The modules are in site-packages/act_utilities. Until ACT-utilities 1.9b
# they were in site-packages/act, which is now the directory of the Python
# package act of ACT (see python/CMakeLists.txt): such an installation has to
# be updated, since its module act.py would hide the package.

if(NOT Python3_Interpreter_FOUND)
    message(FATAL_ERROR "ACT needs Python 3 with the modules of ACT-utilities, but no Python 3 was found")
endif()

set(_act_utilities_site "${CMAKE_INSTALL_PREFIX}/${CMAKE_INSTALL_PYTHON_LIBDIR}")
set(_act_utilities_path "${_act_utilities_site}/act_utilities")
set(_act_utilities_update
    "Update ACT-utilities, either with conda:\n"
    "  conda install -c conda-forge act-utilities\n"
    "or from source, into the same prefix as before, which also removes the old files:\n"
    "  git clone https://github.com/AlexandriaChemistry/ACT-utilities.git\n"
    "  cmake -S ACT-utilities -B ACT-utilities/build -DCMAKE_INSTALL_PREFIX=<prefix>\n"
    "  cmake --install ACT-utilities/build\n")
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
    # ACT-utilities installed for another version of Python?
    file(GLOB _act_utilities_other "${CMAKE_INSTALL_PREFIX}/lib/python*/site-packages/act_utilities/actutils.py")
    if(_act_utilities_other)
        list(GET _act_utilities_other 0 _act_utilities_other)
        message(FATAL_ERROR
            "ACT-utilities is installed in ${_act_utilities_other}, but for another "
            "version of Python than ${Python3_EXECUTABLE} (${Python3_VERSION}), which CMake "
            "found for ACT. Install ACT-utilities with that Python, for instance with "
            "-DPython3_EXECUTABLE=${Python3_EXECUTABLE} or with the same CMAKE_PREFIX_PATH "
            "as ACT, or give ACT the Python of ACT-utilities with -DPython3_EXECUTABLE.")
    endif()
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
# An installation from before the modules moved to site-packages/act_utilities:
# it has the module act.py, which got the name act_legacy.py
if(EXISTS "${_act_utilities_dir}/act.py" AND NOT EXISTS "${_act_utilities_dir}/act_legacy.py")
    message(FATAL_ERROR
        "The ACT-utilities in ${_act_utilities_dir} is too old for this version of ACT: "
        "its module act.py has the name of the Python package act of ACT and would "
        "hide it.\n" ${_act_utilities_update})
endif()
# Files of the old layout that are left where the package act is installed
if(EXISTS "${_act_utilities_site}/act/actutils.py" OR EXISTS "${_act_utilities_site}/act.pth")
    message(FATAL_ERROR
        "There are files of an earlier ACT-utilities, version 1.9b or older, in ${_act_utilities_site}/act, "
        "where the Python package act of ACT is installed.\n" ${_act_utilities_update})
endif()
message(STATUS "Found the Python modules of ACT-utilities in ${_act_utilities_dir}")

# The data files of ACT-utilities, which alexandria reads as well, such as
# atomprops.csv. They are in share/act of the prefix of ACT-utilities, which
# the programs search after that of ACT, such that ACT-utilities may be
# installed in another prefix. ACTDATA is not used here, since it may point
# to an older installation.
execute_process(
    COMMAND ${CMAKE_COMMAND} -E env --unset=ACTDATA "PYTHONPATH=${_act_utilities_path}"
            ${Python3_EXECUTABLE} -c "import actutils; print(actutils.actdata())"
    RESULT_VARIABLE _act_utilities_result
    OUTPUT_VARIABLE ACT_UTILITIES_DATADIR
    ERROR_QUIET
    OUTPUT_STRIP_TRAILING_WHITESPACE)
if(NOT _act_utilities_result EQUAL 0 OR NOT EXISTS "${ACT_UTILITIES_DATADIR}/atomprops.csv")
    message(FATAL_ERROR
        "The data files of ACT-utilities, such as atomprops.csv, were not found next to "
        "its Python modules in ${_act_utilities_dir}. Install a version of ACT-utilities "
        "that has them, see the message above for how.")
endif()
message(STATUS "Found the data files of ACT-utilities in ${ACT_UTILITIES_DATADIR}")

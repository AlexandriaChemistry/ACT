#
# This file is part of the Alexandria Chemistry Toolkit
# https://github.com/AlexandriaChemistry
#
# Copyright (c) 2026, David van der Spoel, Marie-Madeleine Walz, Paul J. van Maaren
#
# The "tidy" target: clean up the code of src/act with clang-tidy and apply
# the fixes in place, in the layout of the ACT code (4 spaces indentation),
# see admin/clang_tidy_fix.py. The checks are in admin/clang-tidy-fix.yaml.
#   make tidy                                # all of src/act
#   make tidy ACT_TIDY_PATHS=src/act/forces  # only these files or directories
#   make tidy ACT_TIDY_ARGS=--dry-run        # only report the fixes
# The script needs Python with PyYAML.

# clang-tidy reads the compiler flags of each file from compile_commands.json
set(CMAKE_EXPORT_COMPILE_COMMANDS ON)

# Looked for when the target is run, so it need not be present when configuring
set(ACT_CLANG_TIDY "clang-tidy" CACHE STRING "clang-tidy executable for make tidy")
mark_as_advanced(ACT_CLANG_TIDY)

add_custom_target(tidy
    COMMAND ${Python3_EXECUTABLE} ${CMAKE_SOURCE_DIR}/admin/clang_tidy_fix.py
            --build-dir ${CMAKE_BINARY_DIR}
            --source-dir ${CMAKE_SOURCE_DIR}
            --clang-tidy ${ACT_CLANG_TIDY}
    WORKING_DIRECTORY ${CMAKE_SOURCE_DIR}
    COMMENT "Cleaning up src/act with clang-tidy"
    USES_TERMINAL
    VERBATIM)

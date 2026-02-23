
set(SOURCE_MD5SUM "unknown" CACHE STRING
    "MD5 sum of the source tarball, normally used only for the pre-release webpage build")
mark_as_advanced(SOURCE_MD5SUM)

# Even if we aren't going to make the full webpage, set up to put all
# the documentation output in the same place, for convenience
set(HTML_OUTPUT_DIR "${CMAKE_CURRENT_BINARY_DIR}/html")
file(MAKE_DIRECTORY ${HTML_OUTPUT_DIR})

option(ACT_BUILD_MANUAL "Whether to try to configure to build the on-line and PDF manual" OFF)
mark_as_advanced(ACT_BUILD_MANUAL)
if (ACT_BUILD_MANUAL)
    add_subdirectory(doxygen)
    add_subdirectory(manual)
endif()

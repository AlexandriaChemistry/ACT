#
# This source file is part of the Alexandria Chemistry Toolkit.
#
# Copyright (C) 2026
#
# Developers:
#             Mohammad Mehdi Ghahremanpour,
#             Paul J. van Maaren,
#             David van der Spoel (Project leader)
#
# This program is free software; you can redistribute it and/or
# modify it under the terms of the GNU General Public License
# as published by the Free Software Foundation; either version 2
# of the License, or (at your option) any later version.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with this program; if not, write to the Free Software
# Foundation, Inc., 51 Franklin Street, Fifth Floor,
# Boston, MA  02110-1301, USA.
#

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

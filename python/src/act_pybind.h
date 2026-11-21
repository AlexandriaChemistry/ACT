/*
 * This source file is part of the Alexandria Chemistry Toolkit.
 *
 * Copyright (C) 2026
 *
 * Developers:
 *             Mohammad Mehdi Ghahremanpour,
 *             Julian Marrades,
 *             Marie-Madeleine Walz,
 *             Paul J. van Maaren,
 *             David van der Spoel (Project leader)
 *
 * This program is free software; you can redistribute it and/or
 * modify it under the terms of the GNU General Public License
 * as published by the Free Software Foundation; either version 2
 * of the License, or (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this program; if not, write to the Free Software
 * Foundation, Inc., 51 Franklin Street, Fifth Floor,
 * Boston, MA  02110-1301, USA.
 */
/*! \internal \brief
 * The headers of pybind11 for the source files of the Python module of ACT.
 *
 * Include this file before any header of ACT. The configuration of Python
 * (pyconfig.h) and that of ACT (config.h) both say what the system has,
 * with macros of the same names and not always the same values. The
 * headers of Python have used theirs by the end of this file, so the names
 * that ACT defines are freed here for the headers of ACT.
 * \author David van der Spoel <david.vanderspoel@icm.uu.se>
 */
#ifndef ACT_PYTHON_PYBIND_H
#define ACT_PYTHON_PYBIND_H

#include <pybind11/numpy.h>
#include <pybind11/pybind11.h>
#include <pybind11/stl.h>
#include <pybind11/stl/filesystem.h>

// The macros of src/config.h.cmakein that tell what the system has
#undef HAVE_BITSCANREVERSE
#undef HAVE_BITSCANREVERSE64
#undef HAVE_BUILTIN_CLZ
#undef HAVE_BUILTIN_CLZLL
#undef HAVE_CNTLZ4
#undef HAVE_CNTLZ8
#undef HAVE_CXX17_REGEX
#undef HAVE_DIRENT_H
#undef HAVE_EXTRAE
#undef HAVE_FEDISABLEEXCEPT
#undef HAVE_FEENABLEEXCEPT
#undef HAVE_FILENO
#undef HAVE_FSEEKO
#undef HAVE_FSYNC
#undef HAVE_IO_H
#undef HAVE_LIBCLN
#undef HAVE_LIBSQLITE3
#undef HAVE_MALLOC_H
#undef HAVE_MEMALIGN
#undef HAVE_MM_MALLOC_H
#undef HAVE_PIPES
#undef HAVE_POSIX_MEMALIGN
#undef HAVE_POSIX_REGEX
#undef HAVE_PWD_H
#undef HAVE_RDTSCP
#undef HAVE_SCHED_AFFINITY
#undef HAVE_SCHED_H
#undef HAVE_SYSCONF
#undef HAVE_SYS_TIME_H
#undef HAVE_UNISTD_H
#undef HAVE_XMMINTRIN_H
#undef HAVE__ALIGNED_MALLOC
#undef HAVE__COMMIT
#undef HAVE__FILENO
#undef HAVE__FSEEKI64
#undef HAVE__MM_MALLOC

#endif

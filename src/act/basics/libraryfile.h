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
 * Implements part of the alexandria program.
 * \author David van der Spoel <david.vanderspoel@icm.uu.se>
 */

#ifndef ACT_LIBFILE_H
#define ACT_LIBFILE_H

#include <string>

namespace alexandria
{

/*! \brief Return the location of a library file.
 * \param[in] filename The file looked after
 * \param[in] useCWD   Whether to search in the current working directory
 * \return either the unmodified filename if present in the current working
 *         directory (if useCWD is true) or the filename in the ACT library
 *         directory (if present there). If the file is not present, an
 *         empty string is returned.
 */
std::string findLibrary(const std::string &filename,
                        bool               useCWD);
                            
} // namespace

#endif

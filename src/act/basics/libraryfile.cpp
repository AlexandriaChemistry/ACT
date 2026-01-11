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

#include "libraryfile.h"

#include <cstdlib>
#include <string>

#include "gromacs/utility/futil.h"
#include "gromacs/utility/stringutil.h"

namespace alexandria
{

std::string findLibrary(const std::string &filename,
                        bool               useCWD)
{
    if (useCWD && gmx_fexist(filename))
    {
        return filename;
    }
    std::string libFile;
    std::string actdir("ACTDATA");
    auto actdata = std::getenv(actdir.c_str());
    if (nullptr != actdata)
    {
        libFile = gmx::formatString("%s/%s", actdata, filename.c_str());
        if (!gmx_fexist(libFile))
        {
            libFile.clear();
        }
    }
    else
    {
        fprintf(stderr, "Please set the environment variable '%s' to point to the correct directory. Maybe you should source the ACTRC file first.", actdir.c_str());
    }
    return libFile;
}

} // namespace

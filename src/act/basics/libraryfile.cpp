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
#include <vector>

#include "buildinfo.h"
#include "gromacs/utility/futil.h"
#include "gromacs/utility/path.h"
#include "gromacs/utility/programcontext.h"
#include "gromacs/utility/stringutil.h"

namespace alexandria
{

/*! \brief The data directory of the installation or source tree of the program
 * The program context finds it from the location of the binary: the
 * share/act directory of the installation, or the share directory of the
 * source tree for a binary in the build tree.
 * \return the directory, empty if not known
 */
static std::string installedDataDirectory()
{
    const gmx::InstallationPrefixInfo prefix = gmx::getProgramContext().installationPrefix();
    if (gmx::isNullOrEmpty(prefix.path))
    {
        return std::string();
    }
    return gmx::Path::join(prefix.path, prefix.bSourceLayout ? "share" : ACT_INSTALL_ACTDATADIR);
}

std::string findLibrary(const std::string &filename,
                        bool               useCWD)
{
    return findLibraryFirst({ filename }, useCWD);
}

std::string findLibraryFirst(const std::vector<std::string> &filenames,
                             bool                            useCWD)
{
    if (useCWD)
    {
        for (const auto &filename : filenames)
        {
            if (gmx_fexist(filename))
            {
                return filename;
            }
        }
    }
    // The directory in the environment variable ACTDATA, if set, comes
    // before that of the installation.
    std::vector<std::string> dirs;
    auto actdata = std::getenv("ACTDATA");
    if (!gmx::isNullOrEmpty(actdata))
    {
        dirs.emplace_back(actdata);
    }
    dirs.push_back(installedDataDirectory());
    for (const auto &dir : dirs)
    {
        if (dir.empty())
        {
            continue;
        }
        for (const auto &filename : filenames)
        {
            std::string libFile = gmx::Path::join(dir, filename);
            if (gmx_fexist(libFile))
            {
                return libFile;
            }
        }
    }
    return std::string();
}

} // namespace

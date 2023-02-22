/*
 * This source file is part of the Alexandria Chemistry Toolkit.
 *
 * Copyright (C) 2023
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
#include "atomprops.h"

#include "act/utility/stringutil.h"
#include "gromacs/utility/exceptions.h"
#include "gromacs/utility/stringutil.h"
#include "gromacs/utility/textreader.h"

namespace alexandria
{

std::map<std::string, AtomProp> readAtomProps()
{
    const char *actdata = getenv("ACTDATA");
    if (nullptr == actdata)
    {
        GMX_THROW(gmx::InvalidInputError("Environment variable ACTDATA is not set"));
    }
    std::string     props = gmx::formatString("%s/atomprops.csv", actdata);
    gmx::TextReader tr(props);
    std::map<std::string, AtomProp> table;
    std::string              tmp;
    while (tr.readLine(&tmp))
    {
        if (tmp.find("#") == std::string::npos)
        {
            auto ptr = split(tmp, ',');
            if (ptr.size() == 4)
            {
                AtomProp ap(ptr[1], 
                            my_atoi(ptr[2].c_str(), "atomnumber"),
                            my_atof(ptr[3].c_str(), "mass"));
                table.insert({ ptr[0], ap });
            }
        }
    }
    return table;
}

} // namespace

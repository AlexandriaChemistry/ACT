/*
 * This source file is part of the Alexandria Chemistry Toolkit.
 *
 * Copyright (C) 2014-2020
 *
 * Developers:
 *             Mohammad Mehdi Ghahremanpour,
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
 * \author Mohammad Mehdi Ghahremanpour <mohammad.ghahremanpour@icm.uu.se>
 * \author David van der Spoel <david.vanderspoel@icm.uu.se>
 */

#include "interactiontype.h"

#include <map>

#include "gromacs/utility/exceptions.h"
#include "gromacs/utility/stringutil.h"

namespace alexandria
{

std::map<InteractionType, std::string> eitNames = {
    { eitBONDS,              "BONDS"              },
    { eitANGLES,             "ANGLES"             },
    { eitLINEAR_ANGLES,      "LINEAR_ANGLES"      },
    { eitPROPER_DIHEDRALS,   "PROPER_DIHEDRALS"   },
    { eitIMPROPER_DIHEDRALS, "IMPROPER_DIHEDRALS" },
    { eitVDW,                "VANDERWAALS"        },
    { eitLJ14,               "LJ14"               },
    { eitPOLARIZATION,       "POLARIZATION"       },
    { eitCONSTR,             "CONSTR"             },
    { eitVSITE2,             "VSITE2"             },
    { eitVSITE3FAD,          "VSITE3FAD"          },
    { eitVSITE3OUT,          "VSITE3OUT"          },
    { eitBONDCORRECTIONS,    "BONDCORRECTIONS"    },
    { eitELECTRONEGATIVITYEQUALIZATION, "ELECTRONEGATIVITYEQUALIZATION" }
};

const std::string &interactionTypeToString(InteractionType iType)
{
    auto en = eitNames.find(iType);
    if (en == eitNames.end())
    {
        GMX_THROW(gmx::InternalError(gmx::formatString("No string corresponding to interaction type %d",
                                                       static_cast<int>(iType)).c_str()));
    }
    return en->second;
}

InteractionType stringToInteractionType(const std::string &name)
{
    for (auto &eit : eitNames)
    {
        if (name == eit.second)
        {
            return eit.first;
        }
    }
    GMX_THROW(gmx::InvalidInputError(gmx::formatString("No such interaction type %s",
                                                       name.c_str()).c_str()));
    return eitBONDS;
}

} // namespace alexandria

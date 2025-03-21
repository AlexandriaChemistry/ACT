/*
 * This source file is part of the Alexandria Chemistry Toolkit.
 *
 * Copyright (C) 2022
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

#include "dataset.h"

    
std::map<iMolSelect, const char *> MolSelect_Names = {
    { iMolSelect::Train,   "Train"   },
    { iMolSelect::Test,    "Test"    },
    { iMolSelect::Ignore,  "Ignore"  }
};

const std::map<iMolSelect, const char *> &iMolSelectNames()
{
    return MolSelect_Names;
}

const char *iMolSelectName(iMolSelect ims)
{
    return MolSelect_Names[ims];
}

bool name2molselect(const std::string &name, iMolSelect *ims)
{
    for (auto iter = MolSelect_Names.begin(); iter != MolSelect_Names.end(); ++iter)
    {
        if (name.compare(iter->second) == 0)
        {
            *ims = iter->first;
            return true;
        }
    }
    return false;
}

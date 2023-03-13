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
#include "allmols.h"

#include "act/utility/stringutil.h"
#include "gromacs/utility/stringutil.h"
#include "gromacs/utility/textreader.h"

namespace alexandria
{

AlexandriaMol::AlexandriaMol(const std::vector<std::string> &line)
{
    if (line.size() >= 20)
    {
        iupac    = line[2];
        formula  = line[3];
        charge   = my_atoi(line[4].c_str(), "charge");
        mult     = my_atoi(line[5].c_str(), "multiplicity");
        mass     = my_atof(line[6].c_str(), "mass");
        cas      = line[7];
        csid     = my_atoi(line[8].c_str(), "chemspider");
        pubid    = my_atoi(line[9].c_str(), "pubchem");
        inchi    = line[10];
        inchikey = line[11];
    }
}

AlexandriaMols::AlexandriaMols()
{
    std::string alexandria;
    auto        actdata = getenv("ACTDATA");
    if (nullptr == actdata)
    {
        alexandria = "alexandria.csv";
    }
    else
    {
        alexandria = gmx::formatString("%s/alexandria.csv", actdata);
    }
    gmx::TextReader tr(alexandria);
    std::string line;
    while (tr.readLine(&line))
    {
        auto words = gmx::splitDelimitedString(line, '|');
        if (words.size() >= 20)
        {
            AlexandriaMol am(words);
            if (!am.inchi.empty())
            {
                mols_.insert({am.inchi, am});
            }
        }
    }
}
      
const AlexandriaMol *AlexandriaMols::find(const std::string &inchi) const
{
    auto mptr = mols_.find(inchi);
    if (mols_.end() == mptr)
    {
        return nullptr;
    }
    return &mptr->second;
}

}

/*
 * This source file is part of the Alexandria Chemistry Toolkit.
 *
 * Copyright (C) 2023-2025
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

#include <algorithm>
#include <string>

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
        for(auto ss : split(line[12], ';'))
        {
            classid.insert(ss);
            std::replace(ss.begin(), ss.end(), ' ', '-');
            classid.insert(ss);
        }
        // Synonyms must include the iupac as well in case we are
        // looking that one up later.
        synonyms.insert(iupac);
        for(auto ss : split(line[13], ';'))
        {
            if (ss.empty())
            {
                continue;
            }
            synonyms.insert(ss);
            std::replace(ss.begin(), ss.end(), ' ', '-');
            synonyms.insert(ss);
        }
        // Add the filename as well
        auto fn = line[14].substr(0, line[14].size()-4);
        if (fn.size() > 0)
        {
            synonyms.insert(fn);
        }
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
                mols_.insert( { am.inchi, am } );
                for(const auto &syn : am.synonyms)
                {
                    nameToInChi_.insert( { syn, am.inchi } );
                }
            }
        }
    }
}

void AlexandriaMols::dump(FILE *fp)
{
    if (fp)
    {
        fprintf(fp, "There are %zu molecules and %zu synonyms\n", mols_.size(), nameToInChi_.size());
        for(const auto &n2i : nameToInChi_)
        {
            fprintf(fp, "%s %s\n", n2i.first.c_str(), n2i.second.c_str());
        }
    }
}
      
const AlexandriaMol *AlexandriaMols::findInChi(const std::string &inchi) const
{
    if (inchi.size() == 0)
    {
        return nullptr;
    }
    auto mptr = mols_.find(inchi);
    if (mols_.end() == mptr)
    {
        return nullptr;
    }
    return &mptr->second;
}

const AlexandriaMol *AlexandriaMols::findMol(const std::string &name) const
{
    auto n2i = nameToInChi_.find(name);
    if (n2i == nameToInChi_.end())
    {
        return nullptr;
    }
    auto mptr = mols_.find(n2i->second);
    if (mols_.end() == mptr)
    {
        return nullptr;
    }
    return &mptr->second;
}

}

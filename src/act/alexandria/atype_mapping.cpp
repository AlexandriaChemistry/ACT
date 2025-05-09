/*
 * This source file is part of the Alexandria Chemistry Toolkit.
 *
 * Copyright (C) 2021,2024
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

#include "actpre.h"

#include "atype_mapping.h"

#include <map>

#include "gromacs/utility/futil.h"
#include "gromacs/utility/stringutil.h"
#include "gromacs/utility/textreader.h"
    
void gaffToAlexandria(const std::string                  &filenm,
                      std::map<std::string, std::string> *g2a,
                      bool                                oneH)
{
    g2a->clear();
    std::string myfilenm(filenm);
    if (myfilenm.empty())
    {
        myfilenm.assign(gmx::findLibraryFile("atomtype_mapping.csv", false, true));
    }
    gmx::TextReader  ttt(myfilenm);
    std::string      line;
    while (ttt.readLine(&line))
    {
        auto words = split(line, '|');
        if (words.size() == 3)
        {
            if (g2a->find(words[0]) == g2a->end())
            {
                std::string hh("h");
                std::string w1 = words[1];
                if (oneH && w1.find(hh) == 0)
                {
                    w1.assign("h");
                }
                g2a->insert(std::pair<std::string, std::string>(words[0], w1));
            }
            else
            {
                gmx_fatal(FARGS, "Duplicate mapping entries in %s for %s",
                          filenm.c_str(), words[0].c_str());
            }
        }
        else
        {
            gmx_fatal(FARGS, "Invalid line in %s: %s", filenm.c_str(), line.c_str());
        }
    }
}

bool renameAtomTypes(alexandria::MolProp                      *mp,
                     const std::map<std::string, std::string> &g2a)
{
    auto expers = mp->experiment();
    for(auto myexp = expers->begin(); myexp < expers->end(); ++myexp)
    {
        auto catoms = myexp->calcAtom();
        for(auto catom = catoms->begin(); catom < catoms->end(); ++catom)
        {
            auto obt = catom->getObtype();
            if (g2a.find(obt) == g2a.end())
            {
                fprintf(stderr, "obType %s cannot be mapped to alexandria", obt.c_str());
                return false;
            }
            else
            {
                if (debug)
                {
                    fprintf(debug, "Replacing %s by %s\n", obt.c_str(),
                            g2a.find(obt)->second.c_str());
                }
                catom->setObtype(g2a.find(obt)->second);
            }
        }
    }
    return true;
}


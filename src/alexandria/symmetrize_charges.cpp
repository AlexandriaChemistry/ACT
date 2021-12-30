/*
 * This source file is part of the Alexandria Chemistry Toolkit.
 *
 * Copyright (C) 2014-2021
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

#include "symmetrize_charges.h"

#include <string>

namespace alexandria
{

void symmetrize_charges(bool                                bQsym, 
                        const t_atoms                      *atoms,
                        const std::vector<TopologyEntry *> &bonds,
                        const Poldata                      *pd,
                        const char                         *symm_string,
                        std::vector<int>                   *sym_charges)
{
    std::string  central, attached;
    int          nrq;
    double       qaver, qsum;

    sym_charges->clear();
    for (int i = 0; i < atoms->nr; i++)
    {
        sym_charges->push_back(i);
    }
    if (bQsym)
    {
        if ((nullptr != symm_string) && (strlen(symm_string) > 0))
        {
            std::vector<std::string> ss = gmx::splitString(symm_string);
            if (static_cast<int>(ss.size()) != atoms->nr)
            {
                gmx_fatal(FARGS, "Wrong number (%d) of atom-numbers in symm_string: expected %d",
                          static_cast<int>(ss.size()), atoms->nr);
            }
            int ii = 0;
            for (auto is = ss.begin();
                 (is < ss.end()); ++is)
            {
                (*sym_charges)[ii] = atoi(is->c_str());
                ii++;
            }
        }
        else
        {
            for (auto symcharges = pd->getSymchargesBegin();
                 symcharges != pd->getSymchargesEnd(); symcharges++)
            {
                for (int i = 0; i < atoms->nr; i++)
                {
                    if (symcharges->getCentral().compare(atoms->atom[i].elem) == 0)
                    {
                        int              hsmin = -1;
                        std::vector<int> hs;
                        for (auto &jj : bonds)
                        {
                            auto j  = static_cast<const Bond *>(jj);
                            auto ai = j->aI();
                            auto aj = j->aJ();
                            
                            if (ai == i && 
                                symcharges->getAttached().compare(atoms->atom[aj].elem) == 0)
                            {
                                hs.push_back(aj);
                            }
                            else if (aj == i &&
                                     symcharges->getAttached().compare(atoms->atom[ai].elem) == 0)
                            {
                                hs.push_back(ai);
                            }
                            if ((hs.size() > 0) && (hsmin == -1 || hs.back() < hsmin))
                            {
                                hsmin = hs.back();
                            }
                        }
                        if ((static_cast<int>(hs.size()) == symcharges->getNumattach()) &&
                            (hsmin != -1))
                        {
                            for (int k = 0; k < symcharges->getNumattach(); k++)
                            {
                                (*sym_charges)[hs[k]] = hsmin;
                            }
                        }
                    }
                }
            }
        }

        for (int i = 0; i < atoms->nr; i++)
        { 
            qsum = 0;
            nrq  = 0;
            for (int j = i+1; j < atoms->nr; j++)
            {
                if ((*sym_charges)[j] == (*sym_charges)[i])
                {
                    qsum += atoms->atom[j].q;
                    nrq++;
                }
            }
            if (0 < nrq)
            {
                qaver = qsum/nrq;
                for (int j = 0; j < atoms->nr; j++)
                {
                    if ((*sym_charges)[j] == (*sym_charges)[i])
                    {
                        atoms->atom[j].q = atoms->atom[j].qB = qaver;
                    }
                }
            }
        }
    }
}

} // namespace alexandria

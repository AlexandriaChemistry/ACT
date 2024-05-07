/*
 * This source file is part of the Alexandria Chemistry Toolkit.
 *
 * Copyright (C) 2014-2020
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
 * \author Mohammad Mehdi Ghahremanpour <mohammad.ghahremanpour@icm.uu.se>
 * \author David van der Spoel <david.vanderspoel@icm.uu.se>
 */

#include "molselect.h"

#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <random>
#include <vector>
#include <strings.h>

#include "gromacs/commandline/filenm.h"
#include "gromacs/commandline/pargs.h"
#include "gromacs/gmxlib/network.h"
#include "gromacs/mdrunutility/mdmodules.h"
#include "gromacs/utility/arraysize.h"
#include "gromacs/utility/coolstuff.h"
#include "gromacs/utility/cstringutil.h"
#include "gromacs/utility/futil.h"
#include "gromacs/utility/smalloc.h"
#include "gromacs/utility/textreader.h"

#include "act/alexandria/actmol.h"
#include "act/alexandria/alex_modules.h"
#include "act/alexandria/molgen.h"
#include "act/molprop/molprop.h"
#include "act/molprop/molprop_xml.h"
#include "act/utility/stringutil.h"

namespace alexandria
{

void MolSelect::addOne(const std::string &iupac,
                       int                index,
                       iMolSelect         ims)
{
    ims_.push_back(IMolSelect(iupac, ims, index));
}

void MolSelect::read(const char *fn)
{
    gmx::TextReader tr(fn);
    std::string     tmp;
    int             index = 0;

    while (tr.readLine(&tmp))
    {
        while (!tmp.empty() && tmp[tmp.length()-1] == '\n')
        {
            tmp.erase(tmp.length()-1);
        }
        auto ptr = split(tmp, '|');
        if ((ptr.size() == 2) && (ptr[0].length() > 1)
            && (ptr[1].length() > 1))
        {
            iMolSelect status;
            if (name2molselect(ptr[1], &status))
            {
                addOne(ptr[0], index++, status);
            }
            else
            {
                fprintf(stderr, "Unknown status '%s' for molecule %s on line %d in file %s\n",
                        ptr[1].c_str(), ptr[0].c_str(), index, fn);
            }
        }
        else
        {
            fprintf(stderr, "Invalid line '%s' in selection file\n",
                    tmp.c_str());
        }
    }
}

bool MolSelect::status(const std::string &iupac, iMolSelect *ims) const
{
    auto imi = std::find_if(ims_.begin(), ims_.end(),
                            [iupac](IMolSelect const &i)
                            {
                                return i.iupac().compare(iupac) == 0;
                            });

    if (imi == ims_.end())
    {
        return false;
    }
    *ims = imi->status();
    return true;
}

bool MolSelect::index(const std::string &iupac, int *index) const
{
    auto imi = std::find_if(ims_.begin(), ims_.end(),
                            [iupac](IMolSelect const &i)
                            {
                                return i.iupac().compare(iupac) == 0;
                            });
    if (imi == ims_.end())
    {
        return false;
    }
    *index = imi->index();
    return true;
}

void MolSelect::bcast(const CommunicationRecord *cr)
{
    std::map<iMolSelect, int> mmm = {
        { iMolSelect::Train,  13 },
        { iMolSelect::Test,   23 },
        { iMolSelect::Ignore, 37 }
    };
    std::map<int, iMolSelect> rmm;
    for(const auto &m : mmm)
    {
        rmm[m.second] = m.first;
    }

    if (cr->isMaster())
    {
        int isize = ims_.size();
        cr->bcast(&isize, cr->comm_world());
        for(int i = 0; i < isize; i++)
        {
            std::string iupac = ims_[i].iupac();
            cr->bcast(&iupac, cr->comm_world());
            int index = ims_[i].index();
            cr->bcast(&index, cr->comm_world());
            int ims = mmm[ims_[i].status()];
            cr->bcast(&ims, cr->comm_world());
        }
    }
    else
    {
        int isize;
        cr->bcast(&isize, cr->comm_world());
        for(int i = 0; i < isize; i++)
        {
            std::string iupac;
            cr->bcast(&iupac, cr->comm_world());
            int         ims, index;
            cr->bcast(&index, cr->comm_world());
            cr->bcast(&ims, cr->comm_world());
            iMolSelect  status = rmm[ims];
            addOne(iupac, index, status);
        }
    }
}

}

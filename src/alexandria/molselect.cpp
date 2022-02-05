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

#include "alex_modules.h"
#include "molgen.h"
#include "molprop/molprop.h"
#include "molprop/molprop_xml.h"
#include "mymol.h"
#include "poldata/poldata.h"
#include "poldata/poldata_xml.h"
#include "utility/stringutil.h"

namespace alexandria
{

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
                ims_.push_back(IMolSelect(ptr[0], status, index++));
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

}


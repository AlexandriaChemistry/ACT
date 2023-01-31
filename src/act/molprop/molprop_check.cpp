/*
 * This source file is part of the Alexandria Chemistry Toolkit.
 *
 * Copyright (C) 2014-2023
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

#include <cstdio>
#include <cstdlib>

#include <map>
#include <string>
#include <vector>

#include "act/alexandria/alex_modules.h"
#include "act/alexandria/fill_inputrec.h"
#include "act/alexandria/actmol.h"
#include "act/molprop/molprop.h"
#include "act/molprop/molprop_util.h"
#include "act/molprop/molprop_xml.h"
#include "act/forcefield/forcefield_xml.h"
#include "act/utility/units.h"
#include "gromacs/commandline/pargs.h"
#include "gromacs/fileio/gmxfio.h"
#include "gromacs/math/functions.h"
#include "gromacs/mdtypes/inputrec.h"
#include "gromacs/utility/strconvert.h"

namespace alexandria
{

int molprop_check(int argc, char*argv[])
{
    static const char               *desc[] = {
    };
    t_filenm                         fnm[] =
    {
        { efXML, "-ff",    "aff",           ffREAD },
        { efXML, "-mp",    "allmols",       ffREAD },
        { efXML, "-mpout", "newmols",       ffOPTWR },
        { efLOG, "-g",     "molprop_check", ffWRITE }
    };
    int NFILE = (sizeof(fnm)/sizeof(fnm[0]));
    bool compress = false;
    
    std::vector<alexandria::MolProp> mp;
    gmx_output_env_t                *oenv;
    
    if (!parse_common_args(&argc, argv, PCA_NOEXIT_ON_ARGS, NFILE, fnm,
                           0, nullptr,
                           sizeof(desc)/sizeof(desc[0]), desc,
                           0, nullptr, &oenv))
    {
        return 0;
    }
    MolPropRead(opt2fn("-mp", NFILE, fnm), &mp);
    
    return 0;
}

} // namespace alexandria

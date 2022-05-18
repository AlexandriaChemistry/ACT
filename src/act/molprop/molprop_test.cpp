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

#include <stdio.h>
#include <stdlib.h>

#include "gromacs/commandline/pargs.h"

#include "alexandria/alex_modules.h"
#include "act/molprop/molprop.h"
#include "act/molprop/molprop_xml.h"
#include "alexandria/mymol.h"
#include "act/poldata/poldata_xml.h"

namespace alexandria
{

int molprop_test(int argc, char*argv[])
{
    static const char               *desc[] = {
        "molprop_test reads a molprop file and writes a new one."
    };
    gmx_output_env_t                *oenv;
    std::vector<alexandria::MolProp> mpt;
    t_filenm                         fnm[] = {
        { efXML, "-mp", "molin",     ffREAD },
        { efXML, "-o", "molout",    ffWRITE }
    };
#define NFILE sizeof(fnm)/sizeof(fnm[0])
    int compress = 1;
    std::vector<t_pargs> pa  = {
        { "-compress", FALSE, etBOOL, {&compress},
          "Compress the output file" }
    };
    if (!parse_common_args(&argc, argv, 0, NFILE, fnm, pa.size(), pa.data(),
                           sizeof(desc)/sizeof(desc[0]), desc, 0, nullptr, &oenv))
    {
        return 0;
    }

    MolPropRead(opt2fn("-mp", NFILE, fnm), &mpt);
    printf("Read %d molecules from %s\n", (int)mpt.size(), opt2fn("-mp", NFILE, fnm));
    MolPropWrite(opt2fn("-o", NFILE, fnm), mpt, compress);

    return 0;
}

} // namespace alexandria

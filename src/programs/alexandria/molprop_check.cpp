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

#include <cstdio>
#include <cstdlib>

#include <map>
#include <string>
#include <vector>

#include "gromacs/commandline/pargs.h"
#include "gromacs/utility/strconvert.h"

#include "alex_modules.h"
#include "molprop.h"
#include "molprop_xml.h"
#include "poldata_xml.h"

int alex_molprop_check(int argc, char*argv[])
{
    static const char               *desc[] = {
        "molprop_check checks calculations for missing hydrogens",
        "and inconsistent dipoles."
    };
    t_filenm                         fnm[] =
    {
        { efDAT, "-f",  "allmols",  ffREAD }
    };
    int                              NFILE   = (sizeof(fnm)/sizeof(fnm[0]));
    std::vector<alexandria::MolProp> mp;
    gmx_output_env_t                *oenv;

    if (!parse_common_args(&argc, argv, PCA_NOEXIT_ON_ARGS, NFILE, fnm,
                           0, nullptr,
                           sizeof(desc)/sizeof(desc[0]), desc,
                           0, nullptr, &oenv))
    {
        return 0;
    }
    MolPropRead(opt2fn("-f", NFILE, fnm), &mp);

    for (auto &m : mp)
    {
        typedef struct
        {
            std::string name;
            rvec        mu;
        } name_mu;
        
        std::vector<name_mu> mus;
        for (auto &ci : m.experimentConst())
        {
            int nH = 0, nC = 0;
            for (auto &cai : ci.calcAtomConst())
            {
                std::string name = cai.getName();
                if (name.compare("H") == 0)
                {
                    nH++;
                }
                else if (name.compare("C") == 0)
                {
                    nC++;
                }
            }
            if (nC > 0 && nH == 0)
            {
                printf("%s #C %d #H %d\n",
                       ci.getDatafile().c_str(),
                       nC, nH);
            }
            rvec   mu;
            tensor Q;
            double value, error, T;
            std::string type;
            if (ci.getVal(type, MPO_DIPOLE, &value, &error,
                          &T, mu, Q))
            {
                name_mu nmu = { ci.getDatafile(), { mu[XX], mu[YY], mu[ZZ] } };
                mus.push_back(nmu);
            }
        }
        // Check dipoles
        for(const auto &mi : mus)
        {
            printf("%s %s %.2f %.2f %.2f\n", m.getMolname().c_str(),
                   mi.name.c_str(), mi.mu[XX], mi.mu[YY], mi.mu[ZZ]);
        }
    }

    return 0;
}

/*
 * This source file is part of the Alexandria Chemistry Toolkit.
 *
 * Copyright (C) 2014-2024
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

#include "actpre.h"

#include <cstdlib>

#include "act/alexandria/alex_modules.h"
#include "act/basics/version.h"
#include "gromacs/commandline/cmdlineinit.h"
#include "gromacs/commandline/cmdlinemodulemanager.h"
#include "gromacs/gmxlib/network.h"
#include "gromacs/mdtypes/commrec.h"
#include "gromacs/utility/binaryinformation.h"
#include "gromacs/utility/exceptions.h"
#include "gromacs/utility/init.h"
#include "gromacs/utility/smalloc.h"

int
main(int argc, char *argv[])
{
    gmx::CommandLineProgramContext &context = gmx::initForCommandLine(&argc, &argv);
    try
    {
        t_commrec *cr = init_commrec();
        gmx::CommandLineModuleManager manager("alexandria", &context);
        registerAlexandriaModules(&manager);
        manager.setQuiet(true);
        setenv("GMX_NB_GENERIC", "1", 1);
        if (MASTER(cr))
        {
            printf("%s", act_welcome().c_str());
        }
        int rc = manager.run(argc, argv);
        gmx::finalizeForCommandLine();
        if (MASTER(cr))
        {
            printf("%s", act_goodbye().c_str());
        }
        else
        {
            // Cleanup
            for(int i = 0; i < argc; i++)
            {
                if (argv[i])
                {
                    sfree(argv[i]);
                }
            }
        }
        done_commrec(cr);
        return rc;
    }
    catch (const std::exception &ex)
    {
        gmx::printFatalErrorMessage(stderr, ex);
        return gmx::processExceptionAtExit(ex);
    }
}

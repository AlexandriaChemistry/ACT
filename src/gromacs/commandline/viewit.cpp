/*
 * This file is part of the GROMACS molecular simulation package.
 *
 * Copyright (c) 1991-2000, University of Groningen, The Netherlands.
 * Copyright (c) 2001-2004, The GROMACS development team.
 * Copyright (c) 2013,2014,2015,2017,2018, by the GROMACS development team, led by
 * Mark Abraham, David van der Spoel, Berk Hess, and Erik Lindahl,
 * and including many others, as listed in the AUTHORS file in the
 * top-level source directory and at http://www.gromacs.org.
 *
 * GROMACS is free software; you can redistribute it and/or
 * modify it under the terms of the GNU Lesser General Public License
 * as published by the Free Software Foundation; either version 2.1
 * of the License, or (at your option) any later version.
 *
 * GROMACS is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
 * Lesser General Public License for more details.
 *
 * You should have received a copy of the GNU Lesser General Public
 * License along with GROMACS; if not, see
 * http://www.gnu.org/licenses, or write to the Free Software Foundation,
 * Inc., 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301  USA.
 *
 * If you want to redistribute modifications to GROMACS, please
 * consider that scientific software is very special. Version
 * control is crucial - bugs must be traceable. We will be happy to
 * consider code for inclusion in the official distribution, but
 * derived work must not be called official GROMACS. Details are found
 * in the README & COPYING files - if they are missing, get the
 * official version at http://www.gromacs.org.
 *
 * To help us fund GROMACS development, we humbly ask that you cite
 * the research papers on the package. Check out http://www.gromacs.org.
 */
/* This file is completely threadsafe - keep it that way! */
#include "actpre.h"

#include "viewit.h"

#include <cstdlib>
#include <cstring>

#include <algorithm>

#include "gromacs/commandline/filenm.h"
#include "gromacs/fileio/oenv.h"
#include "gromacs/utility/arraysize.h"
#include "gromacs/utility/cstringutil.h"
#include "gromacs/utility/fatalerror.h"

static const int   can_view_ftp[] = {
    0,
    efXVG,          efPDB
};
#define NVIEW asize(can_view_ftp)
static const char* view_program[] = {
    nullptr, 
    nullptr,           "xterm -e rasmol"
};

static int can_view(int ftp)
{
    int i;

    for (i = 1; i < NVIEW; i++)
    {
        if (ftp == can_view_ftp[i])
        {
            return i;
        }
    }

    return 0;
}

void do_view(const gmx_output_env_t *oenv, const char *fn, const char *opts)
{
    const char *cmd;
    int         ftp, n;

    if (output_env_get_view(oenv) && fn)
    {
        if (getenv("DISPLAY") == nullptr)
        {
            fprintf(stderr, "Can not view %s, no DISPLAY environment variable.\n", fn);
        }
        else
        {
            ftp = fn2ftp(fn);
            auto env = gmx::formatString("GMX_VIEW_%s", ftp2ext(ftp));
            // Convert to upper case
            std::transform(env.begin(), env.end(), env.begin(), ::toupper);
            switch (ftp)
            {
                case efXVG:
                    if (!(cmd = getenv(env.c_str())) )
                    {
                        if (getenv("GMX_USE_XMGR") )
                        {
                            cmd = "xmgr";
                        }
                        else
                        {
                            cmd = "xmgrace";
                        }
                    }
                    break;
                default:
                    if ( (n = can_view(ftp)) )
                    {
                        if (!(cmd = getenv(env.c_str())) )
                        {
                            cmd = view_program[n];
                        }
                    }
                    else
                    {
                        fprintf(stderr, "Don't know how to view file %s", fn);
                        return;
                    }
            }
            if (strlen(cmd) )
            {
                auto buf = gmx::formatString("%s %s %s &", cmd, opts ? opts : "", fn);
                fprintf(stderr, "Executing '%s'\n", buf.c_str());
                if (0 != system(buf.c_str()) )
                {
                    gmx_fatal(FARGS, "Failed executing command: %s", buf.c_str());
                }
            }
        }
    }
}

void view_all(const gmx_output_env_t *oenv, int nf, t_filenm fnm[])
{
    int i;

    for (i = 0; i < nf; i++)
    {
        if (can_view(fnm[i].ftp) && is_output(&(fnm[i])) &&
            ( !is_optional(&(fnm[i])) || is_set(&(fnm[i])) ) )
        {
            do_view(oenv, fnm[i].filenames[0].c_str(), nullptr);
        }
    }
}

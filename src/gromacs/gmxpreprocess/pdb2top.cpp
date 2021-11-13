/*
 * This file is part of the GROMACS molecular simulation package.
 *
 * Copyright (c) 1991-2000, University of Groningen, The Netherlands.
 * Copyright (c) 2001-2004, The GROMACS development team.
 * Copyright (c) 2013,2014,2015,2016,2017,2018, by the GROMACS development team, led by
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
#include "actpre.h"

#include "pdb2top.h"

#include <cctype>
#include <cmath>
#include <cstdio>
#include <cstring>

#include <algorithm>
#include <string>
#include <vector>

#include "gromacs/fileio/pdbio.h"
#include "gromacs/gmxpreprocess/fflibutil.h"
#include "gromacs/gmxpreprocess/gen_ad.h"
#include "gromacs/gmxpreprocess/gpp_nextnb.h"
#include "gromacs/gmxpreprocess/notset.h"
#include "gromacs/gmxpreprocess/pgutil.h"
#include "gromacs/gmxpreprocess/resall.h"
#include "gromacs/gmxpreprocess/topdirs.h"
#include "gromacs/gmxpreprocess/toputil.h"
#include "gromacs/math/functions.h"
#include "gromacs/math/vec.h"
#include "gromacs/topology/residuetypes.h"
#include "gromacs/topology/symtab.h"
#include "gromacs/utility/binaryinformation.h"
#include "gromacs/utility/cstringutil.h"
#include "gromacs/utility/dir_separator.h"
#include "gromacs/utility/exceptions.h"
#include "gromacs/utility/fatalerror.h"
#include "gromacs/utility/futil.h"
#include "gromacs/utility/niceheader.h"
#include "gromacs/utility/path.h"
#include "gromacs/utility/programcontext.h"
#include "gromacs/utility/smalloc.h"
#include "gromacs/utility/strconvert.h"
#include "gromacs/utility/strdb.h"
#include "gromacs/utility/stringutil.h"
#include "gromacs/utility/textwriter.h"

bool is_int(double x)
{
    const double tol = 1e-4;
    int          ix;

    if (x < 0)
    {
        x = -x;
    }
    ix = gmx::roundToInt(x);

    return (fabs(x-ix) < tol);
}

static void print_top_heavy_H(FILE *out, real mHmult)
{
    if (mHmult == 2.0)
    {
        fprintf(out, "; Using deuterium instead of hydrogen\n\n");
    }
    else if (mHmult == 4.0)
    {
        fprintf(out, "#define HEAVY_H\n\n");
    }
    else if (mHmult != 1.0)
    {
        fprintf(stderr, "WARNING: unsupported proton mass multiplier (%g) "
                "in pdb2top\n", mHmult);
    }
}

void print_top_comment(FILE       *out,
                       const char *filename,
                       const char *ffdir,
                       bool        bITP,
                       const char *remark)
{
    char  ffdir_parent[STRLEN];
    char *p;

    try
    {
        gmx::TextWriter writer(out);
        gmx::niceHeader(&writer, filename, ';');
        writer.writeLine(gmx::formatString(";\tThis is a %s topology file", bITP ? "include" : "standalone"));
        writer.writeLine(";");

        gmx::BinaryInformationSettings settings;
        settings.generatedByHeader(true);
        settings.linePrefix(";\t");
        gmx::printBinaryInformation(&writer, gmx::getProgramContext(), settings);
    }
    GMX_CATCH_ALL_AND_EXIT_WITH_FATAL_ERROR;

    if ((NULL != remark) && (strlen(remark) > 0))
    {
        fprintf(out, "; %s\n", remark);
    }

    if (strchr(ffdir, '/') == nullptr)
    {
        fprintf(out, ";\tForce field was read from the standard GROMACS share directory.\n;\n\n");
    }
    else if (ffdir[0] == '.')
    {
        fprintf(out, ";\tForce field was read from current directory or a relative path - path added.\n;\n\n");
    }
    else
    {
        strncpy(ffdir_parent, ffdir, STRLEN-1);
        ffdir_parent[STRLEN-1] = '\0'; /*make sure it is 0-terminated even for long string*/
        p = strrchr(ffdir_parent, '/');

        *p = '\0';

        fprintf(out,
                ";\tForce field data was read from:\n"
                ";\t%s\n"
                ";\n"
                ";\tNote:\n"
                ";\tThis might be a non-standard force field location. When you use this topology, the\n"
                ";\tforce field must either be present in the current directory, or the location\n"
                ";\tspecified in the GMXLIB path variable or with the 'include' mdp file option.\n;\n\n",
                ffdir_parent);
    }
}

void print_top_header(FILE *out, const char *filename,
                      bool bITP, const char *ffdir, real mHmult,
                      const char *remark)
{
    const char *p;

    print_top_comment(out, filename, ffdir, bITP, remark);

    print_top_heavy_H(out, mHmult);
    fprintf(out, "; Include forcefield parameters\n");

    p = strrchr(ffdir, '/');
    p = (ffdir[0] == '.' || p == nullptr) ? ffdir : p+1;

    fprintf(out, "#include \"%s/%s\"\n\n", p, fflib_forcefield_itp());
}

static void print_top_posre(FILE *out, const char *pr)
{
    fprintf(out, "; Include Position restraint file\n");
    fprintf(out, "#ifdef POSRES\n");
    fprintf(out, "#include \"%s\"\n", pr);
    fprintf(out, "#endif\n\n");
}

static void print_top_water(FILE *out, const char *ffdir, const char *water)
{
    const char *p;
    char        buf[STRLEN];

    fprintf(out, "; Include water topology\n");

    p = strrchr(ffdir, '/');
    p = (ffdir[0] == '.' || p == nullptr) ? ffdir : p+1;
    fprintf(out, "#include \"%s/%s.itp\"\n", p, water);

    fprintf(out, "\n");
    fprintf(out, "#ifdef POSRES_WATER\n");
    fprintf(out, "; Position restraint for each water oxygen\n");
    fprintf(out, "[ position_restraints ]\n");
    fprintf(out, ";%3s %5s %9s %10s %10s\n", "i", "funct", "fcx", "fcy", "fcz");
    fprintf(out, "%4d %4d %10g %10g %10g\n", 1, 1, 1000.0, 1000.0, 1000.0);
    fprintf(out, "#endif\n");
    fprintf(out, "\n");

    sprintf(buf, "%s/ions.itp", p);

    if (fflib_fexist(buf))
    {
        fprintf(out, "; Include topology for ions\n");
        fprintf(out, "#include \"%s\"\n", buf);
        fprintf(out, "\n");
    }
}

static void print_top_system(FILE *out, const char *title)
{
    fprintf(out, "[ %s ]\n", dir2str(d_system));
    fprintf(out, "; Name\n");
    fprintf(out, "%s\n\n", title[0] ? title : "Protein");
}

void print_top_mols(FILE *out,
                    const char *title, const char *ffdir, const char *water,
                    int nincl, char **incls, int nmol, t_mols *mols)
{

    if (nincl > 0)
    {
        fprintf(out, "; Include chain topologies\n");
        for (int i = 0; i < nincl; i++)
        {
            fprintf(out, "#include \"%s\"\n", gmx::Path::getFilename(incls[i]).c_str());
        }
        fprintf(out, "\n");
    }

    if (water)
    {
        print_top_water(out, ffdir, water);
    }
    print_top_system(out, title);

    if (nmol)
    {
        fprintf(out, "[ %s ]\n", dir2str(d_molecules));
        fprintf(out, "; %-15s %5s\n", "Compound", "#mols");
        for (int i = 0; i < nmol; i++)
        {
            fprintf(out, "%-15s %5d\n", mols[i].name, mols[i].nr);
        }
    }
}

void write_top(FILE *out, const char *pr, const char *molname,
               t_atoms *at, bool bRTPresname,
               int bts[], t_params plist[], t_excls excls[],
               gpp_atomtype_t atype, int *cgnr, int nrexcl)
/* NOTE: nrexcl is not the size of *excl! */
{
    if (at && atype && cgnr)
    {
        fprintf(out, "[ %s ]\n", dir2str(d_moleculetype));
        fprintf(out, "; %-15s %5s\n", "Name", "nrexcl");
        fprintf(out, "%-15s %5d\n\n", molname ? molname : "Protein", nrexcl);

        print_atoms(out, atype, at, cgnr, bRTPresname);
        print_bondeds(out, at->nr, d_bonds,      F_BONDS,    bts[ebtsBONDS], plist);
        print_bondeds(out, at->nr, d_constraints, F_CONSTR,   0,              plist);
        print_bondeds(out, at->nr, d_constraints, F_CONSTRNC, 0,              plist);
        print_bondeds(out, at->nr, d_pairs,      F_LJ14,     0,              plist);
        print_excl(out, at->nr, excls);
        print_bondeds(out, at->nr, d_angles,     F_ANGLES,   bts[ebtsANGLES], plist);
        print_bondeds(out, at->nr, d_dihedrals,  F_PDIHS,    bts[ebtsPDIHS], plist);
        print_bondeds(out, at->nr, d_dihedrals,  F_IDIHS,    bts[ebtsIDIHS], plist);
        print_bondeds(out, at->nr, d_cmap,       F_CMAP,     bts[ebtsCMAP],  plist);
        print_bondeds(out, at->nr, d_polarization, F_POLARIZATION,   0,       plist);
        print_bondeds(out, at->nr, d_thole_polarization, F_THOLE_POL, 0,       plist);
        print_bondeds(out, at->nr, d_vsites2,    F_VSITE2,   0,              plist);
        print_bondeds(out, at->nr, d_vsites3,    F_VSITE3,   0,              plist);
        print_bondeds(out, at->nr, d_vsites3,    F_VSITE3FD, 0,              plist);
        print_bondeds(out, at->nr, d_vsites3,    F_VSITE3FAD, 0,              plist);
        print_bondeds(out, at->nr, d_vsites3,    F_VSITE3OUT, 0,              plist);
        print_bondeds(out, at->nr, d_vsites4,    F_VSITE4FD, 0,              plist);
        print_bondeds(out, at->nr, d_vsites4,    F_VSITE4FDN, 0,             plist);

        if (pr)
        {
            print_top_posre(out, pr);
        }
    }
}




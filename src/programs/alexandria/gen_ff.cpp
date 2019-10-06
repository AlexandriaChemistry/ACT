/*! \internal \brief
 * Implements part of the alexandria program.
 * \author David van der Spoel <david.vanderspoel@icm.uu.se>
 *//*
 * This source file is part of the Alexandria Chemistry Toolkit.
 *
 * Copyright (C) 2014-2019 
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
 
 
#include "gmxpre.h"

#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "gromacs/legacyheaders/macros.h"
#include "gromacs/legacyheaders/oenv.h"

#include "gromacs/commandline/pargs.h"
#include "gromacs/gmxpreprocess/pdb2top.h"
#include "gromacs/gmxpreprocess/topdirs.h"
#include "gromacs/math/units.h"
#include "gromacs/topology/atomprop.h"
#include "gromacs/utility/cstringutil.h"

#include "poldata.h"
#include "poldata_xml.h"
#include "stringutil.h"

void print_atypes_tex(FILE *tp, gmx_poldata_t pd, gmx_atomprop_t aps)
{
    char *elem, *desc, *gt_type, *ptype, *btype, *gt_old;
    char *vdwparams;
    int   atomnumber, nline, npage, nr;
    real  mass;

    fprintf(tp, "\\begin{longtable}{cclccc}\n\\hline\n");
    fprintf(tp, "\\caption{Atom types defined by the Alexandria force field}\n", caption);
    fprintf(tp, "\\label{atypes}\\\\\n\hline\n", label);
    fprintf(tp, "Nr. & Type & Description & Elem & $\\alpha$ & Van der Waals\\\\\n");
    fprintf(tp, "\\hline\n");

    gt_old = nullptr;
    nline  = 2;
    npage  = 0;
    nr     = 1;
    while (1 == gmx_poldata_get_atype(pd,
                                      &elem,
                                      &desc,
                                      &gt_type,
                                      &ptype,
                                      &btype,
                                      &vdwparams))
    {
        if (gmx_atomprop_query(aps, epropMass, "", elem, &mass))
        {
            if ((nullptr == gt_old) || (strcmp(gt_old, gt_type) != 0))
            {
                fprintf(tp, "%d & %s & %s & %s & %s & %s & %s\\\\\n",
                        nr++, gt_type, desc, elem, ptype, btype, vdwparams);
                nline++;
            }
        }
        gt_old = gt_type;
    }
    fprintf(tp, "\\hline\n\\end{longtable}\n\n");
}

static void do_brule(FILE *tp, gmx_poldata_t pd)
{
    int    nline, npage, numbonds, nr, bAromatic;
    char   hbuf[1024], colinfo[1024], pbuf[1024];
    double valence;
    char  *rule, *gt_type, *geometry, *neighbors;

    /* Bondtypes */
    strcpy(colinfo, "cccccccc");
    strcpy(hbuf, "Nr. & Rule & Type & Geometry & Valence & Aromatic & \\# Bonds & Neighbors");
    begin_table(tp, "Bonding rules defined in the Alexandria force field to determine atom types.",
                "brules", colinfo, hbuf);

    nline = 2;
    npage = 0;
    nr    = 1;
    while (0 != gmx_poldata_get_bonding_rule(pd, &rule, &gt_type, &geometry, &numbonds,
                                             &valence, &bAromatic, &neighbors))
    {
        if (0 == (nline % maxline))
        {
            end_table(tp);
            sprintf(pbuf, "brule%d", ++npage);
            fprintf(tp, "\\addtocounter{table}{-1}\n");
            begin_table(tp, "Bonding rules, continued", pbuf, colinfo, hbuf);
            nline = 1;
        }
        fprintf(tp, "%d & %s & %s & %s & %g & %d & %d & %s\\\\\n",
                nr++, rule, gt_type, geometry,
                valence, bAromatic, numbonds, neighbors);
        nline++;
    }
    end_table(tp);
}

static int inverse_ifunc_index(directive d, int ftype)
{
    int i;

    for (i = 1; (i < F_NRE); i++)
    {
        int ifi = ifunc_index(d, i);
        if ((ifi == ftype) || (ifi == -1))
        {
            break;
        }
    }
    if ((i < F_NRE) && (i != -1))
    {
        return i;
    }
    else
    {
        fprintf(stderr, "Could not find the right %s type, setting to 1\n",
                dir2str(d));
    }
    return 1;
}

static void do_bad(FILE *fp, FILE *tp, gmx_poldata_t pd)
{
    int    bts[ebtsNR];
    int    nline, npage, nr, ntrain;
    char  *ai, *aj, *ak, *al, *params, *lu;
    char   hbuf[1024], colinfo[1024], pbuf[1024];
    double length, ang, sigma, bondorder;

    lu = gmx_poldata_get_length_unit(pd);
    /* Bondtypes */
    strcpy(colinfo, "ccccccl");
    sprintf(hbuf, "Nr. & i & j & Length (%s) & Ntrain & Bond order & Params", lu);
    begin_table(tp, "Bonds defined in the Alexandria force field. Bond length (standard deviation in parentheses), bond order and Morse potential~\\protect\\cite{Morse29} parameters.",
                (char *)"btypes", colinfo, hbuf);

    fprintf(fp, "\n[ bondtypes ]\n");
    fprintf(fp, "; ; i    j  func       parameters\n");
    nline          = 2;
    npage          = 0;
    nr             = 1;
    bts[ebtsBONDS] = inverse_ifunc_index(d_bonds, gmx_poldata_get_bond_ftype(pd));
    while (0 < gmx_poldata_get_bond(pd, &ai, &aj, &length, &sigma, &ntrain,
                                    &bondorder, &params))
    {
        if (ntrain > 0)
        {
            fprintf(fp, "%-5s  %-5s   %d  %g  %s\n", ai, aj, bts[ebtsBONDS],
                    convert2gmx(length, string2unit(lu)), params);
            if (0 == (nline % maxline))
            {
                end_table(tp);
                sprintf(pbuf, "btype%d", ++npage);
                fprintf(tp, "\\addtocounter{table}{-1}\n");
                begin_table(tp, "Bonds, continued", pbuf, colinfo, hbuf);
                nline = 1;
            }
            fprintf(tp, "%d & %s & %s & %.1f(%.1f) & %d & %g & %s\\\\\n", nr++, ai, aj, length,
                    sigma, ntrain, bondorder, params);
            nline++;
        }
    }
    end_table(tp);

    /* Angletypes */
    strcpy(colinfo, "ccccccl");
    strcpy(hbuf, "Nr. & i & j & k & Angle & Ntrain & Params");
    begin_table(tp, "Angles defined in the Alexandria force field.",
                "angtypes", colinfo, hbuf);

    fprintf(fp, "\n[ angletypes ]\n");
    fprintf(fp, "; ; i    j   k  func       parameters\n");
    nline           = 2;
    npage           = 0;
    nr              = 1;
    bts[ebtsANGLES] = inverse_ifunc_index(d_angles, gmx_poldata_get_angle_ftype(pd));
    while (0 < gmx_poldata_get_angle(pd, &ai, &aj, &ak, &ang, &sigma, &ntrain, &params))
    {
        if (ntrain > 0)
        {
            fprintf(fp, "%-5s  %-5s  %-5s  %d  %g  %s\n", ai, aj, ak, bts[ebtsANGLES], length, params);
            if (0 == (nline % maxline))
            {
                end_table(tp);
                sprintf(pbuf, "angtype%d", ++npage);
                fprintf(tp, "\\addtocounter{table}{-1}\n");
                begin_table(tp, "Angles, continued", pbuf, colinfo, hbuf);
                nline = 1;
            }
            fprintf(tp, "%d & %s & %s & %s & %.2f(%.2f) & %d & %s\\\\\n",
                    nr++, ai, aj, ak, ang, sigma, ntrain, params);
            nline++;
        }
    }
    end_table(tp);

    /* Dihedraltypes */
    strcpy(colinfo, "cccccccl");
    strcpy(hbuf, "Nr. & i & j & k & l & Angle & Ntrain & Params");
    begin_table(tp, "Dihedrals defined in the Alexandria force field.",
                "dihtypes", colinfo, hbuf);

    fprintf(fp, "\n[ dihedraltypes ]\n");
    fprintf(fp, "; ; i    j   k    l  func       parameters\n");
    nline          = 2;
    npage          = 0;
    nr             = 1;
    bts[ebtsPDIHS] = inverse_ifunc_index(d_dihedrals,
                                         gmx_poldata_get_dihedral_ftype(pd, egdPDIHS));
    while (0 < gmx_poldata_get_dihedral(pd, egdPDIHS, &ai, &aj, &ak, &al, &ang, &sigma, &ntrain, &params))
    {
        if (ntrain > 0)
        {
            fprintf(fp, "%-5s  %-5s  %-5s  %-5s  %d  %.1f  %s\n",
                    ai, aj, ak, al, bts[ebtsPDIHS], ang, params);
            if (0 == (nline % maxline))
            {
                end_table(tp);
                sprintf(pbuf, "dihtype%d", ++npage);
                fprintf(tp, "\\addtocounter{table}{-1}\n");
                begin_table(tp, "Dihedrals, continued", pbuf, colinfo, hbuf);
                nline = 1;
            }
            fprintf(tp, "%d & %s & %s & %s & %s & %.1f(%.1f) & %d & %s\\\\\n",
                    nr++, ai, aj, ak, al, ang, sigma, ntrain, params);
            nline++;
        }
    }
    end_table(tp);

    /* Impropertypes */
    strcpy(colinfo, "cccccccl");
    strcpy(hbuf, "Nr. & i & j & k & l & Angle & N & Params");
    begin_table(tp, "Impropers defined in the Alexandria force field.",
                "idihtypes", colinfo, hbuf);

    fprintf(fp, "\n[ dihedraltypes ]\n");
    fprintf(fp, "; ; i    j   k    l  func       parameters\n");
    nline          = 2;
    npage          = 0;
    nr             = 1;
    bts[ebtsIDIHS] = inverse_ifunc_index(d_dihedrals,
                                         gmx_poldata_get_dihedral_ftype(pd, egdIDIHS));
    while (0 < gmx_poldata_get_dihedral(pd, egdIDIHS, &ai, &aj, &ak, &al, &ang, &sigma, &ntrain, &params))
    {
        if (ntrain > 0)
        {
            fprintf(fp, "%-5s  %-5s  %-5s  %-5s  %d  %.1f  %s\n",
                    ai, aj, ak, al, bts[ebtsIDIHS], ang, params);
            if (0 == (nline % maxline))
            {
                end_table(tp);
                sprintf(pbuf, "idihtype%d", ++npage);
                fprintf(tp, "\\addtocounter{table}{-1}\n");
                begin_table(tp, "Impropers, continued", pbuf, colinfo, hbuf);
                nline = 1;
            }
            fprintf(tp, "%d & %s & %s & %s & %s & %.1f(%.1f) & %d & %s\\\\\n",
                    nr++, ai, aj, ak, al, ang, sigma, ntrain, params);
            nline++;
        }
    }
    end_table(tp);
}

int alex_gen_ff(int argc, char *argv[])
{
    static const char               *desc[] = {
        "gen_ff read a force field file (gentop.dat)",
        "and writes a GROMACS topology include file."
    };
    output_env_t                     oenv;

    if (!parse_common_args(&argc, argv, 0, 0, nullptr, 0, nullptr,
                           asize(desc), desc, 0, nullptr, &oenv))
    {
        return 0;
    }

    gmx_atomprop_t aps  = gmx_atomprop_init();
    gmx_poldata_t  pd   = gmx_poldata_read("gentop.dat", aps);
    FILE          *tp   = fopen("forcefield.tex", "w");
    fprintf(tp, "%% Generated by gen_ff\n");
    fprintf(tp, "%% Copyright 2013 David van der Spoel\n");

    FILE *fp = fopen("forcefield.itp", "w");
    fprintf(fp, "; This file is generated from gentop.dat by %s\n", argv[0]);
    fprintf(fp, "; Do not edit this file!\n");
    fprintf(fp, "; This is the force field file for the Alexandria FF\n");
    fprintf(fp, "; Paul J. van Maaren and David van der Spoel\n\n");

    do_atypes(fp, tp, pd, aps);
    do_brule(tp, pd);
    do_bad(fp, tp, pd);

    fclose(fp);
    fclose(tp);
    return 0;
}

/*
 * This file is part of the GROMACS molecular simulation package.
 *
 * Copyright (c) 1991-2000, University of Groningen, The Netherlands.
 * Copyright (c) 2001-2004, The GROMACS development team.
 * Copyright (c) 2013,2014,2015,2018, by the GROMACS development team, led by
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

#include <cstdio>
#include <cstring>

#include <algorithm>

#include "gromacs/domdec/domdec.h"
#include "gromacs/domdec/domdec_struct.h"
#include "gromacs/fileio/gmxfio.h"
#include "gromacs/gmxlib/nrnb.h"
#include "gromacs/mdlib/ns.h"
#include "gromacs/mdtypes/commrec.h"
#include "gromacs/mdtypes/forcerec.h"
#include "gromacs/mdtypes/md_enums.h"
#include "gromacs/mdtypes/nblist.h"
#include "gromacs/utility/fatalerror.h"
#include "gromacs/utility/futil.h"
#include "gromacs/utility/smalloc.h"

static void write_nblist(FILE *out, gmx_domdec_t *dd, t_nblist *nblist, int nDNL)
{
    int                 i, nii, ii, j, zi, zj0, zj1, aj, zj, nj;
    int                 ca1[DD_MAXZONE], np[DD_MAXZONE];
    gmx_domdec_zones_t *dd_zones;

    if (nblist->nri > 0)
    {
        fprintf(out, "elec: %s, vdw: %s, type: %s, Solvent opt: %s\n",
                gmx_nbkernel_elec_names[nblist->ielec],
                gmx_nbkernel_vdw_names[nblist->ivdw],
                gmx_nblist_interaction_names[nblist->type],
                gmx_nblist_geometry_names[nblist->igeometry]);
        fprintf(out, "nri: %d  npair: %d\n", nblist->nri, nblist->nrj);
        if (dd)
        {
            dd_zones = domdec_zones(dd);

            for (zi = 0; zi < dd_zones->n; zi++)
            {
                ca1[zi] = dd->atomGrouping().block(dd_zones->cg_range[zi + 1]).begin();
            }
            i = 0;
            for (zi = 0; zi < dd_zones->nizone && zi < dd_zones->n; zi++)
            {
                zj0 = dd_zones->izone[zi].j0;
                zj1 = dd_zones->izone[zi].j1;
                for (zj = zj0; zj < zj1; zj++)
                {
                    np[zj] = 0;
                }
                while (i < nblist->nri && nblist->iinr[i] < ca1[zi])
                {
                    for (j = nblist->jindex[i]; (j < nblist->jindex[i+1]); j++)
                    {
                        aj = nblist->jjnr[j];
                        zj = zj0;
                        while (aj >= ca1[zj])
                        {
                            zj++;
                        }
                        np[zj]++;
                    }
                    i++;
                }
                fprintf(out, "DD zone %d:", zi);
                for (zj = zj0; zj < zj1; zj++)
                {
                    fprintf(out, " %d %d", zj, np[zj]);
                }
                fprintf(out, "\n");
            }
        }
        if (nDNL >= 2)
        {
            for (i = 0; i < nblist->nri; i++)
            {
                nii = 1;
                if (nDNL >= 3 && nblist->igeometry != GMX_NBLIST_GEOMETRY_PARTICLE_PARTICLE)
                {
                    nii = 3;
                }
                nj = nblist->jindex[i+1] - nblist->jindex[i];
                fprintf(out, "i: %d shift: %d gid: %d nj: %d\n",
                        ddglatnr(dd, nblist->iinr[i]),
                        nblist->shift[i], nblist->gid[i], nj);
                for (ii = 0; ii < nii; ii++)
                {
                    for (j = nblist->jindex[i]; (j < nblist->jindex[i+1]); j++)
                    {
                        fprintf(out, "  i: %5d  j: %5d\n",
                                ddglatnr(dd, nblist->iinr[i]+ii),
                                ddglatnr(dd, nblist->jjnr[j]));
                    }
                }
            }
        }
        fflush(out);
    }
}



void dump_nblist(FILE *out, const t_commrec *cr, t_forcerec *fr, int nDNL)
{
    int  n, i;

    fprintf(out, "Neighborlist:\n");

    for (n = 0; (n < fr->nnblists); n++)
    {
        for (i = 0; (i < eNL_NR); i++)
        {
            write_nblist(out, cr->dd, &fr->nblists[n].nlist_sr[i], nDNL);
        }
    }

}

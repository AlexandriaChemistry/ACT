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
/* Adapted to ACT by DvdS, 2023-06-03 */
#include "actpre.h"

#include "pdbwriter.h"

#include <cctype>
#include <cmath>
#include <cstdlib>
#include <cstring>

#include <map>
#include <string>

#include "gromacs/fileio/gmxfio.h"
#include "gromacs/math/units.h"
#include "gromacs/math/vec.h"
#include "gromacs/pbcutil/pbc.h"
#include "gromacs/topology/atomprop.h"
#include "gromacs/utility/cstringutil.h"
#include "gromacs/utility/fatalerror.h"
#include "gromacs/utility/futil.h"

typedef struct {
    int ai, aj;
} gmx_conection_t;

typedef struct gmx_conect_t {
    int              nconect;
    gmx_bool         bSorted;
    gmx_conection_t *conect;
} gmx_conect_t;

static void gmx_write_pdb_box2(FILE *out, int ePBC, const matrix box)
{
    real alpha, beta, gamma;

    if (ePBC == -1)
    {
        ePBC = guess_ePBC(box);
    }

    if (ePBC == epbcNONE)
    {
        return;
    }

    if (norm2(box[YY])*norm2(box[ZZ]) != 0)
    {
        alpha = RAD2DEG*gmx_angle(box[YY], box[ZZ]);
    }
    else
    {
        alpha = 90;
    }
    if (norm2(box[XX])*norm2(box[ZZ]) != 0)
    {
        beta  = RAD2DEG*gmx_angle(box[XX], box[ZZ]);
    }
    else
    {
        beta  = 90;
    }
    if (norm2(box[XX])*norm2(box[YY]) != 0)
    {
        gamma = RAD2DEG*gmx_angle(box[XX], box[YY]);
    }
    else
    {
        gamma = 90;
    }
    fprintf(out, "REMARK    THIS IS A SIMULATION BOX\n");
    if (ePBC != epbcSCREW)
    {
        fprintf(out, "CRYST1%9.3f%9.3f%9.3f%7.2f%7.2f%7.2f %-11s%4d\n",
                10*norm(box[XX]), 10*norm(box[YY]), 10*norm(box[ZZ]),
                alpha, beta, gamma, "P 1", 1);
    }
    else
    {
        /* Double the a-vector length and write the correct space group */
        fprintf(out, "CRYST1%9.3f%9.3f%9.3f%7.2f%7.2f%7.2f %-11s%4d\n",
                20*norm(box[XX]), 10*norm(box[YY]), 10*norm(box[ZZ]),
                alpha, beta, gamma, "P 21 1 1", 1);

    }
}

namespace alexandria
{

void pdbWriter(FILE                           *out, 
               const char                     *title,
               const std::vector<ActAtom>     &atoms,
               const std::vector<gmx::RVec>   &x,
               const std::vector<std::string> &residueNames,
               int                             ePBC,
               const matrix                    box,
               char                            chain,
               int                             model_nr,
               const std::vector<int>         &index,
               gmx_conect                      conect,
               bool                            renumberAtoms)
{
    gmx_conect_t     *gc = static_cast<gmx_conect_t *>(conect);
    int               i;
    int               resind, resnr;
    enum PDB_record   type;
    unsigned char     ch;
    char              altloc;
    real              occup, bfac;
    bool              bOccup = false;
    const char       *p_restype;
    const char       *p_lastrestype;

    if (title)
    {
        fprintf(out, "TITLE     %s\n", title);
    }
    
    if (box && ( (norm2(box[XX]) != 0.0f) || (norm2(box[YY]) != 0.0f) || (norm2(box[ZZ]) != 0.0f) ) )
    {
        gmx_write_pdb_box2(out, ePBC, box);
    }

    fprintf(out, "MODEL %8d\n", model_nr > 0 ? model_nr : 1);

    p_restype         = nullptr;

    std::map<int, int> reverseIndex;
    if (renumberAtoms)
    {
        for (size_t ii = 0; ii < index.size(); ii++)
        {
            reverseIndex.insert({index[ii], ii});
        }
    }
    unsigned char resic = ' ';
    size_t natoms = atoms.size();
    if (!index.empty())
    {
        natoms = index.size();
    }
    for (size_t ii = 0; ii < natoms; ii++)
    {
        size_t i = ii;
        if (!index.empty())
        {
            i = index[ii];
        }
        resind        = atoms[i].residueNumber();
        resnr         = resind;
        p_lastrestype = p_restype;

        ch = 0;
        if (resnr >= 10000)
        {
            resnr = resnr % 10000;
        }
        t_pdbinfo pdbinfo;
        gmx_pdbinfo_init_default(&pdbinfo);
        
        type   = static_cast<enum PDB_record>(pdbinfo.type);
        altloc = pdbinfo.altloc;
        if (!isalnum(altloc))
        {
            altloc = ' ';
        }
        occup = bOccup ? 1.0 : pdbinfo.occup;
        bfac  = pdbinfo.bfac;

        int atomnr = i+1;
        if (renumberAtoms)
        {
            atomnr = ii+1;
        }
        gmx_fprintf_pdb_atomline(out,
                                 type,
                                 atomnr,
                                 atoms[i].name().c_str(),
                                 altloc,
                                 residueNames[resind].c_str(),
                                 chain,
                                 resnr,
                                 resic,
                                 10*x[i][XX], 10*x[i][YY], 10*x[i][ZZ],
                                 occup,
                                 bfac,
                                 atoms[i].element().c_str());
    }

    fprintf(out, "TER\n");

    if (nullptr != gc)
    {
        /* Write conect records */
        for (i = 0; (i < gc->nconect); i++)
        {
            int ai = gc->conect[i].ai;
            int aj = gc->conect[i].aj;
            if (renumberAtoms)
            {
                ai = reverseIndex.find(ai)->second;
                aj = reverseIndex.find(aj)->second;
            }
            fprintf(out, "CONECT%5d%5d\n", ai+1, aj+1);
        }
    }
    fprintf(out, "ENDMDL\n");
}

} // namespace

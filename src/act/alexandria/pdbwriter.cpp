/*
 * This file is part of the GROMACS molecular simulation package.
 *
 * Copyright (c) 1991-2000, University of Groningen, The Netherlands.
 * Copyright (c) 2001-2004, The GROMACS development team.
 * Copyright (c) 2013,2014,2015,2016,2017,2018, by the GROMACS development team, led by
 * Mark Abraham, David van der Spoel, Berk Hess, and Erik Lindahl,
 * and including many others, as listed in the AUTHORS file in the
 * top-level source directory and at http://www.gromacs.org.
 * Copyright (C) 2026, the Alexandria Chemistry Toolkit developers.
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

#include <algorithm>
#include <cctype>
#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <map>
#include <string>
#include <vector>

#include "act/basics/act_particle.h"
#include "act/topology/actatom.h"
#include "gromacs/fileio/pdbio.h"
#include "gromacs/math/units.h"
#include "gromacs/math/vec.h"
#include "gromacs/math/vectypes.h"
#include "gromacs/pbcutil/pbc.h"
#include "gromacs/topology/atoms.h"
#include "gromacs/utility/basedefinitions.h"
#include "gromacs/utility/real.h"
#include "gromacs/utility/stringutil.h"
#include "pdbwriter.h"

/*! \brief Write the simulation box as a CRYST1 record to a pdb file,
 * nothing is written if there is no periodicity
 * \param[in] out  The file to write to
 * \param[in] ePBC The type of periodic boundary conditions, if -1 it is
 *                 guessed from the box
 * \param[in] box  The box vectors (nm), written in Angstrom
 */
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

std::vector<std::string> particleNames(const std::vector<ActAtom> &atoms)
{
    std::vector<std::string> names;
    int number  = 0;
    int residue = -1;
    for (size_t i = 0; i < atoms.size(); i++)
    {
        if (i == 0 || atoms[i].residueNumber() != residue)
        {
            residue = atoms[i].residueNumber();
            number  = 0;
        }
        if (atoms[i].pType() != ActParticle::Shell)
        {
            number += 1;
        }
        names.push_back(gmx::formatString("%s%d", atoms[i].name().c_str(),
                                          std::max(number, 1)));
    }
    return names;
}

std::vector<std::string> ionCharges(const std::vector<ActAtom> &atoms)
{
    // The number of atoms and the total charge of each residue
    std::map<int, int>    natom;
    std::map<int, double> qtot;
    for (const auto &atom : atoms)
    {
        if (atom.pType() == ActParticle::Atom)
        {
            natom[atom.residueNumber()] += 1;
        }
        // The second charge is on the atom unless it is a virtual site
        qtot[atom.residueNumber()] += atom.charge() + atom.chargeVS();
    }
    std::vector<std::string> charges;
    for (const auto &atom : atoms)
    {
        std::string q;
        auto        res = atom.residueNumber();
        long        iq  = std::lround(qtot[res]);
        if (atom.pType() == ActParticle::Atom && natom[res] == 1 &&
            iq != 0 && std::abs(qtot[res] - iq) < 1e-3)
        {
            q = gmx::formatString("%ld%c", std::abs(iq), iq > 0 ? '+' : '-');
        }
        charges.push_back(q);
    }
    return charges;
}

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
               const std::vector<PdbBond>     &bonds,
               bool                            renumberAtoms)
{
    int               resind, resnr = -1;
    enum PDB_record   type;
    //unsigned char     ch;
    char              altloc;
    real              occup, bfac;
    bool              bOccup = false;
    //const char       *p_restype;
    //const char       *p_lastrestype;

    if (title)
    {
        fprintf(out, "TITLE     %s\n", title);
    }
    
    if (box && ( (norm2(box[XX]) != 0.0f) || (norm2(box[YY]) != 0.0f) || (norm2(box[ZZ]) != 0.0f) ) )
    {
        gmx_write_pdb_box2(out, ePBC, box);
    }

    fprintf(out, "MODEL %8d\n", model_nr > 0 ? model_nr : 1);

    //p_restype         = nullptr;

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
    int minres = 0;
    for (size_t ii = 0; ii < natoms; ii++)
    {
        size_t i = ii;
        if (!index.empty())
        {
            i = index[ii];
        }
        minres = std::min(minres, atoms[i].residueNumber());
    }
    // The names are made from all particles, also those that are not written
    auto anames = particleNames(atoms);
    auto qions  = ionCharges(atoms);
    for (size_t ii = 0; ii < natoms; ii++)
    {
        size_t i = ii;
        if (!index.empty())
        {
            i = index[ii];
        }
        resind        = atoms[i].residueNumber();
        resnr         = resind;
        if (minres == 0)
        {
            resnr += 1;
        }

        if (resnr >= 10000)
        {
            resnr = 1 + resnr % 10000;
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
                                 anames[i].c_str(),
                                 altloc,
                                 residueNames[resind].c_str(),
                                 chain,
                                 resnr,
                                 resic,
                                 10*x[i][XX], 10*x[i][YY], 10*x[i][ZZ],
                                 occup,
                                 bfac,
                                 atoms[i].element().c_str(),
                                 qions[i].empty() ? nullptr : qions[i].c_str());
    }

    // A double bond is listed twice and a triple bond three times, as
    // RDKit and Open Babel read them
    for (const auto &bond : bonds)
    {
        int ai = bond.ai;
        int aj = bond.aj;
        if (renumberAtoms)
        {
            ai = reverseIndex.find(ai)->second;
            aj = reverseIndex.find(aj)->second;
        }
        fprintf(out, "CONECT%5d", ai+1);
        for (int k = 0; k < std::max(1, bond.order); k++)
        {
            fprintf(out, "%5d", aj+1);
        }
        fprintf(out, "\n");
    }
    fprintf(out, "ENDMDL\n");
    fprintf(out, "TER\n");
}

} // namespace

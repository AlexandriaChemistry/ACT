/*
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

#include "gentop_core.h"

#include <ctype.h>
#include <stdlib.h>
#include <string.h>

#include "gromacs/gmxpreprocess/toputil.h"
#include "gromacs/listed-forces/bonded.h"
#include "gromacs/math/units.h"
#include "gromacs/math/vec.h"
#include "gromacs/math/vecdump.h"
#include "gromacs/pbcutil/pbc.h"
#include "gromacs/topology/atomprop.h"
#include "gromacs/utility/futil.h"
#include "gromacs/utility/smalloc.h"

using namespace alexandria;

void calc_angles_dihs(t_params *ang, 
                      t_params *dih, 
                      rvec      x[], 
                      gmx_bool  bPBC,
                      matrix    box)
{
    int           i, ai, aj, ak, al, t1, t2, t3;
    rvec          r_ij, r_kj, r_kl, m, n;
    real          th, costh, ph;
    struct t_pbc  pbc;

    if (bPBC)
    {
        set_pbc(&pbc, -1, box);
    }
    if (debug)
    {
        pr_rvecs(debug, 0, "GENTOP", box, DIM);
    }
    for (i = 0; (i < ang->nr); i++)
    {
        ai = ang->param[i].a[0];
        aj = ang->param[i].a[1];
        ak = ang->param[i].a[2];
        th = RAD2DEG*bond_angle(x[ai], x[aj], x[ak], bPBC ? &pbc : nullptr,
                                r_ij, r_kj, &costh, &t1, &t2);
        if (debug)
        {
            fprintf(debug, "GENTOP: ai=%3d aj=%3d ak=%3d r_ij=%8.3f r_kj=%8.3f th=%8.3f\n",
                    ai, aj, ak, gmx::norm(r_ij), gmx::norm(r_kj), th);
        }
        ang->param[i].c[0] = th;
    }
    for (i = 0; (i < dih->nr); i++)
    {
        ai = dih->param[i].a[0];
        aj = dih->param[i].a[1];
        ak = dih->param[i].a[2];
        al = dih->param[i].a[3];
        ph = RAD2DEG*dih_angle(x[ai], x[aj], x[ak], x[al], bPBC ? &pbc : nullptr,
                               r_ij, r_kj, r_kl, m, n, &t1, &t2, &t3);
        if (debug)
        {
            fprintf(debug, "GENTOP: ai=%3d aj=%3d ak=%3d al=%3d r_ij=%8.3f r_kj=%8.3f r_kl=%8.3f ph=%8.3f\n",
                    ai, aj, ak, al, gmx::norm(r_ij), gmx::norm(r_kj), gmx::norm(r_kl), ph);
        }
        dih->param[i].c[0] = ph;
    }
}

void dump_hybridization(FILE *fp, t_atoms *atoms, int nbonds[])
{
    int i;

    for (i = 0; (i < atoms->nr); i++)
    {
        fprintf(fp, "Atom %5s has %d bonds\n", *atoms->atomname[i], nbonds[i]);
    }
}

static void print_pl(FILE       *fp, 
                     t_params    plist[], 
                     int         ftp, 
                     const char *name,
                     char       ***atomname)
{
    int i, j, nral, nrfp;

    if (plist[ftp].nr > 0)
    {
        fprintf(fp, "\n");
        fprintf(fp, "[ %s ]\n", name);
        nral = interaction_function[ftp].nratoms;
        nrfp = interaction_function[ftp].nrfpA;
        for (i = 0; (i < plist[ftp].nr); i++)
        {
            for (j = 0; (j < nral); j++)
            {
                fprintf(fp, "  %5s", *atomname[plist[ftp].param[i].a[j]]);
            }
            for (j = 0; (j < nrfp); j++)
            {
                fprintf(fp, "  %10.3e", plist[ftp].param[i].c[j]);
            }
            fprintf(fp, "\n");
        }
    }
}

void print_rtp(const char *filenm, const char *title, t_atoms *atoms,
               t_params plist[], int cgnr[], int nbts, int bts[])
{
    FILE *fp;
    int   i;

    fp = gmx_ffopen(filenm, "w");
    fprintf(fp, "; %s\n", title);
    fprintf(fp, "\n");
    fprintf(fp, "[ %s ]\n", *atoms->resinfo[0].name);
    fprintf(fp, "\n");
    fprintf(fp, "[ atoms ]\n");
    for (i = 0; (i < atoms->nr); i++)
    {
        fprintf(fp, "%-8s  %12s  %8.4f  %5d\n",
                *atoms->atomname[i], *atoms->atomtype[i],
                atoms->atom[i].q, cgnr[i]);
    }
    for (i = 0; (i < nbts); i++)
    {
        print_pl(fp, plist, bts[i], interaction_function[bts[i]].name, atoms->atomname);
    }
    fclose(fp);
}

void reset_q(t_atoms *atoms)
{
    int i;

    /* Use values from file */
    for (i = 0; (i < atoms->nr); i++)
    {
        atoms->atom[i].qB = atoms->atom[i].q;
    }
}

void symmetrize_charges(gmx_bool                   bQsym, 
                        t_atoms                   *atoms,
                        ConstPlistWrapperIterator  bonds,
                        const Poldata             *pd,
                        gmx_atomprop_t             aps, 
                        const char                *symm_string,
                        std::vector<int>          &sym_charges)
{
    std::string  central, attached;
    int          ai, aj, anri, anrj;
    int          anr_central, anr_attached, nrq;
    double       qaver, qsum;

    sym_charges.clear();
    for (int i = 0; i < atoms->nr; i++)
    {
        sym_charges.push_back(i);
    }
    if (bQsym)
    {
        if ((nullptr != symm_string) && (strlen(symm_string) > 0))
        {
            std::vector<std::string> ss = gmx::splitString(symm_string);
            if (static_cast<int>(ss.size()) != atoms->nr)
            {
                gmx_fatal(FARGS, "Wrong number (%d) of atom-numbers in symm_string: expected %d",
                          static_cast<int>(ss.size()), atoms->nr);
            }
            int ii = 0;
            for (auto is = ss.begin();
                 (is < ss.end()); ++is)
            {
                sym_charges[ii] = atoi(is->c_str());
                ii++;
            }
        }
        else
        {
            for (auto symcharges = pd->getSymchargesBegin();
                 symcharges != pd->getSymchargesEnd(); symcharges++)
            {
                anr_central  = gmx_atomprop_atomnumber(aps, symcharges->getCentral().c_str());
                anr_attached = gmx_atomprop_atomnumber(aps, symcharges->getAttached().c_str());
                for (int i = 0; i < atoms->nr; i++)
                {
                    if (atoms->atom[i].atomnumber == anr_central)
                    {
                        int              hsmin = -1;
                        std::vector<int> hs;
                        for (auto j = bonds->beginParam(); j < bonds->endParam(); ++j)
                        {
                            ai   = j->a[0];
                            aj   = j->a[1];
                            anri = atoms->atom[ai].atomnumber;
                            anrj = atoms->atom[aj].atomnumber;
                            
                            if ((ai == i) && (anrj == anr_attached))
                            {
                                hs.push_back(aj);
                            }
                            else if ((aj == i) && (anri == anr_attached))
                            {
                                hs.push_back(ai);
                            }
                            if ((hs.size() > 0) && (hsmin == -1 || hs.back() < hsmin))
                            {
                                hsmin = hs.back();
                            }
                        }
                        if ((static_cast<int>(hs.size()) == symcharges->getNumattach()) &&
                            (hsmin != -1))
                        {
                            for (int j = 0; j < symcharges->getNumattach(); j++)
                            {
                                sym_charges[hs[j]] = hsmin;
                            }
                        }
                    }
                }
            }
        }

        for (int i = 0; i < atoms->nr; i++)
        { 
            qsum = 0;
            nrq  = 0;
            for (int j = i+1; j < atoms->nr; j++)
            {
                if (sym_charges[j] == sym_charges[i])
                {
                    qsum += atoms->atom[j].q;
                    nrq++;
                }
            }
            if (0 < nrq)
            {
                qaver = qsum/nrq;
                for (int j = 0; j < atoms->nr; j++)
                {
                    if (sym_charges[j] == sym_charges[i])
                    {
                        atoms->atom[j].q = atoms->atom[j].qB = qaver;
                    }
                }
            }
        }
    }
}

static int *generate_cg_neutral(t_atoms *atoms, gmx_bool bUsePDBcharge)
{
    int     i, n = 1;
    int    *cgnr;
    double  qt = 0;

    snew(cgnr, atoms->nr);
    for (i = 0; (i < atoms->nr); i++)
    {
        if (atoms->pdbinfo && bUsePDBcharge)
        {
            atoms->atom[i].q = atoms->pdbinfo[i].bfac;
        }
        qt     += atoms->atom[i].q;
        cgnr[i] = n;
        if (is_int(qt))
        {
            n++;
            qt = 0;
        }
    }
    return cgnr;
}

static int *generate_cg_group(t_atoms                               *atoms,
                              const std::vector<alexandria::PlistWrapper> &plist)
{
    int        i, j, k, atn, ai, aj, ncg = 1;
    int       *cgnr;
    gmx_bool   bMV;
    int        monovalent[] = { 0, 1, 9, 17, 35, 53, 85 };
    int        nmv          = sizeof(monovalent)/sizeof(monovalent[0]);

    /* Assume that shells and masses have atomnumber 0 */
    snew(cgnr, atoms->nr);
    for (i = 0; (i < atoms->nr); i++)
    {
        cgnr[i] = -1;
    }

    for (i = 0; (i < atoms->nr); i++)
    {
        atn = atoms->atom[i].atomnumber;
        bMV = FALSE;
        for (j = 0; (j < nmv) && !bMV; j++)
        {
            bMV = (atn == monovalent[j]);
        }
        if (!bMV)
        {
            cgnr[i] = ncg++;
        }
    }
    /* Rely on the notion that all H and other monovalent
       atoms are bound to something */
    auto bonds = SearchPlist(plist, F_BONDS);
    if (plist.end() != bonds)
    {
        for (auto j = bonds->beginParam();
             (j < bonds->endParam()); ++j)
        {
            ai  = j->a[0];
            aj  = j->a[1];
            bMV = FALSE;
            atn = atoms->atom[ai].atomnumber;
            for (k = 0; (k < nmv) && !bMV; k++)
            {
                bMV = (atn == monovalent[k]);
            }
            if (bMV)
            {
                if (cgnr[aj] != -1)
                {
                    cgnr[ai] = cgnr[aj];
                }
                else
                {
                    cgnr[ai] = cgnr[aj] = 1;
                }
            }
            else
            {
                bMV = FALSE;
                atn = atoms->atom[aj].atomnumber;
                for (k = 0; (k < nmv) && !bMV; k++)
                {
                    bMV = (atn == monovalent[k]);
                }
                if (bMV)
                {
                    cgnr[aj] = cgnr[ai];
                }
            }
        }
    }
    /* Rely on the notion that all shells are bound to something */
    auto pols = SearchPlist(plist, F_POLARIZATION);
    if (plist.end() != pols)
    {
        for (auto j = pols->beginParam();
             (j < pols->endParam()); ++j)
        {
            ai       = j->a[0];
            aj       = j->a[1];
            cgnr[aj] = cgnr[ai];
        }
        for (i = 0; (i < atoms->nr); i++)
        {
            if (cgnr[i] == -1)
            {
                cgnr[i] = ncg++;
            }
        }
    }
    if (debug)
    {
        fprintf(debug, "There are %d charge groups\n", ncg-1);
    }
    return cgnr;
}

static int *generate_cg_atom(int natom)
{
    int i, *cgnr;

    snew(cgnr, natom);
    for (i = 0; (i < natom); i++)
    {
        cgnr[i] = i+1;
    }

    return cgnr;
}

int *generate_charge_groups(eChargeGroup                           cgtp, 
                            t_atoms                               *atoms,
                            const std::vector<alexandria::PlistWrapper> &plist,
                            bool                                   bUsePDBcharge,
                            real                                  *qtot, 
                            real                                  *mtot)
{
    int i, *cgnr = nullptr;
    //auto pb = alexandria::SearchPlist(plist, F_BONDS);
    //auto ps = alexandria::SearchPlist(plist, F_POLARIZATION);

    switch (cgtp)
    {
        case ecgNeutral:
            cgnr = generate_cg_neutral(atoms, bUsePDBcharge);
            break;
        case ecgGroup:
            cgnr = generate_cg_group(atoms, plist);
            break;
        case ecgAtom:
            cgnr = generate_cg_atom(atoms->nr);
            break;
        default:
            gmx_fatal(FARGS, "Invalid charge group generation type %d", cgtp);
    }
    *qtot = *mtot = 0;
    for (i = 0; (i < atoms->nr); i++)
    {
        *qtot += atoms->atom[i].q;
        *mtot += atoms->atom[i].m;
    }
    return cgnr;
}

static int    *cgnr_copy;
static double *atomnumber;
static int cg_comp(const void *a, const void *b)
{
    int   *aa = (int *)a;
    int   *bb = (int *)b;
    double c;

    int    d = cgnr_copy[*aa] - cgnr_copy[*bb];
    if (d == 0)
    {
        c = atomnumber[*aa] - atomnumber[*bb];
        if (c < 0)
        {
            return -1;
        }
        else if (c > 0)
        {
            return 1;
        }
        else
        {
            return 0;
        }
    }
    else
    {
        return d;
    }
}

void sort_on_charge_groups(int                                   *cgnr, 
                           t_atoms                               *atoms,
                           std::vector<alexandria::PlistWrapper> *pw,
                           rvec                                  x[], 
                           t_excls                               excls[],
                           const char                            *ndxout, 
                           int                                    nmol)
{
    FILE      *fp;
    int        i, j, j0, k, newi, ri, *cg_renum, *ccgg, *inv_renum;
    rvec      *rx;
    t_atom    *ra;
    t_excls   *newexcls;
    char    ***an, ***smn;

    snew(cg_renum, atoms->nr);
    snew(atomnumber, atoms->nr);
    snew(rx, atoms->nr);
    snew(ra, atoms->nr);
    snew(an, atoms->nr);
    for (i = 0; (i < atoms->nr); i++)
    {
        cg_renum[i]   = i;
        atomnumber[i] = 1+i; /*atoms->atom[i].atomnumber;*/
        if ((atoms->atom[i].ptype == eptShell) && (i > 0))
        {
            atomnumber[i] = atomnumber[i-1]+0.1;
        }
    }
    cgnr_copy = cgnr;
    qsort(cg_renum, atoms->nr, sizeof(cg_renum[0]), cg_comp);
    if (debug)
    {
        for (i = 0; (i < atoms->nr); i++)
        {
            fprintf(debug, "cg_renum[%d] = %d\n", i, cg_renum[i]);
        }
    }
    snew(ccgg, atoms->nr);
    for (i = 0; (i < atoms->nr); i++)
    {
        ri = cg_renum[i];
        copy_rvec(x[ri], rx[i]);
        memcpy(&(ra[i]), &(atoms->atom[ri]), sizeof(t_atom));
        an[i]   = atoms->atomname[ri];
        ccgg[i] = cgnr[ri];
    }
    snew(inv_renum, atoms->nr);
    snew(smn, atoms->nr);
    for (i = 0; (i < atoms->nr); i++)
    {
        copy_rvec(rx[i], x[i]);
        memcpy(&(atoms->atom[i]), &(ra[i]), sizeof(t_atom));
        atoms->atomname[i]     = an[i];
        cgnr[i]                = ccgg[i];
        inv_renum[cg_renum[i]] = i;
        smn[i]                 = atoms->atomtype[i];
    }
    for (i = 0; (i < atoms->nr); i++)
    {
        newi               = cg_renum[i];
        atoms->atomtype[i] = smn[newi];
    }
    for (auto i = pw->begin(); i < pw->end(); ++i)
    {
        for (auto j = i->beginParam(); j < i->endParam(); j++)
        {
            for (k = 0; (k < NRAL(i->getFtype())); k++)
            {
                j->a[k] = inv_renum[j->a[k]];
            }
        }
    }
    snew(newexcls, atoms->nr);
    for (i = 0; (i < atoms->nr); i++)
    {
        snew(newexcls[i].e, excls[i].nr);
        newexcls[i].nr = excls[i].nr;
        for (j = 0; (j < excls[i].nr); j++)
        {
            newexcls[i].e[j] = excls[i].e[j];
        }
    }
    for (i = 0; (i < atoms->nr); i++)
    {
        newi = inv_renum[i];
        if (newexcls[i].nr > excls[newi].nr)
        {
            srenew(excls[newi].e, newexcls[i].nr);
        }
        for (j = 0; (j < newexcls[i].nr); j++)
        {
            excls[newi].e[j] = inv_renum[newexcls[i].e[j]];
        }
        excls[newi].nr = newexcls[i].nr;
    }
    if (nullptr != ndxout)
    {
        fp = fopen(ndxout, "w");
        fprintf(fp, "[ number_backward ]\n");
        for (j = 0; (j < nmol); j++)
        {
            j0 = j*atoms->nr;
            for (i = 0; (i < atoms->nr); i++)
            {
                if (atoms->atom[inv_renum[i]].ptype == eptShell)
                {
                    k = j0+inv_renum[i-1]+1;
                }
                else
                {
                    k = j0+inv_renum[i]+1;
                }
                fprintf(fp, " %d", k);
                if (j == 0)
                {
                    cg_renum[inv_renum[i]] = i;
                }
            }
            fprintf(fp, "\n");
        }
        for (j = 0; (j < nmol); j++)
        {
            j0 = j*atoms->nr;
            fprintf(fp, "[ number_forward ]\n");
            for (i = 0; (i < atoms->nr); i++)
            {
                if (atoms->atom[cg_renum[i]].ptype == eptShell)
                {
                    k = j0+cg_renum[i-1]+1;
                }
                else
                {
                    k = j0+cg_renum[i]+1;
                }
                fprintf(fp, " %d", k);
            }
            fprintf(fp, "\n");
        }
        fclose(fp);
    }
    sfree(rx);
    sfree(ra);
    sfree(an);
    sfree(cg_renum);
    sfree(inv_renum);
    sfree(ccgg);
}

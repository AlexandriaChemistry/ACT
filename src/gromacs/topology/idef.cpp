/*
 * This file is part of the GROMACS molecular simulation package.
 *
 * Copyright (c) 1991-2000, University of Groningen, The Netherlands.
 * Copyright (c) 2001-2004, The GROMACS development team.
 * Copyright (c) 2013,2014,2015,2016,2018, by the GROMACS development team, led by
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

#include "idef.h"

#include <cstdio>

#include "gromacs/topology/ifunc.h"
#include "gromacs/utility/fatalerror.h"
#include "gromacs/utility/smalloc.h"
#include "gromacs/utility/txtdump.h"

static void pr_harm(FILE *fp, const t_iparams *iparams, const char *r, const char *kr)
{
    fprintf(fp, "%sA=%12.5e, %sA=%12.5e, %sB=%12.5e, %sB=%12.5e\n",
            r, iparams->harmonic.rA, kr, iparams->harmonic.krA,
            r, iparams->harmonic.rB, kr, iparams->harmonic.krB);
}

void pr_iparams(FILE *fp, t_functype ftype, const t_iparams *iparams)
{
    switch (ftype)
    {
        case F_ANGLES:
        case F_G96ANGLES:
            pr_harm(fp, iparams, "th", "ct");
            break;
        case F_CROSS_BOND_BONDS:
            fprintf(fp, "r1e=%15.8e, r2e=%15.8e, krr=%15.8e\n",
                    iparams->cross_bb.r1e, iparams->cross_bb.r2e,
                    iparams->cross_bb.krr);
            break;
        case F_CROSS_BOND_ANGLES:
            fprintf(fp, "r1e=%15.8e, r1e=%15.8e, r3e=%15.8e, krt=%15.8e\n",
                    iparams->cross_ba.r1e, iparams->cross_ba.r2e,
                    iparams->cross_ba.r3e, iparams->cross_ba.krt);
            break;
        case F_LINEAR_ANGLES:
            fprintf(fp, "klinA=%15.8e, aA=%15.8e, r13A=%15.8e, kUBA=%15.8e, klinB=%15.8e, aB=%15.8e, r13B=%15.8e, kUBB=%15.8e\n",
                    iparams->linangle.klinA, iparams->linangle.aA, iparams->linangle.r13A, iparams->linangle.kUBA,
                    iparams->linangle.klinB, iparams->linangle.aB, iparams->linangle.r13B, iparams->linangle.kUBB);
            break;
        case F_UREY_BRADLEY:
            fprintf(fp, "thetaA=%15.8e, kthetaA=%15.8e, r13A=%15.8e, kUBA=%15.8e, thetaB=%15.8e, kthetaB=%15.8e, r13B=%15.8e, kUBB=%15.8e\n", iparams->u_b.thetaA, iparams->u_b.kthetaA, iparams->u_b.r13A, iparams->u_b.kUBA, iparams->u_b.thetaB, iparams->u_b.kthetaB, iparams->u_b.r13B, iparams->u_b.kUBB);
            break;
        case F_QUARTIC_ANGLES:
            fprintf(fp, "theta=%15.8e", iparams->qangle.theta);
            for (int i = 0; i < 5; i++)
            {
                fprintf(fp, ", c%c=%15.8e", '0'+i, iparams->qangle.c[i]);
            }
            fprintf(fp, "\n");
            break;
        case F_WBHAM:
            fprintf(fp, "a=%15.8e, b=%15.8e, c=%15.8e\n",
                    iparams->bham.a, iparams->bham.b, iparams->bham.c);
            break;
        case F_GBHAM:
            fprintf(fp, "rmin=%15.8e, epsilon=%15.8e, gamma=%15.8e, delta=%15.8e\n",
                    iparams->gbham.rmin, iparams->gbham.epsilon, iparams->gbham.gamma, iparams->gbham.delta);
            break;
        case F_BONDS:
        case F_G96BONDS:
        case F_HARMONIC:
            pr_harm(fp, iparams, "b0", "cb");
            break;
        case F_IDIHS:
            pr_harm(fp, iparams, "xi", "cx");
            break;
        case F_MORSE:
            fprintf(fp, "b0A=%15.8e, cbA=%15.8e, betaA=%15.8e, D0A=%15.8e b0B=%15.8e, cbB=%15.8e, betaB=%15.8e, D0B=%15.8e\n",
                    iparams->morse.b0A, iparams->morse.cbA,
                    iparams->morse.betaA, iparams->morse.D0A,
                    iparams->morse.b0B, iparams->morse.cbB,
                    iparams->morse.betaB, iparams->morse.D0B);
            break;
        case F_CUBICBONDS:
            fprintf(fp, "b0=%15.8e, kb=%15.8e, kcub=%15.8e\n",
                    iparams->cubic.b0, iparams->cubic.kb, iparams->cubic.kcub);
            break;
        case F_CONNBONDS:
            fprintf(fp, "\n");
            break;
        case F_FENEBONDS:
            fprintf(fp, "bm=%15.8e, kb=%15.8e\n", iparams->fene.bm, iparams->fene.kb);
            break;
        case F_TABBONDS:
        case F_TABBONDSNC:
        case F_TABANGLES:
        case F_TABDIHS:
            fprintf(fp, "tab=%d, kA=%15.8e, kB=%15.8e\n",
                    iparams->tab.table, iparams->tab.kA, iparams->tab.kB);
            break;
        case F_POLARIZATION:
            fprintf(fp, "alpha=%15.8e\n", iparams->polarize.alpha);
            break;
            /*case F_HYPER_POL:
            fprintf(fp, "alpha=%15.8e k3=%15.8e k4=%15.8e\n",
                    iparams->hyper_polarize.alpha,
                    iparams->hyper_polarize.k3,
                    iparams->hyper_polarize.k4);
                    break;*/
        case F_ANHARM_POL:
            fprintf(fp, "alpha=%15.8e drcut=%15.8e khyp=%15.8e\n",
                    iparams->anharm_polarize.alpha,
                    iparams->anharm_polarize.drcut,
                    iparams->anharm_polarize.khyp);
            break;
        case F_THOLE_POL:
            fprintf(fp, "a=%15.8e, alpha1=%15.8e, alpha2=%15.8e, rfac=%15.8e\n",
                    iparams->thole.a, iparams->thole.alpha1, iparams->thole.alpha2,
                    iparams->thole.rfac);
            break;
        case F_WATER_POL:
            fprintf(fp, "al_x=%15.8e, al_y=%15.8e, al_z=%15.8e, rOH=%9.6f, rHH=%9.6f, rOD=%9.6f\n",
                    iparams->wpol.al_x, iparams->wpol.al_y, iparams->wpol.al_z,
                    iparams->wpol.rOH, iparams->wpol.rHH, iparams->wpol.rOD);
            break;
        case F_LJ:
            fprintf(fp, "c6=%15.8e, c12=%15.8e\n", iparams->lj.c6, iparams->lj.c12);
            break;
        case F_LJ14:
            fprintf(fp, "c6A=%15.8e, c12A=%15.8e, c6B=%15.8e, c12B=%15.8e\n",
                    iparams->lj14.c6A, iparams->lj14.c12A,
                    iparams->lj14.c6B, iparams->lj14.c12B);
            break;
        case F_LJC14_Q:
            fprintf(fp, "fqq=%15.8e, qi=%15.8e, qj=%15.8e, c6=%15.8e, c12=%15.8e\n",
                    iparams->ljc14.fqq,
                    iparams->ljc14.qi, iparams->ljc14.qj,
                    iparams->ljc14.c6, iparams->ljc14.c12);
            break;
        case F_LJC_PAIRS_NB:
            fprintf(fp, "qi=%15.8e, qj=%15.8e, c6=%15.8e, c12=%15.8e\n",
                    iparams->ljcnb.qi, iparams->ljcnb.qj,
                    iparams->ljcnb.c6, iparams->ljcnb.c12);
            break;
        case F_PDIHS:
        case F_PIDIHS:
            break;
        case F_RBDIHS:
            for (int i = 0; i < NR_RBDIHS; i++)
            {
                fprintf(fp, "%srbcA[%d]=%15.8e", i == 0 ? "" : ", ", i, iparams->rbdihs.rbcA[i]);
            }
            fprintf(fp, "\n");
            for (int i = 0; i < NR_RBDIHS; i++)
            {
                fprintf(fp, "%srbcB[%d]=%15.8e", i == 0 ? "" : ", ", i, iparams->rbdihs.rbcB[i]);
            }
            fprintf(fp, "\n");
            break;
        case F_FOURDIHS:
        {
            /* Use the OPLS -> Ryckaert-Bellemans formula backwards to get
             * the OPLS potential constants back.
             */
            const real *rbcA = iparams->rbdihs.rbcA;
            const real *rbcB = iparams->rbdihs.rbcB;
            real        VA[4], VB[4];

            VA[3] = -0.25*rbcA[4];
            VA[2] = -0.5*rbcA[3];
            VA[1] = 4.0*VA[3]-rbcA[2];
            VA[0] = 3.0*VA[2]-2.0*rbcA[1];

            VB[3] = -0.25*rbcB[4];
            VB[2] = -0.5*rbcB[3];
            VB[1] = 4.0*VB[3]-rbcB[2];
            VB[0] = 3.0*VB[2]-2.0*rbcB[1];

            for (int i = 0; i < NR_FOURDIHS; i++)
            {
                fprintf(fp, "%sFourA[%d]=%15.8e", i == 0 ? "" : ", ", i, VA[i]);
            }
            fprintf(fp, "\n");
            for (int i = 0; i < NR_FOURDIHS; i++)
            {
                fprintf(fp, "%sFourB[%d]=%15.8e", i == 0 ? "" : ", ", i, VB[i]);
            }
            fprintf(fp, "\n");
            break;
        }

        case F_CONSTR:
        case F_CONSTRNC:
            fprintf(fp, "dA=%15.8e, dB=%15.8e\n", iparams->constr.dA, iparams->constr.dB);
            break;
        case F_SETTLE:
            fprintf(fp, "doh=%15.8e, dhh=%15.8e\n", iparams->settle.doh,
                    iparams->settle.dhh);
            break;
        case F_VSITE2:
            fprintf(fp, "a=%15.8e\n", iparams->vsite.a);
            break;
        case F_VSITE3:
        case F_VSITE3FD:
        case F_VSITE3FAD:
            fprintf(fp, "a=%15.8e, b=%15.8e\n", iparams->vsite.a, iparams->vsite.b);
            break;
        case F_VSITE3OUT:
        case F_VSITE4FD:
        case F_VSITE4FDN:
            fprintf(fp, "a=%15.8e, b=%15.8e, c=%15.8e\n",
                    iparams->vsite.a, iparams->vsite.b, iparams->vsite.c);
            break;
        case F_VSITEN:
            fprintf(fp, "n=%2d, a=%15.8e\n", iparams->vsiten.n, iparams->vsiten.a);
            break;
        case F_GB12_NOLONGERUSED:
        case F_GB13_NOLONGERUSED:
        case F_GB14_NOLONGERUSED:
            // These could only be generated by grompp, not written in
            // a .top file. Now that implicit solvent is not
            // supported, they can't be generated, and the values are
            // ignored if read from an old .tpr file. So there is
            // nothing to print.
            break;
        case F_CMAP:
            fprintf(fp, "cmapA=%1d, cmapB=%1d\n", iparams->cmap.cmapA, iparams->cmap.cmapB);
            break;
        case  F_RESTRANGLES:
            pr_harm(fp, iparams, "ktheta", "costheta0");
            break;
        case  F_RESTRDIHS:
            fprintf(fp, "phiA=%15.8e, cpA=%15.8e",
                    iparams->pdihs.phiA, iparams->pdihs.cpA);
            break;
        case  F_CBTDIHS:
            fprintf(fp, "kphi=%15.8e", iparams->cbtdihs.cbtcA[0]);
            for (int i = 1; i < NR_CBTDIHS; i++)
            {
                fprintf(fp, ", cbtcA[%d]=%15.8e", i-1, iparams->cbtdihs.cbtcA[i]);
            }
            fprintf(fp, "\n");
            break;
        case F_COUL_SR:
        case F_DISPERSION:
        case F_REPULSION:
            break;
        default:
            gmx_fatal(FARGS, "unknown function type %d (%s) in %s line %d",
                      ftype, interaction_function[ftype].name, __FILE__, __LINE__);
    }
}

template <typename T>
static void
printIlist(FILE *fp, int indent, const char *title,
           const t_functype *functype, const T &ilist,
           gmx_bool bShowNumbers,
           gmx_bool bShowParameters, const t_iparams *iparams)
{
    int      i, j, k, type, ftype;

    indent = pr_title(fp, indent, title);
    pr_indent(fp, indent);
    fprintf(fp, "nr: %d\n", ilist.size());
    if (ilist.size() > 0)
    {
        pr_indent(fp, indent);
        fprintf(fp, "iatoms:\n");
        for (i = j = 0; i < ilist.size(); )
        {
            pr_indent(fp, indent+INDENT);
            type  = ilist.iatoms[i];
            ftype = functype[type];
            if (bShowNumbers)
            {
                fprintf(fp, "%d type=%d ", j, type);
            }
            j++;
            printf("(%s)", interaction_function[ftype].name);
            for (k = 0; k < interaction_function[ftype].nratoms; k++)
            {
                fprintf(fp, " %3d", ilist.iatoms[i + 1 + k]);
            }
            if (bShowParameters)
            {
                fprintf(fp, "  ");
                pr_iparams(fp, ftype,  &iparams[type]);
            }
            fprintf(fp, "\n");
            i += 1+interaction_function[ftype].nratoms;
        }
    }
}

void pr_ilist(FILE *fp, int indent, const char *title,
              const t_functype *functype, const InteractionList &ilist,
              gmx_bool bShowNumbers,
              gmx_bool bShowParameters, const t_iparams *iparams)
{
    printIlist(fp, indent, title, functype, ilist,
               bShowNumbers, bShowParameters, iparams);
}

void pr_idef(FILE *fp, int indent, const char *title, const t_idef *idef,
             gmx_bool bShowNumbers, gmx_bool bShowParameters)
{
    int i, j;

    if (available(fp, idef, indent, title))
    {
        indent = pr_title(fp, indent, title);
        pr_indent(fp, indent);
        fprintf(fp, "atnr=%d\n", idef->atnr);
        pr_indent(fp, indent);
        fprintf(fp, "ntypes=%d\n", idef->ntypes);
        for (i = 0; i < idef->ntypes; i++)
        {
            pr_indent(fp, indent+INDENT);
            fprintf(fp, "functype[%d]=%s, ",
                    bShowNumbers ? i : -1,
                    interaction_function[idef->functype[i]].name);
            pr_iparams(fp, idef->functype[i], &idef->iparams[i]);
        }
        pr_real(fp, indent, "fudgeQQ", idef->fudgeQQ);

        for (j = 0; (j < F_NRE); j++)
        {
            printIlist(fp, indent, interaction_function[j].longname,
                       idef->functype, idef->il[j], bShowNumbers,
                       bShowParameters, idef->iparams);
        }
    }
}

void init_idef(t_idef *idef)
{
    idef->ntypes           = 0;
    idef->atnr             = 0;
    idef->functype         = nullptr;
    idef->iparams          = nullptr;
    idef->fudgeQQ          = 0.0;
    idef->iparams_posres   = nullptr;
    idef->iparams_fbposres = nullptr;
    for (int f = 0; f < F_NRE; ++f)
    {
        idef->il[f].iatoms          = nullptr;
        idef->il[f].nalloc          = 0;
        idef->il[f].nr              = 0;
        idef->il[f].nr_nonperturbed = 0;
    }
    idef->cmap_grid               = nullptr;
    idef->iparams_posres_nalloc   = 0;
    idef->iparams_fbposres_nalloc = 0;
    idef->ilsort                  = 0;
}

void done_idef(t_idef *idef)
{
    sfree(idef->functype);
    sfree(idef->iparams);
    sfree(idef->iparams_posres);
    sfree(idef->iparams_fbposres);
    for (int f = 0; f < F_NRE; ++f)
    {
        sfree(idef->il[f].iatoms);
    }

    delete idef->cmap_grid;
    init_idef(idef);
}

void copy_ilist(const t_ilist *src, t_ilist *dst)
{
    dst->nr              = src->nr;
    dst->nr_nonperturbed = src->nr_nonperturbed;
    dst->nalloc          = src->nalloc;

    snew(dst->iatoms, dst->nr);
    for (int i = 0; i < dst->nr; ++i)
    {
        dst->iatoms[i] = src->iatoms[i];
    }
}

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
/* This file is completely threadsafe - keep it that way! */
#include "actpre.h"

#include "broadcaststructs.h"

#include <cstring>

#include "gromacs/compat/make_unique.h"
#include "gromacs/gmxlib/network.h"
#include "gromacs/math/vec.h"
#include "gromacs/mdlib/mdrun.h"
#include "gromacs/mdtypes/commrec.h"
#include "gromacs/mdtypes/inputrec.h"
#include "gromacs/mdtypes/md_enums.h"
#include "gromacs/mdtypes/state.h"
#include "gromacs/topology/mtop_util.h"
#include "gromacs/topology/symtab.h"
#include "gromacs/topology/topology.h"
#include "gromacs/utility/fatalerror.h"
#include "gromacs/utility/inmemoryserializer.h"
#include "gromacs/utility/keyvaluetree.h"
#include "gromacs/utility/keyvaluetreeserializer.h"
#include "gromacs/utility/smalloc.h"

static void bc_string(const t_commrec *cr, t_symtab *symtab, char ***s)
{
    int handle;

    if (MASTER(cr))
    {
        handle = lookup_symtab(symtab, *s);
    }
    block_bc(cr, handle);
    if (!MASTER(cr))
    {
        *s = get_symtab_handle(symtab, handle);
    }
}

static void bc_strings(const t_commrec *cr, t_symtab *symtab, int nr, char ****nm)
{
    int     i;
    int    *handle;

    snew(handle, nr);
    if (MASTER(cr))
    {
        for (i = 0; (i < nr); i++)
        {
            handle[i] = lookup_symtab(symtab, (*nm)[i]);
        }
    }
    nblock_bc(cr, nr, handle);

    if (!MASTER(cr))
    {
        snew_bc(cr, *nm, nr);
        for (i = 0; (i < nr); i++)
        {
            (*nm)[i] = get_symtab_handle(symtab, handle[i]);
        }
    }
    sfree(handle);
}

static void bc_strings_resinfo(const t_commrec *cr, t_symtab *symtab,
                               int nr, t_resinfo *resinfo)
{
    int   i;
    int  *handle;

    snew(handle, nr);
    if (MASTER(cr))
    {
        for (i = 0; (i < nr); i++)
        {
            handle[i] = lookup_symtab(symtab, resinfo[i].name);
        }
    }
    nblock_bc(cr, nr, handle);

    if (!MASTER(cr))
    {
        for (i = 0; (i < nr); i++)
        {
            resinfo[i].name = get_symtab_handle(symtab, handle[i]);
        }
    }
    sfree(handle);
}

static void bc_symtab(const t_commrec *cr, t_symtab *symtab)
{
    int       i, nr, len;
    t_symbuf *symbuf;

    block_bc(cr, symtab->nr);
    nr = symtab->nr;
    snew_bc(cr, symtab->symbuf, 1);
    symbuf          = symtab->symbuf;
    symbuf->bufsize = nr;
    snew_bc(cr, symbuf->buf, nr);
    for (i = 0; i < nr; i++)
    {
        if (MASTER(cr))
        {
            len = strlen(symbuf->buf[i]) + 1;
        }
        block_bc(cr, len);
        snew_bc(cr, symbuf->buf[i], len);
        nblock_bc(cr, len, symbuf->buf[i]);
    }
}

static void bc_block(const t_commrec *cr, t_block *block)
{
    block_bc(cr, block->nr);
    snew_bc(cr, block->index, block->nr+1);
    nblock_bc(cr, block->nr+1, block->index);
}

static void bc_blocka(const t_commrec *cr, t_blocka *block)
{
    block_bc(cr, block->nr);
    snew_bc(cr, block->index, block->nr+1);
    nblock_bc(cr, block->nr+1, block->index);
    block_bc(cr, block->nra);
    if (block->nra)
    {
        snew_bc(cr, block->a, block->nra);
        nblock_bc(cr, block->nra, block->a);
    }
}

static void bc_grps(const t_commrec *cr, t_grps grps[])
{
    int i;

    for (i = 0; (i < egcNR); i++)
    {
        block_bc(cr, grps[i].nr);
        snew_bc(cr, grps[i].nm_ind, grps[i].nr);
        nblock_bc(cr, grps[i].nr, grps[i].nm_ind);
    }
}

static void bc_atoms(const t_commrec *cr, t_symtab *symtab, t_atoms *atoms)
{
    block_bc(cr, atoms->nr);
    snew_bc(cr, atoms->atom, atoms->nr);
    nblock_bc(cr, atoms->nr, atoms->atom);
    bc_strings(cr, symtab, atoms->nr, &atoms->atomname);
    block_bc(cr, atoms->nres);
    snew_bc(cr, atoms->resinfo, atoms->nres);
    nblock_bc(cr, atoms->nres, atoms->resinfo);
    bc_strings_resinfo(cr, symtab, atoms->nres, atoms->resinfo);
    /* QMMM requires atomtypes to be known on all nodes as well */
    bc_strings(cr, symtab, atoms->nr, &atoms->atomtype);
    bc_strings(cr, symtab, atoms->nr, &atoms->atomtypeB);
}

static void bc_groups(const t_commrec *cr, t_symtab *symtab,
                      int natoms, gmx_groups_t *groups)
{
    int g, n;

    bc_grps(cr, groups->grps);
    block_bc(cr, groups->ngrpname);
    bc_strings(cr, symtab, groups->ngrpname, &groups->grpname);
    for (g = 0; g < egcNR; g++)
    {
        if (MASTER(cr))
        {
            if (groups->grpnr[g])
            {
                n = natoms;
            }
            else
            {
                n = 0;
            }
        }
        block_bc(cr, n);
        if (n == 0)
        {
            groups->grpnr[g] = nullptr;
        }
        else
        {
            snew_bc(cr, groups->grpnr[g], n);
            nblock_bc(cr, n, groups->grpnr[g]);
        }
    }
    if (debug)
    {
        fprintf(debug, "after bc_groups\n");
    }
}

template <typename AllocatorType>
static void bcastPaddedRVecVector(const t_commrec *cr, gmx::PaddedVector<gmx::RVec, AllocatorType> *v, int numAtoms)
{
    v->resizeWithPadding(numAtoms);
    nblock_bc(cr, makeArrayRef(*v));
}

void broadcastStateWithoutDynamics(const t_commrec *cr, t_state *state)
{
    if (!PAR(cr))
    {
        return;
    }

    /* Broadcasts the state sizes and flags from the master to all ranks
     * in cr->mpi_comm_mygroup.
     */
    block_bc(cr, state->natoms);
    block_bc(cr, state->flags);

    for (int i = 0; i < estNR; i++)
    {
        if (state->flags & (1 << i))
        {
            switch (i)
            {
                case estLAMBDA:
                    nblock_bc(cr, efptNR, state->lambda.data());
                    break;
                case estFEPSTATE:
                    block_bc(cr, state->fep_state);
                    break;
                case estBOX:
                    block_bc(cr, state->box);
                    break;
                case estX:
                    bcastPaddedRVecVector(cr, &state->x, state->natoms);
                    break;
                default:
                    GMX_RELEASE_ASSERT(false, "The state has a dynamic entry, while no dynamic entries should be present");
                    break;
            }
        }
    }
}

static void bc_ilists(const t_commrec *cr, InteractionLists *ilist)
{
    int ftype;

    /* Here we only communicate the non-zero length ilists */
    if (MASTER(cr))
    {
        for (ftype = 0; ftype < F_NRE; ftype++)
        {
            if ((*ilist)[ftype].size() > 0)
            {
                block_bc(cr, ftype);
                int nr = (*ilist)[ftype].size();
                block_bc(cr, nr);
                nblock_bc(cr, nr, (*ilist)[ftype].iatoms.data());
            }
        }
        ftype = -1;
        block_bc(cr, ftype);
    }
    else
    {
        for (ftype = 0; ftype < F_NRE; ftype++)
        {
            (*ilist)[ftype].iatoms.clear();
        }
        do
        {
            block_bc(cr, ftype);
            if (ftype >= 0)
            {
                int nr;
                block_bc(cr, nr);
                (*ilist)[ftype].iatoms.resize(nr);
                nblock_bc(cr, nr, (*ilist)[ftype].iatoms.data());
            }
        }
        while (ftype >= 0);
    }

    if (debug)
    {
        fprintf(debug, "after bc_ilists\n");
    }
}

static void bc_cmap(const t_commrec *cr, gmx_cmap_t *cmap_grid)
{
    int ngrid = cmap_grid->cmapdata.size();
    block_bc(cr, ngrid);
    block_bc(cr, cmap_grid->grid_spacing);

    int nelem = cmap_grid->grid_spacing * cmap_grid->grid_spacing;

    if (ngrid > 0)
    {
        if (!MASTER(cr))
        {
            cmap_grid->cmapdata.resize(ngrid);
        }

        for (int i = 0; i < ngrid; i++)
        {
            nblock_abc(cr, 4*nelem, &cmap_grid->cmapdata[i].cmap);
        }
    }
}

static void bc_ffparams(const t_commrec *cr, gmx_ffparams_t *ffp)
{
    int numTypes = ffp->numTypes();
    block_bc(cr, numTypes);
    block_bc(cr, ffp->atnr);
    nblock_abc(cr, numTypes, &ffp->functype);
    nblock_abc(cr, numTypes, &ffp->iparams);
    block_bc(cr, ffp->reppow);
    block_bc(cr, ffp->fudgeQQ);
    bc_cmap(cr, &ffp->cmap_grid);
}

static void bc_grpopts(const t_commrec *cr, t_grpopts *g)
{
    int i, n;

    block_bc(cr, g->ngtc);
    block_bc(cr, g->ngacc);
    block_bc(cr, g->ngfrz);
    block_bc(cr, g->ngener);
    snew_bc(cr, g->nrdf, g->ngtc);
    snew_bc(cr, g->tau_t, g->ngtc);
    snew_bc(cr, g->ref_t, g->ngtc);
    snew_bc(cr, g->acc, g->ngacc);
    snew_bc(cr, g->nFreeze, g->ngfrz);
    snew_bc(cr, g->egp_flags, g->ngener*g->ngener);

    nblock_bc(cr, g->ngtc, g->nrdf);
    nblock_bc(cr, g->ngtc, g->tau_t);
    nblock_bc(cr, g->ngtc, g->ref_t);
    nblock_bc(cr, g->ngacc, g->acc);
    nblock_bc(cr, g->ngfrz, g->nFreeze);
    nblock_bc(cr, g->ngener*g->ngener, g->egp_flags);
    snew_bc(cr, g->annealing, g->ngtc);
    snew_bc(cr, g->anneal_npoints, g->ngtc);
    snew_bc(cr, g->anneal_time, g->ngtc);
    snew_bc(cr, g->anneal_temp, g->ngtc);
    nblock_bc(cr, g->ngtc, g->annealing);
    nblock_bc(cr, g->ngtc, g->anneal_npoints);
    for (i = 0; (i < g->ngtc); i++)
    {
        n = g->anneal_npoints[i];
        if (n > 0)
        {
            snew_bc(cr, g->anneal_time[i], n);
            snew_bc(cr, g->anneal_temp[i], n);
            nblock_bc(cr, n, g->anneal_time[i]);
            nblock_bc(cr, n, g->anneal_temp[i]);
        }
    }
}

static void bc_fepvals(const t_commrec *cr, t_lambda *fep)
{
    int      i;

    block_bc(cr, fep->nstdhdl);
    block_bc(cr, fep->init_lambda);
    block_bc(cr, fep->init_fep_state);
    block_bc(cr, fep->delta_lambda);
    block_bc(cr, fep->edHdLPrintEnergy);
    block_bc(cr, fep->n_lambda);
    if (fep->n_lambda > 0)
    {
        snew_bc(cr, fep->all_lambda, efptNR);
        nblock_bc(cr, efptNR, fep->all_lambda);
        for (i = 0; i < efptNR; i++)
        {
            snew_bc(cr, fep->all_lambda[i], fep->n_lambda);
            nblock_bc(cr, fep->n_lambda, fep->all_lambda[i]);
        }
    }
    block_bc(cr, fep->sc_alpha);
    block_bc(cr, fep->sc_power);
    block_bc(cr, fep->sc_r_power);
    block_bc(cr, fep->sc_sigma);
    block_bc(cr, fep->sc_sigma_min);
    block_bc(cr, fep->bScCoul);
    nblock_bc(cr, efptNR, &(fep->separate_dvdl[0]));
    block_bc(cr, fep->dhdl_derivatives);
    block_bc(cr, fep->dh_hist_size);
    block_bc(cr, fep->dh_hist_spacing);
    if (debug)
    {
        fprintf(debug, "after bc_fepvals\n");
    }
}

static void bc_expandedvals(const t_commrec *cr, t_expanded *expand, int n_lambda)
{
    block_bc(cr, expand->nstexpanded);
    block_bc(cr, expand->elamstats);
    block_bc(cr, expand->elmcmove);
    block_bc(cr, expand->elmceq);
    block_bc(cr, expand->equil_n_at_lam);
    block_bc(cr, expand->equil_wl_delta);
    block_bc(cr, expand->equil_ratio);
    block_bc(cr, expand->equil_steps);
    block_bc(cr, expand->equil_samples);
    block_bc(cr, expand->lmc_seed);
    block_bc(cr, expand->minvar);
    block_bc(cr, expand->minvar_const);
    block_bc(cr, expand->c_range);
    block_bc(cr, expand->bSymmetrizedTMatrix);
    block_bc(cr, expand->nstTij);
    block_bc(cr, expand->lmc_repeats);
    block_bc(cr, expand->lmc_forced_nstart);
    block_bc(cr, expand->gibbsdeltalam);
    block_bc(cr, expand->wl_scale);
    block_bc(cr, expand->wl_ratio);
    block_bc(cr, expand->init_wl_delta);
    block_bc(cr, expand->bInit_weights);
    snew_bc(cr, expand->init_lambda_weights, n_lambda);
    nblock_bc(cr, n_lambda, expand->init_lambda_weights);
    block_bc(cr, expand->mc_temp);
    if (debug)
    {
        fprintf(debug, "after bc_expandedvals\n");
    }
}

static void bc_simtempvals(const t_commrec *cr, t_simtemp *simtemp, int n_lambda)
{
    block_bc(cr, simtemp->simtemp_low);
    block_bc(cr, simtemp->simtemp_high);
    block_bc(cr, simtemp->eSimTempScale);
    snew_bc(cr, simtemp->temperatures, n_lambda);
    nblock_bc(cr, n_lambda, simtemp->temperatures);
    if (debug)
    {
        fprintf(debug, "after bc_simtempvals\n");
    }
}

static void bc_inputrec(const t_commrec *cr, t_inputrec *inputrec)
{
    // Note that this overwrites pointers in inputrec, so all pointer fields
    // Must be initialized separately below.
    block_bc(cr, *inputrec);
    if (SIMMASTER(cr))
    {
        gmx::InMemorySerializer serializer;
        gmx::serializeKeyValueTree(*inputrec->params, &serializer);
        std::vector<char>       buffer = serializer.finishAndGetBuffer();
        size_t                  size   = buffer.size();
        block_bc(cr, size);
        nblock_bc(cr, size, buffer.data());
    }
    else
    {
        // block_bc() above overwrites the old pointer, so set it to a
        // reasonable value in case code below throws.
        inputrec->params = nullptr;
        std::vector<char> buffer;
        size_t            size;
        block_bc(cr, size);
        nblock_abc(cr, size, &buffer);
        gmx::InMemoryDeserializer serializer(buffer);
        inputrec->params = new gmx::KeyValueTreeObject(
                    gmx::deserializeKeyValueTree(&serializer));
    }

    bc_grpopts(cr, &(inputrec->opts));

    /* even if efep is efepNO, we need to initialize to make sure that
     * n_lambda is set to zero */

    snew_bc(cr, inputrec->fepvals, 1);
    if (inputrec->efep != efepNO || inputrec->bSimTemp)
    {
        bc_fepvals(cr, inputrec->fepvals);
    }
    /* need to initialize this as well because of data checked for in the logic */
    snew_bc(cr, inputrec->expandedvals, 1);
    if (inputrec->bExpanded)
    {
        bc_expandedvals(cr, inputrec->expandedvals, inputrec->fepvals->n_lambda);
    }
    snew_bc(cr, inputrec->simtempvals, 1);
    if (inputrec->bSimTemp)
    {
        bc_simtempvals(cr, inputrec->simtempvals, inputrec->fepvals->n_lambda);
    }
}

static void bc_moltype(const t_commrec *cr, t_symtab *symtab,
                       gmx_moltype_t *moltype)
{
    bc_string(cr, symtab, &moltype->name);
    bc_atoms(cr, symtab, &moltype->atoms);
    if (debug)
    {
        fprintf(debug, "after bc_atoms\n");
    }

    bc_ilists(cr, &moltype->ilist);
    bc_block(cr, &moltype->cgs);
    bc_blocka(cr, &moltype->excls);
}

static void bc_vector_of_rvec(const t_commrec *cr, std::vector<gmx::RVec> *vec)
{
    int numElements = vec->size();
    block_bc(cr, numElements);
    if (!MASTER(cr))
    {
        vec->resize(numElements);
    }
    if (numElements > 0)
    {
        nblock_bc(cr, numElements, as_rvec_array(vec->data()));
    }
}

static void bc_molblock(const t_commrec *cr, gmx_molblock_t *molb)
{
    block_bc(cr, molb->type);
    block_bc(cr, molb->nmol);
    bc_vector_of_rvec(cr, &molb->posres_xA);
    bc_vector_of_rvec(cr, &molb->posres_xB);
    if (debug)
    {
        fprintf(debug, "after bc_molblock\n");
    }
}

static void bc_atomtypes(const t_commrec *cr, t_atomtypes *atomtypes)
{
    block_bc(cr, atomtypes->nr);
}

/*! \brief Broadcasts ir and mtop from the master to all nodes in
 * cr->mpi_comm_mygroup. */
static
void bcast_ir_mtop(const t_commrec *cr, t_inputrec *inputrec, gmx_mtop_t *mtop)
{
    if (debug)
    {
        fprintf(debug, "in bc_data\n");
    }
    bc_inputrec(cr, inputrec);
    if (debug)
    {
        fprintf(debug, "after bc_inputrec\n");
    }
    bc_symtab(cr, &mtop->symtab);
    if (debug)
    {
        fprintf(debug, "after bc_symtab\n");
    }
    bc_string(cr, &mtop->symtab, &mtop->name);
    if (debug)
    {
        fprintf(debug, "after bc_name\n");
    }

    bc_ffparams(cr, &mtop->ffparams);

    int nmoltype = mtop->moltype.size();
    block_bc(cr, nmoltype);
    mtop->moltype.resize(nmoltype);
    for (gmx_moltype_t &moltype : mtop->moltype)
    {
        bc_moltype(cr, &mtop->symtab, &moltype);
    }

    block_bc(cr, mtop->bIntermolecularInteractions);
    if (mtop->bIntermolecularInteractions)
    {
        mtop->intermolecular_ilist = gmx::compat::make_unique<InteractionLists>();
        bc_ilists(cr, mtop->intermolecular_ilist.get());
    }

    int nmolblock = mtop->molblock.size();
    block_bc(cr, nmolblock);
    mtop->molblock.resize(nmolblock);
    for (gmx_molblock_t &molblock : mtop->molblock)
    {
        bc_molblock(cr, &molblock);
    }

    block_bc(cr, mtop->natoms);

    bc_atomtypes(cr, &mtop->atomtypes);

    bc_groups(cr, &mtop->symtab, mtop->natoms, &mtop->groups);

    GMX_RELEASE_ASSERT(!MASTER(cr) || mtop->haveMoleculeIndices, "mtop should have valid molecule indices");
    if (!MASTER(cr))
    {
        mtop->haveMoleculeIndices = true;

        gmx_mtop_finalize(mtop);
    }
}

void init_parallel(t_commrec *cr, t_inputrec *inputrec,
                   gmx_mtop_t *mtop)
{
    bcast_ir_mtop(cr, inputrec, mtop);
}

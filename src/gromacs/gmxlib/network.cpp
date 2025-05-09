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

#include "network.h"

#include "config.h"

#include <cctype>
#include <cstdarg>
#include <cstdlib>
#include <cstring>

#include "gromacs/commandline/filenm.h"
#include "gromacs/mdtypes/commrec.h"
#include "gromacs/utility/basenetwork.h"
#include "gromacs/utility/cstringutil.h"
#include "gromacs/utility/fatalerror.h"
#include "gromacs/utility/futil.h"
#include "gromacs/utility/gmxmpi.h"
#include "gromacs/utility/smalloc.h"

/* The source code in this file should be thread-safe.
      Please keep it that way. */

void gmx_fill_commrec_from_mpi(t_commrec *cr)
{
#if !GMX_MPI
    gmx_call("gmx_fill_commrec_from_mpi");
    GMX_UNUSED_VALUE(cr);
#else
    if (!gmx_mpi_initialized())
    {
        gmx_comm("MPI has not been initialized properly");
    }

    cr->nnodes           = gmx_node_num();
    cr->nodeid           = gmx_node_rank();
    cr->sim_nodeid       = cr->nodeid;
    cr->mpi_comm_mysim   = MPI_COMM_WORLD;
    cr->mpi_comm_mygroup = MPI_COMM_WORLD;
#endif
}

t_commrec *init_commrec()
{
    t_commrec    *cr;

    snew(cr, 1);

#if GMX_LIB_MPI
    gmx_fill_commrec_from_mpi(cr);
#else
    cr->mpi_comm_mysim     = MPI_COMM_NULL;
    cr->mpi_comm_mygroup   = MPI_COMM_NULL;
    //    cr->mpi_act_helpers    = MPI_COMM_NULL;
    //cr->mpi_act_not_master = MPI_COMM_NULL;
    cr->nnodes             = 1;
    cr->sim_nodeid         = 0;
    cr->nodeid             = cr->sim_nodeid;
#endif

    // TODO cr->duty should not be initialized here
    cr->duty = (DUTY_PP | DUTY_PME);

#if GMX_MPI && !MPI_IN_PLACE_EXISTS
    /* initialize the MPI_IN_PLACE replacement buffers */
    snew(cr->mpb, 1);
    cr->mpb->ibuf        = NULL;
    cr->mpb->libuf       = NULL;
    cr->mpb->fbuf        = NULL;
    cr->mpb->dbuf        = NULL;
    cr->mpb->ibuf_alloc  = 0;
    cr->mpb->libuf_alloc = 0;
    cr->mpb->fbuf_alloc  = 0;
    cr->mpb->dbuf_alloc  = 0;
#endif

    return cr;
}

void done_mpi_in_place_buf(mpi_in_place_buf_t *buf)
{
    if (nullptr != buf)
    {
        sfree(buf->ibuf);
        sfree(buf->libuf);
        sfree(buf->fbuf);
        sfree(buf->dbuf);
        sfree(buf);
    }
}

void done_commrec(t_commrec *cr)
{
    done_mpi_in_place_buf(cr->mpb);
    sfree(cr);
}

t_commrec *reinitialize_commrec_for_this_thread(const t_commrec *cro)
{
#if GMX_THREAD_MPI
    t_commrec *cr;

    /* make a thread-specific commrec */
    snew(cr, 1);
    /* now copy the whole thing, so settings like the number of PME nodes
       get propagated. */
    *cr = *cro;

    /* and we start setting our own thread-specific values for things */
    gmx_fill_commrec_from_mpi(cr);

    // TODO cr->duty should not be initialized here
    cr->duty             = (DUTY_PP | DUTY_PME);

    return cr;
#else
    GMX_UNUSED_VALUE(cro);
    return nullptr;
#endif
}

void gmx_setup_nodecomm(FILE gmx_unused *fplog, t_commrec *cr)
{
    gmx_nodecomm_t *nc;

    /* Many MPI implementations do not optimize MPI_Allreduce
     * (and probably also other global communication calls)
     * for multi-core nodes connected by a network.
     * We can optimize such communication by using one MPI call
     * within each node and one between the nodes.
     * For MVAPICH2 and Intel MPI this reduces the time for
     * the global_stat communication by 25%
     * for 2x2-core 3 GHz Woodcrest connected by mixed DDR/SDR Infiniband.
     * B. Hess, November 2007
     */

    nc = &cr->nc;

    nc->bUse = FALSE;
#if !GMX_THREAD_MPI
#if GMX_MPI
    int n, rank;

    // TODO PhysicalNodeCommunicator could be extended/used to handle
    // the need for per-node per-group communicators.
    MPI_Comm_size(cr->mpi_comm_mygroup, &n);
    MPI_Comm_rank(cr->mpi_comm_mygroup, &rank);

    int nodehash = gmx_physicalnode_id_hash();

    if (debug)
    {
        fprintf(debug, "In gmx_setup_nodecomm: splitting communicator of size %d\n", n);
    }


    /* The intra-node communicator, split on node number */
    MPI_Comm_split(cr->mpi_comm_mygroup, nodehash, rank, &nc->comm_intra);
    MPI_Comm_rank(nc->comm_intra, &nc->rank_intra);
    if (debug)
    {
        fprintf(debug, "In gmx_setup_nodecomm: node ID %d rank within node %d\n",
                rank, nc->rank_intra);
    }
    /* The inter-node communicator, split on rank_intra.
     * We actually only need the one for rank=0,
     * but it is easier to create them all.
     */
    MPI_Comm_split(cr->mpi_comm_mygroup, nc->rank_intra, rank, &nc->comm_inter);
    /* Check if this really created two step communication */
    int ng, ni;

    MPI_Comm_size(nc->comm_inter, &ng);
    MPI_Comm_size(nc->comm_intra, &ni);
    if (debug)
    {
        fprintf(debug, "In gmx_setup_nodecomm: groups %d, my group size %d\n",
                ng, ni);
    }

    if (getenv("GMX_NO_NODECOMM") == nullptr &&
        ((ng > 1 && ng < n) || (ni > 1 && ni < n)))
    {
        nc->bUse = TRUE;
        if (fplog)
        {
            fprintf(fplog, "Using two step summing over %d groups of on average %.1f ranks\n\n",
                    ng, (real)n/(real)ng);
        }
        if (nc->rank_intra > 0)
        {
            MPI_Comm_free(&nc->comm_inter);
        }
    }
    else
    {
        /* One group or all processes in a separate group, use normal summing */
        MPI_Comm_free(&nc->comm_inter);
        MPI_Comm_free(&nc->comm_intra);
        if (debug)
        {
            fprintf(debug, "In gmx_setup_nodecomm: not unsing separate inter- and intra-node communicators.\n");
        }
    }
#endif
#else
    /* tMPI runs only on a single node so just use the nodeid */
    nc->rank_intra = cr->nodeid;
#endif
}

void gmx_barrier(const t_commrec gmx_unused *cr)
{
#if !GMX_MPI
    gmx_call("gmx_barrier");
#else
    MPI_Barrier(cr->mpi_comm_mygroup);
#endif
}

void gmx_bcast(int gmx_unused nbytes, void gmx_unused *b, const t_commrec gmx_unused *cr)
{
#if !GMX_MPI
    gmx_call("gmx_bast");
#else
    MPI_Bcast(b, nbytes, MPI_BYTE, MASTERRANK(cr), cr->mpi_comm_mygroup);
#endif
}

void gmx_bcast_sim(int gmx_unused nbytes, void gmx_unused *b, const t_commrec gmx_unused *cr)
{
#if !GMX_MPI
    gmx_call("gmx_bast");
#else
    MPI_Bcast(b, nbytes, MPI_BYTE, MASTERRANK(cr), cr->mpi_comm_mysim);
#endif
}

void gmx_sumd(int gmx_unused nr, double gmx_unused r[], const t_commrec gmx_unused *cr)
{
#if !GMX_MPI
    gmx_call("gmx_sumd");
#else
#if MPI_IN_PLACE_EXISTS
    if (cr->nc.bUse)
    {
        if (cr->nc.rank_intra == 0)
        {
            /* Use two step summing. */
            MPI_Reduce(MPI_IN_PLACE, r, nr, MPI_DOUBLE, MPI_SUM, 0,
                       cr->nc.comm_intra);
            /* Sum the roots of the internal (intra) buffers. */
            MPI_Allreduce(MPI_IN_PLACE, r, nr, MPI_DOUBLE, MPI_SUM,
                          cr->nc.comm_inter);
        }
        else
        {
            /* This is here because of the silly MPI specification
                that MPI_IN_PLACE should be put in sendbuf instead of recvbuf */
            MPI_Reduce(r, nullptr, nr, MPI_DOUBLE, MPI_SUM, 0, cr->nc.comm_intra);
        }
        MPI_Bcast(r, nr, MPI_DOUBLE, 0, cr->nc.comm_intra);
    }
    else
    {
        MPI_Allreduce(MPI_IN_PLACE, r, nr, MPI_DOUBLE, MPI_SUM,
                      cr->mpi_comm_mygroup);
    }
#else
    int i;

    if (nr > cr->mpb->dbuf_alloc)
    {
        cr->mpb->dbuf_alloc = nr;
        srenew(cr->mpb->dbuf, cr->mpb->dbuf_alloc);
    }
    if (cr->nc.bUse)
    {
        /* Use two step summing */
        MPI_Allreduce(r, cr->mpb->dbuf, nr, MPI_DOUBLE, MPI_SUM, cr->nc.comm_intra);
        if (cr->nc.rank_intra == 0)
        {
            /* Sum with the buffers reversed */
            MPI_Allreduce(cr->mpb->dbuf, r, nr, MPI_DOUBLE, MPI_SUM,
                          cr->nc.comm_inter);
        }
        MPI_Bcast(r, nr, MPI_DOUBLE, 0, cr->nc.comm_intra);
    }
    else
    {
        MPI_Allreduce(r, cr->mpb->dbuf, nr, MPI_DOUBLE, MPI_SUM,
                      cr->mpi_comm_mygroup);
        for (i = 0; i < nr; i++)
        {
            r[i] = cr->mpb->dbuf[i];
        }
    }
#endif
#endif
}

void gmx_sumf(int gmx_unused nr, float gmx_unused r[], const t_commrec gmx_unused *cr)
{
#if !GMX_MPI
    gmx_call("gmx_sumf");
#else
#if MPI_IN_PLACE_EXISTS
    if (cr->nc.bUse)
    {
        /* Use two step summing.  */
        if (cr->nc.rank_intra == 0)
        {
            MPI_Reduce(MPI_IN_PLACE, r, nr, MPI_FLOAT, MPI_SUM, 0,
                       cr->nc.comm_intra);
            /* Sum the roots of the internal (intra) buffers */
            MPI_Allreduce(MPI_IN_PLACE, r, nr, MPI_FLOAT, MPI_SUM,
                          cr->nc.comm_inter);
        }
        else
        {
            /* This is here because of the silly MPI specification
                that MPI_IN_PLACE should be put in sendbuf instead of recvbuf */
            MPI_Reduce(r, nullptr, nr, MPI_FLOAT, MPI_SUM, 0, cr->nc.comm_intra);
        }
        MPI_Bcast(r, nr, MPI_FLOAT, 0, cr->nc.comm_intra);
    }
    else
    {
        MPI_Allreduce(MPI_IN_PLACE, r, nr, MPI_FLOAT, MPI_SUM, cr->mpi_comm_mygroup);
    }
#else
    int i;

    if (nr > cr->mpb->fbuf_alloc)
    {
        cr->mpb->fbuf_alloc = nr;
        srenew(cr->mpb->fbuf, cr->mpb->fbuf_alloc);
    }
    if (cr->nc.bUse)
    {
        /* Use two step summing */
        MPI_Allreduce(r, cr->mpb->fbuf, nr, MPI_FLOAT, MPI_SUM, cr->nc.comm_intra);
        if (cr->nc.rank_intra == 0)
        {
            /* Sum with the buffers reversed */
            MPI_Allreduce(cr->mpb->fbuf, r, nr, MPI_FLOAT, MPI_SUM,
                          cr->nc.comm_inter);
        }
        MPI_Bcast(r, nr, MPI_FLOAT, 0, cr->nc.comm_intra);
    }
    else
    {
        MPI_Allreduce(r, cr->mpb->fbuf, nr, MPI_FLOAT, MPI_SUM,
                      cr->mpi_comm_mygroup);
        for (i = 0; i < nr; i++)
        {
            r[i] = cr->mpb->fbuf[i];
        }
    }
#endif
#endif
}

void gmx_sumi(int gmx_unused nr, int gmx_unused r[], const t_commrec gmx_unused *cr)
{
#if !GMX_MPI
    gmx_call("gmx_sumi");
#else
#if MPI_IN_PLACE_EXISTS
    if (cr->nc.bUse)
    {
        /* Use two step summing */
        if (cr->nc.rank_intra == 0)
        {
            MPI_Reduce(MPI_IN_PLACE, r, nr, MPI_INT, MPI_SUM, 0, cr->nc.comm_intra);
            /* Sum with the buffers reversed */
            MPI_Allreduce(MPI_IN_PLACE, r, nr, MPI_INT, MPI_SUM, cr->nc.comm_inter);
        }
        else
        {
            /* This is here because of the silly MPI specification
                that MPI_IN_PLACE should be put in sendbuf instead of recvbuf */
            MPI_Reduce(r, nullptr, nr, MPI_INT, MPI_SUM, 0, cr->nc.comm_intra);
        }
        MPI_Bcast(r, nr, MPI_INT, 0, cr->nc.comm_intra);
    }
    else
    {
        MPI_Allreduce(MPI_IN_PLACE, r, nr, MPI_INT, MPI_SUM, cr->mpi_comm_mygroup);
    }
#else
    int i;

    if (nr > cr->mpb->ibuf_alloc)
    {
        cr->mpb->ibuf_alloc = nr;
        srenew(cr->mpb->ibuf, cr->mpb->ibuf_alloc);
    }
    if (cr->nc.bUse)
    {
        /* Use two step summing */
        MPI_Allreduce(r, cr->mpb->ibuf, nr, MPI_INT, MPI_SUM, cr->nc.comm_intra);
        if (cr->nc.rank_intra == 0)
        {
            /* Sum with the buffers reversed */
            MPI_Allreduce(cr->mpb->ibuf, r, nr, MPI_INT, MPI_SUM, cr->nc.comm_inter);
        }
        MPI_Bcast(r, nr, MPI_INT, 0, cr->nc.comm_intra);
    }
    else
    {
        MPI_Allreduce(r, cr->mpb->ibuf, nr, MPI_INT, MPI_SUM, cr->mpi_comm_mygroup);
        for (i = 0; i < nr; i++)
        {
            r[i] = cr->mpb->ibuf[i];
        }
    }
#endif
#endif
}

void gmx_sumli(int gmx_unused nr, int64_t gmx_unused r[], const t_commrec gmx_unused *cr)
{
#if !GMX_MPI
    gmx_call("gmx_sumli");
#else
#if MPI_IN_PLACE_EXISTS
    if (cr->nc.bUse)
    {
        /* Use two step summing */
        if (cr->nc.rank_intra == 0)
        {
            MPI_Reduce(MPI_IN_PLACE, r, nr, MPI_INT64_T, MPI_SUM, 0,
                       cr->nc.comm_intra);
            /* Sum with the buffers reversed */
            MPI_Allreduce(MPI_IN_PLACE, r, nr, MPI_INT64_T, MPI_SUM,
                          cr->nc.comm_inter);
        }
        else
        {
            /* This is here because of the silly MPI specification
                that MPI_IN_PLACE should be put in sendbuf instead of recvbuf */
            MPI_Reduce(r, nullptr, nr, MPI_INT64_T, MPI_SUM, 0, cr->nc.comm_intra);
        }
        MPI_Bcast(r, nr, MPI_INT64_T, 0, cr->nc.comm_intra);
    }
    else
    {
        MPI_Allreduce(MPI_IN_PLACE, r, nr, MPI_INT64_T, MPI_SUM, cr->mpi_comm_mygroup);
    }
#else
    int i;

    if (nr > cr->mpb->libuf_alloc)
    {
        cr->mpb->libuf_alloc = nr;
        srenew(cr->mpb->libuf, cr->mpb->libuf_alloc);
    }
    if (cr->nc.bUse)
    {
        /* Use two step summing */
        MPI_Allreduce(r, cr->mpb->libuf, nr, MPI_INT64_T, MPI_SUM,
                      cr->nc.comm_intra);
        if (cr->nc.rank_intra == 0)
        {
            /* Sum with the buffers reversed */
            MPI_Allreduce(cr->mpb->libuf, r, nr, MPI_INT64_T, MPI_SUM,
                          cr->nc.comm_inter);
        }
        MPI_Bcast(r, nr, MPI_INT64_T, 0, cr->nc.comm_intra);
    }
    else
    {
        MPI_Allreduce(r, cr->mpb->libuf, nr, MPI_INT64_T, MPI_SUM,
                      cr->mpi_comm_mygroup);
        for (i = 0; i < nr; i++)
        {
            r[i] = cr->mpb->libuf[i];
        }
    }
#endif
#endif
}



#if GMX_MPI
static void gmx_sumd_comm(int nr, double r[], MPI_Comm mpi_comm)
{
#if MPI_IN_PLACE_EXISTS
    MPI_Allreduce(MPI_IN_PLACE, r, nr, MPI_DOUBLE, MPI_SUM, mpi_comm);
#else
    /* this function is only used in code that is not performance critical,
       (during setup, when comm_rec is not the appropriate communication
       structure), so this isn't as bad as it looks. */
    double *buf;
    int     i;

    snew(buf, nr);
    MPI_Allreduce(r, buf, nr, MPI_DOUBLE, MPI_SUM, mpi_comm);
    for (i = 0; i < nr; i++)
    {
        r[i] = buf[i];
    }
    sfree(buf);
#endif
}
#endif

#if GMX_MPI
static void gmx_sumf_comm(int nr, float r[], MPI_Comm mpi_comm)
{
#if MPI_IN_PLACE_EXISTS
    MPI_Allreduce(MPI_IN_PLACE, r, nr, MPI_FLOAT, MPI_SUM, mpi_comm);
#else
    /* this function is only used in code that is not performance critical,
       (during setup, when comm_rec is not the appropriate communication
       structure), so this isn't as bad as it looks. */
    float *buf;
    int    i;

    snew(buf, nr);
    MPI_Allreduce(r, buf, nr, MPI_FLOAT, MPI_SUM, mpi_comm);
    for (i = 0; i < nr; i++)
    {
        r[i] = buf[i];
    }
    sfree(buf);
#endif
}
#endif

void gmx_sumd_sim(int gmx_unused nr, double gmx_unused r[], const gmx_multisim_t gmx_unused *ms)
{
#if !GMX_MPI
    gmx_call("gmx_sumd_sim");
#else
    gmx_sumd_comm(nr, r, ms->mpi_comm_masters);
#endif
}

void gmx_sumf_sim(int gmx_unused nr, float gmx_unused r[], const gmx_multisim_t gmx_unused *ms)
{
#if !GMX_MPI
    gmx_call("gmx_sumf_sim");
#else
    gmx_sumf_comm(nr, r, ms->mpi_comm_masters);
#endif
}

void gmx_sumi_sim(int gmx_unused nr, int gmx_unused r[], const gmx_multisim_t gmx_unused *ms)
{
#if !GMX_MPI
    gmx_call("gmx_sumi_sim");
#else
#if MPI_IN_PLACE_EXISTS
    MPI_Allreduce(MPI_IN_PLACE, r, nr, MPI_INT, MPI_SUM, ms->mpi_comm_masters);
#else
    /* this is thread-unsafe, but it will do for now: */
    int i;

    if (nr > ms->mpb->ibuf_alloc)
    {
        ms->mpb->ibuf_alloc = nr;
        srenew(ms->mpb->ibuf, ms->mpb->ibuf_alloc);
    }
    MPI_Allreduce(r, ms->mpb->ibuf, nr, MPI_INT, MPI_SUM, ms->mpi_comm_masters);
    for (i = 0; i < nr; i++)
    {
        r[i] = ms->mpb->ibuf[i];
    }
#endif
#endif
}

void gmx_sumli_sim(int gmx_unused nr, int64_t gmx_unused r[], const gmx_multisim_t gmx_unused *ms)
{
#if !GMX_MPI
    gmx_call("gmx_sumli_sim");
#else
#if MPI_IN_PLACE_EXISTS
    MPI_Allreduce(MPI_IN_PLACE, r, nr, MPI_INT64_T, MPI_SUM,
                  ms->mpi_comm_masters);
#else
    /* this is thread-unsafe, but it will do for now: */
    int i;

    if (nr > ms->mpb->libuf_alloc)
    {
        ms->mpb->libuf_alloc = nr;
        srenew(ms->mpb->libuf, ms->mpb->libuf_alloc);
    }
    MPI_Allreduce(r, ms->mpb->libuf, nr, MPI_INT64_T, MPI_SUM,
                  ms->mpi_comm_masters);
    for (i = 0; i < nr; i++)
    {
        r[i] = ms->mpb->libuf[i];
    }
#endif
#endif
}

void check_multi_int(FILE *log, const gmx_multisim_t *ms, int val,
                     const char *name,
                     gmx_bool bQuiet)
{
    int     *ibuf, p;
    gmx_bool bCompatible;

    if (nullptr != log && !bQuiet)
    {
        fprintf(log, "Multi-checking %s ... ", name);
    }

    if (ms == nullptr)
    {
        gmx_fatal(FARGS,
                  "check_multi_int called with a NULL communication pointer");
    }

    snew(ibuf, ms->nsim);
    ibuf[ms->sim] = val;
    gmx_sumi_sim(ms->nsim, ibuf, ms);

    bCompatible = TRUE;
    for (p = 1; p < ms->nsim; p++)
    {
        bCompatible = bCompatible && (ibuf[p-1] == ibuf[p]);
    }

    if (bCompatible)
    {
        if (nullptr != log && !bQuiet)
        {
            fprintf(log, "OK\n");
        }
    }
    else
    {
        if (nullptr != log)
        {
            fprintf(log, "\n%s is not equal for all subsystems\n", name);
            for (p = 0; p < ms->nsim; p++)
            {
                fprintf(log, "  subsystem %d: %d\n", p, ibuf[p]);
            }
        }
        gmx_fatal(FARGS, "The %d subsystems are not compatible\n", ms->nsim);
    }

    sfree(ibuf);
}

void check_multi_int64(FILE *log, const gmx_multisim_t *ms,
                       int64_t val, const char *name,
                       gmx_bool bQuiet)
{
    int64_t          *ibuf;
    int               p;
    gmx_bool          bCompatible;

    if (nullptr != log && !bQuiet)
    {
        fprintf(log, "Multi-checking %s ... ", name);
    }

    if (ms == nullptr)
    {
        gmx_fatal(FARGS,
                  "check_multi_int called with a NULL communication pointer");
    }

    snew(ibuf, ms->nsim);
    ibuf[ms->sim] = val;
    gmx_sumli_sim(ms->nsim, ibuf, ms);

    bCompatible = TRUE;
    for (p = 1; p < ms->nsim; p++)
    {
        bCompatible = bCompatible && (ibuf[p-1] == ibuf[p]);
    }

    if (bCompatible)
    {
        if (nullptr != log && !bQuiet)
        {
            fprintf(log, "OK\n");
        }
    }
    else
    {
        // TODO Part of this error message would also be good to go to
        // stderr (from one rank of one sim only)
        if (nullptr != log)
        {
            fprintf(log, "\n%s is not equal for all subsystems\n", name);
            for (p = 0; p < ms->nsim; p++)
            {
                char strbuf[255];
                /* first make the format string */
                snprintf(strbuf, 255, "  subsystem %%d: %s\n",
                         "%" PRId64);
                fprintf(log, strbuf, p, ibuf[p]);
            }
        }
        gmx_fatal(FARGS, "The %d subsystems are not compatible\n", ms->nsim);
    }

    sfree(ibuf);
}

const char *opt2fn_master(const char *opt, int nfile, const t_filenm fnm[],
                          t_commrec *cr)
{
    return SIMMASTER(cr) ? opt2fn(opt, nfile, fnm) : nullptr;
}

void gmx_fatal_collective(int f_errno, const char *file, int line,
                          MPI_Comm comm, gmx_bool bMaster,
                          gmx_fmtstr const char *fmt, ...)
{
    va_list  ap;
    gmx_bool bFinalize;
#if GMX_MPI
    int      result;
    /* Check if we are calling on all processes in MPI_COMM_WORLD */
    MPI_Comm_compare(comm, MPI_COMM_WORLD, &result);
    /* Any result except MPI_UNEQUAL allows us to call MPI_Finalize */
    bFinalize = (result != MPI_UNEQUAL);
#else
    GMX_UNUSED_VALUE(comm);
    bFinalize = TRUE;
#endif

    va_start(ap, fmt);
    gmx_fatal_mpi_va(f_errno, file, line, bMaster, bFinalize, fmt, ap);
    va_end(ap);
}

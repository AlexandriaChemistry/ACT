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

#include "force.h"

#include "config.h"

#include <cassert>
#include <cmath>
#include <cstring>

#include "gromacs/gmxlib/network.h"
#include "gromacs/gmxlib/nrnb.h"
#include "gromacs/gmxlib/nonbonded/nonbonded.h"
#include "gromacs/listed-forces/listed-forces.h"
#include "gromacs/math/vec.h"
#include "gromacs/math/vecdump.h"
#include "gromacs/mdlib/force_flags.h"
#include "gromacs/mdlib/forcerec-threading.h"
#include "gromacs/mdlib/mdrun.h"
#include "gromacs/mdlib/ns.h"
#include "gromacs/mdlib/rf_util.h"
#include "gromacs/mdtypes/commrec.h"
#include "gromacs/mdtypes/enerdata.h"
#include "gromacs/mdtypes/forceoutput.h"
#include "gromacs/mdtypes/forcerec.h"
#include "gromacs/mdtypes/inputrec.h"
#include "gromacs/mdtypes/md_enums.h"
#include "gromacs/pbcutil/ishift.h"
#include "gromacs/pbcutil/mshift.h"
#include "gromacs/pbcutil/pbc.h"
#include "gromacs/timing/wallcycle.h"
#include "gromacs/utility/cstringutil.h"
#include "gromacs/utility/exceptions.h"
#include "gromacs/utility/fatalerror.h"
#include "gromacs/utility/smalloc.h"

void ns(FILE               *fp,
        t_forcerec         *fr,
        matrix              box,
        const gmx_groups_t *groups,
        gmx_localtop_t     *top,
        const t_mdatoms    *md,
        t_nrnb             *nrnb,
        gmx_bool            bFillGrid)
{
    int     nsearch;


    if (!fr->ns->nblist_initialized)
    {
        init_neighbor_list(fp, fr, md->homenr);
    }

    nsearch = search_neighbours(fp, fr, box, top, groups, nrnb, md,
                                bFillGrid);
    if (debug)
    {
        fprintf(debug, "nsearch = %d\n", nsearch);
    }

    /* Check whether we have to do dynamic load balancing */
    /*if ((nsb->nstDlb > 0) && (mod(step,nsb->nstDlb) == 0))
       count_nb(cr,nsb,&(top->blocks[ebCGS]),nns,fr->nlr,
       &(top->idef),opts->ngener);
     */
    if (fr->ns->dump_nl > 0)
    {
        dump_nblist(fp, fr, fr->ns->dump_nl);
    }
}

void do_force_lowlevel(t_forcerec           *fr,
                       const t_inputrec     *ir,
                       const t_idef         *idef,
                       const t_commrec      *cr,
                       const gmx_multisim_t *ms,
                       t_nrnb               *nrnb,
                       gmx_wallcycle_t       wcycle,
                       const t_mdatoms      *md,
                       rvec                  x[],
                       history_t            *hist,
                       rvec                 *forceForUseWithShiftForces,
                       gmx::ForceWithVirial *forceWithVirial,
                       gmx_enerdata_t       *enerd,
                       t_fcdata             *fcd,
                       matrix                box,
                       t_lambda             *fepvals,
                       real                 *lambda,
                       const t_graph        *graph,
                       const t_blocka       *excl,
                       int                   flags)
{
    int         i, j;
    int         donb_flags;
    t_pbc       pbc;
    real        dvdl_dum[efptNR], dvdl_nb[efptNR];

    set_pbc(&pbc, fr->ePBC, box);

    /* reset free energy components */
    for (i = 0; i < efptNR; i++)
    {
        dvdl_nb[i]  = 0;
        dvdl_dum[i] = 0;
    }

    /* We only do non-bonded calculation with group scheme here, the verlet
     * calls are done from do_force_cutsVERLET(). */
    if (fr->cutoff_scheme == ecutsGROUP && (flags & GMX_FORCE_NONBONDED))
    {
        
        donb_flags = 0;
        /* Add short-range interactions */
        donb_flags |= GMX_NONBONDED_DO_SR;

        /* Currently all group scheme kernels always calculate (shift-)forces */
        if (flags & GMX_FORCE_FORCES)
        {
            donb_flags |= GMX_NONBONDED_DO_FORCE;
        }
        if (flags & GMX_FORCE_VIRIAL)
        {
            donb_flags |= GMX_NONBONDED_DO_SHIFTFORCE;
        }
        if (flags & GMX_FORCE_ENERGY)
        {
            donb_flags |= GMX_NONBONDED_DO_POTENTIAL;
        }
	
        wallcycle_sub_start(wcycle, ewcsNONBONDED);
        do_nonbonded(fr, x, forceForUseWithShiftForces, md, excl,
                     &enerd->grpp, nrnb,
                     lambda, dvdl_nb, -1, -1, donb_flags);

        /* If we do foreign lambda and we have soft-core interactions
         * we have to recalculate the (non-linear) energies contributions.
         */
        if (fepvals->n_lambda > 0 && (flags & GMX_FORCE_DHDL) && fepvals->sc_alpha != 0)
        {
            for (i = 0; i < enerd->n_lambda; i++)
            {
                real lam_i[efptNR];

                for (j = 0; j < efptNR; j++)
                {
                    lam_i[j] = (i == 0 ? lambda[j] : fepvals->all_lambda[j][i-1]);
                }
                reset_foreign_enerdata(enerd);
                do_nonbonded(fr, x, forceForUseWithShiftForces, md, excl,
                             &(enerd->foreign_grpp), nrnb,
                             lam_i, dvdl_dum, -1, -1,
                             (donb_flags & ~GMX_NONBONDED_DO_FORCE) | GMX_NONBONDED_DO_FOREIGNLAMBDA);
                sum_epot(&(enerd->foreign_grpp), enerd->foreign_term);
                enerd->enerpart_lambda[i] += enerd->foreign_term[F_EPOT];
            }
        }
        wallcycle_sub_stop(wcycle, ewcsNONBONDED);
    }

    if (fepvals->sc_alpha != 0)
    {
        enerd->dvdl_nonlin[efptVDW] += dvdl_nb[efptVDW];
    }
    else
    {
        enerd->dvdl_lin[efptVDW] += dvdl_nb[efptVDW];
    }

    if (fepvals->sc_alpha != 0)

    /* even though coulomb part is linear, we already added it, beacuse we
       need to go through the vdw calculation anyway */
    {
        enerd->dvdl_nonlin[efptCOUL] += dvdl_nb[efptCOUL];
    }
    else
    {
        enerd->dvdl_lin[efptCOUL] += dvdl_nb[efptCOUL];
    }

    if (debug)
    {
        pr_rvecs(debug, 0, "fshift after SR", fr->fshift, SHIFTS);
    }

    /* Shift the coordinates. Must be done before listed forces and PPPM,
     * but is also necessary for SHAKE and update, therefore it can NOT
     * go when no listed forces have to be evaluated.
     *
     * The shifting and PBC code is deliberately not timed, since with
     * the Verlet scheme it only takes non-zero time with triclinic
     * boxes, and even then the time is around a factor of 100 less
     * than the next smallest counter.
     */


    /* Here sometimes we would not need to shift with NBFonly,
     * but we do so anyhow for consistency of the returned coordinates.
     */
    if (graph)
    {
        shift_self(graph, box, x);
        if (TRICLINIC(box))
        {
            inc_nrnb(nrnb, eNR_SHIFTX, 2*graph->nnodes);
        }
        else
        {
            inc_nrnb(nrnb, eNR_SHIFTX, graph->nnodes);
        }
    }
    /* Check whether we need to do listed interactions or correct for exclusions */
    if (fr->bMolPBC &&
        ((flags & GMX_FORCE_LISTED)
         || EEL_RF(fr->ic->eeltype) || EEL_FULL(fr->ic->eeltype) || EVDW_PME(fr->ic->vdwtype)))
    {
        /* TODO There are no electrostatics methods that require this
           transformation, when using the Verlet scheme, so update the
           above conditional. */
        /* Since all atoms are in the rectangular or triclinic unit-cell,
         * only single box vector shifts (2 in x) are required.
         */
        set_pbc_dd(&pbc, fr->ePBC, nullptr, TRUE, box);
    }

    do_force_listed(wcycle, box, ir->fepvals, cr, ms,
                    idef, x, hist,
                    forceForUseWithShiftForces, forceWithVirial,
                    fr, &pbc, graph, enerd, nrnb, lambda, md, fcd,
                    nullptr, flags);


    /* Do long-range electrostatics and/or LJ-PME, including related short-range
     * corrections.
     */
    {
        /* Is there a reaction-field exclusion correction needed?
         * With the Verlet scheme, exclusion forces are calculated
         * in the non-bonded kernel.
         */
        if (ir->cutoff_scheme != ecutsVERLET && EEL_RF(fr->ic->eeltype))
        {
            real dvdl_rf_excl      = 0;
            enerd->term[F_RF_EXCL] =
                RF_excl_correction(fr, graph, md, excl, false,
                                   x, forceForUseWithShiftForces,
                                   fr->fshift, &pbc, lambda[efptCOUL], &dvdl_rf_excl);

            enerd->dvdl_lin[efptCOUL] += dvdl_rf_excl;
        }
    }

    if (debug)
    {
        print_nrnb(debug, nrnb);
    }

    if (debug)
    {
        pr_rvecs(debug, 0, "fshift after bondeds", fr->fshift, SHIFTS);
    }

}

void init_enerdata(int ngener, int n_lambda, gmx_enerdata_t *enerd)
{
    int i, n2;

    for (i = 0; i < F_NRE; i++)
    {
        enerd->term[i]         = 0;
        enerd->foreign_term[i] = 0;
    }


    for (i = 0; i < efptNR; i++)
    {
        enerd->dvdl_lin[i]     = 0;
        enerd->dvdl_nonlin[i]  = 0;
    }

    n2 = ngener*ngener;
    if (debug)
    {
        fprintf(debug, "Creating %d sized group matrix for energies\n", n2);
    }
    enerd->grpp.nener         = n2;
    enerd->foreign_grpp.nener = n2;
    for (i = 0; (i < egNR); i++)
    {
        snew(enerd->grpp.ener[i], n2);
        snew(enerd->foreign_grpp.ener[i], n2);
    }

    if (n_lambda)
    {
        enerd->n_lambda = 1 + n_lambda;
        snew(enerd->enerpart_lambda, enerd->n_lambda);
    }
    else
    {
        enerd->n_lambda = 0;
    }
}

void destroy_enerdata(gmx_enerdata_t *enerd)
{
    int i;

    for (i = 0; (i < egNR); i++)
    {
        sfree(enerd->grpp.ener[i]);
    }

    for (i = 0; (i < egNR); i++)
    {
        sfree(enerd->foreign_grpp.ener[i]);
    }

    if (enerd->n_lambda)
    {
        sfree(enerd->enerpart_lambda);
    }
}

static real sum_v(int n, const real v[])
{
    real t;
    int  i;

    t = 0.0;
    for (i = 0; (i < n); i++)
    {
        t = t + v[i];
    }

    return t;
}

void sum_epot(gmx_grppairener_t *grpp, real *epot)
{
    int i;

    /* Accumulate energies */
    epot[F_COUL_SR]  = sum_v(grpp->nener, grpp->ener[egCOULSR]);
    epot[F_LJ]       = sum_v(grpp->nener, grpp->ener[egLJSR]);
    epot[F_LJ14]     = sum_v(grpp->nener, grpp->ener[egLJ14]);
    epot[F_COUL14]   = sum_v(grpp->nener, grpp->ener[egCOUL14]);

/* lattice part of LR doesnt belong to any group
 * and has been added earlier
 */
    epot[F_BHAM]     = sum_v(grpp->nener, grpp->ener[egBHAMSR]);

    epot[F_EPOT] = 0;
    for (i = 0; (i < F_EPOT); i++)
    {
        if (i != F_DISRESVIOL && i != F_ORIRESDEV)
        {
            epot[F_EPOT] += epot[i];
        }
    }
}

void sum_dhdl(gmx_enerdata_t *enerd, gmx::ArrayRef<const real> lambda, t_lambda *fepvals)
{
    int    index;

    enerd->dvdl_lin[efptVDW] += enerd->term[F_DVDL_VDW];  /* include dispersion correction */
    enerd->term[F_DVDL]       = 0.0;
    for (int i = 0; i < efptNR; i++)
    {
        if (fepvals->separate_dvdl[i])
        {
            /* could this be done more readably/compactly? */
            switch (i)
            {
                case (efptMASS):
                    index = F_DKDL;
                    break;
                case (efptCOUL):
                    index = F_DVDL_COUL;
                    break;
                case (efptVDW):
                    index = F_DVDL_VDW;
                    break;
                case (efptBONDED):
                    index = F_DVDL_BONDED;
                    break;
                case (efptRESTRAINT):
                    index = F_DVDL_RESTRAINT;
                    break;
                default:
                    index = F_DVDL;
                    break;
            }
            enerd->term[index] = enerd->dvdl_lin[i] + enerd->dvdl_nonlin[i];
            if (debug)
            {
                fprintf(debug, "dvdl-%s[%2d]: %f: non-linear %f + linear %f\n",
                        efpt_names[i], i, enerd->term[index], enerd->dvdl_nonlin[i], enerd->dvdl_lin[i]);
            }
        }
        else
        {
            enerd->term[F_DVDL] += enerd->dvdl_lin[i] + enerd->dvdl_nonlin[i];
            if (debug)
            {
                fprintf(debug, "dvd-%sl[%2d]: %f: non-linear %f + linear %f\n",
                        efpt_names[0], i, enerd->term[F_DVDL], enerd->dvdl_nonlin[i], enerd->dvdl_lin[i]);
            }
        }
    }

    if (fepvals->separate_dvdl[efptBONDED])
    {
        enerd->term[F_DVDL_BONDED] += enerd->term[F_DVDL_CONSTR];
    }
    else
    {
        enerd->term[F_DVDL] += enerd->term[F_DVDL_CONSTR];
    }

    for (int i = 0; i < fepvals->n_lambda; i++)
    {
        /* note we are iterating over fepvals here!
           For the current lam, dlam = 0 automatically,
           so we don't need to add anything to the
           enerd->enerpart_lambda[0] */

        /* we don't need to worry about dvdl_lin contributions to dE at
           current lambda, because the contributions to the current
           lambda are automatically zeroed */

        double &enerpart_lambda = enerd->enerpart_lambda[i + 1];

        for (gmx::index j = 0; j < lambda.size(); j++)
        {
            /* Note that this loop is over all dhdl components, not just the separated ones */
            const double dlam  = fepvals->all_lambda[j][i] - lambda[j];

            enerpart_lambda   += dlam*enerd->dvdl_lin[j];

            /* Constraints can not be evaluated at foreign lambdas, so we add
             * a linear extrapolation. This is an approximation, but usually
             * quite accurate since constraints change little between lambdas.
             */
            if ((j == efptBONDED && fepvals->separate_dvdl[efptBONDED]) ||
                (j == efptFEP && !fepvals->separate_dvdl[efptBONDED]))
            {
                enerpart_lambda += dlam*enerd->term[F_DVDL_CONSTR];
            }

            if (j == efptMASS)
            {
                enerpart_lambda += dlam*enerd->term[F_DKDL];
            }

            if (debug)
            {
                fprintf(debug, "enerdiff lam %g: (%15s), non-linear %f linear %f*%f\n",
                        fepvals->all_lambda[j][i], efpt_names[j],
                        enerpart_lambda - enerd->enerpart_lambda[0],
                        dlam, enerd->dvdl_lin[j]);
            }
        }
    }

    /* The constrain contribution is now included in other terms, so clear it */
    enerd->term[F_DVDL_CONSTR] = 0;
}


void reset_foreign_enerdata(gmx_enerdata_t *enerd)
{
    int  i, j;

    /* First reset all foreign energy components.  Foreign energies always called on
       neighbor search steps */
    for (i = 0; (i < egNR); i++)
    {
        for (j = 0; (j < enerd->grpp.nener); j++)
        {
            enerd->foreign_grpp.ener[i][j] = 0.0;
        }
    }

    /* potential energy components */
    for (i = 0; (i <= F_EPOT); i++)
    {
        enerd->foreign_term[i] = 0.0;
    }
}

void reset_enerdata(gmx_enerdata_t *enerd)
{
    int      i, j;

    /* First reset all energy components. */
    for (i = 0; (i < egNR); i++)
    {
        for (j = 0; (j < enerd->grpp.nener); j++)
        {
            enerd->grpp.ener[i][j] = 0.0;
        }
    }
    for (i = 0; i < efptNR; i++)
    {
        enerd->dvdl_lin[i]    = 0.0;
        enerd->dvdl_nonlin[i] = 0.0;
    }

    /* Normal potential energy components */
    for (i = 0; (i <= F_EPOT); i++)
    {
        enerd->term[i] = 0.0;
    }
    enerd->term[F_DVDL]            = 0.0;
    enerd->term[F_DVDL_COUL]       = 0.0;
    enerd->term[F_DVDL_VDW]        = 0.0;
    enerd->term[F_DVDL_BONDED]     = 0.0;
    enerd->term[F_DVDL_RESTRAINT]  = 0.0;
    enerd->term[F_DKDL]            = 0.0;
    if (enerd->n_lambda > 0)
    {
        for (i = 0; i < enerd->n_lambda; i++)
        {
            enerd->enerpart_lambda[i] = 0.0;
        }
    }
    /* reset foreign energy data - separate function since we also call it elsewhere */
    reset_foreign_enerdata(enerd);
}

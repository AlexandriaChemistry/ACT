/*
 * This file is part of the GROMACS molecular simulation package.
 *
 * Copyright (c) 1991-2000, University of Groningen, The Netherlands.
 * Copyright (c) 2001-2010, The GROMACS development team.
 * Copyright (c) 2012,2014,2015,2016,2017,2018, by the GROMACS development team, led by
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

#include "inputrec.h"

#include <cstdio>
#include <cstdlib>
#include <cstring>

#include <algorithm>

#include "gromacs/math/veccompare.h"
#include "gromacs/math/vecdump.h"
#include "gromacs/mdtypes/md_enums.h"
#include "gromacs/pbcutil/pbc.h"
#include "gromacs/utility/compare.h"
#include "gromacs/utility/cstringutil.h"
#include "gromacs/utility/fatalerror.h"
#include "gromacs/utility/keyvaluetree.h"
#include "gromacs/utility/smalloc.h"
#include "gromacs/utility/snprintf.h"
#include "gromacs/utility/strconvert.h"
#include "gromacs/utility/stringutil.h"
#include "gromacs/utility/textwriter.h"
#include "gromacs/utility/txtdump.h"

//! Macro to select a bool name
#define EBOOL(e)       gmx::boolToString(e)

/* The minimum number of integration steps required for reasonably accurate
 * integration of first and second order coupling algorithms.
 */
const int nstmin_berendsen_tcouple =  5;
const int nstmin_berendsen_pcouple = 10;
const int nstmin_harmonic          = 20;

t_inputrec::t_inputrec()
{
    std::memset(this, 0, sizeof(*this)); // NOLINT(bugprone-undefined-memory-manipulation)
    snew(fepvals, 1);
    snew(expandedvals, 1);
    snew(simtempvals, 1);
}

t_inputrec::~t_inputrec()
{
    done_inputrec(this);
}

static int nst_wanted(const t_inputrec *ir)
{
    if (ir->nstlist > 0)
    {
        return ir->nstlist;
    }
    else
    {
        return 10;
    }
}

int ir_optimal_nstcalcenergy(const t_inputrec *ir)
{
    return nst_wanted(ir);
}

int tcouple_min_integration_steps(int etc)
{
    int n;

    switch (etc)
    {
        case etcNO:
            n = 0;
            break;
        case etcBERENDSEN:
        case etcYES:
            n = nstmin_berendsen_tcouple;
            break;
        case etcVRESCALE:
            /* V-rescale supports instantaneous rescaling */
            n = 0;
            break;
        case etcNOSEHOOVER:
            n = nstmin_harmonic;
            break;
        case etcANDERSEN:
        case etcANDERSENMASSIVE:
            n = 1;
            break;
        default:
            gmx_incons("Unknown etc value");
    }

    return n;
}

int ir_optimal_nsttcouple(const t_inputrec *ir)
{
    int  nmin, nwanted, n;
    real tau_min;
    int  g;

    nmin = tcouple_min_integration_steps(ir->etc);

    nwanted = nst_wanted(ir);

    tau_min = 1e20;
    if (ir->etc != etcNO)
    {
        for (g = 0; g < ir->opts.ngtc; g++)
        {
            if (ir->opts.tau_t[g] > 0)
            {
                tau_min = std::min(tau_min, ir->opts.tau_t[g]);
            }
        }
    }

    if (nmin == 0 || ir->delta_t*nwanted <= tau_min)
    {
        n = nwanted;
    }
    else
    {
        n = static_cast<int>(tau_min/(ir->delta_t*nmin) + 0.001);
        if (n < 1)
        {
            n = 1;
        }
        while (nwanted % n != 0)
        {
            n--;
        }
    }

    return n;
}

int pcouple_min_integration_steps(int epc)
{
    int n;

    switch (epc)
    {
        case epcNO:
            n = 0;
            break;
        case etcBERENDSEN:
        case epcISOTROPIC:
            n = nstmin_berendsen_pcouple;
            break;
        case epcPARRINELLORAHMAN:
        case epcMTTK:
            n = nstmin_harmonic;
            break;
        default:
            gmx_incons("Unknown epc value");
    }

    return n;
}

int ir_optimal_nstpcouple(const t_inputrec *ir)
{
    int  nmin, nwanted, n;

    nmin = pcouple_min_integration_steps(ir->epc);

    nwanted = nst_wanted(ir);

    if (nmin == 0 || ir->delta_t*nwanted <= ir->tau_p)
    {
        n = nwanted;
    }
    else
    {
        n = static_cast<int>(ir->tau_p/(ir->delta_t*nmin) + 0.001);
        if (n < 1)
        {
            n = 1;
        }
        while (nwanted % n != 0)
        {
            n--;
        }
    }

    return n;
}

gmx_bool ir_coulomb_switched(const t_inputrec *ir)
{
    return (ir->coulombtype == eelSWITCH ||
            ir->coulombtype == eelSHIFT ||
            ir->coulombtype == eelENCADSHIFT ||
            ir->coulombtype == eelPMESWITCH ||
            ir->coulombtype == eelPMEUSERSWITCH ||
            ir->coulomb_modifier == eintmodPOTSWITCH ||
            ir->coulomb_modifier == eintmodFORCESWITCH);
}

gmx_bool ir_coulomb_is_zero_at_cutoff(const t_inputrec *ir)
{
    return (ir->cutoff_scheme == ecutsVERLET ||
            ir_coulomb_switched(ir) || ir->coulomb_modifier != eintmodNONE ||
            ir->coulombtype == eelRF_ZERO);
}

gmx_bool ir_coulomb_might_be_zero_at_cutoff(const t_inputrec *ir)
{
    return (ir_coulomb_is_zero_at_cutoff(ir) || ir->coulombtype == eelUSER || ir->coulombtype == eelPMEUSER);
}

gmx_bool ir_vdw_switched(const t_inputrec *ir)
{
    return (ir->vdwtype == evdwSWITCH ||
            ir->vdwtype == evdwSHIFT ||
            ir->vdwtype == evdwENCADSHIFT ||
            ir->vdw_modifier == eintmodPOTSWITCH ||
            ir->vdw_modifier == eintmodFORCESWITCH);
}

gmx_bool ir_vdw_is_zero_at_cutoff(const t_inputrec *ir)
{
    return (ir->cutoff_scheme == ecutsVERLET ||
            ir_vdw_switched(ir) || ir->vdw_modifier != eintmodNONE);
}

gmx_bool ir_vdw_might_be_zero_at_cutoff(const t_inputrec *ir)
{
    return (ir_vdw_is_zero_at_cutoff(ir) || ir->vdwtype == evdwUSER);
}

static void done_lambdas(t_lambda *fep)
{
    if (fep->n_lambda > 0)
    {
        for (int i = 0; i < efptNR; i++)
        {
            sfree(fep->all_lambda[i]);
        }
    }
    sfree(fep->all_lambda);
}

void done_inputrec(t_inputrec *ir)
{
    sfree(ir->opts.nrdf);
    sfree(ir->opts.ref_t);
    sfree(ir->opts.annealing);
    sfree(ir->opts.anneal_npoints);
    sfree(ir->opts.anneal_time);
    sfree(ir->opts.anneal_temp);
    sfree(ir->opts.tau_t);
    sfree(ir->opts.acc);
    sfree(ir->opts.nFreeze);
    sfree(ir->opts.egp_flags);
    done_lambdas(ir->fepvals);
    sfree(ir->fepvals);
    sfree(ir->expandedvals);
    sfree(ir->simtempvals);

    delete ir->params;
}

static void pr_grp_opts(FILE *out, int indent, const char *title, const t_grpopts *opts,
                        gmx_bool bMDPformat)
{
    int i, m, j;

    if (!bMDPformat)
    {
        fprintf(out, "%s:\n", title);
    }

    pr_indent(out, indent);
    fprintf(out, "nrdf%s", bMDPformat ? " = " : ":");
    for (i = 0; (i < opts->ngtc); i++)
    {
        fprintf(out, "  %10g", opts->nrdf[i]);
    }
    fprintf(out, "\n");

    pr_indent(out, indent);
    fprintf(out, "ref-t%s", bMDPformat ? " = " : ":");
    for (i = 0; (i < opts->ngtc); i++)
    {
        fprintf(out, "  %10g", opts->ref_t[i]);
    }
    fprintf(out, "\n");

    pr_indent(out, indent);
    fprintf(out, "tau-t%s", bMDPformat ? " = " : ":");
    for (i = 0; (i < opts->ngtc); i++)
    {
        fprintf(out, "  %10g", opts->tau_t[i]);
    }
    fprintf(out, "\n");

    /* Pretty-print the simulated annealing info */
    fprintf(out, "annealing%s", bMDPformat ? " = " : ":");
    for (i = 0; (i < opts->ngtc); i++)
    {
        fprintf(out, "  %10s", EANNEAL(opts->annealing[i]));
    }
    fprintf(out, "\n");

    fprintf(out, "annealing-npoints%s", bMDPformat ? " = " : ":");
    for (i = 0; (i < opts->ngtc); i++)
    {
        fprintf(out, "  %10d", opts->anneal_npoints[i]);
    }
    fprintf(out, "\n");

    for (i = 0; (i < opts->ngtc); i++)
    {
        if (opts->anneal_npoints[i] > 0)
        {
            fprintf(out, "annealing-time [%d]:\t", i);
            for (j = 0; (j < opts->anneal_npoints[i]); j++)
            {
                fprintf(out, "  %10.1f", opts->anneal_time[i][j]);
            }
            fprintf(out, "\n");
            fprintf(out, "annealing-temp [%d]:\t", i);
            for (j = 0; (j < opts->anneal_npoints[i]); j++)
            {
                fprintf(out, "  %10.1f", opts->anneal_temp[i][j]);
            }
            fprintf(out, "\n");
        }
    }

    pr_indent(out, indent);
    fprintf(out, "acc:\t");
    for (i = 0; (i < opts->ngacc); i++)
    {
        for (m = 0; (m < DIM); m++)
        {
            fprintf(out, "  %10g", opts->acc[i][m]);
        }
    }
    fprintf(out, "\n");

    pr_indent(out, indent);
    fprintf(out, "nfreeze:");
    for (i = 0; (i < opts->ngfrz); i++)
    {
        for (m = 0; (m < DIM); m++)
        {
            fprintf(out, "  %10s", opts->nFreeze[i][m] ? "Y" : "N");
        }
    }
    fprintf(out, "\n");


    for (i = 0; (i < opts->ngener); i++)
    {
        pr_indent(out, indent);
        fprintf(out, "energygrp-flags[%3d]:", i);
        for (m = 0; (m < opts->ngener); m++)
        {
            fprintf(out, " %d", opts->egp_flags[opts->ngener*i+m]);
        }
        fprintf(out, "\n");
    }

    fflush(out);
}

static void pr_matrix(FILE *fp, int indent, const char *title, const rvec *m,
                      gmx_bool bMDPformat)
{
    if (bMDPformat)
    {
        fprintf(fp, "%-10s    = %g %g %g %g %g %g\n", title,
                m[XX][XX], m[YY][YY], m[ZZ][ZZ], m[XX][YY], m[XX][ZZ], m[YY][ZZ]);
    }
    else
    {
        pr_rvecs(fp, indent, title, m, DIM);
    }
}

#define PS(t, s) pr_str(fp, indent, t, s)
#define PI(t, s) pr_int(fp, indent, t, s)
#define PSTEP(t, s) pr_int64(fp, indent, t, s)
#define PR(t, s) pr_real(fp, indent, t, s)
#define PD(t, s) pr_double(fp, indent, t, s)

static void pr_simtempvals(FILE *fp, int indent, const t_simtemp *simtemp, int n_lambda)
{
    PS("simulated-tempering-scaling", ESIMTEMP(simtemp->eSimTempScale));
    PR("sim-temp-low", simtemp->simtemp_low);
    PR("sim-temp-high", simtemp->simtemp_high);
    pr_rvec(fp, indent, "simulated tempering temperatures", simtemp->temperatures, n_lambda, TRUE);
}

static void pr_expandedvals(FILE *fp, int indent, const t_expanded *expand, int n_lambda)
{

    PI("nstexpanded", expand->nstexpanded);
    PS("lmc-stats", elamstats_names[expand->elamstats]);
    PS("lmc-move", elmcmove_names[expand->elmcmove]);
    PS("lmc-weights-equil", elmceq_names[expand->elmceq]);
    if (expand->elmceq == elmceqNUMATLAM)
    {
        PI("weight-equil-number-all-lambda", expand->equil_n_at_lam);
    }
    if (expand->elmceq == elmceqSAMPLES)
    {
        PI("weight-equil-number-samples", expand->equil_samples);
    }
    if (expand->elmceq == elmceqSTEPS)
    {
        PI("weight-equil-number-steps", expand->equil_steps);
    }
    if (expand->elmceq == elmceqWLDELTA)
    {
        PR("weight-equil-wl-delta", expand->equil_wl_delta);
    }
    if (expand->elmceq == elmceqRATIO)
    {
        PR("weight-equil-count-ratio", expand->equil_ratio);
    }
    PI("lmc-seed", expand->lmc_seed);
    PR("mc-temperature", expand->mc_temp);
    PI("lmc-repeats", expand->lmc_repeats);
    PI("lmc-gibbsdelta", expand->gibbsdeltalam);
    PI("lmc-forced-nstart", expand->lmc_forced_nstart);
    PS("symmetrized-transition-matrix", EBOOL(expand->bSymmetrizedTMatrix));
    PI("nst-transition-matrix", expand->nstTij);
    PI("mininum-var-min", expand->minvarmin); /*default is reasonable */
    PI("weight-c-range", expand->c_range);    /* default is just C=0 */
    PR("wl-scale", expand->wl_scale);
    PR("wl-ratio", expand->wl_ratio);
    PR("init-wl-delta", expand->init_wl_delta);
    PS("wl-oneovert", EBOOL(expand->bWLoneovert));

    pr_indent(fp, indent);
    pr_rvec(fp, indent, "init-lambda-weights", expand->init_lambda_weights, n_lambda, TRUE);
    PS("init-weights", EBOOL(expand->bInit_weights));
}

static void pr_fepvals(FILE *fp, int indent, const t_lambda *fep, gmx_bool bMDPformat)
{
    int i, j;

    PR("init-lambda", fep->init_lambda);
    PI("init-lambda-state", fep->init_fep_state);
    PR("delta-lambda", fep->delta_lambda);
    PI("nstdhdl", fep->nstdhdl);

    if (!bMDPformat)
    {
        PI("n-lambdas", fep->n_lambda);
    }
    if (fep->n_lambda > 0)
    {
        pr_indent(fp, indent);
        fprintf(fp, "separate-dvdl%s\n", bMDPformat ? " = " : ":");
        for (i = 0; i < efptNR; i++)
        {
            fprintf(fp, "%18s = ", efpt_names[i]);
            if (fep->separate_dvdl[i])
            {
                fprintf(fp, "  TRUE");
            }
            else
            {
                fprintf(fp, "  FALSE");
            }
            fprintf(fp, "\n");
        }
        fprintf(fp, "all-lambdas%s\n", bMDPformat ? " = " : ":");
        for (i = 0; i < efptNR; i++)
        {
            fprintf(fp, "%18s = ", efpt_names[i]);
            for (j = 0; j < fep->n_lambda; j++)
            {
                fprintf(fp, "  %10g", fep->all_lambda[i][j]);
            }
            fprintf(fp, "\n");
        }
    }
    PI("calc-lambda-neighbors", fep->lambda_neighbors);
    PS("dhdl-print-energy", edHdLPrintEnergy_names[fep->edHdLPrintEnergy]);
    PR("sc-alpha", fep->sc_alpha);
    PI("sc-power", fep->sc_power);
    PR("sc-r-power", fep->sc_r_power);
    PR("sc-sigma", fep->sc_sigma);
    PR("sc-sigma-min", fep->sc_sigma_min);
    PS("sc-coul", EBOOL(fep->bScCoul));
    PI("dh-hist-size", fep->dh_hist_size);
    PD("dh-hist-spacing", fep->dh_hist_spacing);
    PS("separate-dhdl-file", SEPDHDLFILETYPE(fep->separate_dhdl_file));
    PS("dhdl-derivatives", DHDLDERIVATIVESTYPE(fep->dhdl_derivatives));
};

void pr_inputrec(FILE *fp, int indent, const char *title, const t_inputrec *ir,
                 gmx_bool bMDPformat)
{
    const char *infbuf = "inf";

    if (available(fp, ir, indent, title))
    {
        if (!bMDPformat)
        {
            indent = pr_title(fp, indent, title);
        }
        /* Try to make this list appear in the same order as the
         * options are written in the default mdout.mdp, and with
         * the same user-exposed names to facilitate debugging.
         */
        PS("integrator", EI(ir->eI));
        PR("tinit", ir->init_t);
        PR("dt", ir->delta_t);
        PSTEP("nsteps", ir->nsteps);
        PSTEP("init-step", ir->init_step);
        PI("simulation-part", ir->simulation_part);
        PS("comm-mode", ECOM(ir->comm_mode));
        PI("nstcomm", ir->nstcomm);

        /* Langevin dynamics */
        PR("bd-fric", ir->bd_fric);
        PSTEP("ld-seed", ir->ld_seed);

        /* Energy minimization */
        PR("emtol", ir->em_tol);
        PR("emstep", ir->em_stepsize);
        PI("niter", ir->niter);
        PR("fcstep", ir->fc_stepsize);
        PI("nstcgsteep", ir->nstcgsteep);
        PI("nbfgscorr", ir->nbfgscorr);

        /* Test particle insertion */
        PR("rtpi", ir->rtpi);

        /* Output control */
        PI("nstxout", ir->nstxout);
        PI("nstvout", ir->nstvout);
        PI("nstfout", ir->nstfout);
        PI("nstlog", ir->nstlog);
        PI("nstcalcenergy", ir->nstcalcenergy);
        PI("nstenergy", ir->nstenergy);

        /* Neighborsearching parameters */
        PS("cutoff-scheme", ECUTSCHEME(ir->cutoff_scheme));
        PI("nstlist", ir->nstlist);
        PS("ns-type", ENS(ir->ns_type));
        PS("pbc", epbc_names[ir->ePBC]);
        PS("periodic-molecules", EBOOL(ir->bPeriodicMols));
        PR("verlet-buffer-tolerance", ir->verletbuf_tol);
        PR("rlist", ir->rlist);

        /* Options for electrostatics and VdW */
        PS("coulombtype", EELTYPE(ir->coulombtype));
        PS("coulomb-modifier", INTMODIFIER(ir->coulomb_modifier));
        PR("rcoulomb-switch", ir->rcoulomb_switch);
        PR("rcoulomb", ir->rcoulomb);
        if (ir->epsilon_r != 0)
        {
            PR("epsilon-r", ir->epsilon_r);
        }
        else
        {
            PS("epsilon-r", infbuf);
        }
        if (ir->epsilon_rf != 0)
        {
            PR("epsilon-rf", ir->epsilon_rf);
        }
        else
        {
            PS("epsilon-rf", infbuf);
        }
        PS("vdw-type", EVDWTYPE(ir->vdwtype));
        PS("vdw-modifier", INTMODIFIER(ir->vdw_modifier));
        PR("rvdw-switch", ir->rvdw_switch);
        PR("rvdw", ir->rvdw);
        PS("DispCorr", EDISPCORR(ir->eDispCorr));
        PR("table-extension", ir->tabext);

        PR("fourierspacing", ir->fourier_spacing);
        PI("fourier-nx", ir->nkx);
        PI("fourier-ny", ir->nky);
        PI("fourier-nz", ir->nkz);
        PI("pme-order", ir->pme_order);
        PR("ewald-rtol", ir->ewald_rtol);
        PR("ewald-rtol-lj", ir->ewald_rtol_lj);
        PS("lj-pme-comb-rule", ELJPMECOMBNAMES(ir->ljpme_combination_rule));
        PR("ewald-geometry", ir->ewald_geometry);
        PR("epsilon-surface", ir->epsilon_surface);

        /* Options for weak coupling algorithms */
        PS("tcoupl", ETCOUPLTYPE(ir->etc));
        PI("nsttcouple", ir->nsttcouple);
        PI("nh-chain-length", ir->opts.nhchainlength);
        PS("print-nose-hoover-chain-variables", EBOOL(ir->bPrintNHChains));

        PS("pcoupl", EPCOUPLTYPE(ir->epc));
        PS("pcoupltype", EPCOUPLTYPETYPE(ir->epct));
        PI("nstpcouple", ir->nstpcouple);
        PR("tau-p", ir->tau_p);
        pr_matrix(fp, indent, "compressibility", ir->compress, bMDPformat);
        pr_matrix(fp, indent, "ref-p", ir->ref_p, bMDPformat);
        PS("refcoord-scaling", EREFSCALINGTYPE(ir->refcoord_scaling));

        if (bMDPformat)
        {
            fprintf(fp, "posres-com  = %g %g %g\n", ir->posres_com[XX],
                    ir->posres_com[YY], ir->posres_com[ZZ]);
            fprintf(fp, "posres-comB = %g %g %g\n", ir->posres_comB[XX],
                    ir->posres_comB[YY], ir->posres_comB[ZZ]);
        }
        else
        {
            pr_rvec(fp, indent, "posres-com", ir->posres_com, DIM, TRUE);
            pr_rvec(fp, indent, "posres-comB", ir->posres_comB, DIM, TRUE);
        }

        /* CONSTRAINT OPTIONS */
        PS("constraint-algorithm", ECONSTRTYPE(ir->eConstrAlg));
        PS("continuation", EBOOL(ir->bContinuation));

        PS("Shake-SOR", EBOOL(ir->bShakeSOR));
        PR("shake-tol", ir->shake_tol);
        PI("lincs-order", ir->nProjOrder);
        PI("lincs-iter", ir->nLincsIter);
        PR("lincs-warnangle", ir->LincsWarnAngle);

        /* NMR refinement stuff */
        PS("disre", EDISRETYPE(ir->eDisre));
        PS("disre-weighting", EDISREWEIGHTING(ir->eDisreWeighting));
        PS("disre-mixed", EBOOL(ir->bDisreMixed));
        PR("dr-fc", ir->dr_fc);
        PR("dr-tau", ir->dr_tau);
        PR("nstdisreout", ir->nstdisreout);

        PR("orire-fc", ir->orires_fc);
        PR("orire-tau", ir->orires_tau);
        PR("nstorireout", ir->nstorireout);

        /* FREE ENERGY VARIABLES */
        PS("free-energy", EFEPTYPE(ir->efep));
        if (ir->efep != efepNO || ir->bSimTemp)
        {
            pr_fepvals(fp, indent, ir->fepvals, bMDPformat);
        }
        if (ir->bExpanded)
        {
            pr_expandedvals(fp, indent, ir->expandedvals, ir->fepvals->n_lambda);
        }

        /* NON-equilibrium MD stuff */
        PR("cos-acceleration", ir->cos_accel);
        pr_matrix(fp, indent, "deform", ir->deform, bMDPformat);

        /* SIMULATED TEMPERING */
        PS("simulated-tempering", EBOOL(ir->bSimTemp));
        if (ir->bSimTemp)
        {
            pr_simtempvals(fp, indent, ir->simtempvals, ir->fepvals->n_lambda);
        }

        /* USER-DEFINED THINGIES */
        PI("userint1", ir->userint1);
        PI("userint2", ir->userint2);
        PI("userint3", ir->userint3);
        PI("userint4", ir->userint4);
        PR("userreal1", ir->userreal1);
        PR("userreal2", ir->userreal2);
        PR("userreal3", ir->userreal3);
        PR("userreal4", ir->userreal4);

        if (!bMDPformat)
        {
            gmx::TextWriter writer(fp);
            writer.wrapperSettings().setIndent(indent);
            gmx::dumpKeyValueTree(&writer, *ir->params);
        }

        pr_grp_opts(fp, indent, "grpopts", &(ir->opts), bMDPformat);
    }
}
#undef PS
#undef PR
#undef PI

static void cmp_grpopts(FILE *fp, const t_grpopts *opt1, const t_grpopts *opt2, real ftol, real abstol)
{
    int  i, j;

    cmp_int(fp, "inputrec->grpopts.ngtc", -1,  opt1->ngtc, opt2->ngtc);
    cmp_int(fp, "inputrec->grpopts.ngacc", -1, opt1->ngacc, opt2->ngacc);
    cmp_int(fp, "inputrec->grpopts.ngfrz", -1, opt1->ngfrz, opt2->ngfrz);
    cmp_int(fp, "inputrec->grpopts.ngener", -1, opt1->ngener, opt2->ngener);
    for (i = 0; (i < std::min(opt1->ngtc, opt2->ngtc)); i++)
    {
        cmp_real(fp, "inputrec->grpopts.nrdf", i, opt1->nrdf[i], opt2->nrdf[i], ftol, abstol);
        cmp_real(fp, "inputrec->grpopts.ref_t", i, opt1->ref_t[i], opt2->ref_t[i], ftol, abstol);
        cmp_real(fp, "inputrec->grpopts.tau_t", i, opt1->tau_t[i], opt2->tau_t[i], ftol, abstol);
        cmp_int(fp, "inputrec->grpopts.annealing", i, opt1->annealing[i], opt2->annealing[i]);
        cmp_int(fp, "inputrec->grpopts.anneal_npoints", i,
                opt1->anneal_npoints[i], opt2->anneal_npoints[i]);
        if (opt1->anneal_npoints[i] == opt2->anneal_npoints[i])
        {
            auto buf1 = gmx::formatString("inputrec->grpopts.anneal_time[%d]", i);
            auto buf2 = gmx::formatString("inputrec->grpopts.anneal_temp[%d]", i);
            for (j = 0; j < opt1->anneal_npoints[i]; j++)
            {
                cmp_real(fp, buf1.c_str(), j, opt1->anneal_time[i][j], opt2->anneal_time[i][j], ftol, abstol);
                cmp_real(fp, buf2.c_str(), j, opt1->anneal_temp[i][j], opt2->anneal_temp[i][j], ftol, abstol);
            }
        }
    }
    if (opt1->ngener == opt2->ngener)
    {
        for (i = 0; i < opt1->ngener; i++)
        {
            for (j = i; j < opt1->ngener; j++)
            {
                auto buf1 = gmx::formatString("inputrec->grpopts.egp_flags[%d]", i);
                cmp_int(fp, buf1.c_str(), j,
                        opt1->egp_flags[opt1->ngener*i+j],
                        opt2->egp_flags[opt1->ngener*i+j]);
            }
        }
    }
    for (i = 0; (i < std::min(opt1->ngacc, opt2->ngacc)); i++)
    {
        cmp_rvec(fp, "inputrec->grpopts.acc", i, opt1->acc[i], opt2->acc[i], ftol, abstol);
    }
    for (i = 0; (i < std::min(opt1->ngfrz, opt2->ngfrz)); i++)
    {
        cmp_ivec(fp, "inputrec->grpopts.nFreeze", i, opt1->nFreeze[i], opt2->nFreeze[i]);
    }
}

static void cmp_simtempvals(FILE *fp, const t_simtemp *simtemp1, const t_simtemp *simtemp2, int n_lambda, real ftol, real abstol)
{
    int i;
    cmp_int(fp, "inputrec->simtempvals->eSimTempScale", -1, simtemp1->eSimTempScale, simtemp2->eSimTempScale);
    cmp_real(fp, "inputrec->simtempvals->simtemp_high", -1, simtemp1->simtemp_high, simtemp2->simtemp_high, ftol, abstol);
    cmp_real(fp, "inputrec->simtempvals->simtemp_low", -1, simtemp1->simtemp_low, simtemp2->simtemp_low, ftol, abstol);
    for (i = 0; i < n_lambda; i++)
    {
        cmp_real(fp, "inputrec->simtempvals->temperatures", -1, simtemp1->temperatures[i], simtemp2->temperatures[i], ftol, abstol);
    }
}

static void cmp_expandedvals(FILE *fp, const t_expanded *expand1, const t_expanded *expand2, int n_lambda, real ftol, real abstol)
{
    int i;

    cmp_bool(fp, "inputrec->fepvals->bInit_weights", -1, expand1->bInit_weights, expand2->bInit_weights);
    cmp_bool(fp, "inputrec->fepvals->bWLoneovert", -1, expand1->bWLoneovert, expand2->bWLoneovert);

    for (i = 0; i < n_lambda; i++)
    {
        cmp_real(fp, "inputrec->expandedvals->init_lambda_weights", -1,
                 expand1->init_lambda_weights[i], expand2->init_lambda_weights[i], ftol, abstol);
    }

    cmp_int(fp, "inputrec->expandedvals->lambda-stats", -1, expand1->elamstats, expand2->elamstats);
    cmp_int(fp, "inputrec->expandedvals->lambda-mc-move", -1, expand1->elmcmove, expand2->elmcmove);
    cmp_int(fp, "inputrec->expandedvals->lmc-repeats", -1, expand1->lmc_repeats, expand2->lmc_repeats);
    cmp_int(fp, "inputrec->expandedvals->lmc-gibbsdelta", -1, expand1->gibbsdeltalam, expand2->gibbsdeltalam);
    cmp_int(fp, "inputrec->expandedvals->lmc-forced-nstart", -1, expand1->lmc_forced_nstart, expand2->lmc_forced_nstart);
    cmp_int(fp, "inputrec->expandedvals->lambda-weights-equil", -1, expand1->elmceq, expand2->elmceq);
    cmp_int(fp, "inputrec->expandedvals->,weight-equil-number-all-lambda", -1, expand1->equil_n_at_lam, expand2->equil_n_at_lam);
    cmp_int(fp, "inputrec->expandedvals->weight-equil-number-samples", -1, expand1->equil_samples, expand2->equil_samples);
    cmp_int(fp, "inputrec->expandedvals->weight-equil-number-steps", -1, expand1->equil_steps, expand2->equil_steps);
    cmp_real(fp, "inputrec->expandedvals->weight-equil-wl-delta", -1, expand1->equil_wl_delta, expand2->equil_wl_delta, ftol, abstol);
    cmp_real(fp, "inputrec->expandedvals->weight-equil-count-ratio", -1, expand1->equil_ratio, expand2->equil_ratio, ftol, abstol);
    cmp_bool(fp, "inputrec->expandedvals->symmetrized-transition-matrix", -1, expand1->bSymmetrizedTMatrix, expand2->bSymmetrizedTMatrix);
    cmp_int(fp, "inputrec->expandedvals->nstTij", -1, expand1->nstTij, expand2->nstTij);
    cmp_int(fp, "inputrec->expandedvals->mininum-var-min", -1, expand1->minvarmin, expand2->minvarmin); /*default is reasonable */
    cmp_int(fp, "inputrec->expandedvals->weight-c-range", -1, expand1->c_range, expand2->c_range);      /* default is just C=0 */
    cmp_real(fp, "inputrec->expandedvals->wl-scale", -1, expand1->wl_scale, expand2->wl_scale, ftol, abstol);
    cmp_real(fp, "inputrec->expandedvals->init-wl-delta", -1, expand1->init_wl_delta, expand2->init_wl_delta, ftol, abstol);
    cmp_real(fp, "inputrec->expandedvals->wl-ratio", -1, expand1->wl_ratio, expand2->wl_ratio, ftol, abstol);
    cmp_int(fp, "inputrec->expandedvals->nstexpanded", -1, expand1->nstexpanded, expand2->nstexpanded);
    cmp_int(fp, "inputrec->expandedvals->lmc-seed", -1, expand1->lmc_seed, expand2->lmc_seed);
    cmp_real(fp, "inputrec->expandedvals->mc-temperature", -1, expand1->mc_temp, expand2->mc_temp, ftol, abstol);
}

static void cmp_fepvals(FILE *fp, const t_lambda *fep1, const t_lambda *fep2, real ftol, real abstol)
{
    int i, j;
    cmp_int(fp, "inputrec->nstdhdl", -1, fep1->nstdhdl, fep2->nstdhdl);
    cmp_double(fp, "inputrec->fepvals->init_fep_state", -1, fep1->init_fep_state, fep2->init_fep_state, ftol, abstol);
    cmp_double(fp, "inputrec->fepvals->delta_lambda", -1, fep1->delta_lambda, fep2->delta_lambda, ftol, abstol);
    cmp_int(fp, "inputrec->fepvals->n_lambda", -1, fep1->n_lambda, fep2->n_lambda);
    for (i = 0; i < efptNR; i++)
    {
        for (j = 0; j < std::min(fep1->n_lambda, fep2->n_lambda); j++)
        {
            cmp_double(fp, "inputrec->fepvals->all_lambda", -1, fep1->all_lambda[i][j], fep2->all_lambda[i][j], ftol, abstol);
        }
    }
    cmp_int(fp, "inputrec->fepvals->lambda_neighbors", 1, fep1->lambda_neighbors,
            fep2->lambda_neighbors);
    cmp_real(fp, "inputrec->fepvals->sc_alpha", -1, fep1->sc_alpha, fep2->sc_alpha, ftol, abstol);
    cmp_int(fp, "inputrec->fepvals->sc_power", -1, fep1->sc_power, fep2->sc_power);
    cmp_real(fp, "inputrec->fepvals->sc_r_power", -1, fep1->sc_r_power, fep2->sc_r_power, ftol, abstol);
    cmp_real(fp, "inputrec->fepvals->sc_sigma", -1, fep1->sc_sigma, fep2->sc_sigma, ftol, abstol);
    cmp_int(fp, "inputrec->fepvals->edHdLPrintEnergy", -1, fep1->edHdLPrintEnergy, fep1->edHdLPrintEnergy);
    cmp_bool(fp, "inputrec->fepvals->bScCoul", -1, fep1->bScCoul, fep1->bScCoul);
    cmp_int(fp, "inputrec->separate_dhdl_file", -1, fep1->separate_dhdl_file, fep2->separate_dhdl_file);
    cmp_int(fp, "inputrec->dhdl_derivatives", -1, fep1->dhdl_derivatives, fep2->dhdl_derivatives);
    cmp_int(fp, "inputrec->dh_hist_size", -1, fep1->dh_hist_size, fep2->dh_hist_size);
    cmp_double(fp, "inputrec->dh_hist_spacing", -1, fep1->dh_hist_spacing, fep2->dh_hist_spacing, ftol, abstol);
}

void cmp_inputrec(FILE *fp, const t_inputrec *ir1, const t_inputrec *ir2, real ftol, real abstol)
{
    fprintf(fp, "comparing inputrec\n");

    /* gcc 2.96 doesnt like these defines at all, but issues a huge list
     * of warnings. Maybe it will change in future versions, but for the
     * moment I've spelled them out instead. /EL 000820
     * #define CIB(s) cmp_int(fp,"inputrec->"#s,0,ir1->##s,ir2->##s)
     * #define CII(s) cmp_int(fp,"inputrec->"#s,0,ir1->##s,ir2->##s)
     * #define CIR(s) cmp_real(fp,"inputrec->"#s,0,ir1->##s,ir2->##s,ftol)
     */
    cmp_int(fp, "inputrec->eI", -1, ir1->eI, ir2->eI);
    cmp_int64(fp, "inputrec->nsteps", ir1->nsteps, ir2->nsteps);
    cmp_int64(fp, "inputrec->init_step", ir1->init_step, ir2->init_step);
    cmp_int(fp, "inputrec->simulation_part", -1, ir1->simulation_part, ir2->simulation_part);
    cmp_int(fp, "inputrec->ePBC", -1, ir1->ePBC, ir2->ePBC);
    cmp_bool(fp, "inputrec->bPeriodicMols", -1, ir1->bPeriodicMols, ir2->bPeriodicMols);
    cmp_int(fp, "inputrec->cutoff_scheme", -1, ir1->cutoff_scheme, ir2->cutoff_scheme);
    cmp_int(fp, "inputrec->ns_type", -1, ir1->ns_type, ir2->ns_type);
    cmp_int(fp, "inputrec->nstlist", -1, ir1->nstlist, ir2->nstlist);
    cmp_int(fp, "inputrec->nstcomm", -1, ir1->nstcomm, ir2->nstcomm);
    cmp_int(fp, "inputrec->comm_mode", -1, ir1->comm_mode, ir2->comm_mode);
    cmp_int(fp, "inputrec->nstlog", -1, ir1->nstlog, ir2->nstlog);
    cmp_int(fp, "inputrec->nstxout", -1, ir1->nstxout, ir2->nstxout);
    cmp_int(fp, "inputrec->nstvout", -1, ir1->nstvout, ir2->nstvout);
    cmp_int(fp, "inputrec->nstfout", -1, ir1->nstfout, ir2->nstfout);
    cmp_int(fp, "inputrec->nstcalcenergy", -1, ir1->nstcalcenergy, ir2->nstcalcenergy);
    cmp_int(fp, "inputrec->nstenergy", -1, ir1->nstenergy, ir2->nstenergy);
    cmp_double(fp, "inputrec->init_t", -1, ir1->init_t, ir2->init_t, ftol, abstol);
    cmp_double(fp, "inputrec->delta_t", -1, ir1->delta_t, ir2->delta_t, ftol, abstol);
    cmp_real(fp, "inputrec->fourierspacing", -1, ir1->fourier_spacing, ir2->fourier_spacing, ftol, abstol);
    cmp_int(fp, "inputrec->nkx", -1, ir1->nkx, ir2->nkx);
    cmp_int(fp, "inputrec->nky", -1, ir1->nky, ir2->nky);
    cmp_int(fp, "inputrec->nkz", -1, ir1->nkz, ir2->nkz);
    cmp_int(fp, "inputrec->pme_order", -1, ir1->pme_order, ir2->pme_order);
    cmp_real(fp, "inputrec->ewald_rtol", -1, ir1->ewald_rtol, ir2->ewald_rtol, ftol, abstol);
    cmp_int(fp, "inputrec->ewald_geometry", -1, ir1->ewald_geometry, ir2->ewald_geometry);
    cmp_real(fp, "inputrec->epsilon_surface", -1, ir1->epsilon_surface, ir2->epsilon_surface, ftol, abstol);
    cmp_int(fp, "inputrec->bContinuation", -1, static_cast<int>(ir1->bContinuation), static_cast<int>(ir2->bContinuation));
    cmp_int(fp, "inputrec->bShakeSOR", -1, static_cast<int>(ir1->bShakeSOR), static_cast<int>(ir2->bShakeSOR));
    cmp_int(fp, "inputrec->etc", -1, ir1->etc, ir2->etc);
    cmp_int(fp, "inputrec->bPrintNHChains", -1, static_cast<int>(ir1->bPrintNHChains), static_cast<int>(ir2->bPrintNHChains));
    cmp_int(fp, "inputrec->epc", -1, ir1->epc, ir2->epc);
    cmp_int(fp, "inputrec->epct", -1, ir1->epct, ir2->epct);
    cmp_real(fp, "inputrec->tau_p", -1, ir1->tau_p, ir2->tau_p, ftol, abstol);
    cmp_rvec(fp, "inputrec->ref_p(x)", -1, ir1->ref_p[XX], ir2->ref_p[XX], ftol, abstol);
    cmp_rvec(fp, "inputrec->ref_p(y)", -1, ir1->ref_p[YY], ir2->ref_p[YY], ftol, abstol);
    cmp_rvec(fp, "inputrec->ref_p(z)", -1, ir1->ref_p[ZZ], ir2->ref_p[ZZ], ftol, abstol);
    cmp_rvec(fp, "inputrec->compress(x)", -1, ir1->compress[XX], ir2->compress[XX], ftol, abstol);
    cmp_rvec(fp, "inputrec->compress(y)", -1, ir1->compress[YY], ir2->compress[YY], ftol, abstol);
    cmp_rvec(fp, "inputrec->compress(z)", -1, ir1->compress[ZZ], ir2->compress[ZZ], ftol, abstol);
    cmp_int(fp, "refcoord_scaling", -1, ir1->refcoord_scaling, ir2->refcoord_scaling);
    cmp_rvec(fp, "inputrec->posres_com", -1, ir1->posres_com, ir2->posres_com, ftol, abstol);
    cmp_rvec(fp, "inputrec->posres_comB", -1, ir1->posres_comB, ir2->posres_comB, ftol, abstol);
    cmp_real(fp, "inputrec->verletbuf_tol", -1, ir1->verletbuf_tol, ir2->verletbuf_tol, ftol, abstol);
    cmp_real(fp, "inputrec->rlist", -1, ir1->rlist, ir2->rlist, ftol, abstol);
    cmp_real(fp, "inputrec->rtpi", -1, ir1->rtpi, ir2->rtpi, ftol, abstol);
    cmp_int(fp, "inputrec->coulombtype", -1, ir1->coulombtype, ir2->coulombtype);
    cmp_int(fp, "inputrec->coulomb_modifier", -1, ir1->coulomb_modifier, ir2->coulomb_modifier);
    cmp_real(fp, "inputrec->rcoulomb_switch", -1, ir1->rcoulomb_switch, ir2->rcoulomb_switch, ftol, abstol);
    cmp_real(fp, "inputrec->rcoulomb", -1, ir1->rcoulomb, ir2->rcoulomb, ftol, abstol);
    cmp_int(fp, "inputrec->vdwtype", -1, ir1->vdwtype, ir2->vdwtype);
    cmp_int(fp, "inputrec->vdw_modifier", -1, ir1->vdw_modifier, ir2->vdw_modifier);  cmp_real(fp, "inputrec->rvdw_switch", -1, ir1->rvdw_switch, ir2->rvdw_switch, ftol, abstol);
    cmp_real(fp, "inputrec->rvdw", -1, ir1->rvdw, ir2->rvdw, ftol, abstol);
    cmp_real(fp, "inputrec->epsilon_r", -1, ir1->epsilon_r, ir2->epsilon_r, ftol, abstol);
    cmp_real(fp, "inputrec->epsilon_rf", -1, ir1->epsilon_rf, ir2->epsilon_rf, ftol, abstol);
    cmp_real(fp, "inputrec->tabext", -1, ir1->tabext, ir2->tabext, ftol, abstol);

    cmp_int(fp, "inputrec->eDispCorr", -1, ir1->eDispCorr, ir2->eDispCorr);
    cmp_real(fp, "inputrec->shake_tol", -1, ir1->shake_tol, ir2->shake_tol, ftol, abstol);
    cmp_int(fp, "inputrec->efep", -1, ir1->efep, ir2->efep);
    cmp_fepvals(fp, ir1->fepvals, ir2->fepvals, ftol, abstol);
    cmp_int(fp, "inputrec->bSimTemp", -1, static_cast<int>(ir1->bSimTemp), static_cast<int>(ir2->bSimTemp));
    if ((ir1->bSimTemp == ir2->bSimTemp) && (ir1->bSimTemp))
    {
        cmp_simtempvals(fp, ir1->simtempvals, ir2->simtempvals, std::min(ir1->fepvals->n_lambda, ir2->fepvals->n_lambda), ftol, abstol);
    }
    cmp_int(fp, "inputrec->bExpanded", -1, static_cast<int>(ir1->bExpanded), static_cast<int>(ir2->bExpanded));
    if ((ir1->bExpanded == ir2->bExpanded) && (ir1->bExpanded))
    {
        cmp_expandedvals(fp, ir1->expandedvals, ir2->expandedvals, std::min(ir1->fepvals->n_lambda, ir2->fepvals->n_lambda), ftol, abstol);
    }

    cmp_int(fp, "inputrec->eDisre", -1, ir1->eDisre, ir2->eDisre);
    cmp_real(fp, "inputrec->dr_fc", -1, ir1->dr_fc, ir2->dr_fc, ftol, abstol);
    cmp_int(fp, "inputrec->eDisreWeighting", -1, ir1->eDisreWeighting, ir2->eDisreWeighting);
    cmp_int(fp, "inputrec->bDisreMixed", -1, static_cast<int>(ir1->bDisreMixed), static_cast<int>(ir2->bDisreMixed));
    cmp_int(fp, "inputrec->nstdisreout", -1, ir1->nstdisreout, ir2->nstdisreout);
    cmp_real(fp, "inputrec->dr_tau", -1, ir1->dr_tau, ir2->dr_tau, ftol, abstol);
    cmp_real(fp, "inputrec->orires_fc", -1, ir1->orires_fc, ir2->orires_fc, ftol, abstol);
    cmp_real(fp, "inputrec->orires_tau", -1, ir1->orires_tau, ir2->orires_tau, ftol, abstol);
    cmp_int(fp, "inputrec->nstorireout", -1, ir1->nstorireout, ir2->nstorireout);
    cmp_real(fp, "inputrec->em_stepsize", -1, ir1->em_stepsize, ir2->em_stepsize, ftol, abstol);
    cmp_real(fp, "inputrec->em_tol", -1, ir1->em_tol, ir2->em_tol, ftol, abstol);
    cmp_int(fp, "inputrec->niter", -1, ir1->niter, ir2->niter);
    cmp_real(fp, "inputrec->fc_stepsize", -1, ir1->fc_stepsize, ir2->fc_stepsize, ftol, abstol);
    cmp_int(fp, "inputrec->nstcgsteep", -1, ir1->nstcgsteep, ir2->nstcgsteep);
    cmp_int(fp, "inputrec->nbfgscorr", 0, ir1->nbfgscorr, ir2->nbfgscorr);
    cmp_int(fp, "inputrec->eConstrAlg", -1, ir1->eConstrAlg, ir2->eConstrAlg);
    cmp_int(fp, "inputrec->nProjOrder", -1, ir1->nProjOrder, ir2->nProjOrder);
    cmp_real(fp, "inputrec->LincsWarnAngle", -1, ir1->LincsWarnAngle, ir2->LincsWarnAngle, ftol, abstol);
    cmp_int(fp, "inputrec->nLincsIter", -1, ir1->nLincsIter, ir2->nLincsIter);
    cmp_real(fp, "inputrec->bd_fric", -1, ir1->bd_fric, ir2->bd_fric, ftol, abstol);
    cmp_int64(fp, "inputrec->ld_seed", ir1->ld_seed, ir2->ld_seed);
    cmp_real(fp, "inputrec->cos_accel", -1, ir1->cos_accel, ir2->cos_accel, ftol, abstol);
    cmp_rvec(fp, "inputrec->deform(a)", -1, ir1->deform[XX], ir2->deform[XX], ftol, abstol);
    cmp_rvec(fp, "inputrec->deform(b)", -1, ir1->deform[YY], ir2->deform[YY], ftol, abstol);
    cmp_rvec(fp, "inputrec->deform(c)", -1, ir1->deform[ZZ], ir2->deform[ZZ], ftol, abstol);


    cmp_int(fp, "inputrec->userint1", -1, ir1->userint1, ir2->userint1);
    cmp_int(fp, "inputrec->userint2", -1, ir1->userint2, ir2->userint2);
    cmp_int(fp, "inputrec->userint3", -1, ir1->userint3, ir2->userint3);
    cmp_int(fp, "inputrec->userint4", -1, ir1->userint4, ir2->userint4);
    cmp_real(fp, "inputrec->userreal1", -1, ir1->userreal1, ir2->userreal1, ftol, abstol);
    cmp_real(fp, "inputrec->userreal2", -1, ir1->userreal2, ir2->userreal2, ftol, abstol);
    cmp_real(fp, "inputrec->userreal3", -1, ir1->userreal3, ir2->userreal3, ftol, abstol);
    cmp_real(fp, "inputrec->userreal4", -1, ir1->userreal4, ir2->userreal4, ftol, abstol);
    cmp_grpopts(fp, &(ir1->opts), &(ir2->opts), ftol, abstol);
    gmx::TextWriter writer(fp);
    gmx::compareKeyValueTrees(&writer, *ir1->params, *ir2->params, ftol, abstol);
}

gmx_bool inputrecDeform(const t_inputrec *ir)
{
    return (ir->deform[XX][XX] != 0 || ir->deform[YY][YY] != 0 || ir->deform[ZZ][ZZ] != 0 ||
            ir->deform[YY][XX] != 0 || ir->deform[ZZ][XX] != 0 || ir->deform[ZZ][YY] != 0);
}

gmx_bool inputrecDynamicBox(const t_inputrec *ir)
{
    return (ir->epc != epcNO || ir->eI == eiTPI || inputrecDeform(ir));
}

gmx_bool inputrecPreserveShape(const t_inputrec *ir)
{
    return  (ir->epc != epcNO && ir->deform[XX][XX] == 0 &&
             (ir->epct == epctISOTROPIC || ir->epct == epctSEMIISOTROPIC));
}

gmx_bool inputrecNeedMutot(const t_inputrec *ir)
{
    return ((ir->coulombtype == eelEWALD || EEL_PME(ir->coulombtype)) &&
            (ir->ewald_geometry == eewg3DC || ir->epsilon_surface != 0));
}

gmx_bool inputrecExclForces(const t_inputrec *ir)
{
    return (EEL_FULL(ir->coulombtype) || (EEL_RF(ir->coulombtype)));
}

gmx_bool inputrecNptTrotter(const t_inputrec *ir)
{
    return ( ( (ir->eI == eiVV) || (ir->eI == eiVVAK) ) &&
             (ir->epc == epcMTTK) && (ir->etc == etcNOSEHOOVER) );
}

gmx_bool inputrecNvtTrotter(const t_inputrec *ir)
{
    return ( ( (ir->eI == eiVV) || (ir->eI == eiVVAK) ) &&
             (ir->epc != epcMTTK) && (ir->etc == etcNOSEHOOVER) );
}

gmx_bool inputrecNphTrotter(const t_inputrec *ir)
{
    return ( ( (ir->eI == eiVV) || (ir->eI == eiVVAK) ) &&
             (ir->epc == epcMTTK) && (ir->etc != etcNOSEHOOVER) );
}

bool integratorHasConservedEnergyQuantity(const t_inputrec *ir)
{
    if (!EI_MD(ir->eI))
    {
        // Energy minimization or stochastic integrator: no conservation
        return false;
    }
    else if (ir->etc == etcNO && ir->epc == epcNO)
    {
        // The total energy is conserved, no additional conserved quanitity
        return false;
    }
    else
    {
        // Shear stress with Parrinello-Rahman is not supported (tedious)
        bool shearWithPR =
            ((ir->epc == epcPARRINELLORAHMAN || ir->epc == epcMTTK) &&
             (ir->ref_p[YY][XX] != 0 || ir->ref_p[ZZ][XX] != 0 || ir->ref_p[ZZ][YY] != 0));

        return !ETC_ANDERSEN(ir->etc) && !shearWithPR;
    }
}

bool integratorHasReferenceTemperature(const t_inputrec *ir)
{
    return ((ir->etc != etcNO) || EI_SD(ir->eI) || (ir->eI == eiBD) || EI_TPI(ir->eI));
}

int inputrec2nboundeddim(const t_inputrec *ir)
{
    return ePBC2npbcdim(ir->ePBC);
}

int ndof_com(const t_inputrec *ir)
{
    int n = 0;

    switch (ir->ePBC)
    {
        case epbcXYZ:
        case epbcNONE:
            n = 3;
            break;
        case epbcXY:
            n = 3;
            break;
        case epbcSCREW:
            n = 1;
            break;
        default:
            gmx_incons("Unknown pbc in calc_nrdf");
    }

    return n;
}

real maxReferenceTemperature(const t_inputrec &ir)
{
    if (EI_ENERGY_MINIMIZATION(ir.eI) || ir.eI == eiNM)
    {
        return 0;
    }

    if (EI_MD(ir.eI) && ir.etc == etcNO)
    {
        return -1;
    }

    /* SD and BD also use ref_t and tau_t for setting the reference temperature.
     * TPI can be treated as MD, since it needs an ensemble temperature.
     */

    real maxTemperature = 0;
    for (int i = 0; i < ir.opts.ngtc; i++)
    {
        if (ir.opts.tau_t[i] >= 0)
        {
            maxTemperature = std::max(maxTemperature, ir.opts.ref_t[i]);
        }
    }

    return maxTemperature;
}

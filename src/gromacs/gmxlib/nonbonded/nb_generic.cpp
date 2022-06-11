/*
 * This file is part of the GROMACS molecular simulation package.
 *
 * Copyright (c) 1991-2000, University of Groningen, The Netherlands.
 * Copyright (c) 2001-2004, The GROMACS development team.
 * Copyright (c) 2012,2014,2015,2017,2018, by the GROMACS development team, led by
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

#include "nb_generic.h"

#include <cmath>

#include "act/coulombintegrals/coulombintegrals.h"
#include "gromacs/gmxlib/nrnb.h"
#include "gromacs/gmxlib/nonbonded/nb_kernel.h"
#include "gromacs/gmxlib/nonbonded/nonbonded.h"
#include "gromacs/math/functions.h"
#include "gromacs/math/vec.h"
#include "gromacs/mdtypes/md_enums.h"
#include "gromacs/utility/fatalerror.h"

void wang_buckingham(real sigma, real epsilon, real gamma, 
                     real rsq, real rinv,
                     real *vvdw, real *fvdw)
{
    /* Modified Buckingham: JCTC  Volume: 9  Page: 452  Year: 2012 */
    real r           = rsq*rinv;
    real r5          = rsq*rsq*r;
    real r6          = r5*r;
    real sigma2      = sigma*sigma;
    real sigma6      = sigma2*sigma2*sigma2;
    real sigma6_r6   = sigma6 + r6;
    real gamma_3     = 3.0/(gamma + 3.0);
    real gamma3_inv  = 1.0/(1.0 - gamma_3);
    real disp_pre    = 2 * epsilon * gamma3_inv;
    real erep_exp    = gamma_3*std::exp(gamma*(1-(r/sigma)));
    
    real vvdw_disp   = - disp_pre * (sigma6 / sigma6_r6);
    real vvdw_rep    = -vvdw_disp*erep_exp;
    *vvdw            = vvdw_rep + vvdw_disp;
    
    real fvdw_disp   = - disp_pre * sigma6 * 6 * r5 / (sigma6_r6 * sigma6_r6);
    real fvdw_rep    = - fvdw_disp*erep_exp - vvdw_disp*erep_exp*(gamma/sigma);
    *fvdw            = fvdw_rep + fvdw_disp;
}

void coulomb_gaussian(real qq, real izeta, real jzeta,
                      real r, real *velec, real *felec)
{
    *velec       = qq*Coulomb_GG(r, izeta, jzeta);
    *felec       = qq*DCoulomb_GG(r, izeta, jzeta);
}

void
gmx_nb_generic_kernel(t_nblist *                nlist,
                      rvec *                    xx,
                      rvec *                    ff,
                      t_forcerec *              fr,
                      t_mdatoms *               mdatoms,
                      nb_kernel_data_t *        kernel_data,
                      t_nrnb *                  nrnb)
{
    int           ntype, ielec, ivdw;
    real          facel;
    int           n, ii, is3, ii3, k, nj0, nj1, jnr, j3, ggid;
    real          shX, shY, shZ;
    real          felec, fvdw, velec, vvdw, tx, ty, tz;
    real          rinvsq;
    real          iq;
    real          qq, vctot;
    int           nti, nvdwparam;
    int           tj;
    real          r;
    real          rinvsix;
    real          vvdwtot;
    real          vvdw_rep, vvdw_disp;
    real          ix, iy, iz, fix, fiy, fiz;
    real          jx, jy, jz;
    real          dx, dy, dz, rsq, rinv;
    real          izeta, jzeta, irow, jrow;
    real          c6, c12;
    real *        charge;
    real *        zeta;
    real *        row;
    real *        shiftvec;
    real *        vdwparam;
    int *         type;
    real *        fshift;
    real *        velecgrp;
    real *        vvdwgrp;
    real *        x;
    real *        f;
    real          rcutoff, rcutoff2;
    gmx_bool      bExactElecCutoff, bExactVdwCutoff, bExactCutoff;

    x                   = xx[0];
    f                   = ff[0];
    ielec               = nlist->ielec;
    ivdw                = nlist->ivdw;

    fshift              = fr->fshift[0];
    velecgrp            = kernel_data->energygrp_elec;
    vvdwgrp             = kernel_data->energygrp_vdw;

    const interaction_const_t *ic = fr->ic;

    bExactElecCutoff    = (ic->coulomb_modifier != eintmodNONE) || ic->eeltype == eelRF_ZERO;
    bExactVdwCutoff     = (ic->vdw_modifier != eintmodNONE);
    bExactCutoff        = bExactElecCutoff && bExactVdwCutoff;

    if (bExactCutoff)
    {
        rcutoff  = ( ic->rcoulomb > ic->rvdw ) ? ic->rcoulomb : ic->rvdw;
        rcutoff2 = rcutoff*rcutoff;
    }
    else
    {
        /* Fix warnings for stupid compilers */
        rcutoff2 = 1e30;
    }

    /* 3 VdW parameters for Buckingham, otherwise 2 */
    nvdwparam           = (ivdw == GMX_NBKERNEL_VDW_BUCKINGHAM) ? 3 : 2;

    charge              = mdatoms->chargeA;
    zeta                = mdatoms->zetaA;
    row                 = mdatoms->row;
    type                = mdatoms->typeA;
    facel               = fr->ic->epsfac;
    shiftvec            = fr->shift_vec[0];
    vdwparam            = fr->nbfp;
    ntype               = fr->ntype;

    for (n = 0; (n < nlist->nri); n++)
    {
        is3              = 3*nlist->shift[n];
        shX              = shiftvec[is3];
        shY              = shiftvec[is3+1];
        shZ              = shiftvec[is3+2];
        nj0              = nlist->jindex[n];
        nj1              = nlist->jindex[n+1];
        ii               = nlist->iinr[n];
        ii3              = 3*ii;
        ix               = shX + x[ii3+0];
        iy               = shY + x[ii3+1];
        iz               = shZ + x[ii3+2];
        iq               = facel*charge[ii];
        izeta            = zeta[ii];
        irow             = row[ii];
        nti              = nvdwparam*ntype*type[ii];
        vctot            = 0;
        vvdwtot          = 0;
        fix              = 0;
        fiy              = 0;
        fiz              = 0;

        for (k = nj0; (k < nj1); k++)
        {
            jnr              = nlist->jjnr[k];
            j3               = 3*jnr;
            jx               = x[j3+0];
            jy               = x[j3+1];
            jz               = x[j3+2];
            dx               = ix - jx;
            dy               = iy - jy;
            dz               = iz - jz;
            rsq              = dx*dx+dy*dy+dz*dz;
            if (rsq > 0)
            {
                rinv   = gmx::invsqrt(rsq);
                rinvsq = rinv*rinv;
                r      = rsq*rinv;
            }
            else
            {
                r      = 0;
                rinv   = 0;
                rinvsq = 0;
            }
            felec            = 0;
            fvdw             = 0;
            velec            = 0;
            vvdw             = 0;

            if (bExactCutoff && rsq >= rcutoff2)
            {
                continue;
            }

            /* Coulomb interaction. ielec==0 means no interaction */
            if (ielec != GMX_NBKERNEL_ELEC_NONE)
            {
                qq            = iq*charge[jnr];
                jzeta         = zeta[jnr];
                jrow          = row[jnr];
                switch (ielec)
                {
                    case GMX_NBKERNEL_ELEC_NONE:
                        break;

                    case GMX_NBKERNEL_ELEC_COULOMB:
                        /* Vanilla or Gaussian cutoff coulomb */
                        if (irow == 0 && jrow == 0)
                        {
                            coulomb_gaussian(qq, izeta, jzeta,
                                             rsq*rinv, &velec, &felec);
                        }
                        else
                        {
                            if (irow < 0 || jrow < 0)
                            {
                                gmx_fatal(FARGS, "Row cannot be negative in Slater wave function!.\n");
                            }
                            else
                            {
                                velec =  qq*Coulomb_SS(r, irow, jrow, izeta, jzeta);
                                felec = -qq*DCoulomb_SS(r, irow, jrow, izeta, jzeta);
                            }
                        }
                        if (debug)
                        {
                            fprintf(debug, "velec: %0.3f r: %0.3f izeta: %0.3f jzeta: %0.3f irow: %0.3f jrow: %0.3f i: %d j: %d iq: %0.3f jq: %0.3f\n", 
                                    velec, r, izeta, jzeta, irow, jrow, ii, jnr, charge[ii], charge[jnr]);
                        }
                        break;

                    case GMX_NBKERNEL_ELEC_REACTIONFIELD:
                    case GMX_NBKERNEL_ELEC_CUBICSPLINETABLE:
                    case GMX_NBKERNEL_ELEC_EWALD:
                    default:
                        gmx_fatal(FARGS, "Death & horror! No generic coulomb interaction for ielec=%d.\n", ielec);
                }
                vctot           += velec;
            } /* End of coulomb interactions */
            if (debug)
            {
                fprintf(debug, "vctot: %0.3f \n", vctot);
            }

            /* VdW interaction. ivdw==0 means no interaction */
            if (ivdw != GMX_NBKERNEL_VDW_NONE)
            {
                tj               = nti+nvdwparam*type[jnr];
                
                switch (ivdw)
                {
                    case GMX_NBKERNEL_VDW_NONE:
                        break;

                    case GMX_NBKERNEL_VDW_LENNARDJONES:
                        /* Vanilla Lennard-Jones cutoff */
                        c6               = vdwparam[tj];
                        c12              = vdwparam[tj+1];
                        rinvsix          = rinvsq*rinvsq*rinvsq;
                        vvdw_disp        = c6*rinvsix;
                        vvdw_rep         = c12*rinvsix*rinvsix;
                        fvdw             = (12*vvdw_rep-6*vvdw_disp)*rinv;
                        vvdw             = vvdw_rep-vvdw_disp;
                        if (debug)
                        {
                            fprintf(debug, "ai %d aj %d vvdw: %10g c6: %10g c12: %10g\n", ii, jnr, vvdw, c6, c12);
                        }
                        break;

                    case GMX_NBKERNEL_VDW_BUCKINGHAM:
                        wang_buckingham(/*sigma*/   vdwparam[tj+1],
                                        /*epsilon*/ vdwparam[tj+2],
                                        /*gamma*/   vdwparam[tj],
                                        rsq, rinv,
                                        &vvdw,
                                        &fvdw);
                        break;
                        
                    case GMX_NBKERNEL_VDW_CUBICSPLINETABLE:
                    case GMX_NBKERNEL_VDW_LJEWALD:
                    default:
                        gmx_fatal(FARGS, "Death & horror! No generic VdW interaction for ivdw=%d.\n", ivdw);
                }
                vvdwtot         += vvdw;
            } /* end VdW interactions */
            if (debug)
            {
                fprintf(debug, "vvdwtot: %10g \n", vvdwtot);
            }
            
            auto fscal_rinv  = rinv*(felec+fvdw);

            tx               = fscal_rinv*dx;
            ty               = fscal_rinv*dy;
            tz               = fscal_rinv*dz;
            fix              = fix + tx;
            fiy              = fiy + ty;
            fiz              = fiz + tz;
            f[j3+0]          = f[j3+0] - tx;
            f[j3+1]          = f[j3+1] - ty;
            f[j3+2]          = f[j3+2] - tz;
        }

        f[ii3+0]         = f[ii3+0] + fix;
        f[ii3+1]         = f[ii3+1] + fiy;
        f[ii3+2]         = f[ii3+2] + fiz;
        fshift[is3]      = fshift[is3]+fix;
        fshift[is3+1]    = fshift[is3+1]+fiy;
        fshift[is3+2]    = fshift[is3+2]+fiz;
        ggid             = nlist->gid[n];
        velecgrp[ggid]  += vctot;
        vvdwgrp[ggid]   += vvdwtot;
	
    }
    /* Estimate flops, average for generic kernel:
     * 12 flops per outer iteration
     * 50 flops per inner iteration
     */
    inc_nrnb(nrnb, eNR_NBKERNEL_GENERIC, nlist->nri*12 + nlist->jindex[n]*50);
}

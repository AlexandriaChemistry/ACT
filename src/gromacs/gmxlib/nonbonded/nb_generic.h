/*
 * This file is part of the GROMACS molecular simulation package.
 *
 * Copyright (c) 1991-2000, University of Groningen, The Netherlands.
 * Copyright (c) 2001-2008, The GROMACS development team.
 * Copyright (c) 2012,2014,2015, by the GROMACS development team, led by
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

#ifndef _nb_generic_h_
#define _nb_generic_h_

#include "gromacs/gmxlib/nrnb.h"
#include "gromacs/gmxlib/nonbonded/nb_kernel.h"
#include "gromacs/math/vectypes.h"
#include "gromacs/mdtypes/mdatom.h"
#include "gromacs/mdtypes/nblist.h"

void
gmx_nb_generic_kernel(t_nblist *                nlist,
                      rvec *                    x,
                      rvec *                    f,
                      t_forcerec *              fr,
                      t_mdatoms *               mdatoms,
                      nb_kernel_data_t *        kernel_data,
                      t_nrnb *                  nrnb);

static inline void gmx_unused wang_buckingham(real sigma, real epsilon, real gamma,
                                              real rsq, real rinv,
                                              real *vrepulsion,
                                              real *vdispersion,
                                              real *fvdw)
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
    
    real vvdw_disp   = -disp_pre * (sigma6 / sigma6_r6);
    *vdispersion     = vvdw_disp;
    *vrepulsion      = -vvdw_disp*erep_exp;
    
    real fvdw_disp   = - disp_pre * sigma6 * 6 * r5 / (sigma6_r6 * sigma6_r6);
    real fvdw_rep    = - fvdw_disp*erep_exp - vvdw_disp*erep_exp*(gamma/sigma);
    *fvdw            = fvdw_rep + fvdw_disp;
}

#endif

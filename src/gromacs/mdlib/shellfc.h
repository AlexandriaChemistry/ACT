/*
 * This file is part of the GROMACS molecular simulation package.
 *
 * Copyright (c) 1991-2000, University of Groningen, The Netherlands.
 * Copyright (c) 2001-2008, The GROMACS development team.
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
#ifndef GMX_MDLIB_SHELLFC_H
#define GMX_MDLIB_SHELLFC_H

#include <cstdio>

#include "gromacs/math/paddedvector.h"
#include "gromacs/mdlib/vsite.h"
#include "gromacs/timing/wallcycle.h"

struct gmx_enerdata_t;
struct gmx_groups_t;
struct gmx_multisim_t;
struct gmx_shellfc_t;
struct gmx_mtop_t;
struct t_forcerec;
struct t_fcdata;
struct t_graph;
struct t_inputrec;
class t_state;

/* Initialization function, also predicts the initial shell postions.
 */
gmx_shellfc_t *init_shell_flexcon(FILE *fplog,
                                  const gmx_mtop_t *mtop, int nflexcon,
                                  int nstcalcenergy,
                                  bool usingDomainDecomposition);

/* Get the local shell with domain decomposition */
void make_local_shells(const t_mdatoms *md,
                       gmx_shellfc_t   *shfc);

/* Optimize shell positions */
real relax_shell_flexcon(FILE                                     *log,
                         const t_commrec                          *cr,
                         const gmx_multisim_t                     *ms,
                         gmx_bool                                  bVerbose,
                         int64_t                                   mdstep,
                         const t_inputrec                         *inputrec,
                         gmx_bool                                  bDoNS,
                         int                                       force_flags,
                         gmx_localtop_t                           *top,
                         gmx_enerdata_t                           *enerd,
                         t_fcdata                                 *fcd,
                         t_state                                  *state,
                         gmx::ArrayRefWithPadding<gmx::RVec>       f,
                         tensor                                    force_vir,
                         const t_mdatoms                          *md,
                         t_nrnb                                   *nrnb,
                         gmx_wallcycle_t                           wcycle,
                         t_graph                                  *graph,
                         const gmx_groups_t                       *groups,
                         gmx_shellfc_t                            *shfc,
                         t_forcerec                               *fr,
                         double                                    t,
                         rvec                                      mu_tot,
                         const gmx_vsite_t                        *vsite);

/* Print some final output */
void done_shellfc(FILE *fplog, gmx_shellfc_t *shellfc, int64_t numSteps);

#endif

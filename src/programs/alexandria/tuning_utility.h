/*
 * This source file is part of the Alexandria Chemistry Toolkit.
 *
 * Copyright (C) 2014-2020
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

#ifndef TUNNING_UTILITY_H
#define TUNNING_UTILITY_H

#include <cctype>
#include <cmath>
#include <cstdio>
#include <cstdlib>

#include <vector>

#include "gromacs/commandline/pargs.h"
#include "gromacs/commandline/viewit.h"
#include "gromacs/fileio/xvgr.h"
#include "gromacs/mdtypes/commrec.h"
#include "gromacs/statistics/statistics.h"
#include "gromacs/topology/mtop_util.h"
#include "gromacs/utility/coolstuff.h"
#include "gromacs/utility/cstringutil.h"
#include "gromacs/utility/fatalerror.h"
#include "gromacs/utility/futil.h"
#include "gromacs/utility/smalloc.h"

#include "gmx_simple_comm.h"
#include "molgen.h"
#include "optparam.h"
#include "poldata.h"
#include "poldata_xml.h"

namespace alexandria
{

void print_stats(FILE        *fp,
                 const char  *prop,
                 gmx_stats_t  lsq,
                 gmx_bool     bHeader,
                 char        *xaxis,
                 char        *yaxis);

void print_lsq_set(FILE *fp, gmx_stats_t lsq);

void xvgr_symbolize(FILE                   *xvgf,
                    int                     nsym,
                    const char             *leg[],
                    const gmx_output_env_t *oenv);

void print_polarizability(FILE              *fp,
                          alexandria::MyMol *mol,
                          char              *calc_name,
                          real               alpha_toler,
                          real               isopol_toler);

void print_dipole(FILE                      *fp,
                  alexandria::MyMol         *mol,
                  char                      *calc_name,
                  real                       toler);

void print_quadrapole(FILE                  *fp,
                      alexandria::MyMol     *mol,
                      char                  *calc_name,
                      real                   toler);

void print_electric_props(FILE                           *fp,
                          std::vector<alexandria::MyMol> &mymol,
                          const Poldata                  *pd,
                          const gmx::MDLogger            &fplog,
                          real                            hfac,
                          const char                     *lot,
                          const char                     *tabfn,
                          gmx_hw_info_t                  *hwinfo,
                          int                             qcycle,
                          real                            qtol,
                          const char                     *qhisto,
                          const char                     *dipcorr,
                          const char                     *mucorr,
                          const char                     *Qcorr,
                          const char                     *espcorr,
                          const char                     *alphacorr,
                          const char                     *isopolCorr,
                          const char                     *anisopolCorr,
                          const char                     *qCorr,
                          real                            esp_toler,
                          real                            dip_toler,
                          real                            quad_toler,
                          real                            alpha_toler,
                          real                            isopol_toler,
                          const gmx_output_env_t         *oenv,
                          bool                            bPolar,
                          bool                            bDipole,
                          bool                            bQuadrupole,
                          bool                            bfullTensor,
                          IndexCount                     *indexCount,
                          t_commrec                      *cr,
                          real                            efield,
                          bool                            useOffset);

/*! \brief Print header and command line arguments
 *
 * \param[in] fp    File pointer, if nullptr the function returns 
 *                  without doing anything
 * \param[in] pargs The command line arguments
 */
void print_header(FILE                       *fp, 
                  const std::vector<t_pargs> &pargs);

}

#endif

/*
 * This source file is part of the Alexandria Chemistry Toolkit.
 *
 * Copyright (C) 2014-2021
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

#ifndef TUNING_UTILITY_H
#define TUNING_UTILITY_H

#include <cstdio>

#include <vector>

#include "gromacs/commandline/pargs.h"
#include "gromacs/mdtypes/commrec.h"
#include "gromacs/statistics/statistics.h"
#include "gromacs/topology/mtop_util.h"

#include "mymol.h"
#include "poldata.h"

namespace alexandria
{

void print_electric_props(FILE                           *fp,
                          std::vector<alexandria::MyMol> *mymol,
                          const Poldata                  *pd,
                          const gmx::MDLogger            &fplog,
                          const char                     *lot,
                          const char                     *tabfn,
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
                          bool                            bfullTensor,
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

} // namespace alexandria

#endif

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

#ifndef POLDATA_TABLES_H
#define POLDATA_TABLES_H

#include <stdio.h>

#include "poldata.h"

namespace alexandria
{
    /*! \brief
     * Generates a LaTeX table containing the zeta values 
     * for gaussian and/or slater charge models
     *
     * \param[out] fp   File pointer to write to
     * \param[in]  pd   Force field data
     */
    void alexandria_charge_table(FILE            *fp, 
                                 const Poldata   *pd);  
                               
    /*! \brief
     * Generates a LaTeX table containing the chi and eta values 
     * for Alexandria Charge Models
     *
     * \param[out] fp   File pointer to write to
     * \param[in]  pd   Force field data
     */
    void alexandria_eemprops_table(FILE            *fp, 
                                   const Poldata   *pd);  
                               
    /*! \brief
     * Generates a LaTeX table containing the delta_chi and 
     * bond hardness values for Alexandria Charge Models
     *
     * \param[out] fp   File pointer to write to
     * \param[in]  pd   Force field data
     */
    void alexandria_eemprops_corr(const Poldata  *pd,
                                  FILE           *fp);                             
    /*! \brief
     * Generates a LaTeX table containing the subtypes
     * for particle types
     *
     * \param[out] fp   File pointer to write to
     * \param[in]  pd   Force field data
     */
    void alexandria_subtype_table(FILE           *fp,
                                  const Poldata  *pd);
             
} //namespace

#endif

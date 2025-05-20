/*
 * This source file is part of the Alexandria Chemistry Toolkit.
 *
 * Copyright (C) 2014-2025
 *
 * Developers:
 *             Mohammad Mehdi Ghahremanpour, 
 *             Julian Marrades,
 *             Marie-Madeleine Walz,
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

#ifndef FORCEFIELD_TABLES_H
#define FORCEFIELD_TABLES_H

#include <stdio.h>

#include "forcefield.h"

namespace alexandria
{

class ForceFieldTable
{
 private:
    //! File to write to
    FILE  *fp_            = nullptr;

    //! Force field
    const ForceField *pd_ = nullptr;

    //! Minimum number of training points to print
    unsigned int ntrain_  = 1;

    //! Whether to print sigma
    bool printSigma_      = true;
 public:
    /*! \brief
     * Generates a LaTeX tables containing force field data.
     *
     * \param[in] fp     File pointer to write to
     * \param[in] pd     Force field data
     * \param[in] ntrain Minimum number of training points to include
     */
    ForceFieldTable(FILE              *fp, 
                    const ForceField  *pd,
                    int                ntrain,
                    bool               printSigma = true) : fp_(fp), pd_(pd), ntrain_(ntrain), printSigma_(printSigma) {}

    /*! \brief
     * Generates a LaTeX table containing the chi and eta values 
     * for Alexandria Charge Models
     *
     * \param[in] info   Text to add to the caption of the table
     */
    void eemprops_table(const std::string &info);  

    /*! \brief
     * Generates a LaTeX table containing the delta_chi and 
     * bond hardness values for Alexandria Charge Models
     *
     * \param[in] info   Text to add to the caption of the table
     */
    void eemprops_corr_table(const std::string &info);

    /*! \brief
     * Generates a LaTeX table containing the subtypes
     * for particle types
     *
     * \param[in]info Text to add to the caption of the table
     */
    void subtype_table(const std::string &info);

    /*! \brief
     * Generate a LaTeX table for the interaction specified.
     *
     * \param[in] info   Text to add to the caption of the table
     */
    void itype_table(InteractionType    itype,
                     const std::string &info);

};

} //namespace

#endif

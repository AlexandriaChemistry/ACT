/*
 * This source file is part of the Alexandria Chemistry Toolkit.
 *
 * Copyright (C) 2023-2025
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
 * \author David van der Spoel <david.vanderspoel@icm.uu.se>
 */
#ifndef COMBRULEUTIL_H
#define COMBRULEUTIL_H

#include <vector>

#include "act/forcefield/forcefield.h"
#include "act/forcefield/forcefield_parameterlist.h"
#include "gromacs/commandline/pargs.h"

namespace alexandria
{

    class CombRuleUtil
    {
    private:
        // Storage for the command line options
        std::vector<const char *> cr_flag_;
        // Storage the the combination rule
        char *rules_;
    public:
        /*! \brief Utility to make command line information about combrules
         * \param[inout] crinfo Array of strings to be edited
         */
        void addInfo(std::vector<const char *> *crinfo);
        
        /*! \brief Utility to add combination rule arguments to command line
         * \param[inout] pa      The list of arguments
         */ 
        void addPargs(std::vector<t_pargs> *pa);

        /*! \brief Utility to convert strings to combination rules in the FF
         * \param[inout] pd  The force field
         * \return the number of rules that were changed.
         */
        int extract(ForceField *pd);

        /*! \brief Utility to convert old-style combination rule to new
         * \param[inout] vdw     The parameter list
         * \return the number of rules that were changed.
         */
        int convert(ForceFieldParameterList *vdw);
    };

} // namespace aleaxndria
#endif

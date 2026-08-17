/*
 * This source file is part of the Alexandria Chemistry Toolkit.
 *
 * Copyright (C) 2024-2026
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
 
#ifndef ACTMOL_UTIL_H
#define ACTMOL_UTIL_H
    
#include <vector>

namespace alexandria
{
    class  ForceField;
    class  ForceComputer;
    class  ACTMol;
    
    /*! \brief Read a file from the test directories and produce a vector of actmols.
     * \param[in]  molname  The name of the compound
     * \param[in]  pd       The force field
     * \param[in]  fcomp    The force computer
     * \param[out] mps      The ACTMol structures
     */
    void initACTMol(const char          *molname, 
                    ForceField          *pd,
                    ForceComputer       *fcomp,
                    std::vector<ACTMol> *mps);

}

#endif

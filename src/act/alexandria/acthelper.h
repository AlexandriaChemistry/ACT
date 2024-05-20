/*
 * This source file is part of the Alexandria Chemistry Toolkit.
 *
 * Copyright (C) 2021-2024
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
 * \author Julian Ramon Marrades Furquet <julian.marrades@hotmail.es>
 */
#ifndef ACT_ACTHELPER_H
#define ACT_ACTHELPER_H

#include "gromacs/mdtypes/commrec.h"
#include "act/forces/forcecomputer.h"

namespace alexandria
{
    // Declare some classes rather than including headers.
    class StaticIndividualInfo;
    class ACMFitnessComputer;
    class MolGen;
   
    /*! \brief Class that runs the deviation calculations only
     */
    class ACTHelper
    {
    private:
        //! Fitness computer
        ACMFitnessComputer *fitComp_;
        //! Force computer
        ForceComputer      *forceComp_;
    public:
        /*! \brief Constructor
         * \param[in] sii          Static information
         * \param[in] mg           Information about this helpers molecules
         * \param[in] shellToler   Tolerance for minimizing shell positions
         * \param[in] shellMaxIter Max # iterations for the same
         */
        ACTHelper(StaticIndividualInfo *sii,
                  MolGen               *mg,
                  double                shellToler,
                  int                   shellMaxIter);
        
        //! \brief Run the helper process
        void run();    
    };
    
} // namespace alexandria

#endif

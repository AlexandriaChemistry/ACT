/*
 * This source file is part of the Alexandria Chemistry Toolkit.
 *
 * Copyright (C) 2021-2026
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
 * \author Julian Ramon Marrades Furquet <julian@marrad.es>
 */

#ifndef ALEXANDRIA_SENSITIVITY_H
#define ALEXANDRIA_SENSITIVITY_H

#include "act/basics/dataset.h"

namespace ga
{
    class Genome;
}

namespace alexandria
{

    class MsgHandler;
    class StaticIndividualInfo;
    class ACMFitnessComputer;
    class JsonTree;
    
    /*!
     * \brief Perform a sensitivity analysis by systematically changing all parameters and
     * re-evaluating the \f$ \chi^2 \f$.
     * \param[in] msghandler The message and status handler
     * \param[in] genome     Pointer to genome
     * \param[in] ims        Dataset to perform sensitivity analysis on
     * \param[in] jtree      For machine readable output
     * \return true if the current genome represents a local minimum in parameter space
     */
    bool SensitivityAnalysis(MsgHandler           *msghandler,
                             StaticIndividualInfo *sii,
                             ACMFitnessComputer   *fitComp,
                             ga::Genome           *genome,
                             iMolSelect            ims,
                             JsonTree             *jtree);
    
} // namespace

#endif

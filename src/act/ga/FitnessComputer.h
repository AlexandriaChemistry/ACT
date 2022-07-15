/*
 * This source file is part of the Alexandria Chemistry Toolkit.
 *
 * Copyright (C) 2014-2022
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
 * \author Julian Ramon Marrades Furquet <julian.marrades@hotmail.es>
 */


#ifndef GA_FITNESSCOMPUTER_H
#define GA_FITNESSCOMPUTER_H

#include "Genome.h"

namespace ga
{

/*!
 * \brief Abstract class for computing the fitness of an individual
 */
class FitnessComputer
{

public:

    /*!
     * \brief Compute the fitness of a genome
     * \param[in] genome  The genome
     * \param[in] trgtFit The target for fitness computation. Either Train or Test
     * \param[in] verbose Whether to print the components of the fitness
     */
    virtual void compute(Genome                    *genome,
                         iMolSelect                 trgtFit,
                         bool                       verbose = false) = 0;
    
};


} //namespace ga


#endif //GA_FITNESSCOMPUTER_H

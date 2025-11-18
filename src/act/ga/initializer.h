/*
 * This source file is part of the Alexandria Chemistry Toolkit.
 *
 * Copyright (C) 2014-2022,2025
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
 * \author Julian Ramon Marrades Furquet <julian@marrad.es>
 */

#ifndef GA_INITIALIZER_H
#define GA_INITIALIZER_H

#include "individual.h"

#include <memory>
#include <random>

namespace ga
{

/*!
* \brief Abstract class for initializing and randomizing Individual/Genome instances
*/
class Initializer
{

public:

    /*!
     * \brief Initialize an Individual
     * \return   pointer to a new individual
     */
    virtual Individual *initialize() = 0;

    //! Default destructor
    virtual ~Initializer() = default;

    /*! \brief Randomize a Genome object
     * \param[in] genome the Genome to randomize
     */
    virtual void randomizeGenome(Genome *genome) = 0;
};


} //namespace ga


#endif //GA_INITIALIZER_H

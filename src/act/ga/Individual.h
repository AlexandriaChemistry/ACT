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
#ifndef GA_INDIVIDUAL_H
#define GA_INDIVIDUAL_H

#include <cstdio>

#include <map>

#include "Genome.h"

namespace ga
{

//! \brief Abstract individual for genetic algorithms
class Individual
{

protected:

    //! \brief Default constructor
    Individual() {}

public:

    /* * * * * * * * * * * * * * * * * * * * * *
    * BEGIN: Cloning                           *
    * * * * * * * * * * * * * * * * * * * * * */

    /*!
     * \brief Copy the genome from another Individual
     * \param[in] other pointer to another Individual
     */
    virtual void copyGenome(const Genome &other) = 0;

    /* * * * * * * * * * * * * * * * * * * * * *
    * END: Cloning                             *
    * * * * * * * * * * * * * * * * * * * * * */
};


} //namespace ga


#endif //GA_INDIVIDUAL_H

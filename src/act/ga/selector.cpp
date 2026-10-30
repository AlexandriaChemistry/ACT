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
 * \author Julian Ramon Marrades Furquet <julian@marrad.es>
 */


#include "selector.h"

#include <algorithm>
#include <random>

#include "gromacs/utility/gmxassert.h"


namespace ga
{


/* * * * * * * * * * * * * * * * * * * * * *
* BEGIN: RouletteSelector                  *
* * * * * * * * * * * * * * * * * * * * * */

int RouletteSelector::select(const std::vector<Genome> *pop)
{
    GMX_RELEASE_ASSERT(!pop->empty(), "Cannot select from an empty population");
    double num = dis(gen);
    size_t i = 0;
    while (num > 0 and i < pop->size())
    {
        num -= (*pop)[i].probability();
        i++;
    }
    // A draw of exactly zero selects the first genome
    return std::max(i, static_cast<size_t>(1)) - 1;
}

/* * * * * * * * * * * * * * * * * * * * * *
* END: RouletteSelector                    *
* * * * * * * * * * * * * * * * * * * * * */


} //namespace ga

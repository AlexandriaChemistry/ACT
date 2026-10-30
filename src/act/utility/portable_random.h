/*
 * This source file is part of the Alexandria Chemistry Toolkit.
 *
 * Copyright (C) 2026
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
 * Random numbers that are the same on all platforms. The C++ standard
 * specifies the output of the random number engines, but not the
 * algorithms of the integer distributions and of std::shuffle, which
 * differ between standard libraries, such that a seed would give
 * different results on different platforms.
 * \author David van der Spoel <david.vanderspoel@icm.uu.se>
 */
#ifndef ACT_PORTABLE_RANDOM_H
#define ACT_PORTABLE_RANDOM_H

#include <cstdint>
#include <random>
#include <utility>
#include <vector>

namespace alexandria
{

/*! \brief Draw a non-negative seed for another random number generator
 * \param[inout] gen The random number engine
 * \return a seed in [0, 2^31)
 */
inline int portableSeed(std::mt19937 *gen)
{
    return static_cast<int>((*gen)() >> 1);
}

/*! \brief Draw an index
 * \param[inout] gen The random number engine
 * \param[in]    n   The number of indices, at most 2^32
 * \return an index in [0, n)
 */
inline size_t portableIndex(std::mt19937 *gen, size_t n)
{
    return static_cast<size_t>((static_cast<uint64_t>((*gen)()) * n) >> 32);
}

/*! \brief Shuffle a vector with the Fisher-Yates algorithm
 * \param[inout] v   The vector
 * \param[inout] gen The random number engine
 */
template <typename T>
void portableShuffle(std::vector<T> *v, std::mt19937 *gen)
{
    for (size_t i = v->size(); i > 1; i--)
    {
        std::swap((*v)[i-1], (*v)[portableIndex(gen, i)]);
    }
}

} // namespace alexandria

#endif

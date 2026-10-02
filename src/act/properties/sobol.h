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
 * Implements part of the alexandria program.
 * \author David van der Spoel <david.vanderspoel@icm.uu.se>
 */
#ifndef ACT_ALEXANDRIA_SOBOL_H
#define ACT_ALEXANDRIA_SOBOL_H

#include <vector>

/*! \brief Generate a Sobol sequence for integrating in max 6D
 * Useful for second virial calculations.
 * Code reimplemented based on the numerical recipes book.
 */
class SobolSequence
{
private:
    //! Index in the Sobol sequence
    unsigned int              index = 0;
    //! Integer vector
    std::vector<unsigned int> sobol_dim;
    //! Integer vector
    std::vector<unsigned int> int_vec;
public:
    //! Constructor
    SobolSequence();
    /*! \brief Extract a sequence of quasirandom numbers
     * \param[in]  ndim The number of dimensions requested, should be <= 6
     * \param[out] q    Pointer to a vector of doubles of length ndim (or larger)
     */
    void seq(int ndim, std::vector<double> *q);
};

#endif

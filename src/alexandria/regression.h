/*
 * This source file is part of the Alexandria Chemistry Toolkit.
 *
 * Copyright (C) 2014-2020 
 *
 * Developers:
 *             Mohammad Mehdi Ghahremanpour, 
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
 */
 
 
#ifndef ALEXANDRIA_REGRESSION_H
#define ALEXANDRIA_REGRESSION_H

#include <vector>

#include "gromacs/math/vectypes.h"

class MatrixWrapper
{
    public:
        /*! \brief Constructor
         *
         * \param[in] ncolumn Number of columns in the matrix
         * \param[in] nrow    Number of rows
         */
        MatrixWrapper(int ncolumn, int nrow);

        //! \brief Destructor in charge of freeing memory
        ~MatrixWrapper();

        /*! \brief Set a value in the matrix
         *
         * \param[in] col    The column
         * \param[in] row    The row
         * \param[in] value  The value
         */
        void set(int col, int row, double value)
        {
            a_[col][row] = value;
        }

        /*! \brief Get a value from the matrix
         *
         * \param[in] col    The column
         * \param[in] row    The row
         * \returns The value
         */
        double get(int col, int row) const
        {
            return a_[col][row];
        }

        /*! \brief Solve a matrix equation A solution = rhs
         *
         * \param[in]  rhs      Vector of right hand side values.
         * \param[out] solution Pointer to vector of solution
         * \return 0 if ok, non-zero if not.
         */
        int solve(std::vector<double> rhs, std::vector<double> *solution);
    private:
        // The internal data structure
        double **a_;
};

void kabsch_rotation(tensor p, tensor q, tensor rotated_p);

#endif

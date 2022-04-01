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
        //! \brief Default constructor
        MatrixWrapper() {}

        /*! \brief Constructor
         *
         * \param[in] ncolumn Number of columns in the matrix
         * \param[in] nrow    Number of rows
         */
        MatrixWrapper(int ncolumn, int nrow);

        /*! \brief Copy constructor
         * \param[in] source the matrix to copy values from
         */
        MatrixWrapper(const MatrixWrapper &source);

        /*! \brief Move constructor
         * \param[in] source the matrix to move values from
         */
        MatrixWrapper(MatrixWrapper &&source);

        /*! \brief Copy assignment operator
         * \param[in] rhs matrix to copy from
         */
        MatrixWrapper &operator=(const MatrixWrapper &rhs);

        /*! \brief Move assignment operator
         * \param[in] rhs matrix to move from
         */
        MatrixWrapper &operator=(MatrixWrapper &&rhs);

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

        //! \return the number of columns
        int nColumn() const { return ncol_; }

        //! \return the number of rows
        int nRow() const { return nrow_; }

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

        /*! \brief Get a flattened version of the matrix
         * \param[in] order 'R' if row-wise ([row1 - row2 - ...]), or 'C' if column-wise
         *                  ([col1^T - col2^T - ...]), where T stands for transpose
         * \return a vector as the flattened version of the matrix
         */
        std::vector<double> flatten(const char order) const;

        /*! \brief Average the lower and upper triangular parts of the matrix
         * Beware! The matrix should be square! Otherwise, the assert condition will fail.
         */
        void averageTriangle();

        /*! \brief Solve a matrix equation A solution = rhs
         *
         * \param[in]  rhs      Vector of right hand side values.
         * \param[out] solution Pointer to vector of solution
         * \return 0 if ok, non-zero if not.
         */
        int solve(std::vector<double> rhs, std::vector<double> *solution);
    private:
        // The internal data structure
        double **a_ = nullptr;
        // The number of columns
        int      ncol_;
        // The number of rows
        int      nrow_;
};

void kabsch_rotation(tensor p, tensor q, tensor rotated_p);

#endif

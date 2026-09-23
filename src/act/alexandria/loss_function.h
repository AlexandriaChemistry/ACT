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
#ifndef ALEXANDRIA_LOSS_FUNCTION_H
#define ALEXANDRIA_LOSS_FUNCTION_H

#include <string>

namespace alexandria
{

enum class LossFunction {
    //! Mean squared error
    MSE,
    //! Mean absolute error
    MAE,
    //! Huber function, https://doi.org/10.1214/aoms/1177703732
    Huber,
    //! Asinh function, https://doi.org/10.1063/5.0280032 
    Asinh
};

//! \return a string corresponding to the loss function
const char *lossFunctionName(LossFunction lf);

//! \return a LossFunction corresponding to a string
LossFunction lossFunction(const std::string &lf);

/*! \brief Compute the loss according to the chose loss function
 * \param[in] lf     The loss function
 * \param[in] a      Energy scaling factor
 * \param[in] deltaE Difference between reference and calculated value
 * \return the loss
 */
static inline double loss(LossFunction lf,
                          double       a,
                          double       deltaE)
{
    switch (lf)
    {
    case LossFunction::MSE:
        return 0.5*deltaE*deltaE;
    case LossFunction::MAE:
        return std::abs(deltaE);
    case LossFunction::Huber:
        {
            double dea = deltaE/a;
            return a*a*(std::sqrt(dea*dea+1) - 1);
        }
    case LossFunction::Asinh:
        {
            double dea  = deltaE/a;
            double sdea = std::sqrt(dea*dea+1);
            return a*a*(1 - sdea + dea*std::log(dea + sdea));
        }
    }
    return 0;
}

}

#endif

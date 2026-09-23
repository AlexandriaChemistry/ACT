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

#include "loss_function.h"

#include <cstdlib>
#include <map>
#include <string>
#include <vector>

#include "act/utility/jsontree.h"
#include "act/utility/regression.h"
#include "act/utility/stringutil.h"
#include "gromacs/utility/textwriter.h"

namespace alexandria
{

//! \brief Map from enum to string
std::map<LossFunction, const char *> lfMap = {
    { LossFunction::MSE, "MSE" },
    { LossFunction::MAE, "MAE" },
    { LossFunction::Huber, "Huber" },
    { LossFunction::Asinh, "Asinh" }
};

//! \return a string corresponding to a LossFunction.
const char *lossFunctionName(LossFunction lf)
{
    return lfMap[lf];
}


}

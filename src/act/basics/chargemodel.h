/*
 * This source file is part of the Alexandria Chemistry Toolkit.
 *
 * Copyright (C) 2014-2020,2023
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
 * \author Mohammad Mehdi Ghahremanpour <mohammad.ghahremanpour@icm.uu.se>
 * \author David van der Spoel <david.vanderspoel@icm.uu.se>
 */
#ifndef CHARGEMODEL_H
#define CHARGEMODEL_H

#include <map>
#include <string>
#include <vector>

#include "gromacs/utility/basedefinitions.h"

namespace alexandria
{

/*! \brief
 * Enumerated type holding the charge models used in ForceField files
 *
 * \inpublicapi
 * \ingroup module_alexandria
 */
enum class ChargeType {
    //! Point charges
    Point,
    //! Gaussian-distributed charges
    Gaussian,
    //! Slater-distributed charges
    Slater
};

//! \return the list of charge types
std::vector<std::string> qdnames();

//! \brief Return the string corresping to ct
const std::string &chargeTypeName(ChargeType ct);

//! \brief Return the ChargeType corresponding to name
ChargeType name2ChargeType(const std::string &name);

/*! \brief
 * Enumerated type holding the charge generation algorithms
 *
 * \inpublicapi
 * \ingroup module_alexandria
 */
enum class ChargeGenerationAlgorithm {
    NONE, EEM, SQE, ESP, Custom, Read
};

//! \brief Return the string corresping to cg
const std::string &chargeGenerationAlgorithmName(ChargeGenerationAlgorithm cg);

//! \brief Return the ChargeGenerationAlgorithm corresponding to name
ChargeGenerationAlgorithm nameToChargeGenerationAlgorithm(const std::string &name);

} // namespace alexandria
#endif

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
 * Checks of the parameters of a force field, see validateForceField.
 * \author David van der Spoel <david.vanderspoel@icm.uu.se>
 */
#ifndef ACT_FORCEFIELD_VALIDATION_H
#define ACT_FORCEFIELD_VALIDATION_H

#include <string>
#include <vector>

namespace alexandria
{

class ForceField;

//! A problem with a parameter of a force field, see validateForceField
struct ForceFieldIssue
{
    //! The kinds of problems
    enum class Kind
    {
        //! The value was outside the bounds, it was set to the nearest bound when it was read
        ResetToBounds,
        //! The value is outside the bounds, but the parameter may not be changed
        OutOfBounds,
        //! The minimum is larger than the maximum
        MinimumAboveMaximum,
        //! The value is negative for a parameter marked as non-negative
        Negative,
        //! The parameter was defined more than once, the first definition was used
        Duplicate
    };
    //! The kind of problem
    Kind        kind;
    //! The interaction type, or "particle type"
    std::string where;
    //! The identifier of the entry
    std::string identifier;
    //! The type of the parameter
    std::string parameter;
    //! What is wrong and what was done
    std::string text;
    //! Whether reading the force field repaired it, such that writing it removes the problem
    bool        repaired;
};

/*! \brief Check the parameters of a force field
 * The values have to be within their bounds, the minimum must not be larger
 * than the maximum, values of parameters marked as non-negative must not be
 * negative, and each parameter may be defined only once. Reading a force field sets
 * values of trained parameters (mutability Bounded) outside their bounds
 * to the nearest bound, and uses the first definition of a parameter
 * that is defined more than once, see ForceFieldParameter::wasReset and
 * ForceFieldParameterList::duplicates. These count as problems as well,
 * which writing the force field removes.
 * \param[in] pd The force field
 * \return the problems, empty if there are none
 */
std::vector<ForceFieldIssue> validateForceField(const ForceField &pd);

/*! \brief A summary of the problems that reading repaired, for a warning
 * \param[in] issues The problems, see validateForceField
 * \return a sentence, empty if reading repaired nothing
 */
std::string repairedIssuesSummary(const std::vector<ForceFieldIssue> &issues);

} // namespace alexandria

#endif

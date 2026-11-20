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
/*! \internal \file
 * \brief
 * Declares functions that describe a force field in a JsonTree, for
 * programs that read the output of the alexandria tools, such as the Python
 * interface and the graphical user interface.
 *
 * \author David van der Spoel <david.vanderspoel@icm.uu.se>
 */
#ifndef ACT_FORCEFIELD_FORCEFIELD_JSON_H
#define ACT_FORCEFIELD_FORCEFIELD_JSON_H

#include <string>
#include <vector>

#include "act/forcefield/forcefield_validation.h"
#include "act/utility/jsontree.h"

namespace alexandria
{

class ForceField;
class ForceFieldParameter;

/*! \brief Describe one parameter
 *
 * The tree has the type of the parameter, such as sigma, its unit, value,
 * uncertainty, minimum and maximum, the number of training compounds
 * (ntrain), the mutability and whether the parameter may not be negative
 * (nonnegative, a boolean).
 *
 * \param[in] type  The type of the parameter
 * \param[in] param The parameter
 * \param[in] key   The key of the tree, which is not written for an
 *                  element of an array
 * \return the tree
 */
JsonTree parameterToJson(const std::string         &type,
                         const ForceFieldParameter &param,
                         const std::string         &key = "parameter");

/*! \brief Describe a force field with everything that is in its file
 *
 * The tree has the name of the file, the checksum and the time stamp,
 * whether the force field is polarizable, the algorithm that generates
 * charges, the number of exclusions and the relative dielectric constant.
 * Then:
 * - "particletypes" with the "options" that hold for all types and an
 *   array "types": for each type the identifier, the type (Atom, Shell,
 *   VSite, ...), the description, its options and an array "parameters".
 * - an array "interactions": for each interaction its type, function,
 *   canswap, options, "combination_rules" (for each parameter the rule and
 *   the exponent as a parameter) and an array "parameterlists" with the
 *   identifier and the parameters of each entry. Entries of which all
 *   parameters are derived from others (mutability Dependent) are left
 *   out, as they are in the force field file.
 * - an array "symmetric_charges".
 * The order is that of the force field file. Options are texts, as in the
 * file.
 *
 * \param[in] pd  The force field
 * \param[in] key The key of the tree
 * \return the tree
 */
JsonTree forceFieldToJson(const ForceField &pd, const std::string &key = "forcefield");

/*! \brief The name of a kind of problem in the JSON output
 * \param[in] kind The kind
 * \return the name, for instance OutOfBounds
 */
const char *forceFieldIssueKindName(ForceFieldIssue::Kind kind);

/*! \brief Describe the problems that validateForceField found
 *
 * The tree has "ok", true if there are no problems, the number of problems
 * that writing the force field removes ("repaired") and an array "issues":
 * for each problem its kind, where it is (the interaction type or
 * "particle type"), the identifier, the parameter, a text that says what
 * is wrong, and whether reading the force field repaired it.
 *
 * \param[in] issues The problems
 * \param[in] key    The key of the tree
 * \return the tree
 */
JsonTree forceFieldIssuesToJson(const std::vector<ForceFieldIssue> &issues,
                                const std::string                  &key = "validation");

} // namespace alexandria

#endif

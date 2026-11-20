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
 * Declares functions that describe the content of molprop files in a
 * JsonTree, for programs that read the output of the alexandria tools,
 * such as the Python interface and the graphical user interface.
 *
 * \author David van der Spoel <david.vanderspoel@icm.uu.se>
 */
#ifndef ACT_MOLPROP_MOLPROP_JSON_H
#define ACT_MOLPROP_MOLPROP_JSON_H

#include <string>
#include <vector>

#include "act/basics/msg_handler.h"
#include "act/molprop/molprop_util.h"
#include "act/utility/jsontree.h"

namespace alexandria
{

class MolProp;

/*! \brief Describe a molecule of a molprop file, without its data
 *
 * The tree has the name, the IUPAC name, the formula, the number of atoms
 * (that of the fragments together) and the total charge (the sum of that
 * of the fragments), an array "fragments"
 * with for each fragment the InChI, the IUPAC name, the formula, the charge,
 * the multiplicity and the number of atoms, the number of experiments (the
 * calculations or measurements, for instance the conformations of a
 * dimer), "jobtypes" with the number of experiments of each kind, the
 * methods as "method/basisset" and "properties" with, for each property
 * that is present, the number of values over all experiments.
 *
 * \param[in] mp  The molecule
 * \param[in] key The key of the tree
 * \return the tree
 */
JsonTree molPropToJson(const MolProp &mp, const std::string &key = "molecule");

/*! \brief The totals of a set of molecules
 *
 * The tree has the number of molecules, of monomers (one fragment), dimers
 * (two) and larger clusters, of experiments, and "properties" with the
 * number of values of each property.
 *
 * \param[in] mps The molecules
 * \param[in] key The key of the tree
 * \return the tree
 */
JsonTree molPropSummaryToJson(const std::vector<MolProp> &mps, const std::string &key = "summary");

/*! \brief Describe the result of validateMolProps
 *
 * The tree has "checked" (true), "ok", the number of problems, an array
 * "invalid" with the names of the molecules that have problems, and an
 * array "messages" with the level (Error, Warning), the kind and the text
 * of each message.
 *
 * \param[in] validation The result
 * \param[in] messages   The messages that were kept during the validation
 * \param[in] key        The key of the tree
 * \return the tree
 */
JsonTree molPropValidationToJson(const MolPropValidation               &validation,
                                 const std::vector<ACTRecordedMessage> &messages,
                                 const std::string                     &key = "validation");

} // namespace alexandria

#endif

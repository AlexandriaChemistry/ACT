/*
 * This source file is part of the Alexandria Chemistry Toolkit.
 *
 * Copyright (C) 2014-2026
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

#ifndef MOLPROP_UTIL_H
#define MOLPROP_UTIL_H

#include <cstddef>
#include <set>
#include <string>
#include <vector>

#include "act/basics/msg_handler.h"
#include "act/molprop/molprop.h"
#include "gromacs/utility/arrayref.h"

namespace alexandria
{
/*! \brief Reset the molprop index number
 * \param[in] mp The molprops
 */
void generate_index(std::vector<MolProp> *mp);

/*! \brief Will merge multiple molprop files.
 * \param[in] msg_handler For info and debugging
 * \param[in] infiles     MolProp filenames
 * \param[in] mp          The molprop vector
 */
void merge_xml(MsgHandler                             *msg_handler,
               const gmx::ArrayRef<const std::string> &infiles,
               std::vector<alexandria::MolProp>       *mp);

//! The result of validateMolProps
struct MolPropValidation
{
    //! The number of problems
    size_t                nproblem = 0;
    //! The names of the molecules with problems
    std::set<std::string> invalid;
};

/*! \brief Check that a set of molecules is consistent, as edit_mp -validate
 * First each molecule is checked with MolProp::validate, which removes the
 * experiments with the atoms in another order than the reference
 * experiment. Then all molecules are checked together:
 * - each fragment has an InChI as its identifier, and
 * - each compound, a fragment identified by its InChI and formula, has its
 *   atoms in the same order (names and types) in all molecules, such that it can be
 *   treated as one compound. The reference is the order of the compound in
 *   the references, if it is there, or else the order that most molecules
 *   have, the first of those if several are equally common, and
 * - the fragments of a dimer A#B are in the order of A and B, if their
 *   names are those of the dimer. Both A#B and B#A may occur.
 * Each problem is a message of msghandler at the given level.
 * \param[in]    msghandler For the problems
 * \param[inout] mps        The molecules
 * \param[inout] references Molecules whose order of the atoms of a compound
 *                          is the reference, for instance the monomers of
 *                          train_ff -charges, or nullptr. They are checked
 *                          with MolProp::validate as well.
 * \param[in]    level      The level of the messages of the problems
 *                          found here, for instance a warning when
 *                          train_ff leaves out the molecules
 * \return The problems and the molecules that have them
 */
MolPropValidation validateMolProps(MsgHandler           *msghandler,
                                   std::vector<MolProp> *mps,
                                   std::vector<MolProp> *references = nullptr,
                                   ACTStatus             level      = ACTStatus::Error);

/*! \brief A molecule made of one fragment of another molecule
 * The atoms of the fragment in the first experiment that has atoms, with
 * the bonds between them, in a new experiment of an optimization, for
 * instance a monomer of a dimer.
 * \param[in] mp The molecule
 * \param[in] f  The fragment
 * \return The new molecule, named after the fragment
 */
MolProp extractFragment(const MolProp &mp, size_t f);

/*! \brief One molecule for each compound, from which its charges are made
 * For each compound in the molecules, identified by its InChI and formula,
 * the molecule is the first monomer of that compound in the references,
 * for instance those of train_ff -charges, or else the first monomer in
 * all molecules, or else the compound as a fragment of the first of the
 * molecules it is in, see extractFragment. Which one is used is reported
 * at Info level.
 * \param[in] msghandler For the messages
 * \param[in] molecules  The molecules whose compounds are needed
 * \param[in] all        All molecules, whose monomers can be used
 * \param[in] references Monomers to use first, or nullptr
 * \return One molecule per compound
 */
std::vector<MolProp> compoundReferences(MsgHandler                         *msghandler,
                                        const std::vector<const MolProp *> &molecules,
                                        const std::vector<MolProp>         &all,
                                        const std::vector<MolProp>         *references);

} // namespace alexandria

#endif

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

#ifndef ACT_IMPORT_IMPORT_H
#define ACT_IMPORT_IMPORT_H

#include <string>
#include <vector>

#include "act/forcefield/forcefield.h"
#include "act/molprop/experiment.h"
#include "act/molprop/molprop.h"

namespace alexandria
{

class MsgHandler;

/*! \brief
 * Read an sdf, pdb or xyz file using RDKit.
 *
 * \param[in]    msg_handler For debugging and info
 * \param[in]    pd          Alexandria force field
 * \param[in]    filenm      The file to read
 * \param[out]   mp          Pointer to a MolProp vector
 * \param[in]    conf        Conformation the molecule is in [ ignored if nullptr ]
 * \param[in]    jobtype     Calculation type for reading QM output
 * \param[in]    userqtot    Whether the user explicitly set the total charge. If set,
 *                           the qtot below will be used instead of what is read from the input file.
 * \param[out]   qtot        Total charge as deduced by OB from the input. 
 * \param[in]    oneH        Remap all hydrogen atom types to "h"
 */
void importFile(MsgHandler           *msg_handler,
                const ForceField     *pd,
                const std::string    &filenm,
                std::vector<MolProp> *mp,
                const char           *conf,
                JobType               jobtype,
                bool                  userqtot,
                double               *qtot,
                bool                  oneH);

/*! \brief
 * Make a molecule from a SMILES string using RDKit.
 *
 * Hydrogen atoms are added, three-dimensional coordinates are generated
 * with the ETKDG method and optimized roughly with the MMFF94 force field,
 * or with UFF if MMFF94 does not have the parameters. The coordinates are
 * a starting point: minimize them with the ACT force field before use.
 * They are the same every time for the same string, since the random
 * numbers start from a fixed seed.
 *
 * A string that RDKit cannot interpret, for which no coordinates can be
 * made, or that is longer than maxSmilesLength gives an error in the
 * message handler and no molecule.
 *
 * \param[in]    msg_handler For errors, debugging and info
 * \param[in]    pd          Alexandria force field
 * \param[in]    smiles      The SMILES string, e.g. CC(=O)C for acetone
 * \param[out]   mp          Pointer to a MolProp vector, to which the molecule is added
 * \param[in]    oneH        Remap all hydrogen atom types to "h"
 */
void importSmiles(MsgHandler           *msg_handler,
                  const ForceField     *pd,
                  const std::string    &smiles,
                  std::vector<MolProp> *mp,
                  bool                  oneH);

//! The longest SMILES string that importSmiles accepts
const size_t maxSmilesLength = 2000;

}

#endif

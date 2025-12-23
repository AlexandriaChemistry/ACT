/*
 * This source file is part of the Alexandria Chemistry Toolkit.
 *
 * Copyright (C) 2014-2024
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

#include <vector>

#include "act/molprop/molprop.h"

class MsgHandler;

/*! \brief
 * Read a Gaussian log file or other file supported by OpenBabel
 *
 * \param[in]    msg_handler For debugging and info
 * \param[in]    pd          Alexandria force field
 * \param[in]    filenm      The file to read
 * \param[out]   mp          Pointer to a MolProp vector
 * \param[in]    molnm       Molecule name to override the one from the filename [ ignored if nullptr ]
 * \param[in]    iupac       IUPAC name to override the one from the filename [ ignored if nullptr ]
 * \param[in]    conf        Conformation the molecule is in [ ignored if nullptr ]
 * \param[inout] method      Theoretical chemistry method detected in the file.
 * \param[inout] basis       Basis set used for the calculation. If not empty, will override
 *                           what is detected in the file. If empty, will get the detected value.
 * \param[in]    maxpot      Maximum number of electrostatic potential data points to store
 * \param[in]    nsymm       Symmetry number for this molecule. If zero it will be detected from
 *                           the input.
 * \param[in]    jobtype     Calculation type for reading QM output
 * \param[in]    userqtot    Whether the user explicitly set the total charge. If set,
 *                           the qtot below will be used instead of what is read from the input file.
 * \param[out]   qtot        Total charge as deduced by OB from the input. 
 * \param[in]    addHydrogen If true, hydrogens will be added
 * \param[out]   box         Unit cell if present, in gromacs format.
 * \param[in]    oneH        Remap all hydrogen atom types to "h"
 */
void importFile(alexandria::MsgHandler           *msg_handler,
                const alexandria::ForceField     *pd,
                const std::string                &filenm,
                std::vector<alexandria::MolProp> *mp,
                const char                       *molnm,
                const char                       *iupac,
                const char                       *conf,
                std::string                      *method,
                std::string                      *basis,
                int                               maxpot,
                int                               nsymm,
                const char                       *jobtype,
                bool                              userqtot,
                double                           *qtot,
                bool                              addHydrogen,
                matrix                            box,
                bool                              oneH);

/*! \brief Add atomtype to a Molprop object
 *
 * The atom names and coordinates in the first Experiment
 * (calculation) are used to determine atom types using
 * routines from either OpenBabel or RDKit.
 * \param[inout] mmm Molprop object
 * \return true if successful, false otherwise
 */
bool SetMolpropAtomTypesAndBonds(alexandria::MolProp *mmm);

#endif

/*
 * This source file is part of the Alexandria Chemistry Toolkit.
 *
 * Copyright (C) 2026
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

#include <cstdio>
#include <string>
#include <vector>

#include "gromacs/fileio/pdbio.h"
#include "gromacs/math/vectypes.h"

namespace alexandria
{

class ActAtom;

/*! \brief Make unique names for the particles in each residue
 * The name is the atom name followed by the number of the atom or
 * virtual site in its residue, such as C1 or H2, the numbering starts
 * at one in every residue. A shell gets the number of the atom
 * preceding it. The OpenMM force field and the PDB file use the same
 * names, such that OpenMM can find the bonds from them.
 * \param[in] atoms The particles
 * \return The names
 */
std::vector<std::string> particleNames(const std::vector<ActAtom> &atoms);

/*! \brief Determine the formal charges of monatomic ions
 * A residue with one atom is an ion if the charges of its particles,
 * that is the atom, its shells and virtual sites and its second charge,
 * add up to a nonzero integer. Then
 * the atom gets that formal charge, such as 1+ or 2-, as in columns
 * 79-80 of a PDB file. Other particles get an empty string, since
 * ACT does not know their formal charges.
 * \param[in] atoms The particles
 * \return The formal charges
 */
std::vector<std::string> ionCharges(const std::vector<ActAtom> &atoms);

//! A bond as a PDB file has it, in its CONECT records
struct PdbBond
{
    //! The indices of the two particles, from zero
    int ai;
    int aj;
    //! The bond order, 1, 2 or 3: the number of times the bond is listed
    int order;
};

/*! \brief Write a pdb file
 * \param[in] out           The file
 * \param[in] title         The title, nullptr for none
 * \param[in] atoms         The particles
 * \param[in] x             Their coordinates (nm)
 * \param[in] residueNames  The names of the residues
 * \param[in] ePBC          The periodic boundary conditions
 * \param[in] box           The box, nullptr for none
 * \param[in] chain         The chain identifier
 * \param[in] model_nr      The number of the model
 * \param[in] index         The particles to write, all if empty
 * \param[in] bonds         The bonds, written as CONECT records with the
 *                          second particle as many times as the bond
 *                          order, as PDB files have double and triple
 *                          bonds. Empty for no CONECT records.
 * \param[in] renumberAtoms Whether the particles that are written are
 *                          numbered from one without gaps
 */
void pdbWriter(FILE                           *out,
               const char                     *title,
               const std::vector<ActAtom>     &atoms,
               const std::vector<gmx::RVec>   &x,
               const std::vector<std::string> &residueNames,
               int                             ePBC,
               const matrix                    box,
               char                            chain,
               int                             model_nr,
               const std::vector<int>         &index,
               const std::vector<PdbBond>     &bonds,
               bool                            renumberAtoms=false);

} // namespace

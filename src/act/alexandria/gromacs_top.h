/*
 * This source file is part of the Alexandria Chemistry Toolkit.
 *
 * Copyright (C) 2021-2026
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

#ifndef ACT_GROMACS_TOP_H
#define ACT_GROMACS_TOP_H

#include <cstdio>
#include <string>
#include <vector>

struct t_blocka;
struct t_excls;
struct t_mols;

namespace alexandria
{

class ActAtom;
class ForceField;
class Topology;
class MsgHandler;

//! \brief Print the gromacs t_mols
void print_top_mols(FILE *out, const char *title,
                    int nmol, const t_mols *mols);

/*! \brief Print the topology header including commercials
 * Without bItp the defaults, the atom types and the van der Waals
 * parameters of the pairs of particle types in the topology are printed.
 * GROMACS only supports Lennard-Jones interactions and point charges.
 * \param[in] msg_handler For warnings about what GROMACS does not support
 * \param[in] fp          The file to print to
 * \param[in] pd          The force field
 * \param[in] topology    The molecule
 * \param[in] commercials Lines for the comment at the top
 * \param[in] bItp        Whether this is an include topology, without
 *                        the force field
 */
void print_top_header(MsgHandler                     *msg_handler,
                      FILE                           *fp,
                      const ForceField               *pd,
                      const Topology                 *topology,
                      const std::vector<std::string> &commercials,
                      bool                            bItp);

/*! \brief Write a molecule type of a gromacs topology
 *
 * The bonded interactions, the polarization of shells, the virtual
 * sites, in the definitions of GROMACS, and the exclusions.
 * \param[in] msg_handler  For warnings etc.
 * \param[in] out          The file to print to
 * \param[in] molname      The molecule name
 * \param[in] topology     ACT structure
 * \param[in] pd           The force field
 * \param[in] atoms        The particles, with their charges and residues,
 *                         if not those of the topology: for the topology
 *                         of one compound of several, the particles of
 *                         that compound in the complete system
 * \param[in] residueNames The residue names that the residue numbers of
 *                         atoms refer to, if not those of the topology
 */
void write_top(MsgHandler                     *msg_handler,
               FILE                           *out,
               const char                     *molname,
               const Topology                 *topology,
               const ForceField               *pd,
               const std::vector<ActAtom>     *atoms = nullptr,
               const std::vector<std::string> *residueNames = nullptr);

/*! \brief Write the molecule types and the molecules of a system of compounds
 *
 * Each compound becomes a molecule of GROMACS. Compounds that are the
 * same, that is with the same particles, charges, interactions and
 * parameters, share a molecule type, and the molecules are listed in
 * the order of the particles in the system.
 * \param[in] msg_handler    For warnings etc.
 * \param[in] out            The file to print to
 * \param[in] title          The name of the system
 * \param[in] names          A name for each compound, the molecule types
 *                           are named after these
 * \param[in] compounds      The topology of each compound, in the order
 *                           of the particles in the system
 * \param[in] topology       The system, of which the charges are taken
 * \param[in] pd             The force field
 * \param[in] printMolecules Whether to print the system and its molecules
 * \return whether the topology was written. It is not, and nothing is
 * printed, if the compounds do not make up the system.
 */
bool write_top_compounds(MsgHandler                     *msg_handler,
                         FILE                           *out,
                         const char                     *title,
                         const std::vector<std::string> &names,
                         const std::vector<Topology>    &compounds,
                         const Topology                 *topology,
                         const ForceField               *pd,
                         bool                            printMolecules);

} // namespace alexandria

#endif

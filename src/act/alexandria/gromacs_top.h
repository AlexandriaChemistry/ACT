/*
 * This source file is part of the Alexandria Chemistry Toolkit.
 *
 * Copyright (C) 2021-2025
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

/*! \brief Write a gromacs topology
 * \param[in] msg_handler For warnings etc.
 * \param[in] out         The file to print to
 * \param[in] molname     The molecule name
 * \param[in] topology    ACT structure
 * \param[in] pd          The force field
 */
void write_top(MsgHandler          *msg_handler,
               FILE                *out,
               const char          *molname,
               const Topology      *topology,
               const ForceField    *pd);

} // namespace alexandria

#endif

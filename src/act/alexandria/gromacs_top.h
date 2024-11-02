/*
 * This source file is part of the Alexandria Chemistry Toolkit.
 *
 * Copyright (C) 2021-2024
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

void print_top_mols(FILE *out, const char *title,
                    int nmol, t_mols *mols);

void print_top_header(FILE                           *fp,
                      const ForceField               *pd,
                      bool                            bPol,
                      const std::vector<std::string> &commercials,
                      bool                            bItp);

void write_top(FILE                *out,
               char                *molname,
               const Topology      *topology,
               const ForceField    *pd);

} // namespace alexandria

#endif

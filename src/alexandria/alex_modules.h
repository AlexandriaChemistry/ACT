/*
 * This source file is part of the Alexandria Chemistry Toolkit.
 *
 * Copyright (C) 2014-2020 
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
 
/*! \internal \brief
 * Implements part of the alexandria program.
 * \author Mohammad Mehdi Ghahremanpour <mohammad.ghahremanpour@icm.uu.se>
 * \author David van der Spoel <david.vanderspoel@icm.uu.se>
 */
 
#ifndef GMX_ALEXMODULES_H
#define GMX_ALEXMODULES_H

namespace alexandria
{

int gentop(int argc, char *argv[]);
int simulate(int argc, char *argv[]);
int tune_ff(int argc, char *argv[]);
int edit(int argc, char *argv[]);
int bastat(int argc, char *argv[]);
int analyze(int argc, char *argv[]);
int merge_mp(int argc, char *argv[]);
int merge_pd(int argc, char *argv[]);
int mp2csv(int argc, char *argv[]);
int molprop_test(int argc, char *argv[]);
int molprop_check(int argc, char *argv[]);
}

namespace gmx
{
class CommandLineModuleManager;
} // namespace gmx

/*! \internal \brief
 * Registers all alex command-line modules.
 *
 * \param[in] manager  Command-line module manager to receive the modules.
 * \throws    std::bad_alloc if out of memory.
 *
 * Registers all modules corresponding to pre-5.0 binaries such that
 * they can be run through \p manager.
 */
void registerAlexandriaModules(gmx::CommandLineModuleManager *manager);

#endif

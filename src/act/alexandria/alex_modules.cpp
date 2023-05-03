/*
 * This source file is part of the Alexandria Chemistry Toolkit.
 *
 * Copyright (C) 2014-2023
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


#include "actpre.h"

#include "alex_modules.h"

#include <cstdio>

#include "gromacs/commandline/cmdlinemodule.h"
#include "gromacs/commandline/cmdlinemodulemanager.h"

//! Initializer for a module that defaults to nice level zero.
static void initSettingsNoNice(gmx::CommandLineModuleSettings *settings)
{
    settings->setDefaultNiceLevel(0);
}
/*! \brief
 * Convenience function for creating and registering a module.
 *
 * \param[in] manager          Module manager to which to register the module.
 * \param[in] mainFunction     Main function to wrap.
 * \param[in] name             Name for the new module.
 * \param[in] shortDescription One-line description for the new module.
 */
static void registerModule(gmx::CommandLineModuleManager                *manager,
                           gmx::CommandLineModuleManager::CMainFunction  mainFunction,
                           const char *name, const char *shortDescription)
{
    manager->addModuleCMainWithSettings(name, shortDescription, mainFunction,
                                        &initSettingsNoNice);
}

void registerAlexandriaModules(gmx::CommandLineModuleManager *manager)
{
    // Modules from this directory
    registerModule(manager, &alexandria::gentop, "gentop",
                   "Generate a molecular topology and coordinates based on structure files or quantum chemistry output from Gaussian. Only inputs for OpenMM can be generated at this point in time.");
    registerModule(manager, &alexandria::simulate, "simulate",
                   "Perform a MD simulation and generate a trajectory.");
    registerModule(manager, &alexandria::min_complex, "min_complex",
                   "Generate inputs for an energy scan.");
    registerModule(manager, &alexandria::b2, "b2",
                   "Compute second virial coefficient as a function of temperature.");
    registerModule(manager, &alexandria::nma, "nma",
                   "Perform normal mode analysis and compute thermochemistry properties.");
    registerModule(manager, &alexandria::tune_ff, "tune_ff",
                   "Optimize force field parameters.");
    registerModule(manager, &alexandria::geometry_ff, "geometry_ff",
                   "Deduce bond/angle/dihedral distributions from a set of strucures and add those to a force field file.");
    registerModule(manager, &alexandria::analyze, "analyze",
                   "Analyze molecular- or force field properties from a database and generate publication quality tables in LaTeX.");
    registerModule(manager, &alexandria::edit_ff, "edit_ff",
                   "Manipulate and compare force field files in various ways and test whether reading and writing works.");
    registerModule(manager, &alexandria::gen_ff, "gen_ff",
                   "Generate a force field file from a user specification.");
    registerModule(manager, &alexandria::edit_mp, "edit_mp",
                   "Utility to merge a number of molecular property files and a SQLite database. Can also test reading and writing the molecular property file. It can also check the molecular property file for missing hydrogens and for whether it is possible to generate topologies for all compounds. Finally it can generate charges for all compounds.");

    registerModule(manager, &alexandria::merge_ff, "merge_ff",
                   "Utility to merge a number of force field files and write a new file with average parameters. Can also write a LaTeX table.");

    {
        gmx::CommandLineModuleGroup group =
            manager->addModuleGroup("Alexandria simulation tools");
        group.addModule("nma");
        group.addModule("simulate");
        group.addModule("b2");
    }
    {
        gmx::CommandLineModuleGroup group =
            manager->addModuleGroup("Generating topologies and other simulation input");
        group.addModule("gentop");
    }
    {
        gmx::CommandLineModuleGroup group =
            manager->addModuleGroup("Utilities");
        group.addModule("geometry_ff");
        group.addModule("edit_ff");
        group.addModule("merge_ff");
        group.addModule("tune_ff");
        group.addModule("analyze");
        group.addModule("edit_mp");
    }
}

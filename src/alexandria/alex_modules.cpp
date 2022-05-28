/*
 * This source file is part of the Alexandria Chemistry Toolkit.
 *
 * Copyright (C) 2014-2022
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
                   "Generate a molecular topology and coordinates based on structure files or quantum chemistry output from Gaussian.");
    registerModule(manager, &alexandria::tune_ff, "tune_ff",
                   "Optimize force field parameters.");
    registerModule(manager, &alexandria::bastat, "bastat",
                   "Deduce bond/angle/dihedral distributions from a set of strucures and create a new force field file.");
    registerModule(manager, &alexandria::analyze, "analyze",
                   "Analyze molecular- or force field properties from a database and generate publication quality tables in LaTeX.");
    registerModule(manager, &alexandria::edit, "edit",
                   "Manipulate force field files in various ways and test whether reading and writing works.");
    registerModule(manager, &alexandria::molprop_test, "molprop_test",
                   "Test reading and writing the molecular property file.");
    registerModule(manager, &alexandria::molprop_check, "molprop_check",
                   "Check the molecular property file for missing hydrogens and for whether it is possible to generate topologies for all compounds.");
    registerModule(manager, &alexandria::mp2csv, "mp2csv",
                   "Utility to dump a molecular property file to a spreadsheet.");
    registerModule(manager, &alexandria::merge_mp, "merge_mp",
                   "Utility to merge a number of molecular property files and a SQLite database.");
    registerModule(manager, &alexandria::merge_pd, "merge_pd",
                   "Utility to merge a number of gentop files and write a new file with average parameters. Can also write a LaTeX table.");

    {
        gmx::CommandLineModuleGroup group =
            manager->addModuleGroup("Alexandria core tools");
        group.addModule("bastat");
        group.addModule("tune_ff");
        group.addModule("molprop_check");
    }
    {
        gmx::CommandLineModuleGroup group =
            manager->addModuleGroup("Generating topologies and other simulation input");
        group.addModule("gentop");
    }
    {
        gmx::CommandLineModuleGroup group =
            manager->addModuleGroup("Poldata utilities");
        group.addModule("merge_pd");
    }
    {
        gmx::CommandLineModuleGroup group =
            manager->addModuleGroup("Molprop utilities");
        group.addModule("analyze");
        group.addModule("merge_mp");
        group.addModule("molprop_test");
        group.addModule("mp2csv");
    }
    {
        gmx::CommandLineModuleGroup group =
            manager->addModuleGroup("Testing stuff and funky utilities");
        group.addModule("edit");
        group.addModule("molprop_test");
    }
}

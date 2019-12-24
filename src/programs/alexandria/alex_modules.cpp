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


#include "gmxpre.h"

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
    registerModule(manager, &alex_gentop, "gentop",
                   "Generate topology for structure files");
    registerModule(manager, &alex_tune_fc, "tune_fc",
                   "Optimize force field parameters");
    registerModule(manager, &alex_tune_eem, "tune_eem",
                   "Optimize parameters of the EEM algorithm");
    registerModule(manager, &alex_tune_pol, "tune_pol",
                   "Optimize atomic polarizabilities");
    //    registerModule(manager, &alex_tune_zeta, "tune_zeta",
    //             "Optimize the distribution of Gaussian and Slater charges");
    registerModule(manager, &alex_bastat, "bastat",
                   "Deduce bond/angle/dihedral distributions from a set of strucures");
    registerModule(manager, &alex_analyze, "analyze",
                   "Analyze molecular- or force field properties from a database and generate tables");
    registerModule(manager, &alex_gen_table, "gen_table",
                   "Generate tables for interaction functions used in mdrun");
    registerModule(manager, &alex_poldata_test, "poldata_test",
                   "Test the force field file I/O");
    registerModule(manager, &alex_fit_qs_zeta, "fit_qs_zeta",
                   "Fit qshell and zeta based on polarizability");
    registerModule(manager, &alex_molprop_test, "molprop_test",
                   "Test the molecular property file I/O");
    registerModule(manager, &alex_molprop_check, "molprop_check",
                   "Check the molecular property file for missing hydrogens");
    registerModule(manager, &alex_qm2molprop, "qm2molprop",
                   "Convert Gaussian or Psi4 output to molecular property file");
    registerModule(manager, &alex_mp2csv, "mp2csv",
                   "Utility to dump a molecular property file to a spreadsheet");
    registerModule(manager, &alex_merge_mp, "merge_mp",
                   "Utility to merge a number of molecular property files and a SQLite database");
    registerModule(manager, &alex_merge_pd, "merge_pd",
                   "Utility to merge a number of gentop files");
    registerModule(manager, &alex_molselect, "molselect",
                   "Utility to generate random samples from molprop database");

    {
        gmx::CommandLineModuleGroup group =
            manager->addModuleGroup("Alexandria core tools");
        group.addModule("bastat");
        group.addModule("tune_pol");
        //group.addModule("tune_zeta");
        group.addModule("tune_eem");
        group.addModule("tune_fc");
        group.addModule("qm2molprop");
        group.addModule("molprop_check");
    }
    {
        gmx::CommandLineModuleGroup group =
            manager->addModuleGroup("Generating topologies and other simulation input");
        group.addModule("gentop");
        group.addModule("gen_table");
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
        group.addModule("mp2csv");
        group.addModule("molselect");
    }
    {
        gmx::CommandLineModuleGroup group =
            manager->addModuleGroup("Testing stuff and funky utilities");
        group.addModule("poldata_test");
        group.addModule("molprop_test");
    }
}

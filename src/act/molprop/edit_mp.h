/*
 * This source file is part of the Alexandria Chemistry Toolkit.
 *
 * Copyright (C) 2026
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
/*! \internal \file
 * \brief
 * The editing of molprop files of alexandria edit_mp, without its command line.
 * \author David van der Spoel <david.vanderspoel@icm.uu.se>
 */
#ifndef ACT_EDIT_MP_H
#define ACT_EDIT_MP_H

#include <string>
#include <vector>

struct gmx_output_env_t;

namespace alexandria
{

class MsgHandler;

/*! \brief What alexandria edit_mp does, as its command line gives it
 *
 * The program fills it from the command line, a program without one sets
 * the members. runEditMP checks it, see check.
 */
struct EditMPOptions
{
    //! The molprop files to merge (-mp)
    std::vector<std::string> inputFiles;
    //! The molprop file to write (-o), none if empty
    std::string              outputFile;
    //! A force field to check the compounds with, which leaves out those that fail (-ff), none if empty
    std::string              forceFieldFile;
    //! A database with experimental data to add (-db), none if empty
    std::string              databaseFile;
    //! The description of the compounds that are written in JSON (-json), none if empty
    std::string              jsonFile;
    //! Whether the molprop file is compressed (-compress)
    bool                     compress        = false;
    //! The temperature in K of the properties to take from the database (-temperature)
    double                   temperature     = 298.15;
    //! Whether the names of the compounds come from the compound table of ACT (-molinfo)
    bool                     molinfo         = true;
    //! Whether to write a histogram of the energies of each compound (-ehisto)
    bool                     energyHistogram = false;
    //! A warning for an energy below this, in kJ/mol (-ewarnLow)
    double                   energyWarnLow   = -20;
    //! A warning for an energy above this, in kJ/mol (-ewarnHi)
    double                   energyWarnHigh  = 100;
    //! Whether to check the order of the atoms and compounds (-validate)
    bool                     validate        = false;

    /*! \brief Check that the options make sense
     * \throws gmx::InvalidInputError that says which option is wrong
     */
    void check() const;
};

/*! \brief Merge, check and write molprop files, as alexandria edit_mp does
 * \param[in] msghandler For output and errors
 * \param[in] options    What to do
 * \param[in] oenv       The output environment for the histograms of the
 *                       energies; one with the defaults is made if it is nullptr
 * \throws gmx::InvalidInputError for options that make no sense, see
 * EditMPOptions::check, gmx::FileIOError if a file cannot be read or written
 */
void runEditMP(MsgHandler *msghandler, const EditMPOptions &options, gmx_output_env_t *oenv = nullptr);

} // namespace alexandria

#endif

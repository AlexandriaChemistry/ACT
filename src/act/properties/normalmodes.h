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
 * The normal mode analysis of alexandria nma, without its command line.
 * \author David van der Spoel <david.vanderspoel@icm.uu.se>
 */
#ifndef ACT_NORMALMODES_H
#define ACT_NORMALMODES_H

#include <string>

#include "act/alexandria/confighandler.h"
#include "act/import/compound_reader.h"

struct gmx_output_env_t;

namespace alexandria
{

class MsgHandler;

/*! \brief What alexandria nma does, as its command line gives it
 *
 * The program fills it from the command line, a program without one sets
 * the members. runNma checks it, see check.
 */
struct NmaOptions
{
    //! The force field file (-ff)
    std::string             forceFieldFile;
    //! The minimized structure to write (-c), none if empty
    std::string             confout;
    //! The infrared spectrum to write (-ir), none if empty
    std::string             spectrumFile;
    //! The file for the results in JSON (-json); if empty they go to the log
    std::string             jsonFile;
    //! The width in cm^-1 of the lines of the spectrum (-linewidth)
    double                  lineWidth      = 24;
    //! The tolerance for relaxing the shells (-shelltoler), see SimulationConfigHandler::shellTolerance
    double                  shellTolerance = 1e-6;
    //! The settings of the minimization that comes first (-minimize, -minalg, -toler and so on)
    SimulationConfigHandler simulation;
    //! The molecule and its charges (-f, -smiles, -db, -qalg and so on)
    CompoundReader          compounds;

    //! The energy is minimized first, unless the settings say otherwise
    NmaOptions() { simulation.setMinimize(true); }

    /*! \brief Check that the options make sense
     * \throws gmx::InvalidInputError that says which option is wrong
     */
    void check() const;

    /*! \brief Check the settings, without the files and the input
     *
     * For a molecule in memory, such as one of Python. check() does this
     * as well.
     * \throws gmx::InvalidInputError that says which option is wrong
     */
    void checkSettings() const;
};

/*! \brief The normal modes of a molecule, as alexandria nma computes them
 *
 * Reads the force field and the molecule, minimizes the energy unless the
 * settings say otherwise, computes the frequencies, intensities and the
 * thermochemistry, and writes the files that the options name.
 * \param[in] msghandler For output and errors
 * \param[in] options    What to do
 * \param[in] oenv       The output environment for the spectrum; one with
 *                       the defaults is made if it is nullptr
 * \return 0 if the analysis was done, 1 if the minimization or reading the
 * input failed
 * \throws gmx::InvalidInputError for options that make no sense, see
 * NmaOptions::check
 */
int runNma(MsgHandler *msghandler, const NmaOptions &options, gmx_output_env_t *oenv = nullptr);

} // namespace alexandria

#endif

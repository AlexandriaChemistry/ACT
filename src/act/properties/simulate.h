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
 * The simulation of alexandria simulate, without its command line.
 * \author David van der Spoel <david.vanderspoel@icm.uu.se>
 */
#ifndef ACT_SIMULATE_H
#define ACT_SIMULATE_H

#include <memory>
#include <string>

#include "act/alexandria/confighandler.h"
#include "act/import/compound_reader.h"

struct gmx_output_env_t;

namespace alexandria
{

class CommunicationRecord;
class MsgHandler;
class ReRunner;

/*! \brief What alexandria simulate does, as its command line gives it
 *
 * The program fills it from the command line, a program without one sets
 * the members. runSimulate checks it, see check.
 */
struct SimulateOptions
{
    //! The force field file (-ff)
    std::string               forceFieldFile;
    //! The trajectory to write (-o), none if empty
    std::string               trajectoryFile;
    //! The energies to write (-e), none if empty
    std::string               energyFile;
    //! The final structure to write (-c), none if empty
    std::string               confout;
    //! Index file with the atoms to freeze in the minimization (-freeze), none if empty
    std::string               freezeFile;
    //! The molprop file to write (-omp), none if empty
    std::string               molpropFile;
    //! The file for part of the output in JSON (-json); if empty it goes to the log
    std::string               jsonFile;
    //! The file for the progress, a line of JSON for each energy (-progress), none if empty
    std::string               progressFile;
    //! The tolerance for relaxing the shells (-shelltoler), see SimulationConfigHandler::shellTolerance
    double                    shellTolerance = 1e-6;
    //! The settings of the minimization and the simulation (-minimize, -nsteps and so on)
    SimulationConfigHandler   simulation;
    //! The molecule and its charges (-f, -smiles, -db, -qalg and so on)
    CompoundReader            compounds;
    /*! \brief The energies of the structures of a trajectory instead of a
     * simulation (-traj and the options for dimers)
     *
     * Only the program sets it, from its command line, since the
     * generator of dimers cannot be set up otherwise yet. Without it, or
     * without a trajectory in it, the molecule is simulated.
     */
    std::shared_ptr<ReRunner> rerun;

    /*! \brief Check that the options make sense
     * \throws gmx::InvalidInputError that says which option is wrong
     */
    void check() const;
};

/*! \brief Minimize and simulate a molecule, as alexandria simulate does
 *
 * Reads the force field and the molecule, computes a single point, or
 * minimizes, simulates, and writes the files that the options name. Or,
 * with a trajectory in options.rerun, computes the energies of its
 * structures.
 * \param[in] msghandler For output and errors
 * \param[in] options    What to do
 * \param[in] oenv       The output environment for the energies; one with
 *                       the defaults is made if it is nullptr
 * \param[in] cr         The processes, for the timings and the statistics
 *                       of the shells; without, these are not given
 * \return 0 if it went well, 1 if the minimization or reading the input failed
 * \throws gmx::InvalidInputError for options that make no sense, see
 * SimulateOptions::check
 */
int runSimulate(MsgHandler                *msghandler,
                const SimulateOptions     &options,
                gmx_output_env_t          *oenv = nullptr,
                const CommunicationRecord *cr   = nullptr);

} // namespace alexandria

#endif

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

/*! \internal \brief
 * Implements part of the alexandria program.
 * \author David van der Spoel <david.vanderspoel@icm.uu.se>
 */

#ifndef ACT_MINIMIZE_H
#define ACT_MINIMIZE_H

#include <string>
#include <vector>

#include "act/alexandria/confighandler.h"
#include "act/alexandria/molhandler.h"
#include "act/import/compound_reader.h"
#include "gromacs/math/vectypes.h"

namespace alexandria
{

class ACTMol;
class ForceComputer;
class ForceField;
class JsonTree;
class MsgHandler;

/*! \brief Minimize the energy of a molecule and report the results
 *
 * The energies and multipoles are written before and after the
 * minimization, and for a dimer the interaction energies after the
 * minimization. This is used by alexandria minimize and by
 * alexandria simulate -minimize.
 * \param[in]    msghandler For output and errors
 * \param[in]    pd         The force field
 * \param[in]    forceComp  The force computer
 * \param[in]    sch        The minimization settings
 * \param[in]    freezeFile Index file with atoms to freeze, may be nullptr
 * \param[inout] actmol     The molecule
 * \param[inout] xmin       The starting coordinates, on return the minimized ones
 * \param[inout] jtree      The energies before and after are added here
 * \param[out]   epot       The potential energy after the minimization, if not nullptr
 * \return the status of the minimization, see minimumUsable
 */
eMinimizeStatus minimizeMolecule(MsgHandler                    *msghandler,
                                 const ForceField              *pd,
                                 const ForceComputer           *forceComp,
                                 const SimulationConfigHandler &sch,
                                 const char                    *freezeFile,
                                 ACTMol                        *actmol,
                                 std::vector<gmx::RVec>        *xmin,
                                 JsonTree                      *jtree,
                                 double                        *epot = nullptr);

/*! \brief What alexandria minimize does, as its command line gives it
 *
 * The program fills it from the command line, a program without one sets
 * the members. runMinimize checks it, see check.
 */
struct MinimizeOptions
{
    //! The force field file (-ff)
    std::string             forceFieldFile;
    //! The minimized structure to write (-c), none if empty
    std::string             confout;
    //! Index file with the atoms to freeze (-freeze), none if empty
    std::string             freezeFile;
    //! The molprop file to write with the minimized structure (-omp), none if empty
    std::string             molpropFile;
    //! The file for the results in JSON (-json); if empty they go to the log
    std::string             jsonFile;
    //! The tolerance for relaxing the shells (-shelltoler), see SimulationConfigHandler::shellTolerance
    double                  shellTolerance = 1e-6;
    //! The settings of the minimization (-minalg, -toler, -maxiter and so on)
    SimulationConfigHandler simulation;
    //! The molecule and its charges (-f, -smiles, -db, -qalg and so on)
    CompoundReader          compounds;

    /*! \brief Check that the options make sense
     * \throws gmx::InvalidInputError that says which option is wrong
     */
    void check() const;
};

/*! \brief Minimize the energy of a molecule as alexandria minimize does
 *
 * Reads the force field and the molecule, minimizes, and writes the files
 * that the options name.
 * \param[in] msghandler For output and errors
 * \param[in] options    What to do
 * \return 0 if the minimum or the structure with the lowest energy was
 * written, 1 if the minimization or reading the input failed
 * \throws gmx::InvalidInputError for options that make no sense, see
 * MinimizeOptions::check
 */
int runMinimize(MsgHandler *msghandler, const MinimizeOptions &options);

} // namespace alexandria

#endif

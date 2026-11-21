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
 * The topologies of alexandria gentop, without its command line.
 * \author David van der Spoel <david.vanderspoel@icm.uu.se>
 */
#ifndef ACT_GENTOP_H
#define ACT_GENTOP_H

#include <array>
#include <optional>
#include <string>

#include "act/import/compound_reader.h"

struct gmx_output_env_t;

namespace alexandria
{

class MsgHandler;

/*! \brief What alexandria gentop does, as its command line gives it
 *
 * The program fills it from the command line, a program without one sets
 * the members. runGentop checks it, see check. With several compounds
 * the names of the topology and structure files get the number of the
 * compound in front, 1_, 2_ and so on.
 */
struct GentopOptions
{
    //! The force field file (-ff)
    std::string                          forceFieldFile;
    //! The GROMACS topology to write (-p, or -oi for an include topology), none if empty
    std::string                          topologyFile;
    //! Whether the topology is an include topology (-oi)
    bool                                 includeTopology  = false;
    //! The OpenMM force field file to write (-openmm), none if empty
    std::string                          openmmFile;
    //! The structure to write (-c), none if empty
    std::string                          structureFile;
    //! The compounds with their particles in JSON (-json), none if empty
    std::string                          jsonFile;
    //! The box of the structure file in nm (-box); a cube of 5 nm if not given
    std::optional<std::array<double, 3>> box;
    //! Whether the structure file has the shells (-ws)
    bool                                 writeShells      = false;
    //! The mass of a shell in the OpenMM file, zero for that of the force field (-mDrude)
    double                               drudeMass        = 0.1;
    //! The least number of data points that a parameter in the OpenMM file was trained on (-ntrain)
    int                                  minTraining      = 1;
    //! Whether each particle type in the OpenMM file gets a number (-numberAtypes)
    bool                                 numberTypes      = true;
    //! The residue names of the compounds (-resnames), see openMMResidueNames
    std::string                          residueNames;
    //! The spacing in nm of the grid for the electrostatic potential (-spacing)
    double                               spacing          = 0.01;
    //! The space in nm around the compound for the potential (-border)
    double                               border           = 0.2;
    //! A reference potential to read (-ref), none if empty
    std::string                          referencePotential;
    //! The files about the potential to write (-pc, -pdbdiff, -pot, -rho, -his, -diff, -diffhist), none if empty
    std::string                          potentialComparison, pdbDifference, potential, density,
                                         potentialHistogram, potentialDifference, differenceHistogram;
    //! The molecules and their charges (-f, -smiles, -db, -sel, -qalg and so on)
    CompoundReader                       compounds;

    /*! \brief Check that the options make sense
     * \throws gmx::InvalidInputError that says which option is wrong
     */
    void check() const;
};

/*! \brief Make the topologies of compounds, as alexandria gentop does
 *
 * Reads the force field and the compounds, makes their charges, and
 * writes the files that the options name.
 * \param[in] msghandler For output and errors
 * \param[in] options    What to do
 * \param[in] oenv       The output environment for the histograms of the
 *                       potential; one with the defaults is made if it is
 *                       nullptr
 * \return 0 if it went well, 1 if no compound got a topology or another
 * error was reported to the message handler
 * \throws gmx::InvalidInputError for options that make no sense, see
 * GentopOptions::check
 */
int runGentop(MsgHandler *msghandler, const GentopOptions &options, gmx_output_env_t *oenv = nullptr);

} // namespace alexandria

#endif

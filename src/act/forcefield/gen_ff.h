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
 * The generation of a force field of alexandria gen_ff, without its command line.
 * \author David van der Spoel <david.vanderspoel@icm.uu.se>
 */
#ifndef ACT_GEN_FF_H
#define ACT_GEN_FF_H

#include <string>
#include <vector>

#include "act/basics/chargemodel.h"
#include "act/forcefield/potential.h"

namespace alexandria
{

/*! \brief What alexandria gen_ff does, as its command line gives it
 *
 * The program fills it from the command line, a program without one sets
 * the members. runGenFF checks it, see check.
 */
struct GenFFOptions
{
    //! The table of the particle types (-f)
    std::string              atomTypesFile;
    //! The files with virtual sites (-vs), none if empty
    std::vector<std::string> virtualSiteFiles;
    //! The force field to write (-o)
    std::string              outputFile;
    //! The force field in JSON with the result of the checks (-json), none if empty
    std::string              jsonFile;
    //! The number of bonds within which nonbonded interactions are excluded (-nexcl)
    int                      nexcl                 = 2;
    //! The relative dielectric constant (-epsilonr)
    double                   epsilonr              = 1;
    //! The largest distance in nm between a shell and its core, 0 for none (-shellmaxdistance)
    double                   shellMaxDistance      = 0.02;
    //! Whether charges of particles in a symmetric group are the same (-qsymm)
    bool                     symmetricCharges      = true;
    //! The distribution of the charges (-qdist)
    ChargeDistributionType   chargeDistribution    = ChargeDistributionType::Gaussian;
    //! The potential of the bonds (-bondfn)
    Potential                bondPotential         = Potential::HARMONIC_BONDS;
    //! The potential of the angles (-anglefn)
    Potential                anglePotential        = Potential::HARMONIC_ANGLES;
    //! The potential of the proper dihedrals (-dihfn)
    Potential                dihedralPotential     = Potential::FOURIER_DIHEDRALS;
    //! The potential of the van der Waals interactions (-vdwfn)
    Potential                vdwPotential          = Potential::LJ12_6;
    //! Whether the van der Waals correction gets the Born-Mayer potential (-vdwcorr)
    bool                     vdwCorrection         = false;
    //! The potential of the induction correction (-icfn)
    Potential                inductionCorrection   = Potential::MACDANIEL_SCHMIDT;
    //! The combination rules, separated by spaces (-cr), see CombRuleUtil
    std::string              combinationRules;

    /*! \brief Check that the options make sense
     * \throws gmx::InvalidInputError that says which option is wrong
     */
    void check() const;
};

//! \return the potentials that gen_ff -bondfn offers
const std::vector<Potential> &genFFBondPotentials();
//! \return the potentials that gen_ff -anglefn offers
const std::vector<Potential> &genFFAnglePotentials();
//! \return the potentials that gen_ff -dihfn offers
const std::vector<Potential> &genFFDihedralPotentials();
//! \return the potentials that gen_ff -vdwfn offers
const std::vector<Potential> &genFFVdwPotentials();
//! \return the potentials that gen_ff -icfn offers
const std::vector<Potential> &genFFInductionCorrectionPotentials();

/*! \brief Make a force field from a table of particle types, as alexandria gen_ff does
 * \param[in] options What to do
 * \throws gmx::InvalidInputError for options that make no sense, see
 * GenFFOptions::check, or a table that cannot be used, gmx::FileIOError
 * if a file cannot be read or written
 */
void runGenFF(const GenFFOptions &options);

} // namespace alexandria

#endif

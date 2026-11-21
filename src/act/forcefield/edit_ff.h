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
 * The editing of a force field of alexandria edit_ff, without its command line.
 * \author David van der Spoel <david.vanderspoel@icm.uu.se>
 */
#ifndef ACT_EDIT_FF_H
#define ACT_EDIT_FF_H

#include <optional>
#include <string>

namespace alexandria
{

class MsgHandler;

/*! \brief What alexandria edit_ff does, as its command line gives it
 *
 * The program fills it from the command line, a program without one sets
 * the members. runEditFF checks it, see check. One thing is done, the
 * first of: with a second force field, copying, replacing or implanting
 * its parameters, or else comparing the two; checking the bounds;
 * validating; analysing or plotting; removing particle types; changing
 * parameters. The combination rules are applied after any of these.
 */
struct EditFFOptions
{
    //! The force field to read (-ff)
    std::string           forceFieldFile;
    //! A second force field, the source of parameters or for a comparison (-ff2), none if empty
    std::string           secondForceFieldFile;
    //! The force field to write if it changed or if forceWrite (-o), none if empty
    std::string           outputFile;
    //! The file for the parameters of the particle types with -ana (-dump), none if empty
    std::string           dumpFile;
    //! The force field in JSON after the changes (-json), none if empty
    std::string           jsonFile;
    //! The kind of parameter, for instance zeta (-p)
    std::string           parameter;
    //! The particle or interaction types, separated by spaces, all if empty (-a)
    std::string           particles;
    //! The new minimum (-min), none if not set
    std::optional<double> minimum;
    //! The new value (-val), none if not set
    std::optional<double> value;
    //! The new maximum (-max), none if not set
    std::optional<double> maximum;
    //! The factor to scale the values by (-scale), none if not set
    std::optional<double> scale;
    //! The new mutability, Fixed or Bounded (-mut), none if not set
    std::optional<std::string> mutability;
    //! Set (true) or clear (false) that the parameters may not be negative (-nonneg), none if not set
    std::optional<bool>   nonNegative;
    //! The factor for new bounds around the value, see -limits; a value that is not above zero is ignored
    std::optional<double> limits;
    //! Remove the particle types and the interactions that have them (-del)
    bool                  remove     = false;
    //! Change also parameters that are not mutable (-force)
    bool                  force      = false;
    //! Stretch the bounds to the values (-stretch)
    bool                  stretch    = false;
    //! Check which parameters are at their bounds (-bounds)
    bool                  bounds     = false;
    //! Check the parameters (-validate)
    bool                  validate   = false;
    //! Plot interactions (-plot), with analyze
    bool                  plot       = false;
    //! Write the force field also without changes (-write)
    bool                  forceWrite = false;
    //! What to analyse: EEM, BONDS or OTHER (-ana), none if empty
    std::string           analyze;
    //! Copy the parameters of this kind that are missing from the second force field (-copy_missing)
    std::string           copyMissing;
    //! Replace the parameters of this kind by those of the second force field (-replace)
    std::string           replace;
    //! Implant the trained parameters of this kind of the second force field (-implant)
    std::string           implant;
    //! The combination rules, separated by spaces (-cr), see CombRuleUtil
    std::string           combinationRules;

    /*! \brief Check that the options make sense
     * \throws gmx::InvalidInputError that says which option is wrong
     */
    void check() const;
};

/*! \brief Edit, check or analyse a force field, as alexandria edit_ff does
 * \param[in] msghandler For output and errors
 * \param[in] options    What to do
 * \throws gmx::InvalidInputError for options that make no sense, see
 * EditFFOptions::check, gmx::FileIOError if a file cannot be read or written
 */
void runEditFF(MsgHandler *msghandler, const EditFFOptions &options);

} // namespace alexandria

#endif

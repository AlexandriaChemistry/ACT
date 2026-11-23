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
 * Helpers for tests that keep what a tool writes as reference data: the
 * JSON output, text files with numbers, and what an edit changed in a
 * force field. Such tests fix the results of the tools as they are, such
 * that a change of the code that handles the options of a tool (or of
 * anything else) cannot change a result unnoticed.
 *
 * \author David van der Spoel <david.vanderspoel@icm.uu.se>
 */
#ifndef ACT_TESTUTILS_OUTPUT_CHECKS_H
#define ACT_TESTUTILS_OUTPUT_CHECKS_H

#include <nlohmann/json.hpp>
#include <set>
#include <string>
#include <vector>

#include "act/forcefield/forcefield.h"
#include "testutils/refdata.h"
#include "testutils/testasserts.h"

namespace alexandria
{

/*! \brief A tolerance for results that may come from another build
 *
 * The reference data are made once, and checked in builds in single and
 * double precision, with and without explicit SIMD, on several platforms.
 * The precision of the build chooses the absolute tolerance, not the type
 * of the value, since the values are doubles in both.
 * \param[in] singleAbsolute The difference allowed in single precision
 * \param[in] doubleAbsolute The difference allowed in double precision
 * \param[in] relative       The relative difference allowed in both
 * \return The tolerance
 */
gmx::test::FloatingPointTolerance buildTolerance(double singleAbsolute,
                                                 double doubleAbsolute,
                                                 double relative);

/*! \brief Compare a JSON value with the reference data
 *
 * Objects become compounds with a part per key, in the order of the
 * keys, arrays become sequences. All numbers are compared as floating
 * point numbers, with the default tolerance of the checker, such that an
 * integer in one build and a real number in another compare equal.
 * Strings and booleans must be the same.
 * \param[in] checker The checker to add the value to
 * \param[in] value   The value
 * \param[in] id      The name of the value, nullptr in a sequence
 * \param[in] skip    Keys of objects, at any depth, that are not compared,
 *                    for instance file names and times
 */
void checkJson(gmx::test::TestReferenceChecker *checker,
               const nlohmann::json            &value,
               const char                      *id,
               const std::set<std::string>     &skip = {});

/*! \brief Compare a text file with the reference data, numbers within the
 * default tolerance of the checker
 *
 * \param[in] checker  The checker
 * \param[in] fileName The file
 * \param[in] skip     Lines that contain one of these texts are left out,
 *                     for instance those with a time stamp. A negative
 *                     zero, such as -0.000, counts as a zero, and a run
 *                     of spaces as one space.
 * \param[in] id       The name of the contents in the reference data
 */
void checkTextFile(gmx::test::TestReferenceChecker *checker,
                   const std::string               &fileName,
                   const std::vector<std::string>  &skip,
                   const char                      *id);

/*! \brief What differs between two force fields, as JSON
 *
 * For each particle type and each interaction: what was added, what was
 * removed, and for each parameter that differs its value, bounds,
 * mutability, number of training points, unit and whether it must not be
 * negative, before and after. Also the options of the force field, of the
 * particle types and of the interactions. Equal force fields give an
 * empty object.
 * \param[in] before The force field before
 * \param[in] after  The force field after
 * \return The differences
 */
nlohmann::json forceFieldChanges(const ForceField &before,
                                 const ForceField &after);

/*! \brief Read a text file without its comments
 * \param[in] fileName The file
 * \return The lines that do not start with #, for instance those of an
 * xvg file without the time and the command line it was written with
 */
std::string textWithoutComments(const std::string &fileName);

/*! \brief Read a JSON file
 * \param[in] fileName The file
 * \return The contents, a discarded value if the file is no JSON
 */
nlohmann::json readJsonOutput(const std::string &fileName);

} // namespace alexandria

#endif

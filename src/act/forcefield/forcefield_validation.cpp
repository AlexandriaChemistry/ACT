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
 * Checks of the parameters of a force field, see forcefield_validation.h.
 * \author David van der Spoel <david.vanderspoel@icm.uu.se>
 */
#include "forcefield_validation.h"

#include "act/basics/interactiontype.h"
#include "act/basics/mutability.h"
#include "act/forcefield/forcefield.h"
#include "gromacs/utility/stringutil.h"

namespace alexandria
{

/*! \brief Check one parameter
 * \param[in]    where      The interaction type, or "particle type"
 * \param[in]    identifier The identifier of the entry
 * \param[in]    type       The type of the parameter
 * \param[in]    p          The parameter
 * \param[inout] issues     To add the problems to
 */
static void checkParameter(const std::string            &where,
                           const std::string            &identifier,
                           const std::string            &type,
                           const ForceFieldParameter    &p,
                           std::vector<ForceFieldIssue> *issues)
{
    auto add = [&](ForceFieldIssue::Kind kind, const std::string &text, bool repaired)
    {
        issues->push_back({ kind, where, identifier, type, text, repaired });
    };
    if (p.minimum() > p.maximum())
    {
        add(ForceFieldIssue::Kind::MinimumAboveMaximum,
            gmx::formatString("minimum %g is larger than maximum %g", p.minimum(), p.maximum()), false);
    }
    if (p.wasReset())
    {
        add(ForceFieldIssue::Kind::ResetToBounds,
            gmx::formatString("value %g was outside [%g, %g] and was set to %g",
                              p.valueBeforeReset(), p.minimum(), p.maximum(), p.value()), true);
    }
    else if (p.value() < p.minimum() || p.value() > p.maximum())
    {
        add(ForceFieldIssue::Kind::OutOfBounds,
            gmx::formatString("value %g is outside [%g, %g], not changed since the mutability is %s",
                              p.value(), p.minimum(), p.maximum(), mutabilityName(p.mutability()).c_str()), false);
    }
    // Only a negative value counts, many force fields have parameters that
    // are marked non-negative with a negative minimum, which only prevents
    // setting a negative minimum, see ForceFieldParameter::setMinimum
    if (p.nonNegative() && p.value() < 0)
    {
        add(ForceFieldIssue::Kind::Negative,
            gmx::formatString("marked non-negative, but the value is %g; "
                              "the value or the attribute nonnegative is wrong",
                              p.value()), false);
    }
}

std::vector<ForceFieldIssue> validateForceField(const ForceField &pd)
{
    std::vector<ForceFieldIssue> issues;
    const std::string particleType("particle type");
    for (const auto &pt : pd.particleTypesConst())
    {
        for (const auto &p : pt.second.parametersConst())
        {
            checkParameter(particleType, pt.first.id(), p.first, p.second, &issues);
        }
        for (const auto &type : pt.second.duplicates())
        {
            issues.push_back({ ForceFieldIssue::Kind::Duplicate, particleType, pt.first.id(), type,
                               "defined more than once, the first definition was used", true });
        }
    }
    for (const auto &fs : pd.forcesConst())
    {
        const std::string &where = interactionTypeToString(fs.first);
        for (const auto &entry : fs.second.parametersConst())
        {
            for (const auto &p : entry.second)
            {
                checkParameter(where, entry.first.id(), p.first, p.second, &issues);
            }
        }
        for (const auto &dup : fs.second.duplicates())
        {
            issues.push_back({ ForceFieldIssue::Kind::Duplicate, where, dup.first.id(), dup.second,
                               "defined more than once, the first definition was used", true });
        }
    }
    return issues;
}

std::string repairedIssuesSummary(const std::vector<ForceFieldIssue> &issues)
{
    size_t nreset = 0, nduplicate = 0;
    for (const auto &issue : issues)
    {
        nreset     += issue.kind == ForceFieldIssue::Kind::ResetToBounds ? 1 : 0;
        nduplicate += issue.kind == ForceFieldIssue::Kind::Duplicate ? 1 : 0;
    }
    if (nreset + nduplicate == 0)
    {
        return std::string();
    }
    return gmx::formatString("%zu parameter(s) were outside their bounds and were set to the nearest bound, "
                             "and %zu were defined more than once, of which the first definition was used. "
                             "alexandria edit_ff -validate lists them and writes a corrected force field with -o.",
                             nreset, nduplicate);
}

} // namespace alexandria

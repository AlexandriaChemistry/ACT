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
 * Implements the helpers of output_checks.h.
 *
 * \author David van der Spoel <david.vanderspoel@icm.uu.se>
 */
#include "actpre.h"

#include <cstdint>
#include <fstream>
#include <map>
#include <regex>
#include <sstream>

#include "act/basics/interactiontype.h"
#include "act/basics/mutability.h"
#include "act/forcefield/potential.h"
#include "output_checks.h"

namespace alexandria
{

gmx::test::FloatingPointTolerance buildTolerance(double singleAbsolute,
                                                 double doubleAbsolute,
                                                 double relative)
{
    // The checker takes the tolerance for the type of the value, and the
    // values are doubles also in single precision, so the precision of
    // the build chooses
    const double absolute = GMX_DOUBLE ? doubleAbsolute : singleAbsolute;
    return gmx::test::FloatingPointTolerance(absolute, absolute,
                                             relative, relative,
                                             UINT64_MAX, UINT64_MAX, false);
}

void checkJson(gmx::test::TestReferenceChecker *checker,
               const nlohmann::json            &value,
               const char                      *id,
               const std::set<std::string>     &skip)
{
    if (value.is_object())
    {
        auto compound = checker->checkCompound("Object", id);
        for (const auto &item : value.items())
        {
            if (skip.count(item.key()) == 0)
            {
                checkJson(&compound, item.value(), item.key().c_str(), skip);
            }
        }
    }
    else if (value.is_array())
    {
        auto sequence = checker->checkSequenceCompound(id, value.size());
        for (const auto &item : value)
        {
            checkJson(&sequence, item, nullptr, skip);
        }
    }
    else if (value.is_number())
    {
        checker->checkDouble(value.get<double>(), id);
    }
    else if (value.is_boolean())
    {
        checker->checkBoolean(value.get<bool>(), id);
    }
    else if (value.is_string())
    {
        checker->checkString(value.get<std::string>(), id);
    }
    else
    {
        checker->checkString("null", id);
    }
}

void checkTextFile(gmx::test::TestReferenceChecker *checker,
                   const std::string               &fileName,
                   const std::vector<std::string>  &skip,
                   const char                      *id)
{
    std::ifstream in(fileName);
    if (!in.good())
    {
        ADD_FAILURE() << "Cannot read " << fileName;
        return;
    }
    std::ostringstream kept;
    std::string        line;
    // A minus before a number that is zero, such as -0 or -0.000
    const std::regex   negativeZero("(^|[^0-9.eE])-(0*\\.?0+)(?![0-9.eE])");
    // Spaces and tabs
    const std::regex   spaces("[ \\t]+");
    while (std::getline(in, line))
    {
        bool leaveOut = false;
        for (const auto &text : skip)
        {
            leaveOut = leaveOut || line.find(text) != std::string::npos;
        }
        if (!leaveOut)
        {
            // A zero is written as -0.000 on one platform and as 0.000 on
            // another, which the comparison of the numbers takes for
            // different texts. Every space is a token of its own, so that
            // the width of a column would count as well: runs of spaces
            // are one space.
            const auto zero = std::regex_replace(line, negativeZero, "$1$2");
            kept << std::regex_replace(zero, spaces, " ") << '\n';
        }
    }
    checker->checkTextBlockWithNumbers(kept.str(), id);
}

namespace
{

/*! \brief A force field parameter as JSON
 * \param[in] parameter The parameter
 * \return Its properties
 */
nlohmann::json parameterJson(const ForceFieldParameter &parameter)
{
    return {
        { "value", parameter.value() },
        { "minimum", parameter.minimum() },
        { "maximum", parameter.maximum() },
        { "mutability", mutabilityName(parameter.mutability()) },
        { "ntrain", parameter.ntrain() },
        { "unit", parameter.unit() },
        { "nonnegative", parameter.nonNegative() }
    };
}

/*! \brief Whether two parameters are the same in all that is written to a force field file
 * \param[in] a One parameter
 * \param[in] b Another parameter
 * \return Whether they are the same
 */
bool sameParameter(const ForceFieldParameter &a, const ForceFieldParameter &b)
{
    return parameterJson(a) == parameterJson(b);
}

/*! \brief The differences between two sets of parameters
 * \param[in] before The parameters before, by name
 * \param[in] after  The parameters after, by name
 * \return For each name that differs, the parameter before and after;
 * null where there is none
 */
template <class Map>
nlohmann::json parameterChanges(const Map &before, const Map &after)
{
    auto changes = nlohmann::json::object();
    for (const auto &b : before)
    {
        const auto a = after.find(b.first);
        if (a == after.end())
        {
            changes[b.first] = { { "before", parameterJson(b.second) }, { "after", nullptr } };
        }
        else if (!sameParameter(b.second, a->second))
        {
            changes[b.first] = { { "before", parameterJson(b.second) }, { "after", parameterJson(a->second) } };
        }
    }
    for (const auto &a : after)
    {
        if (before.find(a.first) == before.end())
        {
            changes[a.first] = { { "before", nullptr }, { "after", parameterJson(a.second) } };
        }
    }
    return changes;
}

/*! \brief The differences between two maps of options
 * \param[in] before The options before
 * \param[in] after  The options after
 * \return For each option that differs, its value before and after
 */
nlohmann::json optionChanges(const std::map<std::string, std::string> &before,
                             const std::map<std::string, std::string> &after)
{
    auto changes = nlohmann::json::object();
    for (const auto &b : before)
    {
        const auto a = after.find(b.first);
        if (a == after.end() || a->second != b.second)
        {
            changes[b.first] = { { "before", b.second },
                                 { "after", a == after.end() ? nlohmann::json(nullptr) : nlohmann::json(a->second) } };
        }
    }
    for (const auto &a : after)
    {
        if (before.find(a.first) == before.end())
        {
            changes[a.first] = { { "before", nullptr }, { "after", a.second } };
        }
    }
    return changes;
}

/*! \brief Add the changes of one thing to a list of changes, if there are any
 * \param[in,out] changes The list, an object
 * \param[in]     key     The name of the thing
 * \param[in]     change  Its changes
 */
void addChange(nlohmann::json *changes, const std::string &key, const nlohmann::json &change)
{
    if (!change.empty())
    {
        (*changes)[key] = change;
    }
}

} // namespace

nlohmann::json forceFieldChanges(const ForceField &before,
                                 const ForceField &after)
{
    auto changes = nlohmann::json::object();

    // The properties of the whole force field
    auto general = nlohmann::json::object();
    if (before.polarizable() != after.polarizable())
    {
        general["polarizable"] = { { "before", before.polarizable() }, { "after", after.polarizable() } };
    }
    if (before.nexcl() != after.nexcl())
    {
        general["nexcl"] = { { "before", before.nexcl() }, { "after", after.nexcl() } };
    }
    if (before.epsilonr() != after.epsilonr())
    {
        general["epsilonr"] = { { "before", before.epsilonr() }, { "after", after.epsilonr() } };
    }
    if (before.chargeGenerationAlgorithm() != after.chargeGenerationAlgorithm())
    {
        general["chargegeneration"] = { { "before", static_cast<int>(before.chargeGenerationAlgorithm()) },
                                        { "after", static_cast<int>(after.chargeGenerationAlgorithm()) } };
    }
    addChange(&general, "options", optionChanges(before.forceFieldOptions(), after.forceFieldOptions()));
    addChange(&changes, "forcefield", general);

    // The particle types
    auto types = nlohmann::json::object();
    for (const auto &b : before.particleTypesConst())
    {
        const auto a = after.particleTypesConst().find(b.first);
        if (a == after.particleTypesConst().end())
        {
            types[b.first.id()] = "removed";
            continue;
        }
        auto change = nlohmann::json::object();
        addChange(&change, "options", optionChanges(b.second.optionsConst(), a->second.optionsConst()));
        addChange(&change, "parameters", parameterChanges(b.second.parametersConst(), a->second.parametersConst()));
        addChange(&types, b.first.id(), change);
    }
    for (const auto &a : after.particleTypesConst())
    {
        if (before.particleTypesConst().find(a.first) == before.particleTypesConst().end())
        {
            types[a.first.id()] = { { "added", parameterChanges(ForceFieldParameterMap(),
                                                                a.second.parametersConst()) } };
        }
    }
    addChange(&changes, "particletypes", types);

    // The interactions
    auto interactions = nlohmann::json::object();
    for (const auto &b : before.forcesConst())
    {
        const auto &name = interactionTypeToString(b.first);
        const auto  a    = after.forcesConst().find(b.first);
        if (a == after.forcesConst().end())
        {
            interactions[name] = "removed";
            continue;
        }
        auto change = nlohmann::json::object();
        if (b.second.potential() != a->second.potential())
        {
            change["potential"] = { { "before", potentialToString(b.second.potential()) },
                                    { "after", potentialToString(a->second.potential()) } };
        }
        addChange(&change, "options", optionChanges(b.second.option(), a->second.option()));
        auto identifiers = nlohmann::json::object();
        const auto &pb = b.second.parametersConst();
        const auto &pa = a->second.parametersConst();
        for (const auto &ib : pb)
        {
            const auto ia = pa.find(ib.first);
            if (ia == pa.end())
            {
                identifiers[ib.first.id()] = "removed";
            }
            else
            {
                addChange(&identifiers, ib.first.id(), parameterChanges(ib.second, ia->second));
            }
        }
        for (const auto &ia : pa)
        {
            if (pb.find(ia.first) == pb.end())
            {
                identifiers[ia.first.id()] = { { "added", parameterChanges(ForceFieldParameterMap(),
                                                                           ia.second) } };
            }
        }
        addChange(&change, "parameters", identifiers);
        addChange(&interactions, name, change);
    }
    for (const auto &a : after.forcesConst())
    {
        if (before.forcesConst().find(a.first) == before.forcesConst().end())
        {
            // The whole interaction, with all its parameters
            auto added = nlohmann::json::object();
            added["potential"] = potentialToString(a.second.potential());
            added["options"]   = optionChanges({}, a.second.option());
            auto identifiers   = nlohmann::json::object();
            for (const auto &ia : a.second.parametersConst())
            {
                identifiers[ia.first.id()] = parameterChanges(ForceFieldParameterMap(), ia.second);
            }
            added["parameters"] = identifiers;
            interactions[interactionTypeToString(a.first)] = { { "added", added } };
        }
    }
    addChange(&changes, "interactions", interactions);
    return changes;
}

std::string textWithoutComments(const std::string &fileName)
{
    std::ifstream in(fileName);
    std::string   text;
    for (std::string line; std::getline(in, line); )
    {
        if (line.empty() || line[0] != '#')
        {
            text += line + "\n";
        }
    }
    return text;
}

nlohmann::json readJsonOutput(const std::string &fileName)
{
    std::ifstream in(fileName);
    return nlohmann::json::parse(in, nullptr, false);
}

} // namespace alexandria

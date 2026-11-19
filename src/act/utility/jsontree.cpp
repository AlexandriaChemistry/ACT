/*
 * This source file is part of the Alexandria Chemistry Toolkit.
 *
 * Copyright (C) 2022,2025,2026
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
 * JSON serialization uses the nlohmann/json library.
 * \author David van der Spoel <david.vanderspoel@icm.uu.se>
 */

#include <cerrno>
#include <cmath>
#include <cstddef>
#include <cstdio>
#include <cstdlib>
#include <nlohmann/json.hpp>
#include <nlohmann/json_fwd.hpp>
#include <string>

#include "gromacs/utility/futil.h"
#include "gromacs/utility/stringutil.h"
#include "jsontree.h"

namespace alexandria
{

/*! \brief Remove spaces at the start and the end of a text
 * \param[in] text The text
 * \return the text without the spaces around it
 */
static std::string stripSpaces(const std::string &text)
{
    const char  *space = " \t\n\r";
    const size_t first = text.find_first_not_of(space);
    if (first == std::string::npos)
    {
        return "";
    }
    return text.substr(first, text.find_last_not_of(space) - first + 1);
}

/*! \brief Interpret a text as a number
 * \param[in]  text    The text, without spaces around it
 * \param[out] integer Set if the text is an integer that fits in a long long
 * \param[out] real    Set if the text is a number, also for integers
 * \param[out] isInt   Whether the text is an integer
 * \return whether the whole text is a number. That includes nan and inf.
 */
static bool textToNumber(const std::string &text,
                         long long         *integer,
                         double            *real,
                         bool              *isInt)
{
    *isInt = false;
    if (text.empty())
    {
        return false;
    }
    char *end = nullptr;
    errno     = 0;
    *integer  = std::strtoll(text.c_str(), &end, 10);
    if (*end == '\0' && errno == 0)
    {
        *isInt = true;
        *real  = static_cast<double>(*integer);
        return true;
    }
    end   = nullptr;
    *real = std::strtod(text.c_str(), &end);
    return *end == '\0';
}

std::string jsonNumberText(double value)
{
    if (!std::isfinite(value))
    {
        return std::isnan(value) ? "nan" : "inf";
    }
    // Fifteen digits are always kept by a double, seventeen always suffice.
    for (int digits = 15; digits < 17; digits++)
    {
        const std::string text = gmx::formatString("%.*g", digits, value);
        if (std::strtod(text.c_str(), nullptr) == value)
        {
            return text;
        }
    }
    return gmx::formatString("%.17g", value);
}

std::string jsonNumberText(float value)
{
    if (!std::isfinite(value))
    {
        return std::isnan(value) ? "nan" : "inf";
    }
    // Six digits are always kept by a float, nine always suffice.
    for (int digits = 6; digits < 9; digits++)
    {
        const std::string text = gmx::formatString("%.*g", digits, static_cast<double>(value));
        if (std::strtof(text.c_str(), nullptr) == value)
        {
            return text;
        }
    }
    return gmx::formatString("%.9g", static_cast<double>(value));
}

void JsonTree::addValueUnit(const std::string &key,
                            const std::string &value,
                            const std::string &unit)
{
    JsonTree  jt(key);
    long long integer = 0;
    double    real    = 0;
    bool      isInt   = false;
    if (textToNumber(stripSpaces(value), &integer, &real, &isInt))
    {
        jt.addNumber("value", value);
    }
    else
    {
        jt.addObject("value", value);
    }
    jt.addObject("unit", unit);
    objects_.push_back(jt);
}

/*! \brief Convert the text of a number to a JSON value
 * \param[in] value The text
 * \return a JSON integer or real number, null for a number that is not
 * finite and a string if the text is not a number after all.
 */
static nlohmann::ordered_json numberToNlohmann(const std::string &value)
{
    long long integer = 0;
    double    real    = 0;
    bool      isInt   = false;
    if (!textToNumber(stripSpaces(value), &integer, &real, &isInt))
    {
        return value;
    }
    if (isInt)
    {
        return integer;
    }
    if (!std::isfinite(real))
    {
        return nullptr;
    }
    return real;
}

/*! \brief Convert a JsonTree node's content to an nlohmann::ordered_json value.
 * If the node has a value, it becomes a JSON string or a JSON number.
 * If the node has children, they become members of a JSON object
 * preserving insertion order, or the elements of a JSON array.
 * If both value and children are present, a warning is printed.
 * \param[in]  tree  The JsonTree node to convert
 * \return An nlohmann::ordered_json value representing this node's content
 */
static nlohmann::ordered_json toNlohmannValue(const JsonTree &tree)
{
    if (!tree.value().empty())
    {
        if (!tree.objects().empty())
        {
            fprintf(stderr, "JsonTree object %s has both a value %s and %zu branches. Please fix code.\n",
                    tree.key().c_str(), tree.value().c_str(), tree.objects().size());
        }
        if (tree.isNumber())
        {
            return numberToNlohmann(tree.value());
        }
        if (tree.isBool())
        {
            return tree.value() == "true";
        }
        return tree.value();
    }
    else if (tree.isArray())
    {
        // Also when it is empty
        nlohmann::ordered_json arr = nlohmann::ordered_json::array();
        for (const auto &child : tree.objects())
        {
            arr.push_back(toNlohmannValue(child));
        }
        return arr;
    }
    else if (!tree.objects().empty())
    {
        nlohmann::ordered_json obj = nlohmann::ordered_json::object();
        for (const auto &child : tree.objects())
        {
            obj[child.key()] = toNlohmannValue(child);
        }
        return obj;
    }
    else
    {
        return nullptr;
    }
}

//! Add n spaces to a string.
static void add_space(std::string *str, int n)
{
    for(int i = 0; i < n; i++)
    {
        *str += " ";
    }
}

const std::string JsonTree::writeString(bool json, int *indent) const
{
    if (json)
    {
        nlohmann::ordered_json doc = nlohmann::ordered_json::object();
        doc[key_] = toNlohmannValue(*this);
        return doc.dump(*indent);
    }
    else
    {
        std::string str;
        add_space(&str, *indent);
        // The elements of an array have no key
        if (!key_.empty())
        {
            str += gmx::formatString("%s:", key_.c_str());
        }
        *indent += 2;
        if (!value_.empty())
        {
            str += gmx::formatString("%s%s", key_.empty() ? "" : " ", value_.c_str());
        }
        else
        {
            size_t count = 0;
            bool leaf = true;
            for (const auto &obj : objects_)
            {
                leaf = leaf && !obj.branched();
            }
            if (!leaf)
            {
                str += "\n";
            }
            for (const auto &obj : objects_)
            {
                str += obj.writeString(json, indent);
                if (count < objects_.size()-1)
                {
                    str += ",";
                }
                if (leaf)
                {
                    if (count < objects_.size()-1)
                    {
                        str += " ";
                    }
                }
                else
                {
                    str += "\n";
                }
                count++;
            }
            if (!leaf)
            {
                add_space(&str, *indent);
            }
        }
        *indent -= 2;
        return str;
    }
}
        
void JsonTree::write(const std::string &fileName,
                     bool               json,
                     int                indent)
{
    FILE *fp   = gmx_ffopen(fileName, "w");
    fprintf(fp, "%s", writeString(json, &indent).c_str());
    gmx_ffclose(fp);
}
        
} // namespace alexandria

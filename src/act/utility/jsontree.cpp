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

#include "jsontree.h"

#include <nlohmann/json.hpp>

#include "gromacs/utility/futil.h"
#include "gromacs/utility/stringutil.h"

namespace alexandria
{

/*! \brief Convert a JsonTree node's content to an nlohmann::ordered_json value.
 * If the node has a value, it becomes a JSON string.
 * If the node has children, they become members of a JSON object
 * preserving insertion order.
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
        return tree.value();
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
        return doc.dump();
    }
    else
    {
        std::string str;
        add_space(&str, *indent);
        str += gmx::formatString("%s:", key_.c_str());
        *indent += 2;
        if (!value_.empty())
        {
            str += gmx::formatString(" %s", value_.c_str());
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
                     bool               json)
{
    FILE *fp   = gmx_ffopen(fileName.c_str(), "w");
    int indent = 0;
    fprintf(fp, "%s", writeString(json, &indent).c_str());
    gmx_ffclose(fp);
}
        
} // namespace alexandria

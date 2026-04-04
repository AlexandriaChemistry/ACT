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
 * JSON serialization uses the RapidJSON library.
 * \author David van der Spoel <david.vanderspoel@icm.uu.se>
 */

#include "jsontree.h"

#include <rapidjson/document.h>
#include <rapidjson/stringbuffer.h>
#include <rapidjson/writer.h>

#include "gromacs/utility/futil.h"
#include "gromacs/utility/stringutil.h"

namespace alexandria
{

/*! \brief Convert a JsonTree node's content to a rapidjson::Value.
 * If the node has a value, it becomes a JSON string.
 * If the node has children, they become members of a JSON object.
 * If both value and children are present, a warning is printed.
 * \param[in]  tree  The JsonTree node to convert
 * \param[in]  alloc The RapidJSON allocator
 * \return A rapidjson::Value representing this node's content
 */
static rapidjson::Value toRapidJsonValue(const JsonTree                       &tree,
                                         rapidjson::Document::AllocatorType   &alloc)
{
    if (!tree.value().empty())
    {
        if (!tree.objects().empty())
        {
            fprintf(stderr, "JsonTree object %s has both a value %s and %zu branches. Please fix code.\n",
                    tree.key().c_str(), tree.value().c_str(), tree.objects().size());
        }
        rapidjson::Value v;
        v.SetString(tree.value().c_str(),
                    static_cast<rapidjson::SizeType>(tree.value().size()),
                    alloc);
        return v;
    }
    else if (!tree.objects().empty())
    {
        rapidjson::Value obj(rapidjson::kObjectType);
        for (const auto &child : tree.objects())
        {
            rapidjson::Value key;
            key.SetString(child.key().c_str(),
                          static_cast<rapidjson::SizeType>(child.key().size()),
                          alloc);
            rapidjson::Value val = toRapidJsonValue(child, alloc);
            obj.AddMember(key, val, alloc);
        }
        return obj;
    }
    else
    {
        return rapidjson::Value(rapidjson::kNullType);
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
        rapidjson::Document doc;
        doc.SetObject();
        auto &alloc = doc.GetAllocator();

        rapidjson::Value key;
        key.SetString(key_.c_str(),
                      static_cast<rapidjson::SizeType>(key_.size()),
                      alloc);
        rapidjson::Value val = toRapidJsonValue(*this, alloc);
        doc.AddMember(key, val, alloc);

        rapidjson::StringBuffer                          buffer;
        rapidjson::Writer<rapidjson::StringBuffer>       writer(buffer);
        doc.Accept(writer);
        return buffer.GetString();
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
